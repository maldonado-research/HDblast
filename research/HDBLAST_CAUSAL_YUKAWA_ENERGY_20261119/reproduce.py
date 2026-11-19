"""Authenticate the archived experiment and rerun its independent calculations."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
NATIVE = REPO / "research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native"


def load(path):
    return json.loads(path.read_bytes())


def require(test, message):
    if not test:
        raise ValueError(message)


def authenticate():
    paths = {}
    for line in (HERE / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split("  ", 1)
        path = (HERE / name).resolve()
        require(path.is_relative_to(HERE) and name not in paths, name)
        paths[name] = digest
        require(hashlib.sha256(path.read_bytes()).hexdigest() == digest, name)
    files = {p.relative_to(HERE).as_posix() for p in HERE.rglob("*")
             if p.is_file() and "__pycache__" not in p.parts and p.name != "SHA256SUMS"}
    require(files == paths.keys(), "Study manifest coverage")
    sources = load(HERE / "SOURCES.json")
    require(sources["base"] == "b7da8ef48935d5ac34892d7f692954408133576d", "Historical revision")
    require({(REPO / r["path"]).resolve() for r in sources["sources"]} ==
            {NATIVE / "FORK.json", NATIVE / "WALLS.zip"}, "Native input inventory")
    for row in sources["sources"]:
        blob = (REPO / row["path"]).read_bytes()
        require(len(blob) == row["bytes"] and hashlib.sha256(blob).hexdigest() == row["sha256"],
                row["path"])
    print(f"Verified {len(paths)} study files and both native inputs.", flush=True)


def compare(actual, expected, field=""):
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and expected.keys() == actual.keys(), field)
        for key in expected:
            compare(actual[key], expected[key], field + "/" + key)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), field)
        for i, (a, b) in enumerate(zip(actual, expected, strict=True)):
            compare(a, b, field + "/" + str(i))
    elif isinstance(expected, float):
        # Finite-coupling extrapolations difference large vacuum energies. Their
        # independent convergence and symplectic tests run inside check/verify.py.
        require(math.isclose(actual, expected, rel_tol=3e-6, abs_tol=3e-7),
                f"{field}: {actual} vs {expected}")
    else:
        require(type(actual) is type(expected) and actual == expected, field)


def science():
    checks = load(HERE / "solver/CHECKS.json")
    require(max(checks["max_sampled_absolute_balance_residuals"].values()) < 4e-7, "Sector balance")
    budgets = load(HERE / "solver/BUDGETS.json")["runs"]
    tail = []
    for cutoff, key in [(240, "fine_wide_240_N160_T4096"), (480, "fine_wide_480_N240_T8192")]:
        x = budgets[key]["retain"]
        b, c = [x[k]["checkpoints"]["cut"] for k in ("bare", "static")]
        energy = b["energy"]["daughter"]
        require(math.isclose(energy, c["energy"]["daughter"], abs_tol=1e-12), str(cutoff))
        require(abs(b["energy"]["parent"] + energy) > 1 and b["work"]["switch"] < 0, key)
        tail.append((energy, c["work"]["switch"] + c["work"]["switch_counterterm"]))
    require(abs(tail[1][0]-tail[0][0]) < 5e-7 and tail[1][1]-tail[0][1] > 1.9, "Cutoff test")
    pulses = load(HERE / "pulses/PULSES.json")["runs"]
    for run in pulses.values():
        for arm, signs in run["pulse_minus_quiet"].items():
            for scheme, row in signs.items():
                require((row["daughter"] < 0) == (arm == "plus"), "Pulse signs")
                require(abs(row["independent_switch"] - row["switch"]) < 6e-7, "Pulse work")
                require(abs(row["closure"]) < 6e-7, scheme)
    clocks = load(HERE / "clock/CLOCK.json")["grids"]
    for name, grid in clocks.items():
        arms = grid["arms"]
        flux = pulses[name + "_120"]
        baseline = flux["quiet"]["bare"]["daughter"]
        for arm, row in arms.items():
            delta = 0. if arm == "quiet" else flux["pulse_minus_quiet"][arm]["bare"]["daughter"]
            require(abs(row["energy"]["120"]["proper_time"] - baseline - delta) < 2e-7,
                    "Independent proper-time arm energy")
        for cutoff, quiet in arms["quiet"]["energy"].items():
            require(abs(quiet["proper_time"] - quiet["scalar"]) < 1e-9, "Quiet clocks")
            changes = [arms[arm]["energy"][cutoff]["change_of_pulse_contrast_pp"]
                       for arm in ("plus", "minus", "far")]
            require(changes[0] < -.20 and changes[1] > .55 and changes[2] < -.08,
                    "Fixed scalar-clock predictions")
    state = load(HERE / "memory/STATE_DIFFERENCE.json")
    require(state["checks"]["largest_absolute_difference_causal_vs_direct"] < 1e-7,
            "Continuum state-energy comparison")
    check = load(HERE / "check/CHECKS.json")
    require(check["native_energy_comparisons"] == 30 and
            check["max_native_energy_relative_difference"] < 1e-7, "Direct spectrum")


def replay():
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")

    def run(script, *args, emit=None):
        print(f"Recomputing {script.name}", flush=True)
        cmd = [sys.executable, "-I", "-B", str(script), *map(str, args)]
        if emit is None:
            subprocess.run(cmd, env=env, check=True)
        else:
            with emit.open("wb") as out:
                subprocess.run(cmd, env=env, stdout=out, check=True)

    with tempfile.TemporaryDirectory() as folder:
        temp = Path(folder)
        for sub in ("solver", "check", "memory", "pulses", "theory", "clock"):
            (temp / sub).mkdir()
            for code in (HERE / sub).glob("*.py"):
                shutil.copyfile(code, temp / sub / code.name)
        work = temp / "solver"
        run(work / "budget.py", NATIVE, "--output", work / "BUDGETS.json")
        run(work / "check.py", work / "BUDGETS.json", "--output", work / "CHECKS.json")
        run(temp / "check/verify.py", "--native", NATIVE, "--full")
        run(temp / "memory/retarded.py", "--native", NATIVE, emit=temp / "memory/MEMORY.json")
        run(temp / "memory/state_difference.py", "--native", NATIVE, emit=temp / "memory/STATE_DIFFERENCE.json")
        run(temp / "theory/dispersion.py", NATIVE, emit=temp / "theory/DISPERSION.json")
        run(temp / "pulses/pulses.py", NATIVE, "--output", temp / "pulses/PULSES.json")
        run(temp / "clock/clock_test.py", NATIVE)
        for sub in ("solver", "check", "memory", "pulses", "theory", "clock"):
            for baseline in sorted((HERE / sub).glob("*.json")):
                relative = baseline.relative_to(HERE)
                compare(load(temp / relative), load(baseline), relative.as_posix())
    print("All numerical calculations reproduced and compared.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Recompute every numerical report")
    options = parser.parse_args()
    authenticate()
    science()
    if options.verify:
        replay()


if __name__ == "__main__":
    main()
