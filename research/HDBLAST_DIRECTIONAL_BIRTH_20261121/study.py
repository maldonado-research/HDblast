"""Verify the frozen HDblast directional-birth study and reproduce its numerical products."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
PRIOR = ROOT / "prior"
PARENT_COMMIT = "dc99ff5c5894abcff3e11e02e4b97e19da31de00"
PRODUCTS = ("science/DIRECTIONAL","audit/AUDIT","theory/HISTORY",
            "scattering/THERMAL","scattering/RESONANCE","scattering/LIVE","feedback/FEEDBACK","cosmology/BUDGET",
            "boost/SAMPLING","boost/BOOST_CHECK","boost/EVOLUTION","boost/REFINEMENT","boost/FLAVOR","scattering/CONVERGENCE")


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def authenticate(external_prior=None):
    manifest = ROOT / "SHA256SUMS"
    records = {}
    for line in manifest.read_text().splitlines():
        checksum, name = line.split("  ", 1)
        if name in records:
            raise ValueError(f"Duplicate manifest entry {name}")
        records[name] = checksum
    files = {p.relative_to(ROOT).as_posix(): p for p in ROOT.rglob("*") if p.is_file()
             and p != manifest and "__pycache__" not in p.parts}
    if files.keys() != records.keys():
        raise ValueError(f"Some files have changed: {files.keys() ^ records.keys()}")
    for name, path in files.items():
        if path.is_symlink() or digest(path) != records[name]:
            raise ValueError(f"File differs from manifest: {name}")
    previous = json.loads((PRIOR / "source.json").read_text())
    if previous["commit"] != PARENT_COMMIT:
        raise ValueError("Unexpected predecessor revision")
    for name, checksum in previous["sha256"].items():
        if digest(PRIOR / name) != checksum:
            raise ValueError(f"Study input has changed: {name}")
        if external_prior is not None and digest(external_prior / name) != checksum:
            raise ValueError(f"Preceding repository study differs: {name}")
    print(f"Verified {len(records)} files including {len(previous['sha256'])} pinned inputs.", flush=True)


def compare(value, ref, path):
    if type(value) is not type(ref):
        raise ValueError(f"Changed data type at {path}")
    if isinstance(ref, dict):
        if value.keys() != ref.keys():
            raise ValueError(f"Changed keys at {path}")
        for k, item in ref.items():
            compare(value[k], item, path + "/" + k)
    elif isinstance(ref, list):
        if len(value) != len(ref):
            raise ValueError(f"Changed list length at {path}")
        for i, (x, y) in enumerate(zip(value, ref)):
            compare(x, y, path + "/" + str(i))
    elif isinstance(ref, float):
        rel_tol, abs_tol = 1e-6, 2e-12
        field = path.rsplit("/", 1)[-1]
        if path.startswith("theory/HISTORY/runs/"):
            rel_tol = 5e-5
            if "/mass_4e12/" in path and (field.startswith(("global_hard_", "history_extra_", "one_way_"))
                                          or field == "hard_with_reset_18"):
                # The independent direct integral samples secondary grazing folds.
                # The largest compiler difference across five grids was 0.51 percent.
                rel_tol = .01
        if ((path.startswith("boost/SAMPLING/") and (field == "mean" or "/scrambles/" in path))
                or (path.startswith("boost/BOOST_CHECK/sobol/") and field in ("exact", "sobol_mean"))):
            abs_tol = 1e-30
        if field.endswith("_balance_error"):
            abs_tol = 5e-11
        if not math.isfinite(value) or not math.isclose(value, ref, rel_tol=rel_tol, abs_tol=abs_tol):
            raise ValueError(f"{path}: regenerated {value!r}, archived {ref!r}")
    elif value != ref:
        raise ValueError(f"{path}: regenerated {value!r}, archived {ref!r}")


def check(products):
    directional = products["science/DIRECTIONAL"]["runs"]
    audit = products["audit/AUDIT"]["histories"]
    for mass in ("mass_4e10", "mass_4e12"):
        for grid in ("native", "half_step"):
            result = audit[mass][grid]
            feedback_grid = "half" if grid == "half_step" else "native"
            if not (2 < directional[mass]["runs"][grid]["hot_factor"] < 2.4
                    and abs(result["coverage"] - .99998) < 1e-10):
                raise ValueError(f"Growth factor or measured birth interval changed: {mass}, {grid}")

            local = directional[mass]["runs"][grid]["action"]["extended"]
            required = local["denied_over_all_requested"]
            if abs(result["loss_per_total_requested"]-required) > 3e-8:
                raise ValueError(f"Independent directional integral differs: {mass}, {grid}")
            if not (0 <= local["pooled_denied_over_all_requested"] <= required):
                raise ValueError(f"Total capacity or directional loss changed: {mass}, {grid}")
            if mass == "mass_4e12" and not (
                    2.8e-5 < required < 3.e-5 and .24 < local["peak_denied_over_instantaneous_request"] < .27
                    and local["pooled_denied_over_all_requested"] == 0):
                raise ValueError("The integrated and directional high-mass tests no longer separate")
            repeated = products["theory/HISTORY"]["runs"][mass][grid]
            if not (repeated["global_hard_decomposed"] > required and
                    abs(repeated["total_weight"]-1) < .001):
                raise ValueError(f"Repeated-label check changed: {mass}, {grid}")
            for law, meshes in products["feedback/FEEDBACK"][mass].items():
                if feedback_grid not in meshes:
                    continue
                live = meshes[feedback_grid]
                if max(live["surface_balance_error"], live["gas_balance_error"], live["oscillator_balance_error"]) > 1e-5:
                    raise ValueError(f"Energy transfer check failed: {mass}, {grid}, {law}")
                if live["stop_turnaround"] != 1 or live["hot_radius_factor"] < 2:
                    raise ValueError(f"Hot growth interval changed: {mass}, {grid}, {law}")
    resonance = products["scattering/RESONANCE"]
    if not (abs(resonance["kernel_moments"][0]-1) < 1e-8 and abs(resonance["kernel_moments"][1]-20) < 1e-7
            and resonance["scalar"]["max_rate_over_required_at_sharp_sweep_upper"] < .7
            and all(bin_run["minimum_ratio_lower"]>1 for bin_run in resonance["vector"]["runs"].values())):
        raise ValueError("Thermal transport or resonance capacity changed")
    for case, live in products["scattering/LIVE"]["runs"].items():
        if not (live["stop_turnaround"] and
                abs(live["end_rest"]+live["lost_rest"]-live["tried_rest"]) < 1e-10):
            raise ValueError(f"Hot-interval collision continuation changed: {case}")
        if case.startswith("early_pole_alone"):
            if (live["peak_hot_fraction"] >= .60 or live["hot_radius_factor"] >= 2):
                raise ValueError(f"Early elastic channel comparison changed: {case}")
        elif live["hot_radius_factor"] < 2:
            raise ValueError(f"Hot doubling changed: {case}")
    boosted = products["boost/EVOLUTION"]["runs"]
    for case, live in boosted.items():
        if (live["stop_turnaround"] != 1 or
                max(live[k] for k in ("surface_balance_error", "gas_balance_error", "oscillator_balance_error")) > 1e-5):
            raise ValueError(f"Full-history collision or conservation check failed for {case}")
        if case.startswith("pole_") and not (2.14 < live["hot_radius_factor"] < 2.17 and
                                             .97 < live["peak_hot_fraction"] < .99):
            raise ValueError(f"Fixed-pole hot interval changed for {case}")
    if (boosted["zero_drag"]["peak_hot_fraction"] > .60 or
            boosted["zero_drag"]["hot_radius_factor"] >= 2):
        raise ValueError("The zero-drag comparison changed")
    for name, row in products["boost/FLAVOR"].items():
        lower, upper = (2.05, 2.09) if name.startswith("inclusive_") else (1.94, 1.98)
        if not (lower < row["hot_radius_factor"] < upper and row["stop_turnaround"] == 1):
            raise ValueError(f"Shared three-flavor pole changed for {name}")
    for name, row in products["boost/REFINEMENT"]["runs"].items():
        if abs(row["hot_radius_factor"]["relative_change"]) > 1e-4:
            raise ValueError(f"Dense rate-table refinement changed hot-radius growth in {name}")
    cosmology = products["cosmology/BUDGET"]
    benchmark = cosmology["runs"]["mass_4e12"]["directional_cap_native"]
    if not (1.5e88 < cosmology["photon_entropy_Hubble_sphere_in_kB"] < 1.7e88
            and 2.4e44 < benchmark["radiation_entropy_upper_in_kB"] < 2.7e44
            and benchmark["photon_entropy_over_upper"] > 6e43):
        raise ValueError("Closed-wall equilibrium radiation entropy comparison changed")


def compare_table(actual, frozen):
    # Entries are natural logs. Their absolute differences bound relative rate errors.
    with actual.open() as a, frozen.open() as b:
        header = a.readline().split()
        if header != b.readline().split():
            raise ValueError("The boosted collision table axes changed")
        rows, cols = map(int, header[:2])
        worst = 0.
        for i in range(rows):
            u, v = [float(x) for x in a.readline().split()], [float(x) for x in b.readline().split()]
            if len(u) != cols or len(v) != cols or not all(map(math.isfinite, u + v)):
                raise ValueError(f"Unexpected boosted table data in row {i}")
            worst = max(worst, max(abs(x - y) for x, y in zip(u, v)))
        if a.read().strip() or b.read().strip() or worst > 1e-7:
            raise ValueError(f"Boosted collision table mismatch, maximum log difference {worst}")
        print(f"Verified {rows * cols} boosted collision rates.", flush=True)


def replay(destination):
    destination = destination.resolve()
    if destination.is_relative_to(ROOT) or (destination.exists() and
            (not destination.is_dir() or any(destination.iterdir()))):
        raise ValueError("Choose a new or empty output directory outside the study.")
    destination.mkdir(parents=True, exist_ok=True)
    outputs = {n: destination / (n + ".json") for n in PRODUCTS}
    for path in outputs.values():
        path.parent.mkdir(exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def run(name, *args):
        print("Replaying " + name, flush=True)
        subprocess.run([sys.executable, "-I", "-B", str(ROOT / name), *map(str,args)],
                       check=True, env=env)

    def scattering():
        run("scattering/history.py","--base",PRIOR,"--output",outputs["scattering/THERMAL"])
        run("scattering/resonance.py","--history",outputs["scattering/THERMAL"],"--sweep",PRIOR/"kinetics/SWEEP.json",
            "--output",outputs["scattering/RESONANCE"])
        run("scattering/live.py","--base",PRIOR,"--jobs","2","--output",outputs["scattering/LIVE"])

    def boosted():
        grid = destination / "boost/BOOST_TABLE.dat"
        run("boost/thermal.py", "--output", grid, "--jobs", "2")
        compare_table(grid, ROOT / "boost/BOOST_TABLE.dat")
        run("boost/sampling.py", "--output", outputs["boost/SAMPLING"])
        run("boost/check.py", "--table", grid, "--sampling", outputs["boost/SAMPLING"], "--output", outputs["boost/BOOST_CHECK"])
        run("boost/evolve.py", "--base", PRIOR, "--table", grid, "--jobs", "2", "--output", outputs["boost/EVOLUTION"])
        run("scattering/convergence.py", "--output", outputs["scattering/CONVERGENCE"])
        run("boost/flavor.py", "--base", PRIOR, "--table", grid, "--output", outputs["boost/FLAVOR"], "--jobs", "2")
        run("boost/refine.py", "--base", PRIOR, "--coarse", outputs["boost/EVOLUTION"],
            "--output", outputs["boost/REFINEMENT"], "--jobs", "2")

    with tempfile.TemporaryDirectory(prefix="hdblast-directional-history-") as work:
        native = Path(work)
        with ThreadPoolExecutor(max_workers=3) as pool:
            rates = pool.submit(scattering)
            boost = pool.submit(boosted)
            feedback = pool.submit(run,"feedback/feedback.py","--base",PRIOR,"--jobs","2",
                                   "--output",outputs["feedback/FEEDBACK"])
            run("science/directional.py","--prior",PRIOR,"--jobs","2","--native-output",native,
                "--output",outputs["science/DIRECTIONAL"])
            angular_audit = pool.submit(run,"audit/audit.py","--native",native,"--prior",PRIOR,
                                        "--output",outputs["audit/AUDIT"])
            run("theory/history.py","--traces",native / "traces.json","--output",outputs["theory/HISTORY"])
            angular_audit.result()
            feedback.result()
            rates.result()
            boost.result()
    run("cosmology/budget.py","--base",PRIOR,"--feedback",outputs["feedback/FEEDBACK"], "--evolution",outputs["boost/EVOLUTION"],
        "--output",outputs["cosmology/BUDGET"])
    values = {name: json.loads(path.read_bytes()) for name, path in outputs.items()}
    check(values)
    for name, obj in values.items():
        compare(obj, json.loads((ROOT / (name + ".json")).read_bytes()), name)
    print(f"All {len(values)} numerical products and physical checks passed.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-prior", type=Path, help="Compare the vendored files with the preceding repository study")
    parser.add_argument("--verify", action="store_true", help="Recompute the numerical outputs")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if args.output_dir and not args.verify:
        parser.error("--output-dir requires --verify")
    authenticate(args.check_prior)
    if args.verify and args.output_dir:
        replay(args.output_dir)
    elif args.verify:
        with tempfile.TemporaryDirectory(prefix="hdblast-directional-replay-") as temp:
            replay(Path(temp))
