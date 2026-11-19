"""Verify the archived two-pulse study and optionally replay the numerical experiments."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
BASE = "44370269720415893dcd8d5bc802697e629d7655"
PINS_SHA = "07c490004870506ab38fa68cb61c3f48735cf70476ef195fc7fabadc2b12dcc9"
ARMS = {"quiet", "plus", "minus", "far"}


def require(condition, description):
    if not condition:
        raise ValueError(description)


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def authenticate(repo):
    manifest = HERE / "SHA256SUMS"
    files = set()
    for line in manifest.read_text().splitlines():
        expected, name = line.split("  ", 1)
        path = (HERE / name).resolve()
        require(path.is_relative_to(HERE) and path != manifest and path not in files,
                f"Invalid or duplicate manifest entry {name}")
        require(digest(path) == expected, f"Study file checksum mismatch at {name}")
        files.add(path)
    actual = {p.resolve() for p in HERE.rglob("*") if p.is_file() and p != manifest
              and "__pycache__" not in p.relative_to(HERE).parts}
    require(actual == files, f"Capsule inventory differs: {actual ^ files}")
    pins = HERE / "SOURCES.json"
    require(digest(pins) == PINS_SHA, "Historical manifest checksum")
    manifest_data = read(pins)
    require(manifest_data["base"] == BASE, "Historical revision")
    inputs = set()
    for entry in manifest_data["sources"]:
        path = (repo / entry["path"]).resolve()
        require(path.is_relative_to(repo) and path not in inputs, f"Invalid input {path}")
        require(path.stat().st_size == entry["bytes"] and digest(path) == entry["sha256"],
                f"Historical input size or checksum mismatch at {entry['path']}")
        inputs.add(path)
    print(f"Authenticated {len(files)} study files and {len(inputs)} historical inputs.", flush=True)


def tolerances(path):
    group = path[0]
    if group == "family":
        if "coefficients" in path:
            # The objective is nearly flat near some maxima. The achieved gain
            # and its independently generated certificate are both checked below.
            return 1e-7, 2e-7
        return 0., 2e-12
    if group in ("budget", "controls"):
        if any(k in path for k in ("direct_opposite_sign_maximum_normalized_residuals", "max_pointwise_ward")):
            return 1e-6, 2e-9
        return 2e-7, 2e-6
    if "archived_A_max_relative_error" in path:
        # Residual against a fixed spectrum, obtained with a different solver.
        return 0., 1e-9
    if group == "emission" and any(k in path for k in ("real", "imag")):
        # Off-diagonal overlaps can nearly cancel. Population and energy
        # integrals, and the optimized small gain, are compared separately.
        return 5e-8, 1e-11
    return 3e-8, 7e-13


def compare(actual, expected, path):
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and expected.keys() == actual.keys(), f"Field names at {path}")
        for k, value in expected.items():
            compare(actual[k], value, path + (k,))
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), f"Array size at {path}")
        for i, (new, old) in enumerate(zip(actual, expected)):
            compare(new, old, path + (str(i),))
    elif isinstance(expected, float):
        rel, absolute = tolerances(path)
        require(type(actual) in (int, float) and math.isfinite(actual)
                and math.isclose(actual, expected, rel_tol=rel, abs_tol=absolute),
                f"At {'/'.join(path)} computed {actual}, expected {expected}")
    elif path[0] == "family" and path[-1] == "boxes":
        # This is the work used by the exact rational subdivision certificate.
        # The subdivision count can depend on the floating-point proposed maximum.
        require(type(actual) is int and actual > 0, f"Certificate subdivisions at {path}")
    else:
        require(type(actual) is type(expected) and actual == expected, f"Value at {path}")


def spectra(path):
    totals = {}
    n_rows = 0
    with path.open(newline="") as rebuilt, (HERE / "emission/SPECTRA.csv").open(newline="") as saved:
        new, old = csv.DictReader(rebuilt), csv.DictReader(saved)
        require(new.fieldnames == old.fieldnames, "Daughter CSV column names")
        for a, b in zip(new, old, strict=True):
            key = (a["arm"], a["source"])
            require(key == (b["arm"], b["source"]), f"CSV identity at row {n_rows}")
            for field in new.fieldnames[2:]:
                x, y = float(a[field]), float(b[field])
                require(math.isfinite(x) and math.isclose(x, y, rel_tol=3e-8, abs_tol=3e-14),
                        f"CSV row {n_rows}/{field}: computed {x}, expected {y}")
            moments = totals.setdefault(key, [0., 0.])
            moments[0] += float(a["one_sign_comoving_number"])
            moments[1] += float(a["one_sign_conformal_energy"])
            n_rows += 1
    require(n_rows == 4 * 13 * 64, "Daughter histogram coverage")
    return totals


def quadratic(matrix, vector):
    return math.fsum(x * matrix[i][j] * y for i, x in enumerate(vector)
                     for j, y in enumerate(vector))


def physical_checks(emission, criteria, stationary, family, budget, controls, gravity, sums):
    require(set(emission["arms"]) == set(criteria["histories"]) == set(stationary) == ARMS,
            "Coverage of the four shell histories")
    require(gravity["native_runs"] == 18 and gravity["upstream_revision"] == BASE, "Native gravity runs")
    for arm in sorted(ARMS):
        row = emission["arms"][arm]
        a = row["scale_factor"]
        cutoff120 = row["cutoffs"]["120"]
        bound = criteria["histories"][arm]
        inv3 = cutoff120["inverse_cubed_pair_mean_energy_matrix"]["real"][0][0]
        require(inv3 > 0 and math.isclose(bound["initial_pair_mean_energy_inverse_third"], inv3,
                                           rel_tol=5e-8), f"{arm} inverse-third-moment cross-check")
        cont = bound["bounds"]["whole_post_gate"]["continuum"]
        require(0 < bound["num_to_chemical_initial"] < cont["lowest_upper"] < .223,
                f"{arm} continuum gate bound")
        profiles = cutoff120["profiles"]
        for label, report in profiles.items():
            n, e = sums[arm, label]
            require(math.isclose(n, report["pair_count"], rel_tol=3e-8)
                    and math.isclose(2*e, report["both_sign_conformal_energy"], rel_tol=3e-8),
                    f"{arm}/{label} binned number and energy cross-check")
        for label, switch in (("first", "A"), ("second", "B"), ("both", "A+B")):
            native = budget["cases"]["fine_" + arm]["static"][label]["full"]["proper_energy"]["daughter"]
            expected = profiles[switch]["both_sign_conformal_energy"]/a
            require(abs(native - expected) < 2e-6, f"{arm}/{switch} native/spectral energy cross-check")
        fixed = stationary[arm]["B"]["cutoffs"]["120"]
        require(fixed["vacuum_switch_work"] > 0 and math.isclose(
            fixed["vacuum_switch_work"], fixed["vacuum_daughter_energy"] + fixed["vacuum_parent_energy"],
            rel_tol=2e-12), f"{arm} stationary vacuum work account")
        for cutoff, arms in family.items():
            matrices = row["cutoffs"][cutoff]
            number = matrices["pair_number_matrix"]["real"]
            energy = matrices["both_sign_conformal_energy_matrix"]["real"]
            for label, result in arms[arm].items():
                early, late = result["coefficients"]
                low = -1. if label == "signed" else 0.
                require(low <= early <= 1 and low <= late <= 1, f"Gate domain for {arm}/{label}")
                gates = dict(A=1., B=0., first=early, last=late)
                vec = [gates[name] for name in emission["matrix_basis_order"]]
                score = (quadratic(number, vec)/number[0][0]) / (quadratic(energy, vec)/energy[0][0])**.75
                require(abs(score-result["ratio_to_A"]) < 2e-13 and
                        0 <= result["certified_upper"] - score < 1.1e-11, f"{arm}/{cutoff} optimality certificate")
    quiet = emission["arms"]["quiet"]
    for cutoff, data in controls["quiet_late_increments_per_cutoff"].items():
        profiles = quiet["cutoffs"][cutoff]["profiles"]
        gain = (profiles["A+B"]["both_sign_conformal_energy"] -
                profiles["A"]["both_sign_conformal_energy"])/quiet["scale_factor"]
        require(abs(gain-data["static"]["change_energy"]["daughter"]) < 2e-6,
                f"Native and spectral added energy at cutoff {cutoff}")
    for arm in ("quiet", "minus"):
        wall = gravity["results"]["fine_" + arm]["increment_0.005"]["readouts"]["2.16"]
        added = budget["cases"]["fine_" + arm]["static"]["added"]["late"]["change_energy"]["daughter"]
        require(abs(wall["pair"]-added) < 7e-6, f"{arm} gravity and Hamiltonian daughter source")
        require(abs(wall["h"]-wall["h_trace"]) < 8e-9, f"{arm} independent expansion readout")


def replay(repo, jobs):
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with tempfile.TemporaryDirectory(prefix="hdblast-renewed-") as temp:
        root = Path(temp)
        def run(script, *args):
            print(f"Recomputing {script}", flush=True)
            subprocess.run([sys.executable, "-I", "-B", str(HERE / script), *map(str, args)],
                           env=env, cwd=root, check=True)

        def result(section, name):
            new = read(root / name)
            compare(new, read(HERE / section / name), (("controls" if name == "CONTROLS.json" else section),))
            return new

        run("emission/emission.py", "--repo", repo, "--output", root)
        emission = result("emission", "EMISSION.json")
        sums = spectra(root / "SPECTRA.csv")
        run("criteria/criteria.py", "--repo", repo, "--output", root / "CRITERIA.json")
        criteria = result("criteria", "CRITERIA.json")
        run("stationary/stationary.py", "--repo", repo, "--output", root / "STATIONARY.json")
        stationary = result("stationary", "STATIONARY.json")
        run("family/family.py", "--emission", root / "EMISSION.json", "--output", root / "FAMILY.json")
        family = result("family", "FAMILY.json")
        tables = root / "tables"
        run("budget/budget.py", "--repo", repo, "--tables", tables, "--jobs", jobs,
            "--output", root / "BUDGET.json")
        budget = result("budget", "BUDGET.json")
        run("budget/controls.py", "--repo", repo, "--tables", tables, "--work", root / "cutoff-controls",
            "--jobs", jobs, "--output", root / "CONTROLS.json")
        controls = result("budget", "CONTROLS.json")
        run("gravity/gravity.py", "--repo", repo, "--tables", tables, "--jobs", jobs,
            "--output", root / "GRAVITY.json", "--compare", HERE / "gravity/GRAVITY.json")
        gravity = read(root / "GRAVITY.json")
        for name, sha in gravity["source_sha256"].items():
            require(digest(tables / name) == sha, f"Input to native solve: {name}")
        physical_checks(emission, criteria, stationary, family, budget, controls, gravity, sums)
    print("Reproduced the spectra, full-frequency and gate-family bounds, stationary work, "
          "finite-cutoff Hamiltonian budgets and 18 native gravity evolutions.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=HERE.parents[1],
                        help="Repository root with the pinned research/ files")
    parser.add_argument("--verify", action="store_true", help="Regenerate the numerical experiments")
    parser.add_argument("--jobs", type=int, default=2, help="Number of parallel numerical workers")
    options = parser.parse_args()
    require(options.jobs > 0, "Number of workers must be positive")
    repository = options.repo.resolve()
    authenticate(repository)
    if options.verify:
        replay(repository, options.jobs)
