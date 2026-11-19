"""Authenticate the self-contained modulus capsule and optionally recompute its physics."""
import argparse
import csv
import hashlib
import json
import math
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
BASE = "c2e8ccb69ae4e48ca7860a59e31a8bcf2969bfad"
PINS_SHA = "ec2f1682c23c90b504a6d8d54ddcb04c6f67e09e7b245e2398e4524ea5f251f1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def authenticate():
    inventory = HERE / "SHA256SUMS"
    listed = set()
    for line in inventory.read_text().splitlines():
        sha, name = line.split("  ", 1)
        path = (HERE / name).resolve()
        require(path.is_relative_to(HERE) and path != inventory and path not in listed, f"Invalid entry {name}")
        require(digest(path) == sha, f"Study checksum mismatch for {name}")
        listed.add(path)
    present = {p.resolve() for p in HERE.rglob("*") if p.is_file() and p != inventory
               and "__pycache__" not in p.relative_to(HERE).parts}
    require(present == listed, f"Capsule inventory differs: {present ^ listed}")
    manifest = HERE / "SOURCES.json"
    require(digest(manifest) == PINS_SHA, "Historical manifest checksum mismatch")
    pins = read(manifest)
    require(pins["base"] == BASE, "Historical revision mismatch")
    for source in pins["sources"]:
        local = HERE / source["bundled"]
        require(local.resolve() in listed and local.stat().st_size == source["bytes"]
                and digest(local) == source["sha256"], f"Historical input mismatch for {source['path']}")
    print(f"Authenticated {len(listed)} study files including {len(pins['sources'])} bundled historical inputs.", flush=True)


def tolerances(path):
    group = path[0]
    key = path[-1].lower()
    if group == "hilltop_check":
        return 4e-5, 3e-6 if any("margin" in p for p in path) else 2e-8
    if group == "hilltop":
        if "residual" in key or "error" in key or key in {"v", "log_margin", "stationary_theta"}:
            return 3e-5, 3e-9
        if "estimated_log10_seed_shifts" in path:
            return 3e-5, 1e-7
        return 3e-5, 1e-40
    if group == "widths":
        return (0., 2e-13) if "error" in key else (2e-10, 1e-40)
    if group == "modes":
        if any("residual" in str(p).lower() or "max_abs" in str(p)
               or "error" in str(p) or "refinement_relative" in str(p) for p in path):
            return 1e-5, 2e-11
        return 1e-5, 1e-16
    if group == "spectra":
        if "moment_relative_errors" in path or "maximum_two_solver_relative_difference" in path:
            return 0., 4e-8
        return 3e-7, 1e-35
    if group == "curvature":
        if "residual" in key or key in {"velocity", "h_at_max", "h_over_gamma"}:
            return 3e-7, 2e-10
        return 2e-6, 1e-35
    # Oscillator phases and pressure-work cancellations depend on integration order.
    if group == "checks" or "residual" in key or "checks" in path:
        return 3e-4, 3e-7
    if key in {"gamma", "radiation", "radiation_fraction", "radiation_source_comoving4"}:
        return 4e-5, 2e-12
    return 3e-5, 2e-6


def compare(actual, expected, path):
    if isinstance(expected, dict):
        require(type(actual) is dict and actual.keys() == expected.keys(), f"Object fields at {path}")
        for k in expected:
            compare(actual[k], expected[k], (*path, k))
    elif isinstance(expected, list):
        require(type(actual) is list and len(actual) == len(expected), f"List length at {path}")
        for i, (new, old) in enumerate(zip(actual, expected)):
            compare(new, old, (*path, str(i)))
    elif isinstance(expected, float):
        rel, absolute = tolerances(path)
        require(type(actual) in (int, float) and math.isfinite(actual)
                and math.isclose(actual, expected, rel_tol=rel, abs_tol=absolute),
                f"{'/'.join(path)}: computed {actual}, expected {expected}")
    else:
        require(type(expected) is type(actual) and expected == actual, f"Value at {'/'.join(path)}")


def csv_check(new, old):
    lines = old.name == "LINES.csv"
    n_rows = 0
    with new.open(newline="") as a, old.open(newline="") as b:
        fresh, saved = csv.DictReader(a), csv.DictReader(b)
        require(saved.fieldnames == fresh.fieldnames, f"Column names in {old.name}")
        for n_rows, (r, s) in enumerate(zip(fresh, saved, strict=True), 1):
            for field in saved.fieldnames:
                if field in {"case", "epoch", "anchor", "phi_min", "detuning"}:
                    require(r[field] == s[field], f"Row {n_rows} {field} in {old.name}")
                else:
                    x, y = float(r[field]), float(s[field])
                    rel, absolute = (1e-4, 1e-13) if lines else (3e-7, 2e-15)
                    require(math.isfinite(x) and math.isclose(x, y, rel_tol=rel, abs_tol=absolute),
                            f"Row {n_rows}, {field} of {old.name}: {x} vs {y}")
    require(n_rows > 0, f"Empty {old.name}")
    return n_rows


def cross_checks(widths, background, modes, curvature, spectra, hilltop):
    require(abs(widths["scalar_to_original_ledger_ratio"] * 32 * math.pi - 1) < 1e-13,
            "Real scalar trace normalization")
    universal = curvature["universal"]
    require(abs(universal["g_critical"] - universal["g_from_angle"]) < 3e-10,
            "Independent parametrizations at the equality/turnaround boundary")
    require(abs(universal["g_critical"] ** (-2/3) - universal["k_over_GammaM_to_2thirds"]) < 3e-10,
            "Dimensionless curvature normalization")
    cases = modes["cases"]
    require(cases["conformal"]["energy_density"] == cases["flat"]["energy_density"] == 0,
            "Conformal and static nulls")
    for label in ("short", "reference", "long"):
        case = cases[label]
        require(case["number_density"] > 0 and math.isclose(
            case["energy_density"]/case["number_density"], case["mean_momentum"], rel_tol=1e-10),
            f"Mode energy per daughter, {label}")
    for pm in ("0.5", "-0.5"):
        name = "phi_m=" + pm
        w = widths["branches"][name]
        b = background["branches"][name]
        c = curvature["branches"][name]
        expected = w["alpha"]**2 * w["m_over_Htop"]**3 / (32*math.pi)
        require(math.isclose(b["gamma"], expected*b["Htop_over_Mpl"]**2, rel_tol=1e-11),
                f"Common trace width on branch {name}")
        require(modes["archived_A"][name]["archived_derivative_magnitude"] == w["alpha"],
                f"Finite-pulse brane coupling on {name}")
        for d, anchor in c["flat_m_over_50"].items():
            flat = spectra["branches"][pm][d]
            eq = flat["epochs"]["equality"]
            reference = anchor["flat_equality"]
            require(math.isclose(reference["H_over_Gamma"], eq["H_over_gamma"], rel_tol=3e-7)
                    and math.isclose(reference["gamma_elapsed"], eq["gamma_elapsed_time"], rel_tol=3e-7),
                    f"Independent equality time and expansion on {name} at {d}")
            require(0 < eq["bose_count_fraction"] <= flat["all_time_bose_count_upper"] < 1,
                    f"Spectral occupancy bound on {name} at {d}")
        turn = b["closed"]["harmonic"]["turnaround"]
        fraction = turn["radiation"]/(turn["scalar"]+turn["radiation"])
        require(abs(fraction-c["closed_archived_launch_at_t_0p001"]["harmonic"]["radiation_fraction"]) < 1e-12
                and fraction < .5, f"Actual closed-trajectory endpoint on {name}")
        for detuning, fit in hilltop["branches"][name]["detunings"].items():
            primary = fit["primary"]
            root = primary["anchors"]["24_cross"]
            coordinate = fit["gamma"]*root["a"]**3*root["scalar"]
            require(math.isclose(coordinate, root["critical_rate"], rel_tol=4e-6),
                    f"Closed launch at curvature threshold on {name} at {detuning}")
            time_brackets = fit["brackets_all_handoff_points"]
            require(time_brackets[0]["max_log_margin"] < 0 < time_brackets[1]["min_log_margin"]
                    and primary["log10_bounce_seed_over_archive"] < 0,
                    f"Opposite first-turn outcomes across the hilltop bracket on {name} at {detuning}")



def replay(jobs):
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    archive = HERE / "archive"
    ledger = archive / "BLAST_REHEATING_RESULTS.json"
    with tempfile.TemporaryDirectory(prefix="hdblast-modulus-") as temp:
        fresh = Path(temp)
        def run(name, *args, script=None):
            source = HERE/name/((script or name)+".py")
            print(f"Recomputing {source.relative_to(HERE)}", flush=True)
            subprocess.run([sys.executable, "-I", "-B", str(source), *map(str, args)],
                           cwd=fresh, env=env, check=True)

        with ThreadPoolExecutor(max_workers=jobs) as workers:
            widths_job = workers.submit(run, "widths", "--ledger", ledger, "--output", fresh/"WIDTHS.json")
            background_job = workers.submit(run, "background", "--archive", archive,
                                            "--output", fresh/"BACKGROUND.json", "--checks-output", fresh/"CHECKS.json")
            modes_job = workers.submit(run, "modes", "--archive", archive, "--output", fresh)
            hilltop_job = workers.submit(run, "hilltop", "--archive", archive, "--output", fresh/"HILLTOP.json")
            widths_job.result()
            background_job.result()
            curvature_job = workers.submit(run, "curvature", "--widths", fresh/"WIDTHS.json",
                                          "--background", fresh/"BACKGROUND.json", "--output", fresh/"CURVATURE.json")
            spectra_job = workers.submit(run, "spectra", "--ledger", ledger,
                                        "--background", fresh/"BACKGROUND.json", "--output", fresh)
            hilltop_job.result()
            independent = workers.submit(run, "hilltop", "--archive", archive, "--input", fresh/"HILLTOP.json",
                                         "--output", fresh/"HILLTOP_CHECK.json", script="check")
            for job in (modes_job, independent, curvature_job, spectra_job):
                job.result()
        results = {}
        for name, section in (("widths","widths"),("background","background"),("checks","background"),
                              ("modes","modes"),("curvature","curvature"),("spectra","spectra"),("hilltop","hilltop")):
            file = name.upper() + ".json"
            data = read(fresh/file)
            compare(data, read(HERE/section/file), (name,))
            results[name] = data
        check = read(fresh/"HILLTOP_CHECK.json")
        compare(check, read(HERE/"hilltop/CHECK.json"), ("hilltop_check",))
        for branch in check.values():
            for row in branch.values():
                if "actual_bounce_bracket_margins" in row:
                    lo, hi = row["actual_bounce_bracket_margins"]
                    require(lo < 0 < hi, "Independent nonlinear launch bracket")
        row_counts = [csv_check(fresh/file, HERE/section/file)
                      for section, file in (("modes", "LINES.csv"), ("spectra", "SPECTRA.csv"))]
        cross_checks(*(results[name] for name in ("widths","background","modes","curvature","spectra","hilltop")))
    print(f"Verified the widths, homogeneous trajectories, finite-pulse modes, hilltop launches, curvature and spectra "
          f"including {sum(row_counts)} CSV rows.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Regenerate the numerical study in a temporary directory")
    parser.add_argument("--jobs", type=int, default=2, help="Number of concurrent numerical workers")
    args = parser.parse_args()
    require(args.jobs > 0, "--jobs must be positive")
    authenticate()
    if args.verify:
        replay(args.jobs)
