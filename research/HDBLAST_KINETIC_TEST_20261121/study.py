"""Regenerate the HDblast particle relaxation calculations using shared wall histories."""
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
PARENT_COMMIT = "7fe36012da64e98add3915acd0356c733e717670"
PRODUCTS = ("moments/MOMENTS", "box/BOX", "history/HISTORY", "ballistic/BALLISTIC", "closure/FLOW", "box/QUICK", "history/CHECK", "box/AZIMUTH", "driven/DRIVEN", "ballistic/TUBE", "closure/SPATIAL", "driven/ORIENTED", "closure/CONVERGENCE", "driven/LARGE", "ballistic/STRESS", "driven/ORIENTED_CHECK", "closure/CHECK", "live/LIVE", "ballistic/REVERSAL", "ballistic/KICK", "live/ANGULAR_NOISE", "live/EXTENDED")


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


def compare_ensemble(actual, reference, path, differences):
    """Check a finite-ensemble result at a time marker shared with the wall solver."""
    if actual.keys() != reference.keys():
        differences.append(f"{path}: ensemble fields changed")
        return
    mean = actual["mean"]
    ref = reference["mean"]
    if mean.keys() != ref.keys() or not all(math.isfinite(x) for x in mean.values()):
        differences.append(f"{path}: invalid ensemble means")
        return
    centered = actual.get("standard_error", {})
    ref_centered = reference.get("standard_error", {})
    if centered.keys() != ref_centered.keys() or any(not math.isfinite(x) or x < 0 for x in centered.values()):
        differences.append(f"{path}: invalid standard errors")
        return
    ranges = actual["seed_range"]
    if ranges.keys() != reference["seed_range"].keys():
        differences.append(f"{path}: changed seed-range fields")
        return
    for name, (lower, upper) in ranges.items():
        if (not math.isfinite(lower) or not math.isfinite(upper) or
                lower > mean[name] + 1.e-12 * max(1., abs(lower)) or
                upper < mean[name] - 1.e-12 * max(1., abs(upper))):
            differences.append(f"{path}/{name}: seed range excludes its mean")

    for name in ("t", "n", "injection"):
        compare(mean[name], ref[name], path + "/mean/" + name, differences)
    if mean["pair_error"] > 1.e-7 or mean["energy_residual"] > 2.e-6:
        differences.append(f"{path}: energy ledger or elastic collision balance failed")
    if min(mean[k] for k in ("gas_comoving", "mean_lab_energy", "replaced_total_density")) <= 0:
        differences.append(f"{path}: nonpositive particle energy or wall density")

    fraction_fields = {"never_hit", "landau_fraction", "averaged_landau_fraction",
                       "longitudinal_share", "averaged_longitudinal_share", "rest_longitudinal_share",
                       "energy_cdf_distance", "tube_fraction", "thermal_tube_fraction"}
    for name in fraction_fields & mean.keys():
        if not -1.e-10 <= mean[name] <= 1.0000000001:
            differences.append(f"{path}/mean/{name}: fraction outside physical range")
    for name in ("charge_velocity", "landau_velocity"):
        if abs(mean[name]) > 1.:
            differences.append(f"{path}/mean/{name}: superluminal velocity")

    free = "_no_scattering/" in path
    if free and (mean["hits_per_particle"] != 0 or mean["never_hit"] != 1):
        differences.append(f"{path}: particle scattered in collisionless control")
    if "fixed_field/" in path and abs(mean["field_change"]) > 1.e-12:
        differences.append(f"{path}: prescribed-field control changed the field")

    # Late self-consistent electric oscillations dephase across compilers. For the
    # free-streaming control test the integrated diagnostics and the absence of hits.
    late = path.endswith(("/radius_doubling", "/turnaround", "/reference_first_hot"))
    for name in (fraction_fields - {"thermal_tube_fraction"}) & mean.keys():
        if free and late and name not in {"never_hit", "averaged_landau_fraction", "averaged_longitudinal_share"}:
            continue
        width = lambda row: (row[name][1] - row[name][0]) / 2 if name in row else 0.
        total_error = math.hypot(centered.get(name, width(ranges)),
                                 ref_centered.get(name, width(reference["seed_range"])))
        # Each run samples a finite set of particles on a fixed wall-time grid.
        allowance = 3. / math.sqrt(mean["n"]) + 4. * total_error
        if abs(mean[name] - ref[name]) > allowance:
            differences.append(f"{path}/mean/{name}: {mean[name]} versus {ref[name]}, allowance {allowance}")


def compare_live(actual, reference, path, differences):
    """Check configuration, conservation and heating in a finite particle realization."""
    if actual.keys() != reference.keys() or not all(math.isfinite(x) for x in actual.values()):
        differences.append(f"{path}: incomplete or nonfinite particle trajectory")
        return
    for name in ("particles", "particle_count", "seed", "branch", "cfl", "stop_turnaround"):
        compare(actual[name], reference[name], path + "/" + name, differences)

    # The two seeds and numerical resolutions have distinct sampled electric
    # phases. Compare growth and the bulk evolution and test thresholds directly.
    for name in ("maximum_radius", "maximum_tau", "end_gas_density", "end_electric_work",
                 "end_gas_pressure_work", "hot_radius_growth"):
        if not math.isclose(actual[name], reference[name], rel_tol=.05):
            differences.append(f"{path}/{name}: {actual[name]} versus {reference[name]}")
    for name in ("end_rest", "integrated_birth_energy", "end_sampled_rest", "min_areal_radius"):
        compare(actual[name], reference[name], path + "/" + name, differences)
    if not math.isclose(actual["end_sampled_rest"], actual["end_rest"], rel_tol=1.e-6):
        differences.append(f"{path}: final injected and sampled rest masses differ")
    if actual["particle_energy_error_over_peak"] > .005 or actual["max_event_conservation_error"] > 1.e-10:
        differences.append(f"{path}: energy or scattering conservation failed")
    for name, value in actual.items():
        if name.endswith(("_empirical_tube_share", "_angular_trace_share", "_lab_dispersion_fraction", "_shape_ks")):
            if not -1.e-10 <= value <= 1.0000000001:
                differences.append(f"{path}/{name}: invalid fraction")
        if name.endswith("_shape_pressure_check") and value > 1.e-8:
            differences.append(f"{path}/{name}: independently evaluated pressure differs")
    if actual["branch"] == 0:
        if actual["scatter_events"] != 0 or actual["peak_lab_dispersion"] >= .90:
            differences.append(f"{path}: collisionless dispersion control failed")
    elif actual["scatter_events"] <= 0 or not .95 < actual["peak_lab_dispersion"] <= 1.0000000001:
        differences.append(f"{path}: late particle heating comparison changed")
    if actual["particles"] == 2048:
        if not (.90 < actual["doubling_lab_dispersion_fraction"] <= 1
                and 1 < actual["lab_hot_growth"] < 2
                and actual["doubling_shape_ks"] > .1):
            differences.append(f"{path}: high-resolution heating or spectral comparison changed")
    for name, value in actual.items():
        if name.endswith("_start_radius"):
            prefix = name.removesuffix("_start_radius")
            growth = actual[prefix + "_growth"]
            if value <= 0 or not math.isclose(actual[prefix + "_end_radius"] / value, growth, rel_tol=1.e-7):
                differences.append(f"{path}/{prefix}: interval endpoints and growth disagree")


def compare(actual, reference, path, differences):
    if type(actual) is not type(reference):
        differences.append(f"{path}: changed data type")
    elif isinstance(reference, dict):
        if path.startswith(("live/LIVE/", "live/EXTENDED/")) and reference.get("particles", 0) > 0:
            compare_live(actual, reference, path, differences)
            return
        if path.startswith("driven/") and "/groups/" in path and "mean" in reference:
            compare_ensemble(actual, reference, path, differences)
            return
        if actual.keys() != reference.keys():
            differences.append(f"{path}: different fields {actual.keys() ^ reference.keys()}")
        for key in sorted(reference.keys() & actual.keys()):
            compare(actual[key], reference[key], path + "/" + key, differences)
    elif isinstance(reference, list):
        if len(actual) != len(reference):
            differences.append(f"{path}: changed array length")
        for i, (left, right) in enumerate(zip(actual, reference)):
            compare(left, right, f"{path}/{i}", differences)
    elif isinstance(reference, float):
        abs_tol = 3.e-12
        # Balance residuals are differences between two accumulated energy ledgers.
        if path.startswith(("closure/", "live/")) and path.endswith(("_balance_error", "/friedmann_error")):
            abs_tol = 1.e-10
        # Preserve relative precision even for the extremely small dilute rates
        # and probabilities that enter the thermal exclusion.
        if (path.startswith("closure/CHECK/") and path.rsplit("/", 1)[-1] in {"angular_reference", "kernel"}
                or path.startswith("ballistic/TUBE/") and path.endswith("/maxwell_juttner_tube_fraction")
                or path.startswith("moments/MOMENTS/") and (path.endswith("/C") or "/TM/C/" in path)):
            abs_tol = 1.e-28
        if not math.isfinite(actual) or not math.isclose(actual, reference, rel_tol=2.e-6, abs_tol=abs_tol):
            differences.append(f"{path}: regenerated {actual!r}, archived {reference!r}")
    elif actual != reference:
        differences.append(f"{path}: regenerated {actual!r}, archived {reference!r}")


def check_particles(products):
    # Two smooth positive spectra have the same first moments and different drag.
    witnesses = products["moments/MOMENTS"]["startup"]["finite_band_Synge"]["witnesses"]
    high, low = witnesses["high_band"], witnesses["low_band"]
    if min(high["C"], low["C"]) <= 0 or not 3.7 < high["C"] / low["C"] < 4.0:
        raise ValueError("Equal-moment finite-band comparison changed")
    if not all(0 < state["per_spin_f_upper"] < 1 for state in (high, low)):
        raise ValueError("Trial spectrum exceeds the vacancy limit")

    box = products["box/BOX"]["runs"]
    for name, row in box.items():
        if max(row["maximum_residual"].values()) > 1e-8:
            raise ValueError(f"Elastic conservation failed in {name}")
    narrow, thermal = box["drift_scale_mono"], box["drift_scale_thermal"]
    index = narrow["times"].index(8.)
    gap = narrow["mean"]["ks_energy_pooled"][index]
    if not (abs(narrow["mean"]["charge_odd_momentum"][index]) < .03 and
            gap > .1 and thermal["mean"]["ks_energy_pooled"][index] < .02):
        raise ValueError("Loss of drift and relaxation of the energy distribution changed")

    for name, row in products["ballistic/BALLISTIC"]["runs"].items():
        mid = row["half_work"]["fine"]["mean_survival"]
        first_hot = row["first_hot"]["fine"]["mean_survival"]
        if not (0 < first_hot < .15 and .25 < mid < .70):
            raise ValueError(f"Collisionless survival estimate changed for {name}")
        stage = products["history/HISTORY"]["runs"][name]["first_hot_time"]
        if not (0 < stage["fraction_zero"] < .01 and stage["mean_events_per_particle"] > 1):
            raise ValueError(f"Thermal event clock changed for {name}")

    tubes = products["ballistic/TUBE"]["runs"]
    for name in ("shared_native", "shared_half"):
        early = tubes[name]["early_field_work"]["0.01"]
        if not (early["mean_survival_floor"] > .54 and early["maxwell_juttner_tube_fraction"] < 7e-8):
            raise ValueError(f"Early transverse-momentum exclusion changed for {name}")
    for name in ("dedicated_native", "dedicated_half"):
        half_work = products["ballistic/STRESS"]["runs"][name]["absolute_work_50_percent"]["cell_enclosure"]
        if not half_work["fluid_current_velocity_high"] < half_work["allowed_current_velocity_low"]:
            raise ValueError(f"Current constraint changed for {name}")
    for group, stages in products["driven/ORIENTED"]["groups"].items():
        gyro_group = "_".join(group.split("_")[:-1] + ["kinetic"])
        if gyro_group not in products["driven/DRIVEN"]["groups"]:
            continue
        for stage, row in stages.items():
            if "gyrotropic_landau_fraction" in row:
                same_run = products["driven/DRIVEN"]["groups"][gyro_group][stage]["mean"]["landau_fraction"]
                if abs(row["gyrotropic_landau_fraction"] - same_run) > 1e-12:
                    raise ValueError(f"Different underlying gyrotropic history in {group}/{stage}")


def check_physics(products):
    if "moments/MOMENTS" in products:
        check_particles(products)
    if "closure/FLOW" not in products:
        return
    flows = products["closure/FLOW"]
    for name, row in flows.items():
        if not row["stop_turnaround"] or max(row[k] for k in (
                "surface_balance_error", "gas_balance_error", "gas_primitive_error")) > 1e-5:
            raise ValueError(f"Fluid closure conservation or turnaround failed in {name}")


def replay(destination, suite="all"):
    destination = destination.resolve()
    if destination.is_relative_to(ROOT) or (destination.exists() and
            (not destination.is_dir() or any(destination.iterdir()))):
        raise ValueError("Choose a new or empty output directory outside the study.")
    destination.mkdir(parents=True, exist_ok=True)
    def family(key):
        prefix = key.split("/", 1)[0]
        return prefix if prefix in {"closure", "live"} else "particles"

    outputs = {key: destination / (key + ".json") for key in PRODUCTS if suite == "all" or family(key) == suite}
    for path in outputs.values():
        path.parent.mkdir(exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def run(script, *args):
        print("Replaying " + script, flush=True)
        subprocess.run([sys.executable, "-I", "-B", str(ROOT / script), *map(str, args)],
                       check=True, env=env)

    def closure():
        with ThreadPoolExecutor(max_workers=2) as pool:
            tasks = [pool.submit(run, "closure/flow.py", "--prior", PRIOR, "--output",
                                 outputs["closure/" + name], "--jobs", "2", *flags)
                     for name, flags in (("FLOW", ()), ("SPATIAL", ("--spatial",)),
                                         ("CONVERGENCE", ("--convergence",)))]
            run("closure/check.py", "--output", outputs["closure/CHECK"])
            for task in tasks:
                task.result()

    def live():
        run("live/angular_noise.py", "--output", outputs["live/ANGULAR_NOISE"])
        with ThreadPoolExecutor(max_workers=2) as pool:
            extended = pool.submit(run, "live/run.py", "--prior", PRIOR, "--extended", "--output",
                                   outputs["live/EXTENDED"], "--jobs", "1")
            run("live/run.py", "--prior", PRIOR, "--output", outputs["live/LIVE"], "--jobs", "2")
            extended.result()

    def particles():
        with tempfile.TemporaryDirectory(prefix="hdblast-kinetic-history-") as work:
            native = Path(work)
            with ThreadPoolExecutor(max_workers=3) as pool:
                box = pool.submit(run, "box/run.py", "--prior", PRIOR, "--jobs", "2", "--output", outputs["box/BOX"])
                run("history/export.py", "--prior", PRIOR, "--outdir", native, "--jobs", "2")
                large = pool.submit(run, "driven/run.py", "--large", "--native", native,
                                    "--jobs", "2", "--output", outputs["driven/LARGE"])
                driven = pool.submit(run, "driven/run.py", "--native", native,
                                     "--jobs", "2", "--output", outputs["driven/DRIVEN"])
                optical = pool.submit(run, "history/optical.py", "--native", native,
                                      "--jobs", "2", "--output", outputs["history/HISTORY"])
                ballistic = pool.submit(run, "ballistic/ballistic.py", "--traces", native,
                                        "--jobs", "2", "--output", outputs["ballistic/BALLISTIC"])
                run("moments/counterexample.py", "--history", native, "--output", outputs["moments/MOMENTS"])
                run("box/run.py", "--prior", PRIOR, "--quick", "--jobs", "2", "--output", outputs["box/QUICK"])
                box.result()
                run("box/azimuth.py", "--input", outputs["box/BOX"], "--output", outputs["box/AZIMUTH"])
                optical.result()
                for script, name in (("reversal", "REVERSAL"), ("stress", "STRESS"), ("tube", "TUBE")):
                    run("ballistic/" + script + ".py", "--traces", native, "--history", outputs["history/HISTORY"],
                        "--output", outputs["ballistic/" + name])
                run("ballistic/kick.py", "--tube", outputs["ballistic/TUBE"], "--native", native,
                    "--output", outputs["ballistic/KICK"])
                run("history/check.py", "--native", native, "--history", outputs["history/HISTORY"],
                    "--output", outputs["history/CHECK"], "--jobs", "2")
                driven.result()
                run("driven/oriented.py", "--native", native, "--gyro", outputs["driven/DRIVEN"],
                    "--output", outputs["driven/ORIENTED"], "--jobs", "2")
                run("driven/oriented_check.py", "--output", outputs["driven/ORIENTED_CHECK"])
                large.result()
                ballistic.result()

    suites = {"particles": particles, "closure": closure, "live": live}
    if suite == "all":
        with ThreadPoolExecutor(max_workers=3) as pool:
            tasks = [pool.submit(method) for method in suites.values()]
            for task in tasks:
                task.result()
    else:
        suites[suite]()

    results = {key: json.loads(path.read_bytes()) for key, path in outputs.items()}
    check_physics(results)
    differences = []
    for key, actual in results.items():
        reference = json.loads((ROOT / (key + ".json")).read_bytes())
        compare(actual, reference, key, differences)
    if differences:
        count = len(differences)
        raise ValueError(f"{count} differences in regenerated results:\n" + "\n".join(differences[:100]))
    print(f"All {len(outputs)} numerical products and physical checks passed.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-prior", type=Path, help="Compare the vendored files with the preceding repository study")
    parser.add_argument("--verify", action="store_true", help="Recompute the numerical outputs")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--suite", choices=("all", "particles", "closure", "live"), default="all",
                        help="Run one independent part of the complete study")
    args = parser.parse_args()
    if args.suite != "all" and not args.verify:
        parser.error("--suite requires --verify")
    if args.output_dir and not args.verify:
        parser.error("--output-dir requires --verify")
    authenticate(args.check_prior)
    if args.verify and args.output_dir:
        replay(args.output_dir, args.suite)
    elif args.verify:
        with tempfile.TemporaryDirectory(prefix="hdblast-kinetic-replay-") as temp:
            replay(Path(temp), args.suite)
