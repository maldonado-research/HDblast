"""Verify the study files or regenerate its numerical results on fresh wall histories."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
PARENT = ROOT / "parent"
PARENT_MANIFEST = "c38d6fa2f0d49b81e3ebb84bfd676a1a7fe0b4ab65c1641b86592ecbf59d09f6"


def digest(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def manifest(folder):
    records = {}
    for line in (folder / "SHA256SUMS").read_text().splitlines():
        checksum, name = line.split("  ", 1)
        if name in records or name.startswith("/") or ".." in Path(name).parts:
            raise ValueError(f"Invalid manifest entry {name}")
        records[name] = checksum
    return records


def authenticate():
    expected = manifest(ROOT)
    files = {p.relative_to(ROOT).as_posix(): p for p in ROOT.rglob("*")
             if p.is_file() and "__pycache__" not in p.parts and p != ROOT / "SHA256SUMS"}
    if files.keys() != expected.keys():
        raise ValueError(f"Different file set: {files.keys() ^ expected.keys()}")
    for name, path in files.items():
        if path.is_symlink() or digest(path) != expected[name]:
            raise ValueError(f"File differs: {name}")
    if digest(PARENT / "SHA256SUMS") != PARENT_MANIFEST:
        raise ValueError("Predecessor manifest changed")
    original = manifest(PARENT)
    supplied = [name.removeprefix("parent/") for name in files
                if name.startswith("parent/") and name != "parent/SHA256SUMS"]
    for name in supplied:
        if digest(PARENT / name) != original[name]:
            raise ValueError(f"Input differs from preceding study: {name}")
    print(f"Verified {len(files)} study files and {len(supplied)} predecessor inputs.", flush=True)


def compare(actual, expected, name):
    if isinstance(expected, dict) and isinstance(actual, dict):
        if actual.keys() != expected.keys():
            raise ValueError(f"Different fields in {name}: {actual.keys() ^ expected.keys()}")
        for key, value in expected.items():
            compare(actual[key], value, f"{name}/{key}")
    elif isinstance(expected, list) and isinstance(actual, list):
        if len(actual) != len(expected):
            raise ValueError(f"Different number of elements in {name}")
        for i, value in enumerate(expected):
            compare(actual[i], value, f"{name}/{i}")
    elif isinstance(expected, float) and type(actual) in (int, float):
        if not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=2.e-6, abs_tol=1.e-9):
            raise ValueError(f"{name}: regenerated {actual}, archived {expected}")
    elif type(expected) != type(actual) or actual != expected:
        raise ValueError(f"{name}: regenerated {actual!r}, archived {expected!r}")


def run(script, *args):
    print(f"Replaying {script}", flush=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    subprocess.run([sys.executable, "-I", "-B", str(ROOT / script), *map(str, args)],
                   check=True, env=env)


def load_json(path):
    return json.loads(path.read_bytes())


def compare_history(history):
    import numpy as np
    compare(load_json(history / "markers.json"), load_json(ROOT / "history/markers.json"), "history/markers")
    for template in sorted((ROOT / "history").glob("*.npz")):
        with np.load(template, allow_pickle=False) as ref, np.load(history / template.name, allow_pickle=False) as actual:
            if set(ref.files) != set(actual.files):
                raise ValueError(f"Different history fields: {template.name}")
            for key in ref.files:
                a, b = actual[key], ref[key]
                if a.shape != b.shape or not np.isfinite(a).all():
                    raise ValueError(f"Invalid shape or samples in {template.name}/{key}")
                # Native stepping comes from the original solver. Resolve both individual
                # samples and the near-zero residuals from subtracting large impulses.
                floor = 2.e-10 * max(float(np.max(np.abs(b))), 1.e-12)
                np.testing.assert_allclose(a, b, rtol=1.e-6, atol=floor, err_msg=f"{template.name}/{key}")
        print(f"Matched every sample in {template.name}", flush=True)


def check_conclusions(output):
    """Check the two exclusions and the finite-wave example against their physical targets."""
    clustering = load_json(output / "global/CLUMPING.json")
    harmonic = load_json(output / "causal/CURRENT_CAP.json")["runs"]
    incoming = load_json(output / "bulk/INCOMING.json")
    particle = load_json(output / "particles/RESULTS.json")["runs"]
    for channel, early_limit, later_limit in (("dedicated", 15., 1.43), ("shared", 14., 1.29)):
        for mesh in ("native", "half"):
            key = f"{channel}_{mesh}"
            for work, threshold in (("0.001", early_limit), ("0.01", later_limit)):
                limit = clustering[key][work]["general_fields_eps_0.01"]
                if not limit["exclusion_at_unit_density"] or limit["excluded_cap_lower"] <= threshold:
                    raise ValueError(f"All-mode concentration exclusion changed at {key}/{work}")
            late_signal = incoming[key]["markers"]["0.01"]
            required = late_signal["bounds"]["0.01"]
            if required["minimum_incoming_profile_multiple"] <= (300 if channel == "dedicated" else 190):
                raise ValueError(f"Recorded-tail comparison changed at {key}")
            single = harmonic[key]["0.01"]
            if single["enclosing_entire_marker_cell"]["rms_upper"] >= single["previous_conditional_rms_impulse_lower"]:
                raise ValueError(f"A single initially empty Fourier pair now reaches the target at {key}")
        sample = [r for r in particle if r["history"] == channel + "_half" and r["resolution"] == "spectral"
                  and r["k_physical_opening"] == 50 and r["initial_peak_field"] == .04
                  and r["feedback"] and r["polarization"] == "By"]
        if len(sample) != 1:
            raise ValueError(f"Missing nonlinear finite-wave comparison for {channel}")
        final = sample[0]["result"]
        target = harmonic[channel + "_half"]["0.01"]["previous_conditional_rms_impulse_lower"]
        if final["transverse_force_integral_rms"] <= target or final["transverse_cdf_distance"] <= .25:
            raise ValueError(f"Finite-wave impulse and nonthermal-distribution comparison changed for {channel}")
        key = channel + "_half"
        marker = harmonic[key]["0.01"]
        summary = {
            "density_ratio_floor_at_0.1pct_work_for_1pp": clustering[key]["0.001"]["general_fields_eps_0.01"]["excluded_cap_lower"],
            "density_ratio_floor_at_1pct_work_for_1pp": clustering[key]["0.01"]["general_fields_eps_0.01"]["excluded_cap_lower"],
            "incoming_profile_minimum_factor": incoming[key]["markers"]["0.01"]["bounds"]["0.01"]["minimum_incoming_profile_multiple"],
            "single_mode_impulse_cap": marker["enclosing_entire_marker_cell"]["rms_upper"],
            "wave_case_transverse_cdf_gap": final["transverse_cdf_distance"],
        }
        print("Benchmark", channel, json.dumps(summary, sort_keys=True), flush=True)


def replay(dest, jobs):
    dest = dest.resolve()
    if dest.is_relative_to(ROOT) or (dest.exists() and (not dest.is_dir() or any(dest.iterdir()))):
        raise ValueError("Choose a new or empty output directory outside the study")
    dest.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="hdblast-transverse-histories-") as tmp:
        native = Path(tmp) / "native"
        run("parent/history/export.py", "--prior", PARENT / "prior", "--outdir", native, "--jobs", jobs)
        history = Path(tmp) / "history"
        run("history/extract.py", "--native", native, "--output-dir", history)
        compare_history(history)
        (dest / "history").mkdir(exist_ok=True)
        (dest / "history/markers.json").write_bytes((history / "markers.json").read_bytes())
        products = [
            ("history/budget.py", "history/BUDGET.json", ("--history", history)),
            ("history/focus.py", "history/FOCUS.json", ("--history", history, "--prior", PARENT / "prior", "--jobs", jobs)),
            ("dispersion/frozen.py", "dispersion/FROZEN.json", ("--history", history)),
            ("modes/noise_bound.py", "modes/NOISE_BOUND.json", ("--history", history)),
            ("modes/check.py", "modes/MODE_CHECK.json", ("--history", history)),
            ("causal/current_cap.py", "causal/CURRENT_CAP.json", ("--history", history)),
            ("modes/run.py", "modes/MODES.json", ("--history", history, "--jobs", jobs)),
            ("particles/check.py", "particles/CHECK.json", ()),
            ("particles/run.py", "particles/RESULTS.json", ("--history", history)),
            ("bulk/incoming.py", "bulk/INCOMING.json", ("--history", history, "--prior", PARENT / "prior", "--jobs", jobs)),
            ("global/clumping.py", "global/CLUMPING.json", ("--history", history)),
            ("energy/exposure.py", "energy/EXPOSURE.json", ("--history", history)),
        ]
        for script, product, args in products:
            out = dest / product
            out.parent.mkdir(parents=True, exist_ok=True)
            run(script, *args, "--output", out)
            compare(load_json(out), load_json(ROOT / product), product)
            if product == "history/BUDGET.json":
                (history / "BUDGET.json").write_bytes(out.read_bytes())
    check_conclusions(dest)
    print(f"All {len(products)} numerical products, four native histories and benchmark checks passed.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Rebuild histories and numerical results")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    opt = parser.parse_args()
    if opt.jobs < 1 or (opt.output_dir and not opt.verify):
        parser.error("--jobs must be positive and --output-dir requires --verify")
    authenticate()
    if opt.verify and opt.output_dir:
        replay(opt.output_dir, opt.jobs)
    elif opt.verify:
        with tempfile.TemporaryDirectory(prefix="hdblast-transverse-output-") as temp:
            replay(Path(temp), opt.jobs)
