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
import shutil

ROOT = Path(__file__).resolve().parent
PARENT = ROOT / "parent"
PARENT_MANIFEST = "9ce9e06b5574bf7d7c90cfc4fbedb4dd59c70311526e42652e19a31c137e6583"


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
    compare(load_json(history / "markers.json"), load_json(PARENT / "history/markers.json"), "history/markers")
    for template in sorted((PARENT / "history").glob("*.npz")):
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


def rebuild_histories(tmp, jobs):
    """Replace the previous study's stored time series with their independently regenerated values."""
    native = tmp / "native"
    run("parent/parent/history/export.py", "--prior", PARENT / "parent/prior",
        "--outdir", native, "--jobs", jobs)
    scratch = tmp / "prior"
    shutil.copytree(PARENT, scratch)
    run("parent/history/extract.py", "--native", native, "--output-dir", scratch / "history")
    compare_history(scratch / "history")
    shutil.rmtree(native)
    return scratch


def replay(dest, jobs):
    dest = dest.resolve()
    if dest.is_relative_to(ROOT) or (dest.exists() and (not dest.is_dir() or any(dest.iterdir()))):
        raise ValueError("Choose a new or empty output directory outside the study")
    dest.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="hdblast-prebirth-histories-") as tmp:
        prior = rebuild_histories(Path(tmp), jobs)
        transfer = dest / "transfer"
        geometry = dest / "geometry"
        raw = Path(tmp) / "incident"
        run("transfer/capture.py", "--prior", prior, "--output-dir", raw, "--jobs", jobs)
        run("transfer/run.py", "--prior", prior, "--output-dir", transfer, "--raw-dir", raw, "--jobs", jobs)
        run("geometry/probe.py", "--prior", prior / "parent/prior", "--history", prior / "history",
            "--output-dir", geometry, "--jobs", jobs)
        run("geometry/probe.py", "--prior", prior / "parent/prior", "--history", prior / "history",
            "--output-dir", geometry, "--spectrum", "--jobs", jobs)
        run("geometry/probe.py", "--prior", prior / "parent/prior", "--history", prior / "history",
            "--output-dir", geometry, "--gauge-controls", "--jobs", jobs)
        orbits = dest / "orbits"
        orbits.mkdir()
        run("orbits/run.py", "--prior", prior, "--transfer", transfer / "traces",
            "--jobs", jobs, "--output", orbits / "LOCAL.json")
        run("orbits/run.py", "--prior", prior, "--transfer", geometry / "traces",
            "--jobs", jobs, "--output", orbits / "BULK.json")
        run("orbits/calibrate.py", "--prior", prior, "--transfer", geometry / "traces",
            "--orbit-product", orbits / "BULK.json", "--output", orbits / "CALIBRATION.json", "--jobs", jobs)
        run("orbits/axial_budget.py", "--history", prior / "history", "--output", orbits / "AXIAL.json")
        run("orbits/frame.py", "--history", prior / "history", "--calibration", orbits / "CALIBRATION.json",
            "--axial-budget", orbits / "AXIAL.json", "--output", orbits / "FRAME.json")
        control = dest / "control"
        control.mkdir()
        run("control/optimal.py", "--transfer", raw, "--output", control / "OPTIMAL.json")
        sphere = dest / "sphere"
        sphere.mkdir()
        run("sphere/round.py", "--history", prior / "history", "--geometry", geometry / "traces",
            "--output", sphere / "SPHERE.json", "--jobs", 1)
        run("sphere/all_patterns.py", "--history", prior / "history", "--geometry", geometry / "spectrum_traces",
            "--output", sphere / "ALL_PATTERNS.json")
    results = [
        "transfer/transfer.json",
        "geometry/PROBE.json", "geometry/SPECTRUM.json", "geometry/GAUGE.json",
        "orbits/LOCAL.json", "orbits/BULK.json", "orbits/CALIBRATION.json",
        "orbits/AXIAL.json", "orbits/FRAME.json", "control/OPTIMAL.json",
        "sphere/SPHERE.json", "sphere/ALL_PATTERNS.json",
    ]
    for rel in results:
        compare(load_json(dest / rel), load_json(ROOT / rel), rel)
        print("Reproduced", rel, flush=True)
    print_benchmarks(dest)
    print(f"Rebuilt the four wall histories and all {len(results)} prebirth-transfer results.", flush=True)


def print_benchmarks(output):
    """Expose the regenerated benchmark independently of the checked-in reference tables."""
    benchmarks = reported_values(output)
    for name, values in benchmarks.items():
        line = json.dumps(values, sort_keys=True, allow_nan=False)
        print("Benchmark", name, line, flush=True)
        if os.environ.get("GITHUB_ACTIONS") == "true":
            print(f"::notice file=study.py,line=1,title=HDblast {name}::{line}", flush=True)

def reported_values(output):
    particle = load_json(output / 'orbits/BULK.json')
    geometry = load_json(output / 'geometry/PROBE.json')
    spectrum = load_json(output / 'geometry/SPECTRUM.json')
    axial = load_json(output / 'orbits/AXIAL.json')
    calibration = load_json(output / 'orbits/CALIBRATION.json')
    frames = load_json(output / 'orbits/FRAME.json')
    patterns = load_json(output / 'sphere/ALL_PATTERNS.json')
    gauge = load_json(output / 'geometry/GAUGE.json')
    cases = {}
    for channel in ('dedicated', 'shared'):
        label = channel + '_half'
        candidate = particle['best_base'][label]
        rows = [r for r in particle['runs'] if r['history'] == label and r['feedback']
                and r['resolution'] == 'spectral' and r['marker'] == '0.01'
                and r['mode_index'] == candidate['mode_index'] and r['drive_amplitude'] == candidate['amplitude']]
        if len(rows) != 1:
            raise ValueError('Missing refined particle comparison for ' + label)
        row = rows[0]
        result = row['result']
        modes = geometry[label]['modes']
        kappa, mode = min(modes.items(), key=lambda item: abs(item[1]['k_birth'] - row['k_physical_birth']))
        if abs(mode['k_birth'] - row['k_physical_birth']) > 1.e-8:
            raise ValueError('Unmatched five-dimensional wave number')
        if geometry[label]['max_absolute_unforced_field'] != 0:
            raise ValueError('An unforced harmonic acquired a field')
        low = abs(modes['5']['birth']['E'])
        high = max(abs(m['birth']['E']) for n, m in modes.items() if int(n) >= 50)
        if high > 1.e-5*low:
            raise ValueError('The small-scale field exceeds the stated spectral contrast')
        band = spectrum[label]['modes']
        peak = max(band, key=lambda n: math.hypot(band[n]['birth']['E'], band[n]['birth']['B']))
        birth = band[peak]['birth']
        minimum = axial[label]['0.01']['energy_over_thermal_lower']['7']
        if peak != '7' or minimum < (7.2 if channel == 'dedicated' else 12.1):
            raise ValueError('The transfer peak or one-polarization energy exclusion changed for ' + label)
        fine = calibration[label]['best_fine']['result']
        all_modes = patterns[label]['0.01']
        uniform_floor = all_modes['bounds']['uniform_birth_positions_finite_travel']['minimum_launch_rms_to_avoid_excess']
        if all_modes['real_angular_coefficients'] != 44200 or uniform_floor < .15:
            raise ValueError('The spherical angular-band energy floor changed for ' + label)
        held_electric = gauge['dedicated_half_hold_release']['modes'][peak]['birth']['E']
        drifting = frames[label]['comparisons_by_resolution']['fine']
        count_test = all_modes['transverse_cdf_test']['uniform_birth_positions_finite_travel']
        cases[channel] = {
            'uniform_birth_minimum_rms_for_thermal_counts_within_one_point': count_test['necessary_launch_rms_for_one_percentage_point'],
            'minimum_cdf_discrepancy_for_one_percent_launch': count_test['witnesses']['0.01']['guaranteed_excess_cdf'],
            'mean_energy_over_drifting_thermal_prediction': drifting['mean_energy_over_boosted_prediction'],
            'lowest_pressure_over_drifting_thermal_prediction': drifting['least_principal_stress_fraction_over_prediction'],
            'real_angular_coefficients': all_modes['real_angular_coefficients'],
            'uniform_birth_minimum_launch_rms_for_stationary_energy': uniform_floor,
            'birth_electric_amplitude_gain_over_held_final_coupling': abs(birth['E']/held_electric),
            'maximum_birth_response_harmonic_2_to_40': int(peak),
            'peak_mode_birth_field_magnitude_per_unit_launch': math.hypot(birth['E'], birth['B']),
            'peak_mode_comoving_wall_energy': birth['comoving_wall_energy'],
            'peak_mode_net_incoming_work': birth['net_incoming_work'],
            'peak_mode_gauge_work': birth['gauge_work'],
            'linear_polarization_energy_over_stationary_thermal_lower': minimum,
            'refined_amplitude_cdf_gap': fine['cdf_gap'],
            'refined_amplitude_largest_pressure_share': fine['stress_eigen_max'],
            'incoming_coefficient': row['drive_amplitude'],
            'bulk_harmonic': int(kappa),
            'birth_E_per_unit_incoming_coefficient': mode['birth']['E'],
            'birth_B_per_unit_incoming_coefficient': mode['birth']['B'],
            'transverse_cdf_gap': result['cdf_gap'],
            'principal_pressure_max_share': result['stress_eigen_max'],
            'mean_energy_over_stationary_thermal_target': result['energy_over_thermal'],
        }
    return cases


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
        with tempfile.TemporaryDirectory(prefix="hdblast-prebirth-output-") as temp:
            replay(Path(temp), opt.jobs)
