"""Authenticate the collision study or regenerate its spectra, rates, and equilibrium tests."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
BASE = "49756c922760880a9f5da0e4da6f2d8cce36738a"
PINS_SHA = "31bcda3be54ed06584a44f19ff0dc6ad6c8717425f95ba2a931986897751c77b"
NATIVE = "research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native"
QUANTUM = "research/HDBLAST_QUANTUM_MEMORY_20261119/quantum"
INPUTS = {f"{NATIVE}/FORK.json", f"{NATIVE}/WALLS.zip",
          f"{QUANTUM}/spectrum.py", f"{QUANTUM}/SPECTRUM.json"}
OUTPUTS = {
    "daughters": "DISTRIBUTIONS.json",
    "transport": "CONTACT.json",
    "envelope": "ENVELOPE.json",
    "shape": "SHAPE.json",
    "survival": "SURVIVAL.json",
}


def check(ok, message):
    if not ok:
        raise ValueError(message)


def load(path):
    return json.loads(path.read_bytes())


def checksum(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def authenticate(repo):
    manifest = HERE / "SHA256SUMS"
    listed = set()
    for line in manifest.read_text().splitlines():
        digest, name = line.split("  ", 1)
        path = (HERE / name).resolve()
        check(path.is_relative_to(HERE) and path != manifest and path not in listed and
              re.fullmatch(r"[0-9a-f]{64}", digest), f"Invalid manifest record {name}")
        check(checksum(path) == digest, f"Study checksum differs for {name}")
        listed.add(path)
    physical = {p.resolve() for p in HERE.rglob("*") if p.is_file() and p != manifest
                and "__pycache__" not in p.relative_to(HERE).parts}
    check(physical == listed, f"Study inventory mismatch: {physical ^ listed}")
    check(checksum(HERE / "SOURCES.json") == PINS_SHA, "Historical source manifest checksum")
    pins = load(HERE / "SOURCES.json")
    check(pins["base"] == BASE, "Historical revision mismatch")
    sources = pins["sources"]
    check(len(sources) == len(INPUTS) and {s["path"] for s in sources} == INPUTS,
          "Historical source inventory mismatch")
    for row in sources:
        path = (repo / row["path"]).resolve()
        check(path.is_relative_to(repo) and path.stat().st_size == row["bytes"]
              and checksum(path) == row["sha256"], f"Historical input mismatch at {row['path']}")
    print(f"Authenticated {len(listed)} study files and {len(sources)} historical inputs.", flush=True)


def compare(actual, expected, path):
    if isinstance(expected, dict):
        check(isinstance(actual, dict) and expected.keys() == actual.keys(), f"Field names at {path}")
        for key in expected:
            compare(actual[key], expected[key], (*path, key))
    elif isinstance(expected, list):
        check(isinstance(actual, list) and len(actual) == len(expected), f"List size at {path}")
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, (*path, str(i)))
    elif type(expected) in (float, int) and type(actual) in (float, int):
        # The physical rates and fractions include very small positive values.
        # Solver-to-archive residuals and convergence diagnostics measure roundoff
        # from two separate computations and have their own absolute scales.
        atol = 3e-15
        if "relative_archived_energy_residual" in path:
            atol = 1e-9
        elif any(word in key for key in path for word in (
                "refinement", "residual", "discrepancy", "closure", "disagreement", "roundoff",
                "relative_difference", "angular_moment_relative_error", "band_additivity")):
            atol = 2e-12
        check(math.isfinite(actual) and math.isclose(actual, expected, rel_tol=3e-8, abs_tol=atol),
              f"Numeric difference at {'/'.join(path)}: {actual} vs {expected}")
    else:
        check(type(actual) is type(expected) and actual == expected, f"Value at {path}")


def compare_spectra(output):
    # The last bins include small differences of antiderivatives. Compare their
    # absolute errors as well as all positive bin masses and differential values.
    with output.open(newline="") as new, (HERE / "daughters/SPECTRA.csv").open(newline="") as old:
        reader = csv.DictReader(new)
        reference = csv.DictReader(old)
        check(reader.fieldnames == reference.fieldnames, "Spectral CSV fields")
        rows = []
        for i, (a, b) in enumerate(zip(reader, reference, strict=True)):
            check((a["arm"], a["band"]) == (b["arm"], b["band"]), f"Spectral CSV row {i}")
            for key in reader.fieldnames[2:]:
                x, y = float(a[key]), float(b[key])
                check(math.isfinite(x) and math.isclose(x, y, rel_tol=3e-8, abs_tol=3e-14),
                      f"Spectral CSV row {i}/{key}: {x} vs {y}")
            rows.append(a)
    return rows


def physical_checks(data, rows):
    arms = {"quiet", "plus", "minus", "far"}
    for part, result in data.items():
        check(set(result.get("arms", result.get("histories"))) == arms, f"History coverage in {part}")
    y = data["survival"]["illustrative_y"]
    spin = data["survival"]["spin_dof_per_sign"]
    for arm in sorted(arms):
        density = data["daughters"]["arms"][arm]["bands"]["full"]["one_sign_comoving_moments"]
        n, energy, inv = (density[k] for k in ("0", "1", "-1"))
        mean = energy / n
        env = data["envelope"]["histories"][arm]
        contact = data["transport"]["arms"][arm]["readouts"]["2.16"]
        cap = data["survival"]["histories"][arm]
        shape = data["shape"]["arms"][arm]["full_band_chemical_comparison"]["elastic_number_preserving_FD_example"]
        def agree(a, b, description):
            check(math.isclose(a, b, rel_tol=5e-8, abs_tol=2e-15),
                  f"{arm} {description}: {a} vs {b}")
        agree(energy, env["moments_per_species"]["all"]["M1"], "independent emitted energy")
        agree(inv, env["moments_per_species"]["all"]["Mminus1"], "independent inverse-energy moment")
        agree(8/3 * mean * energy * contact["wall"]["J"],
              contact["cases"]["full_on_full"]["mean_optical_depth_per_lambda"], "mean contact depth")
        agree(shape["fugacity"], cap["thermal_fugacity"], "Fermi equilibrium fugacity")
        agree(shape["comoving_temperature_H0"], cap["thermal_comoving_scale"], "Fermi equilibrium scale")
        # Fill the minimum-momentum Fermi sphere and eliminate its radius to
        # evaluate an independent, evolving-target first-event upper bound.
        max_inverse = 1.5 * (spin / (6 * math.pi**2))**(1/3) * (y**2 * n)**(2/3)
        eta = env["readouts"]["2.16"]["conformal_path"]
        coarse_cap = -math.expm1(-8*math.pi * max_inverse * eta * inv / n)
        agree(coarse_cap, cap["upper_first_event_fraction_number"], "Pauli first-event bound")
        # Enclose the integrated feedback inequality by left and right sums
        # of its decreasing reciprocal slope, using the other module's moments.
        evolving = env["readouts"]["2.16"]["evolving_targets"]
        sharp = evolving["first_collision_fraction_number_upper"]
        exposure = 8*math.pi * eta * y**2 * inv**2 / n
        coefficient = max_inverse / (y**2 * inv)
        subdivisions = 2048
        width = sharp/subdivisions
        values = [1/(1 + coefficient*(i*width)**(2/3)) for i in range(subdivisions+1)]
        check(width*math.fsum(values[1:]) <= exposure <= width*math.fsum(values[:-1]),
              f"{arm} exposure lies outside independent quadrature bounds")
        check(0 < sharp < coarse_cap and sharp < 74e-6 and coarse_cap < 0.00273,
              f"{arm} final first-collision fraction bounds")
        # Independent finite-bin Maxwell-Boltzmann CDF and total-variation sum.
        bins = [r for r in rows if r["arm"] == arm and r["band"] == "full"]
        def mb_tail(k):
            b = 3*k/mean
            return math.exp(-b) * (1 + b + b*b/2)
        excess = math.fsum(max(0., float(r["one_sign_comoving_number"])/n
                               - (mb_tail(float(r["k_low"])) - mb_tail(float(r["k_high"]))))
                           for r in bins)
        agree(excess, cap["bin_widths"]["0.5"]["initial_TV_MB"], "classical equilibrium distance")
        bench = cap["bin_widths"]["0.5"]
        check(bench["final_TV_FD_lower"] >= bench["initial_TV_FD"] - coarse_cap - 2e-13,
              f"{arm} evolving-target bound compared with equilibrium")


def verify(repo):
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
    with tempfile.TemporaryDirectory(prefix="hdblast-kinetic-") as folder:
        temp = Path(folder)
        data = {}
        for section, name in OUTPUTS.items():
            file = temp / name
            target = temp if section == "daughters" else file
            args = [sys.executable, "-I", "-B", str(HERE / section / (section + ".py")),
                    "--repo", str(repo), "--output", str(target)]
            if section == "shape":
                args += ["--daughters", str(temp / OUTPUTS["daughters"])]
            print(f"Recomputing {section}/{name}", flush=True)
            subprocess.run(args, env=env, check=True, cwd=temp)
            data[section] = load(file)
            compare(data[section], load(HERE / section / name), (section,))
        rows = compare_spectra(temp / "SPECTRA.csv")
        physical_checks(data, rows)
    print(f"Reproduced five kinetic calculations and {len(rows)} daughter-spectrum bins on four walls.")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--repo", type=Path, default=HERE.parents[1],
                     help="Root containing the pinned research/ inputs")
    cli.add_argument("--verify", action="store_true", help="Regenerate and compare all numerical results")
    args = cli.parse_args()
    repo = args.repo.resolve()
    authenticate(repo)
    if args.verify:
        verify(repo)
