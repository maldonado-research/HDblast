"""Regenerate and verify the finite-bath experiments from the public shell inputs."""
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
REPO = HERE.parent.parent
MODELS = ("physical", "sheared", "sheared_prepared")
ARMS = ("quiet", "plus", "minus", "far")


def compare_numbers(actual, expected, context, absolute=5e-9, relative=2e-8):
    if isinstance(expected, dict) and isinstance(actual, dict) and actual.keys() == expected.keys():
        for key in expected:
            compare_numbers(actual[key], expected[key], context + "/" + key, absolute, relative)
        return
    if isinstance(expected, list) and isinstance(actual, list) and len(actual) == len(expected):
        for i, (left, right) in enumerate(zip(actual, expected)):
            compare_numbers(left, right, f"{context}/{i}", absolute, relative)
        return
    if isinstance(expected, (int, float)) and isinstance(actual, (int, float)):
        if math.isfinite(actual) and math.isclose(actual, expected, abs_tol=absolute, rel_tol=relative):
            return
    elif type(actual) is type(expected) and actual == expected:
        return
    raise ValueError(f"{context}: calculated {actual!r}, archived {expected!r}")


def archived(path):
    return json.loads((HERE / path).read_bytes())


def execute(script, *args):
    flags = ["-O"] if sys.flags.optimize else []
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    raw = subprocess.check_output([sys.executable, "-I", "-B", *flags, HERE / script, *args], env=env)
    return json.loads(raw)


def energy_endpoints(report):
    output = {}
    for arm in ARMS:
        values = report["arms"][arm]
        row = {"contrasts": values["contrasts"]}
        for model in MODELS:
            run = values[model]
            row[model] = {"end": run["end"]["redshifted_density"],
                          "min_a4_H": run["min_signed_a4_radiation"]["radiation_a4"]["H_loss"]}
        output[arm] = row
    return output


def curve(path):
    with path.open(newline="") as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def compare_independent(primary, transfer):
    for arm in ("quiet", "plus"):
        for model, key in zip(MODELS, ("H", "B", "B_prepared")):
            exact_step = transfer["joint_refined"][arm + "_" + key]
            row = primary["arms"][arm][model]
            for left, right in (("H_loss", "R_H"), ("count", "R_count")):
                compare_numbers(row["end"]["redshifted_density"][left], exact_step[right],
                                arm + "/" + model + "/" + left, absolute=5e-5, relative=1e-7)
            compare_numbers(row["min_signed_a4_radiation"]["radiation_a4"]["H_loss"],
                            exact_step["min_a4_R_H"], arm + "/" + model + "/minimum",
                            absolute=1e-5, relative=1e-7)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="also recompute the Fock-space and time evolutions")
    opt = parser.parse_args()
    for line in (HERE/"SHA256SUMS").read_text().splitlines():
        digest, filename = line.split("  ", 1)
        if hashlib.sha256((HERE/filename).read_bytes()).hexdigest() != digest:
            raise ValueError("Input changed: " + filename)
    power = execute("controls/local_power.py", REPO)
    compare_numbers(power, archived("controls/POWER.json"), "ground-state power")
    primary = archived("evolution/EVOLUTION48.json")
    transfer = archived("independent/TRANSFER.json")
    compare_independent(primary, transfer)
    quadratures = [HERE/f"evolution/EVOLUTION{n}.json" for n in (24, 48, 96)]
    summary = execute("evolution/compare.py", *quadratures)
    compare_numbers(summary, archived("evolution/CHECK.json"), "quadrature controls")
    print("Verified local powers and independent endpoint and minimum-energy comparisons.", flush=True)
    if not opt.full:
        return
    fock = execute("collision/collision.py")
    compare_numbers(fock, archived("collision/COLLISION.json"), "finite collision")
    separate = execute("independent/bath_transfer.py", REPO)
    compare_numbers(separate, transfer, "independent transfer", absolute=3e-7)
    with tempfile.TemporaryDirectory() as work:
        generated = Path(work) / "quiet.csv"
        for order, tolerance in ((24, "1e-9"), (48, "2e-10"), (96, "2e-10")):
            args = ["--curves", generated] if order == 48 else []
            fresh = execute("evolution/evolve.py", REPO, "--order", str(order), "--tol", tolerance, *args)
            expected = archived(f"evolution/EVOLUTION{order}.json")
            compare_numbers(energy_endpoints(fresh), energy_endpoints(expected),
                            f"evolution order {order}", absolute=1e-6)
            if order == 48:
                compare_independent(fresh, separate)
        compare_numbers(curve(generated), curve(HERE/"evolution/QUIET48.csv"), "quiet time series", absolute=1e-6)
    print("Reproduced the finite collision, three mode quadratures, independent transfer and quiet curves.")


if __name__ == "__main__":
    main()
