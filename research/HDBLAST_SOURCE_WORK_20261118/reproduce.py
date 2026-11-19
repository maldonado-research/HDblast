"""Run the fixed-history and reciprocal B1 tests from the published capsule."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / "HDBLAST_PULSE_PARTICLE_RESPONSE_20261118"


def read(name):
    return json.loads((HERE / name).read_bytes())


def compare(actual, expected, label, rel, absolute):
    if isinstance(expected, dict) and isinstance(actual, dict) and actual.keys() == expected.keys():
        for k in expected:
            compare(actual[k], expected[k], f"{label}/{k}", rel, absolute)
        return
    if isinstance(expected, list) and isinstance(actual, list) and len(actual) == len(expected):
        for i,(v,w) in enumerate(zip(actual, expected)):
            compare(v, w, f"{label}/{i}", rel, absolute)
        return
    if type(expected) is float and type(actual) in (float, int):
        if math.isfinite(actual) and math.isclose(actual, expected, rel_tol=rel, abs_tol=absolute):
            return
    elif type(actual) is type(expected) and expected == actual:
        return
    raise ValueError(f"Changed {label}: {actual!r}; expected {expected!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="regenerate numerical field, work and radiation reports")
    full = parser.parse_args().full
    opts = ["-O"] if sys.flags.optimize else []
    python = [sys.executable, "-I", "-B", *opts]
    args = ["--full"] if full else []
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    for script in ("replay_science.py", "coupled/verify_coupled.py"):
        subprocess.run([*python, HERE / script, *args], check=True, env=env)
    subprocess.run([*python, HERE / "coupled/bands/verify_bands.py"], check=True, env=env)

    def regenerate(script, report, positional=(), rel=1e-6, absolute=1e-10):
        data = json.loads(subprocess.check_output([*python, HERE/script, *positional], env=env))
        compare(data, read(report), report, rel, absolute)
        print(f"Reproduced {report}", flush=True)
    regenerate("interface/interface.py", "interface/INTERFACE.json", rel=1e-6, absolute=1e-10)
    if full:
        regenerate("review/check_interface.py","review/INTERFACE_CHECK.json", rel=.01, absolute=2e-9)
        regenerate("review/replay_native.py","review/NATIVE_TRANSFER.json",
                   (HERE/"coupled/WALLS.zip", PRIOR), rel=3e-6, absolute=2e-9)


if __name__ == "__main__":
    main()
