"""Check the study files and physical controls; --full regenerates all numerical studies."""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
DEFAULT_ARCHIVE = HERE.parent / "HDBLAST_COHERENT_DECAY_20261120"


def check_manifest():
    entries = (HERE / "SHA256SUMS").read_text().splitlines()
    for entry in entries:
        expected, name = entry.split("  ", 1)
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != expected:
            raise ValueError("File checksum differs: " + name)
    print(f"Checked {len(entries)} file hashes.", flush=True)


def run(script, *args, stdout=None):
    subprocess.run([sys.executable, "-B", str(HERE / script), *map(str, args)],
                   cwd=HERE, check=True, stdout=stdout)


def compare(actual, reference, path="result"):
    if isinstance(reference, dict):
        if actual.keys() != reference.keys():
            raise ValueError(f"Fields differ in {path}")
        for key in reference:
            # Generated native binary traces can differ by floating-point library.
            # Source digests and all reported witness measurements are still compared.
            if key == "trace_sha256" and re.fullmatch(r"[0-9a-f]{64}", actual[key]):
                continue
            compare(actual[key], reference[key], f"{path}.{key}")
    elif isinstance(reference, list):
        if len(actual) != len(reference):
            raise ValueError(f"Length differs in {path}")
        for i, value in enumerate(reference):
            compare(actual[i], value, f"{path}[{i}]")
    elif isinstance(reference, float):
        if not math.isclose(actual, reference, rel_tol=2e-5, abs_tol=2e-10):
            raise ValueError(f"Numerical mismatch in {path}: {actual} vs {reference}")
    elif actual != reference:
        raise ValueError(f"Value differs in {path}")


def matches(result, stored):
    compare(json.loads(result.read_text()), json.loads((HERE / stored).read_text()), stored)
    print("Reproduced " + stored, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    args = parser.parse_args()
    archive = args.archive.resolve()
    check_manifest()
    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        drag = work / "drag.json"
        with drag.open("w") as dest:
            run("transport/drag.py", archive / "kinetics/SWEEP.json", stdout=dest)
        matches(drag, "transport/DRAG.json")
        run("memory/check.py")
        if args.full:
            pulse = work / "pulse.json"
            run("pulse/pulse_modes.py", "--output", pulse)
            matches(pulse, "pulse/PULSE_MODES.json")
            angles = work / "angles.json"
            slices = work / "slices"
            run("replay/replay.py", "--archive", archive, "--output", angles, "--slices", slices)
            matches(angles, "replay/results.json")
            independent = work / "independent.json"
            run("independent/check.py", "--source", archive, "--compare", angles, "--output", independent)
            matches(independent, "independent/CHECK.json")
            memory = work / "memory.json"
            run("memory/reproduce.py", "--slices", slices, "--reference", archive / "kinetics/SWEEP.json",
                "--output", memory)
            matches(memory, "memory/MEMORY.json")
            moment = work / "moment.json"
            run("replay/moment.py", "--slices", slices, "--output", moment)
            matches(moment, "replay/MOMENT.json")


if __name__ == "__main__":
    main()
