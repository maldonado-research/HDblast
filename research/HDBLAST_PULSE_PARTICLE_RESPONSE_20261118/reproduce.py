"""Replay the pulse packet from a repository or an extracted Zenodo directory."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

CAPSULE = Path(__file__).resolve().parent
REPOSITORY = CAPSULE.parents[1]
INPUTS = Path("research/HDBLAST_INCIDENT_PULSE_20261118")


def json_file(name):
    return json.loads((CAPSULE / name).read_bytes())


def check_bytes():
    records = json_file("SOURCES.json")["sources"]
    for record in records:
        path = REPOSITORY / record["path"]
        raw = path.read_bytes()
        if len(raw) != record["bytes"] or hashlib.sha256(raw).hexdigest() != record["sha256"]:
            raise ValueError(f"Historical source differs: {path}")
    print(f"Verified {len(records)} historical input files", flush=True)

    manifest = CAPSULE / "SHA256SUMS"
    if manifest.is_file():
        for line in manifest.read_text().splitlines():
            digest, name = line.split(maxsplit=1)
            path = CAPSULE / name.removeprefix("*")
            if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                raise ValueError(f"Capsule file differs: {name}")
        print("Verified capsule SHA256SUMS", flush=True)


def compare(expected, actual, key, rel, absolute):
    """Check the full structure. Permit rounding for floating diagnostics and observables."""
    if isinstance(expected, dict) and isinstance(actual, dict) and expected.keys() == actual.keys():
        for name in expected:
            compare(expected[name], actual[name], f"{key}/{name}", rel, absolute)
        return
    if isinstance(expected, list) and isinstance(actual, list) and len(expected) == len(actual):
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f"{key}/{index}", rel, absolute)
        return
    if type(expected) is float and type(actual) in (float, int):
        if not (math.isfinite(actual) and math.isfinite(expected)):
            raise ValueError(f"Nonfinite numerical value at {key}")
        # The ratio of two small fourth-order benchmark errors loses several digits.
        if "BENCHMARK.json" in key and "/reduction/" in key:
            rel = 1e-4
        # Cross-platform summation of the wide-band difference of vacuum terms
        # and subtraction of near-equal positive energies lose absolute precision.
        if "/5-80/" in key and "energy_" in key:
            absolute = max(absolute, 4e-10)
        if math.isclose(expected, actual, rel_tol=rel, abs_tol=absolute):
            return
    elif type(expected) is type(actual) and expected == actual:
        return
    raise ValueError(f"Replay differs at {key}: expected {expected!r}, got {actual!r}")


def replay(script, record, *args, rel=2e-8, absolute=3e-12):
    flags = ["-O"] if sys.flags.optimize else []
    command = [sys.executable, "-I", "-B", *flags, CAPSULE / script, *args]
    result = json.loads(subprocess.check_output(command))
    compare(json_file(record), result, record, rel, absolute)
    print(f"Reproduced {record}", flush=True)
    return result


def cross_checks(modes, validity, detector, direct):
    for grid in ("g4", "g8"):
        for case, wall in (("plus", "delayed"), ("minus", "delayed_negative")):
            final = modes["results"][f"{grid}/{case}"]["0.15"]
            for seed, vacuum in (("instant", "instantaneous_vacuum"), ("wkb", "WKB_tangent")):
                band = final[seed]["band"]["5-18/32"]
                for key, arm in (("energy_excitation", wall), ("energy_excitation_quiet", "quiet")):
                    cap = validity[grid][arm]["excitation_envelopes"]["primary"]["upper"][vacuum]
                    if not 0 <= band[key] <= cap["rho_excitation_over_H0_fourth"]:
                        raise ValueError(f"Excitation bound exceeded: {grid}, {arm}, {seed}")
            if validity[grid][wall]["certified_comoving_k_ceiling"] < 18:
                raise ValueError("Primary band exceeds the on-shell frequency screen")

    # Two distinct integrators, clock interpolations and momentum quadratures.
    late = direct["s1200_k64"]
    production = detector["grids"]["g8"]["runs"]
    for name in ("plus", "minus", "half_plus", "half_minus"):
        measured = production[name]["probability_difference"]["band_5_18"][0][1][0]
        if not math.isclose(late[name]["difference"], measured, rel_tol=3e-5, abs_tol=1e-14):
            raise ValueError(f"Independent detector replay disagrees for {name}")
    print("Verified finite-band energy ceilings and four signed detector comparisons", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true",
                        help="recompute all quantum modes, detector controls and DOP853 check (NumPy and SciPy)")
    args = parser.parse_args()
    check_bytes()
    archive = REPOSITORY / INPUTS / "evolution/wall_histories.zip"
    source = REPOSITORY / INPUTS / "B1_summary.json"
    direct = replay("detector/replay_detector.py", "detector/REPLAY.json", archive, source,
                    absolute=1e-14)
    # make_pulse and bounds serialize their reported decimal values as exact strings.
    replay("scattering/make_pulse.py", "scattering/PULSES.json", archive)
    replay("scattering/bounds.py", "scattering/BOUNDS.json")
    replay("scattering/scatter.py", "scattering/SCATTER.json", absolute=1e-13)
    validity = replay("validity/validity.py", "validity/VALIDITY.json", archive, source)
    replay("modes/audit_modes.py", "modes/SUMMARY.json")
    if args.full:
        modes = replay("modes/spectator.py", "modes/MODES.json", archive)
        detector = replay("detector/detector.py", "detector/DETECTOR.json", archive, source)
        replay("modes/benchmark.py", "modes/BENCHMARK.json", absolute=1e-14)
        replay("far/sensitivity.py", "far/SENSITIVITY.json", rel=2e-4, absolute=1e-7)
        replay("scattering/independent.py", "scattering/scatter-independent.json", absolute=3e-13)
    else:
        modes = json_file("modes/MODES.json")
        detector = json_file("detector/DETECTOR.json")
    cross_checks(modes, validity, detector, direct)
    flags = ["-O"] if sys.flags.optimize else []
    full = ["--full"] if args.full else []
    for script in ("replay_rolling.py", "far/replay_far.py"):
        subprocess.run([sys.executable, "-I", "-B", *flags, CAPSULE / script, *full], check=True)


if __name__ == "__main__":
    main()
