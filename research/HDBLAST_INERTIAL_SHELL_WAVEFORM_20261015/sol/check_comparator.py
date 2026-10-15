#!/usr/bin/env python3
"""Fresh normal/-O comparator replay, signed anchors and mathematical faults."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
PRODUCER = HERE / "inertial_comparator.py"
SHELL_PATH = "inputs/previous/archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json"
SHELL_PIN = "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def invoke(source_root, optimized=False, fault=None, count=4096):
    command = [sys.executable, "-I", "-B"]
    if optimized:
        command.append("-O")
    command += [str(PRODUCER), "--source-root", str(source_root), "--mass-cells", str(count)]
    if fault:
        command += ["--fault", fault]
    return subprocess.run(command, text=True, capture_output=True, check=False)


def run(source_root):
    # A relocated minimal input tree proves there is no implicit dependency on
    # other capsule files or on the archived reference-only producer.
    data = (Path(source_root) / SHELL_PATH).read_bytes()
    require(sha256(data).hexdigest() == SHELL_PIN, "input authentication failed")
    controls = []
    with tempfile.TemporaryDirectory(prefix="comparator-minimal-", dir=HERE) as temporary:
        minimal = Path(temporary)
        target = minimal / SHELL_PATH
        target.parent.mkdir(parents=True)
        target.write_bytes(data)
        outputs = []
        for optimized in (False, True):
            result = invoke(minimal, optimized)
            require(result.returncode == 0, "positive replay failed: " + result.stderr)
            output = json.loads(result.stdout)
            require(output["python_optimization"] == int(optimized), "optimization metadata mismatch")
            output.pop("python_optimization")
            outputs.append(output)
            tail = [F(x) for x in output["H0_tail_continuous_2_3"]]
            center = [F(x) for x in output["H0_center"]]
            require(F(3, 10**8) < tail[0] < tail[1] < F(6, 10**8), "independent signed tail gate failed")
            require(center[1] < F(-4, 10**4), "independent center inertia gate failed")
            for fault, expected in (
                ("reverse_kernel_initial", "kernel derivative strictly positive"),
                ("reverse_kernel_sign", "positive postpulse integrand"),
                ("omit_local_inertia", "center retains negative local inertia"),
                ("reverse_q_derivative", "positive incoming lobe"),
            ):
                bad = invoke(minimal, optimized, fault, count=64)
                require(bad.returncode != 0, "fault unexpectedly passed: " + fault)
                require(expected in bad.stderr, "fault failed at unrelated gate: " + bad.stderr)
                controls.append({"fault": fault, "optimized": optimized, "returncode": bad.returncode, "diagnostic": bad.stderr.strip()})
        require(outputs[0] == outputs[1], "normal/-O mathematical outputs differ")
        # Authentication is an executable failure path in both modes.
        target.write_bytes(data + b"\n")
        for optimized in (False, True):
            bad = invoke(minimal, optimized, count=64)
            require(bad.returncode != 0 and "authentication failed" in bad.stderr, "changed input was not rejected")
            controls.append({"fault": "changed_authenticated_input", "optimized": optimized, "returncode": bad.returncode, "diagnostic": bad.stderr.strip()})
    return {
        "status": "PASS_INERTIAL_COMPARATOR_FRESH_REPLAY_AND_FAULTS",
        "python_optimization": sys.flags.optimize,
        "normal_optimized_mathematics_identical": True,
        "minimal_input_file_count": 1,
        "normal_optimized_control_count": len(controls),
        "controls": controls,
        "signed_comparator_tail": outputs[0]["H0_tail_continuous_2_3"],
        "producer_sha256": sha256(PRODUCER.read_bytes()).hexdigest(),
        "scope": "Comparator-only arithmetic audit. Does not certify the analytic coupled PDE assumptions or actual-minus-comparator estimate.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(run(args.source_root), indent=2, sort_keys=True))
    except (RuntimeError, ValueError, OSError, KeyError) as exc:
        print("FAIL: " + str(exc), file=sys.stderr)
        raise SystemExit(1)
