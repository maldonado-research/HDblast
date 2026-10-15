"""Authenticated fresh replay of three exact boundary-algebra checkers.

The six positive processes and twenty-two injected-corruption processes do
not constitute a computer-certified proof of the analytic PDE estimates.
"""

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time


CHECKERS = [
    {
        "name": "astra",
        "script": "astra/verify_uniform_boundary.py",
        "golden": "astra/NORMAL.json",
        "controls": [
            ("--corrupt", "robin_sign", "Robin_to_Neumann_boundary_jet"),
            ("--corrupt", "lapse_factor", "actual_tangent_row_00"),
            ("--corrupt", "claimed_symmetry", "false_full_Robin_symmetry_claim"),
        ],
    },
    {
        "name": "sol",
        "script": "sol/verify_uniform_robin.py",
        "golden": "sol/verification_results.json",
        "controls": [
            ("--inject-corruption", "boundary-sign", "B is minus Jacobian of the exact boundary RHS"),
            ("--inject-corruption", "conjugation-sign", "correct conjugation sign"),
            ("--inject-corruption", "omitted-quadratic-derivative", "all equation (5) metric-wave sources"),
            ("--inject-corruption", "noncommuting-time", "noncommuting Duhamel derivative"),
        ],
    },
    {
        "name": "referee",
        "script": "referee/check_boundary_matrix_and_glancing.py",
        "golden": "referee/BOUNDARY_CHECK_RESULTS.json",
        "controls": [
            ("--inject-corruption", "wrong-linear-scalar", "exact direct-metric linearization disagrees with Robin matrix"),
            ("--inject-corruption", "wrong-trace-factor", "exact direct-metric linearization disagrees with Robin matrix"),
            ("--inject-corruption", "wrong-glancing-rank", "glancing harmonic constraint rank is not four"),
            ("--inject-corruption", "wrong-conjugation", "Robin conjugation leaves a nonzero boundary term"),
        ],
    },
]


def file_record(path, relative_to=None):
    raw = path.read_bytes()
    return {
        "path": path.relative_to(relative_to).as_posix() if relative_to else str(path),
        "bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
    }


def scientific(checker, data):
    # This is the only ignored value: it intentionally records Python -O.
    if checker == "astra":
        return {key: value for key, value in data.items() if key != "optimized_python"}
    return data


def run_checker(root, output, checker, mode, control=None):
    suffix = control[1] if control else "positive"
    work = output / mode / checker["name"] / suffix
    work.mkdir(parents=True)
    source = root / checker["script"]
    script = work / source.name
    shutil.copy2(source, script)
    if script.read_bytes() != source.read_bytes():
        raise RuntimeError("Fresh script copy differs: " + checker["script"])
    python_flags = ["-I", "-B"] + (["-O"] if mode == "optimized" else [])
    arguments = ["--output", "fresh.json"]
    if control:
        arguments += list(control[:2])
    started = time.monotonic()
    with (work / "execution.log").open("w", encoding="utf-8") as log:
        result = subprocess.run(
            [sys.executable, *python_flags, str(script), *arguments],
            cwd=work, stdout=log, stderr=subprocess.STDOUT,
        )
    row = {
        "kind": "actual_corruption" if control else "positive",
        "mode": mode,
        "script": checker["script"],
        "source_sha256": file_record(source)["sha256"],
        "python_flags": python_flags,
        "arguments": arguments,
        "exit_code": result.returncode,
        "wall_seconds": round(time.monotonic() - started, 3),
        "copied_script": file_record(script, output),
        "log": file_record(work / "execution.log", output),
    }
    fresh = work / "fresh.json"
    if fresh.is_file():
        row["output"] = file_record(fresh, output)
    if control:
        diagnostic = (work / "execution.log").read_text(encoding="utf-8")
        if result.returncode != 1 or control[2] not in diagnostic:
            raise RuntimeError("Actual corruption did not reject for the expected reason: "
                               + checker["name"] + "/" + control[1])
        if fresh.is_file() and json.loads(fresh.read_text(encoding="utf-8")).get("status") not in ("FAIL", "fail"):
            raise RuntimeError("Corrupt output did not report failure: " + checker["name"])
        row.update({
            "corruption": control[1],
            "expected_failure_exit_code": 1,
            "expected_reason": control[2],
            "expected_reason_detected": True,
        })
    else:
        if result.returncode != 0 or not fresh.is_file():
            raise RuntimeError("Fresh positive checker failed: " + checker["script"])
        actual = json.loads(fresh.read_text(encoding="utf-8"))
        golden = root / checker["golden"]
        expected = json.loads(golden.read_text(encoding="utf-8"))
        if actual.get("status") not in ("PASS", "pass"):
            raise RuntimeError("Fresh checker did not report success: " + checker["script"])
        if scientific(checker["name"], actual) != scientific(checker["name"], expected):
            raise RuntimeError("Scientific output differs: " + checker["script"])
        if checker["name"] == "astra" and actual.get("optimized_python") is not (mode == "optimized"):
            raise RuntimeError("Astra optimization flag differs from the process mode")
        row.update({
            "golden": file_record(golden, root),
            "scientific_output_matches": True,
            "ignored_top_level_fields": ["optimized_python"] if checker["name"] == "astra" else [],
        })
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-sha256", required=True,
                        help="Externally pinned SHA256 of the frozen MANIFEST.json")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="New or empty directory outside the frozen capsule")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location("angular_contact_integrity", root / "verify_manifest.py")
    integrity = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(integrity)
    authentication = integrity.verify(root, args.manifest_sha256)
    requested_output = args.output_dir.absolute()
    if requested_output.is_symlink():
        raise RuntimeError("A real output directory is required")
    output = requested_output.resolve()
    if output.is_relative_to(root):
        raise RuntimeError("Output directory must be outside the frozen capsule")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise RuntimeError("New or empty output directory required")
    sympy_version = importlib.metadata.version("sympy")
    if sympy_version != "1.14.0":
        raise RuntimeError("Sympy 1.14.0 is required for this frozen replay")
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for mode in ("normal", "optimized"):
        for checker in CHECKERS:
            rows.append(run_checker(root, output, checker, mode))
        for checker in CHECKERS:
            for control in checker["controls"]:
                rows.append(run_checker(root, output, checker, mode, control))
    positive = sum(row["kind"] == "positive" for row in rows)
    corruptions = sum(row["kind"] == "actual_corruption" for row in rows)
    if positive != 6 or corruptions != 22:
        raise RuntimeError("Incomplete fresh replay")
    receipt = {
        "status": "PASS",
        "capsule": "HDBLAST_ANGULAR_UNIFORM_CONTACT_20261015",
        "authentication": authentication,
        "runner": file_record(Path(__file__).resolve(), root),
        "integrity_verifier": file_record(root / "verify_manifest.py", root),
        "requirements": file_record(root / "requirements.txt", root),
        "python": {"version": sys.version, "executable": sys.executable,
                   "implementation": sys.implementation.name},
        "sympy_version": sympy_version,
        "successful_checker_runs": positive,
        "actual_corruption_controls": corruptions,
        "steps": rows,
        "scope": "Fresh exact boundary, conjugation, and symbol algebra replay. "
                 "The written analytic estimates are separately reviewed mathematical arguments; "
                 "these process results do not computer-certify them or establish an exact nonlinear spacetime.",
    }
    (output / "REPRODUCTION_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "successful_checker_runs": positive,
                      "actual_corruption_controls": corruptions, "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
