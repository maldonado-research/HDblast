"""Authenticate the closed nonlinear-constraint-ball capsule before replay.

The external manifest pin binds every local code, proof and input byte.
Finite replay does not prove the retained continuous analytic premises.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import re
import sys

REQUIRED_MEMBERS = {
    "SOURCE_PINS.json", "INPUT_PROVENANCE.json", "DATA.json", "requirements-checks.txt",
    "checks/check_checkpoint.py", "checks/evaluate_nonlinear_ball.py",
    "checks/run_controls.py", "checks/run_math_controls.py", "checks/constraint_algebra.py",
    "checks/weight_degree_ledger.py", "checks/ball_arithmetic.py", "checks/seed_envelope.py",
    "README.md", "CONTINUATION.md", "REPRODUCIBILITY.md", "AUDIT.md", "LITERATURE.md",
    "NONLINEAR_CONSTRAINT_BALL.md", "WEIGHTED_EMBEDDING_PRODUCTS.md",
    "PHYSICAL_CORRECTION_MAPS.md", "SEED_ENVELOPE.md",
}


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def load_json(raw):
    def unique(items):
        result = {}
        for key, value in items:
            need(key not in result, "Duplicate JSON key: " + key)
            result[key] = value
        return result

    def nonfinite(value):
        raise ValueError("Nonfinite JSON constant: " + value)

    return json.loads(raw.decode("utf-8"), object_pairs_hook=unique,
                      parse_float=str, parse_constant=nonfinite)


def regular(path):
    return path.is_file() and all(not item.is_symlink() for item in (path, *path.parents))


def authenticate(root, pin):
    need(re.fullmatch(r"[0-9a-f]{64}", pin) is not None,
         "Manifest pin must be 64 lowercase hexadecimal characters")
    manifest_path = root / "MANIFEST.json"
    need(regular(manifest_path), "Regular manifest required")
    raw = manifest_path.read_bytes()
    need(sha(raw) == pin, "Manifest SHA256 mismatch")
    manifest = load_json(raw)
    need(isinstance(manifest, dict) and set(manifest) == {"schema", "files"},
         "Unexpected manifest keys")
    need(manifest["schema"] == "HDBLAST_NONLINEAR_CONSTRAINT_BALL_MANIFEST_V1",
         "Manifest schema mismatch")
    members = manifest["files"]
    need(isinstance(members, dict) and members, "Nonempty file inventory required")
    for name, row in members.items():
        need(isinstance(name, str), "Manifest member name must be a string")
        path = PurePosixPath(name)
        need("\\" not in name and not path.is_absolute() and bool(path.parts)
             and path.as_posix() == name and ".." not in path.parts
             and name != "MANIFEST.json", "Unsafe manifest member")
        need(isinstance(row, dict) and set(row) == {"bytes", "sha256"},
             "Unexpected member metadata")
        need(type(row["bytes"]) is int and row["bytes"] >= 0
             and isinstance(row["sha256"], str)
             and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is not None,
             "Invalid member metadata")
    need(REQUIRED_MEMBERS <= set(members), "Required capsule core missing")
    actual = set()
    for path in root.rglob("*"):
        need(not path.is_symlink(), "Capsule symlink rejected")
        need(path.is_file() or path.is_dir(), "Nonregular capsule member")
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
    need(actual == set(members) | {"MANIFEST.json"}, "Closed file roster mismatch")
    for name, row in members.items():
        path = root / name
        need(regular(path), "Regular source required: " + name)
        data = path.read_bytes()
        need(len(data) == row["bytes"] and sha(data) == row["sha256"],
             "Source bytes differ: " + name)
        if path.suffix.lower() == ".json":
            load_json(data)
    return len(actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    script = Path(__file__).absolute()
    need(regular(script), "Regular checker path required")
    root = script.resolve(strict=True).parent.parent
    count = authenticate(root, args.manifest_sha256)
    output = args.output.absolute()
    need(".." not in output.parts, "Parent components in output path rejected")
    need(not output.is_relative_to(root), "Output must be outside capsule")
    need(not output.exists() and not output.is_symlink(), "Output already exists")
    need(output.parent.is_dir() and all(not item.is_symlink() for item in output.parents),
         "Regular existing output directory required")
    sys.dont_write_bytecode = True
    module_path = root / "checks/evaluate_nonlinear_ball.py"
    need(regular(module_path), "Regular pinned evaluator required")
    spec = importlib.util.spec_from_file_location("nonlinear_ball_evaluator", module_path)
    need(spec is not None and spec.loader is not None, "Evaluator loader unavailable")
    evaluator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(evaluator)
    mathematics = evaluator.evaluate(root)
    need(isinstance(mathematics, dict) and mathematics.get("status") == "PASS",
         "Mathematical evaluation failed")
    names = mathematics.get("named_checks")
    need(isinstance(names, list) and names
         and all(isinstance(name, str) and name for name in names)
         and len(set(names)) == len(names)
         and type(mathematics.get("named_check_count")) is int
         and mathematics["named_check_count"] == len(names),
         "Mathematical check roster differs")
    need(isinstance(mathematics.get("exact"), dict) and mathematics["exact"],
         "Exact arithmetic report required")
    report = {
        "schema": "HDBLAST_NONLINEAR_CONSTRAINT_BALL_CHECKPOINT_V1",
        "status": "PASS", "manifest_sha256": args.manifest_sha256,
        "closed_files_authenticated": count,
        "named_check_count": len(names), "mathematics": mathematics,
        "scope": "Closed-source authentication and finite exact constraint algebra, "
                 "weight allocations, source-model and numerical-ledger checks. "
                 "Weighted products, dotted closure, inverse and implicit-family "
                 "existence remain written conditional analytic arguments. "
                 "No spacetime evolution or physical confirmation follows.",
    }
    encoded = json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n"
    with output.open("x", encoding="utf-8") as stream:
        stream.write(encoded)
    print(json.dumps({"status": "PASS", "checks": len(names),
                      "manifest_sha256": args.manifest_sha256}, sort_keys=True))


if __name__ == "__main__":
    main()
