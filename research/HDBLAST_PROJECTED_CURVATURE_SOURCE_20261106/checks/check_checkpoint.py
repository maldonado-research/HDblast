"""Authenticate a closed projected-curvature-source capsule and evaluate finite arithmetic.

The continuous proper-cut and PDE premises are not verified by this runner.
An externally supplied manifest pin binds all capsule code, input and
document bytes. Canonical output must stay outside the capsule.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import re
import sys


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
                      parse_constant=nonfinite)


def regular(path):
    return path.is_file() and all(not p.is_symlink() for p in (path, *path.parents))


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
    need(manifest["schema"] == "HDBLAST_PROJECTED_CURVATURE_SOURCE_MANIFEST_V1",
         "Manifest schema mismatch")
    members = manifest["files"]
    need(isinstance(members, dict) and members, "Nonempty file inventory required")
    for name, row in members.items():
        path = PurePosixPath(name)
        need("\\" not in name and not path.is_absolute()
             and path.as_posix() == name and ".." not in path.parts
             and name != "MANIFEST.json", "Unsafe manifest member")
        need(isinstance(row, dict) and set(row) == {"bytes", "sha256"},
             "Unexpected member metadata")
        need(type(row["bytes"]) is int and row["bytes"] >= 0
             and isinstance(row["sha256"], str)
             and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is not None,
             "Invalid member metadata")
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
    need(output.parent.is_dir() and all(not p.is_symlink() for p in
         (output.parent, *output.parents)), "Regular existing output directory required")
    # Keep a closed capsule closed even if the caller omits recommended -B.
    sys.dont_write_bytecode = True
    module_path = root / "evaluate_curvature_source.py"
    need(regular(module_path), "Regular pinned evaluator required")
    spec = importlib.util.spec_from_file_location("curvature_source_evaluator", module_path)
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
        "schema": "HDBLAST_PROJECTED_CURVATURE_SOURCE_CHECKPOINT_V1",
        "status": "PASS", "manifest_sha256": args.manifest_sha256,
        "closed_files_authenticated": count,
        "named_check_count": mathematics["named_check_count"],
        "mathematics": mathematics,
        "scope": "Closed-source authentication and finite exact arithmetic only. "
                 "The retained actual master certificate is consumed, not rerun. "
                 "The conditional Gauss decomposition, actual first-scalar source "
                 "bounds and positive weighted quadratic Ricci-source lower bound "
                 "are evaluated at amplitude zero. The scalar-second-variation "
                 "and Weyl compensation is necessary only if D2 is nonnegative. "
                 "Continuous derivation, state/gauge, reconstruction and PDE "
                 "arguments need separate review. No D2 sign, actual second-Ricci "
                 "variation bound, nonzero amplitude radius or physical energy "
                 "identification follows.",
    }
    encoded = json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n"
    with output.open("x", encoding="utf-8") as stream:
        stream.write(encoded)
    print(json.dumps({"status": "PASS", "checks": mathematics["named_check_count"],
                      "manifest_sha256": args.manifest_sha256}, sort_keys=True))


if __name__ == "__main__":
    main()
