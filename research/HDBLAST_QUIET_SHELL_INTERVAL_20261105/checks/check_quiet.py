"""Authenticate one closed capsule and check its finite exact arithmetic.

This does not verify PDE existence, localized uniqueness, the nonlinear
speed bound, physical embedding, a parameter radius, or originality.
Canonical output must be outside the capsule, preventing self-reference.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path, PurePosixPath
import re


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
    return json.loads(raw.decode("utf-8"), object_pairs_hook=unique)


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
    need(set(manifest) == {"schema", "files"}, "Unexpected manifest keys")
    need(manifest["schema"] == "HDBLAST_QUIET_SHELL_MANIFEST_V1",
         "Manifest schema mismatch")
    members = manifest["files"]
    need(isinstance(members, dict) and members, "Nonempty file inventory required")
    for name, row in members.items():
        path = PurePosixPath(name)
        need(isinstance(name, str) and "\\" not in name and not path.is_absolute()
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
    return len(actual)


def exact_checks(data):
    expected = {
        "ell": "1/1000000", "full_initial_collar_in_ell": "1/2",
        "actual_radial_speed_bound": "2", "quiet_reset_in_ell": "1/4",
        "total_q_horizon_in_ell": "25/16",
        "coordinate_slab_horizon_in_ell": "13/8",
        "active_q_duration_in_ell": "21/16", "duration_ratio": "21/25",
        "first_variation_ratio_bound": "21/25",
        "second_variation_H_term_ratio_bound": "9261/15625",
        "second_variation_C_term_ratio_bound": "441/625",
        "second_variation_P2_term_ratio_bound": "21/25",
        "old_clock_guard": "1/50", "delayed_clock_guard": "1/42",
        "retained_slab_margin_in_ell": "1/32",
    }
    need(set(data) == set(expected) | {"schema"}, "Unexpected input keys")
    need(data["schema"] == "HDBLAST_QUIET_SHELL_INPUTS_V1", "Input schema mismatch")
    need(all(isinstance(data[key], str) for key in expected),
         "Rational inputs must be strings")
    x = {key: Q(data[key]) for key in expected}
    need(x == {key: Q(value) for key, value in expected.items()},
         "Registered arithmetic contract differs")
    passed = []

    def check(name, condition):
        need(condition, "Exact arithmetic failed: " + name)
        passed.append(name)

    ell = x["ell"]
    d = ell * x["full_initial_collar_in_ell"]
    speed = x["actual_radial_speed_bound"]
    reset = ell * x["quiet_reset_in_ell"]
    total = ell * x["total_q_horizon_in_ell"]
    slab = ell * x["coordinate_slab_horizon_in_ell"]
    active = ell * x["active_q_duration_in_ell"]
    theta = x["duration_ratio"]
    check("full_collar_precedes_first_seed", 0 < d == ell / 2 < ell)
    check("protected_wedge_front", d - speed * reset == 0)
    check("positive_active_horizon", 0 < reset < total and active == total - reset)
    check("proper_duration_ratio_cancels_fixed_radius", active / total == theta)
    check("phi_convexity_input", 0 < theta < 1)
    check("first_variation_factor", x["first_variation_ratio_bound"] == theta)
    check("second_variation_factors",
          x["second_variation_H_term_ratio_bound"] == theta ** 3
          and x["second_variation_C_term_ratio_bound"] == theta ** 2
          and x["second_variation_P2_term_ratio_bound"] == theta)
    check("old_and_delayed_clock_thresholds",
          x["old_clock_guard"] == ell / (32 * total)
          and x["delayed_clock_guard"] == ell / (32 * active)
          and x["delayed_clock_guard"] > x["old_clock_guard"])
    margin = ell * x["retained_slab_margin_in_ell"]
    check("strict_guard_leaves_requested_slab_margin",
          slab - total == 2 * margin
          and x["delayed_clock_guard"] * active == margin
          and 1 - x["delayed_clock_guard"] == Q(41, 42))
    # Uniform S3 coordinate moments: E[n_i^4]=3/(4*6),
    # E[n_i^2 n_j^2]=1/(4*6). Volume(S3)=2*pi^2.
    mean_y2 = 2 * Q(3, 24) - 2 * Q(1, 24)
    check("round_quadrupole_denominator", mean_y2 == Q(1, 6)
          and 2 * mean_y2 == Q(1, 3) and mean_y2 > 0)
    check("quiet_actual_area_quotient", Q(0) / mean_y2 == 0)
    return passed, x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    script = Path(__file__).absolute()
    need(regular(script), "Regular checker path required")
    root = script.resolve(strict=True).parent.parent
    count = authenticate(root, args.manifest_sha256)
    passed, x = exact_checks(load_json((root / "INPUTS.json").read_bytes()))
    output = args.output.absolute()
    need(".." not in output.parts, "Parent components in output path rejected")
    need(not output.is_relative_to(root), "Output must be outside capsule")
    need(not output.exists() and not output.is_symlink(), "Output already exists")
    need(output.parent.is_dir() and all(not p.is_symlink() for p in
         (output.parent, *output.parents)), "Regular existing output directory required")
    report = {
        "schema": "HDBLAST_QUIET_SHELL_EXACT_REPORT_V1", "status": "PASS",
        "manifest_sha256": args.manifest_sha256,
        "closed_files_authenticated": count, "named_checks": passed,
        "named_check_count": len(passed),
        "exact": {key: str(value) for key, value in x.items()},
        "sphere_integral_Y2": "pi^2/3",
        "scope": "Finite exact arithmetic and closed-source authentication only. "
                 "The continuous quiet-wedge, endpoint and convexity arguments "
                 "require review of the proof. PDE existence, localized uniqueness, "
                 "actual nonlinear speed/state bounds, embedding, coefficient "
                 "suprema, parameter amplitude and originality are not verified.",
    }
    with output.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(report, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "checks": len(passed),
                      "manifest_sha256": args.manifest_sha256}, sort_keys=True))


if __name__ == "__main__":
    main()
