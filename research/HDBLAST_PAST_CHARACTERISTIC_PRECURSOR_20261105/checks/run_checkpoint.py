"""Authenticate and reproduce the compact past-precursor arithmetic capsule."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def relative(name):
    need(isinstance(name, str) and name and "\\" not in name, "RELATIVE_PATH")
    path = PurePosixPath(name)
    need(not path.is_absolute() and path.as_posix() == name
         and all(part not in ("", ".", "..") for part in name.split("/")),
         "RELATIVE_PATH")
    return Path(*path.parts)


def regular(path):
    return path.is_file() and all(not p.is_symlink() for p in (path, *path.parents))


def authenticate(expected_pin=None):
    manifest = ROOT / "MANIFEST.json"
    need(regular(manifest), "MANIFEST_REGULAR_FILE")
    raw = manifest.read_bytes()
    if expected_pin is not None:
        need(sha(raw) == expected_pin, "MANIFEST_PIN")
    data = json.loads(raw)
    need(data["schema"] == "HDBLAST_PAST_PRECURSOR_INVENTORY_V1", "MANIFEST_SCHEMA")
    rows = data["files"]
    need(isinstance(rows, dict) and rows and "MANIFEST.json" not in rows,
         "MANIFEST_MEMBERS")
    actual = set()
    for path in ROOT.rglob("*"):
        need(not path.is_symlink(), "CAPSULE_SYMLINK")
        if path.is_file():
            actual.add(path.relative_to(ROOT).as_posix())
        else:
            need(path.is_dir(), "CAPSULE_SPECIAL_FILE")
    need(actual == set(rows) | {"MANIFEST.json"}, "CLOSED_CAPSULE_INVENTORY")
    for name, row in rows.items():
        path = ROOT / relative(name)
        need(regular(path), "MEMBER_REGULAR_FILE:" + name)
        content = path.read_bytes()
        need(type(row["bytes"]) is int and len(content) == row["bytes"]
             and sha(content) == row["sha256"], "MEMBER_BYTES:" + name)
    return sha(raw), rows


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--manifest-sha256")
    parser.add_argument("--authenticate-only", action="store_true")
    args = parser.parse_args()
    pin, members = authenticate(args.manifest_sha256)
    if args.authenticate_only:
        need(args.output_dir is None, "AUTHENTICATION_HAS_NO_OUTPUT_DIRECTORY")
        print(json.dumps({"status": "PASS_CLOSED_PAST_PRECURSOR_SOURCE",
                          "manifest_sha256": pin, "files": len(members) + 1}, sort_keys=True))
        return
    need(args.output_dir is not None, "OUTPUT_DIRECTORY_REQUIRED")
    out = args.output_dir.absolute()
    need(all(not p.is_symlink() for p in (out, *out.parents)), "OUTPUT_SYMLINK")
    out = out.resolve()
    need(not out.exists(), "FRESH_OUTPUT_DIRECTORY")
    need(not out.is_relative_to(ROOT), "OUTPUT_OUTSIDE_CAPSULE")
    plan = json.loads((ROOT / "checks/REPLAY_PLAN.json").read_bytes())
    need(plan["schema"] == "HDBLAST_PAST_PRECURSOR_REPLAY_PLAN_V1", "PLAN_SCHEMA")
    components = plan["components"]
    names = [component["name"] for component in components]
    need(names and len(names) == len(set(names)), "UNIQUE_COMPONENT_NAMES")
    out.mkdir(parents=True)
    pairs = []
    for component in components:
        name = component["name"]
        need(relative(name).as_posix() == name and "/" not in name, "COMPONENT_NAME")
        script = ROOT / relative(component["script"])
        need(component["script"] in members and regular(script), "PINNED_COMPONENT")
        outputs = []
        for optimized, mode in ((False, "NORMAL"), (True, "OPTIMIZED")):
            target = out / (name + "_" + mode + ".json")
            argv = []
            for argument in component["arguments"]:
                need(isinstance(argument, str), "COMPONENT_ARGUMENT")
                if argument == "{output}":
                    argument = str(target)
                elif argument.startswith("{capsule}/"):
                    argument = str(ROOT / relative(argument[len("{capsule}/"):]))
                need("{output}" not in argument and "{capsule}" not in argument,
                     "UNKNOWN_ARGUMENT_PLACEHOLDER")
                argv.append(argument)
            command = [sys.executable, "-I", "-S", "-B"]
            if optimized:
                command.append("-O")
            result = subprocess.run(command + [str(script)] + argv,
                                    capture_output=True, timeout=90, check=False)
            need(result.returncode == 0, "COMPONENT_EXIT:" + name + ":" + mode)
            need(result.stderr == b"", "COMPONENT_STDERR:" + name + ":" + mode)
            need(len(result.stdout.splitlines()) == 1,
                 "COMPONENT_SINGLE_SUMMARY:" + name + ":" + mode)
            summary = json.loads(result.stdout)
            need(isinstance(summary, dict) and regular(target), "COMPONENT_REPORT:" + name)
            report = target.read_bytes()
            parsed = json.loads(report)
            reference_name = component["references"][mode]
            need(reference_name in members, "PINNED_REFERENCE:" + name)
            reference = (ROOT / relative(reference_name)).read_bytes()
            need(report == reference, "CANONICAL_REPORT_BYTES:" + name + ":" + mode)
            need(parsed["status"] == summary["status"], "COMPONENT_SUMMARY_STATUS:" + name)
            outputs.append({"mode": mode, "report": target.name,
                            "bytes": len(report), "sha256": sha(report),
                            "status": parsed["status"]})
        need(outputs[0]["sha256"] == outputs[1]["sha256"], "NORMAL_OPTIMIZED_EQUAL:" + name)
        pairs.append({"component": name, "reports": outputs,
                      "normal_optimized_identical": True})
    after_pin, after_members = authenticate(pin)
    need(after_pin == pin and after_members == members, "CAPSULE_UNCHANGED_AFTER_REPLAY")
    receipt = {
        "schema": "HDBLAST_PAST_PRECURSOR_REPLAY_RECEIPT_V1",
        "status": "PASS_PAST_PRECURSOR_EXACT_REPLAY",
        "manifest_sha256": pin, "closed_files_authenticated": len(members) + 1,
        "canonical_normal_optimized_pairs": len(pairs), "pairs": pairs,
        "scientific_children_executed": 2 * len(pairs),
        "scope": "Exact arithmetic and source authentication. Continuous proofs and retained branch/reconstruction premises remain explicit. No numerical PDE evolution or background, bump or source-oracle query is performed by these reviewed standard-library components.",
    }
    with (out / "VERIFICATION.json").open("xb") as stream:
        stream.write(canonical(receipt))
    print(json.dumps({"status": receipt["status"],
                      "canonical_normal_optimized_pairs": len(pairs),
                      "scientific_children_executed": 2 * len(pairs),
                      "manifest_sha256": pin}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError, subprocess.SubprocessError) as error:
        print(json.dumps({"status": "FAIL", "witness": str(error)}, sort_keys=True))
        sys.exit(2)
