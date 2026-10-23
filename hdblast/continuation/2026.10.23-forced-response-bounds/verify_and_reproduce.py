#!/usr/bin/env python3
"""Authenticate concrete external pins and exact source/ZIP closure before replay.

Standard library only. No network or Git operation. Every CLI path/hash is
explicit; a workflow must authenticate this wrapper before executing it.
"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import signal
import stat
import subprocess
import sys
import zipfile

CAPSULE = "HDBLAST_FORCED_RESPONSE_BOUNDS_20261023"
EVIDENCE = "hdblast/continuation/2026.10.23-forced-response-bounds"
DEPENDENCY = "HDBLAST_GEOMETRIC_BULK_SOURCES_20261016"
DEPENDENCY_ZIP_SHA = "a300d732fddc8c6e3d309a7ab955860e0835b3d3c0cafb9ca4621767fbeebe49"
DEPENDENCY_MANIFEST_SHA = "3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2"
DEPENDENCY_ZIP_BYTES = 2251635
DEPENDENCY_FILES = 169
DEPENDENCY_SOURCE_BYTES = 15173440
CONTROLS = {"numeric": 9, "bounds": 5, "total": 14}
ROLES = {
    "actual_cone": ("compute/certify_exterior_response.py", "compute/ACTUAL_CONE_RESPONSE_{mode}.json"),
    "independent_numeric_review": ("numeric-review/check_cone_response.py", "numeric-review/REVIEW_FINAL_{mode}.json"),
    "numeric_rejections": ("numeric-review/check_rejection_controls.py", "numeric-review/CONTROLS_FINAL_{mode}.json"),
    "manuscript_display": ("manuscript/build_causal_section.py", "manuscript/MANUSCRIPT_ARITHMETIC.json"),
    "exact_bounds": ("check_bounds.py", "golden/BOUNDS.json"),
    "bounds_rejections": ("check_controls.py", "golden/CONTROLS.json"),
}
PROCESS_FIELDS = (
    "actual_source_integration_processes", "independent_numerical_review_processes",
    "exact_bound_arithmetic_processes", "manuscript_display_arithmetic_processes",
    "numerical_control_processes", "bound_control_processes",
)
PIN_STEMS = ("manifest", "archive", "runner", "verifier", "source_freeze",
             "official_fresh_reproduction_receipt")
PIN_KEYS = {"capsule", "source_files_including_manifest", "expected_top_level_children",
            "expected_targeted_controls_per_mode", "replay_timeout_seconds"}
PIN_KEYS |= {stem+suffix for stem in PIN_STEMS for suffix in ("_sha256", "_bytes")}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def unique(pairs):
    out = {}
    for key, value in pairs:
        need(key not in out, "DUPLICATE_JSON_KEY:"+key)
        out[key] = value
    return out


def reject_constant(value):
    raise ValueError("NONFINITE_JSON_CONSTANT:"+value)


def decode(raw):
    return json.loads(raw, object_pairs_hook=unique, parse_constant=reject_constant)


def canonical(name):
    need(isinstance(name, str), "PATH_NOT_STRING")
    pure = PurePosixPath(name)
    need(bool(name) and name != "." and not pure.is_absolute() and ".." not in pure.parts
         and "\\" not in name and pure.as_posix() == name
         and all(ord(c) >= 32 for c in name), "NONCANONICAL_PATH:"+name)
    return name


def concrete_hash(value):
    return isinstance(value, str) and re.fullmatch("[0-9a-f]{64}", value) is not None and value != "0"*64


def real_path(path, directory=False):
    path = Path(path)
    need(path.is_absolute() and str(path) == str(path.resolve()), "EXPLICIT_CANONICAL_ABSOLUTE_PATH")
    for parent in path.parents:
        need(parent.is_dir() and not parent.is_symlink(), "NONREAL_PATH_PARENT")
    mode = path.lstat().st_mode
    need(stat.S_ISDIR(mode) if directory else stat.S_ISREG(mode), "NONREGULAR_PATH")
    return path


def read_regular(path):
    path = real_path(path)
    raw = path.read_bytes()
    return raw, {"bytes": len(raw), "sha256": sha(raw), "mode": stat.S_IMODE(path.stat().st_mode)}


def pinned_file(path, pins, stem):
    raw, record = read_regular(path)
    need(record["bytes"] == pins[stem+"_bytes"] and record["sha256"] == pins[stem+"_sha256"],
         "EXTERNAL_FILE_PIN:"+stem)
    return raw, record


def validate_pins(pins):
    need(isinstance(pins, dict) and set(pins) == PIN_KEYS, "EXACT_PIN_SCHEMA_NO_WORKFLOW_CYCLE")
    need(pins["capsule"] == CAPSULE, "PIN_CAPSULE")
    for stem in PIN_STEMS:
        need(concrete_hash(pins[stem+"_sha256"]), "CONCRETE_PIN_REQUIRED:"+stem)
        need(type(pins[stem+"_bytes"]) is int and pins[stem+"_bytes"] > 0,
             "POSITIVE_PIN_SIZE_REQUIRED:"+stem)
    need(type(pins["source_files_including_manifest"]) is int and pins["source_files_including_manifest"] > 1,
         "CONCRETE_SOURCE_COUNT")
    need(type(pins["expected_top_level_children"]) is int and pins["expected_top_level_children"] == 12,
         "EXACT_TWELVE_CHILDREN_REQUIRED")
    need(pins["expected_targeted_controls_per_mode"] == CONTROLS
         and all(type(x) is int for x in pins["expected_targeted_controls_per_mode"].values()),
         "EXACT_CONTROL_COUNTS_REQUIRED")
    need(type(pins["replay_timeout_seconds"]) is int and 0 < pins["replay_timeout_seconds"] <= 1200,
         "AGGREGATE_REPLAY_TIMEOUT_AT_MOST_1200")
    return pins


def manifest_rows(raw, capsule, expected_sha, count, expected_bytes=None):
    need(sha(raw) == expected_sha and (expected_bytes is None or len(raw) == expected_bytes),
         "MANIFEST_EXTERNAL_PIN")
    manifest = decode(raw)
    need(manifest["capsule"] == capsule and isinstance(manifest["files"], list), "MANIFEST_IDENTITY")
    rows = {}
    for row in manifest["files"]:
        name = canonical(row["path"])
        need(name != "MANIFEST.json" and name not in rows, "DUPLICATE_OR_RESERVED_MANIFEST_MEMBER")
        need(type(row["bytes"]) is int and row["bytes"] >= 0 and concrete_hash(row["sha256"]),
             "MANIFEST_MEMBER_PIN_SHAPE")
        rows[name] = row
    need(len(rows)+1 == count, "MANIFEST_SOURCE_COUNT")
    return rows


def validate_payload(payload, rows):
    need(set(payload) == set(rows) | {"MANIFEST.json"}, "EXACT_SOURCE_FILE_CLOSURE")
    for name, row in rows.items():
        raw = payload[name]
        need(len(raw) == row["bytes"] and sha(raw) == row["sha256"], "SOURCE_MEMBER_PIN:"+name)


def source_inventory(root, pins, readonly=False):
    root = real_path(root, directory=True)
    manifest_raw, _ = read_regular(root/"MANIFEST.json")
    rows = manifest_rows(manifest_raw, CAPSULE, pins["manifest_sha256"],
                         pins["source_files_including_manifest"], pins["manifest_bytes"])
    expected_dirs = {"."}
    for name in rows:
        expected_dirs.update(p.as_posix() for p in PurePosixPath(name).parents)
    files, directories, modes = {}, {"."}, {".": stat.S_IMODE(root.stat().st_mode)}
    need(not readonly or not root.stat().st_mode & 0o222, "EXTRACTED_ROOT_NOT_READONLY")
    for path in root.rglob("*"):
        need(not path.is_symlink(), "SOURCE_SYMLINK")
        name = canonical(path.relative_to(root).as_posix())
        mode = path.lstat().st_mode
        need(stat.S_ISREG(mode) or stat.S_ISDIR(mode), "SOURCE_SPECIAL_FILE")
        need(not readonly or not mode & 0o222, "EXTRACTED_MEMBER_NOT_READONLY")
        modes[name] = stat.S_IMODE(mode)
        if stat.S_ISDIR(mode):
            directories.add(name)
        else:
            files[name] = path.read_bytes()
    need(directories == expected_dirs, "EXACT_SOURCE_DIRECTORY_CLOSURE")
    need(files.get("MANIFEST.json") == manifest_raw, "MANIFEST_CHANGED_DURING_INVENTORY")
    validate_payload(files, rows)
    need(len(files["reproduce.py"]) == pins["runner_bytes"]
         and sha(files["reproduce.py"]) == pins["runner_sha256"], "RUNNER_EXTERNAL_PIN")
    record = {"files": {name: {"bytes": len(raw), "sha256": sha(raw)} for name, raw in sorted(files.items())},
              "modes": modes, "files_including_manifest": len(files),
              "source_bytes": sum(len(raw) for raw in files.values())}
    return files, record


def zip_payload(raw, prefix, count, total_bytes, expected=None):
    payload = {}
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        infos = archive.infolist()
        need(len(infos) == count and sum(info.file_size for info in infos) == total_bytes,
             "ZIP_COMPLETE_COUNT_OR_SIZE")
        for info in infos:
            name = canonical(info.filename)
            need(name.startswith(prefix+"/"), "ZIP_CAPSULE_PREFIX")
            relative = canonical(name[len(prefix)+1:])
            need(relative not in payload, "ZIP_DUPLICATE_MEMBER")
            need(not info.is_dir() and not info.flag_bits & 1
                 and stat.S_IFMT(info.external_attr >> 16) == stat.S_IFREG, "ZIP_NONREGULAR_MEMBER")
            if expected is not None:
                need(relative in expected and info.file_size == len(expected[relative]), "ZIP_UNEXPECTED_MEMBER")
            data = archive.read(info)
            need(len(data) == info.file_size, "ZIP_MEMBER_SIZE")
            if expected is not None:
                need(data == expected[relative], "ZIP_SOURCE_BYTES_DIFFER")
            payload[relative] = data
    need(len(payload) == count and (expected is None or set(payload) == set(expected)), "ZIP_EXACT_FILE_CLOSURE")
    return payload


def dependency_payload(raw):
    need(len(raw) == DEPENDENCY_ZIP_BYTES and sha(raw) == DEPENDENCY_ZIP_SHA, "RETAINED_ZIP_EXTERNAL_PIN")
    payload = zip_payload(raw, DEPENDENCY, DEPENDENCY_FILES, DEPENDENCY_SOURCE_BYTES)
    rows = manifest_rows(payload["MANIFEST.json"], DEPENDENCY, DEPENDENCY_MANIFEST_SHA, DEPENDENCY_FILES)
    validate_payload(payload, rows)
    record = {"capsule": DEPENDENCY, "archive_bytes": len(raw), "archive_sha256": sha(raw),
              "manifest_sha256": DEPENDENCY_MANIFEST_SHA, "files_including_manifest": len(payload),
              "uncompressed_bytes": sum(len(raw) for raw in payload.values()),
              "complete_archive_and_manifest_authenticated": True}
    return payload, record


def control_report(report, role):
    if role == "numeric_rejections":
        rows = report["rejections"]
        need(report["status"] == "PASS_TARGETED_FALSE_CERTIFICATE_REJECTIONS"
             and report["baseline_passed"] is True and report["count"] == 9 and len(rows) == 9
             and len({row["case"] for row in rows}) == 9
             and all(row["rejected"] is True and bool(row["error"]) for row in rows), "NUMERIC_NINE_CONTROLS")
    elif role == "bounds_rejections":
        rows = report["controls"]
        need(report["status"] == "PASS_FIVE_TARGETED_REJECTION_CONTROLS" and len(rows) == 5
             and len({row["name"] for row in rows}) == 5
             and all(row["status"].startswith("REJECTED_") for row in rows)
             and report["upstream_scientific_producers_replayed"] is False, "BOUND_FIVE_CONTROLS")


def validate_receipt(report, source, pins, dependency_record, replay_root=None):
    need(report["status"] == "PASS_AUTHENTICATED_FORCED_RESPONSE_REPRODUCTION"
         and report["manifest_sha256"] == pins["manifest_sha256"]
         and report["files_including_manifest"] == pins["source_files_including_manifest"], "REPLAY_RECEIPT_IDENTITY")
    need(report["top_level_children"] == 12 and report["targeted_in_process_controls_per_mode"] == CONTROLS
         and all(report[key] == 2 for key in PROCESS_FIELDS), "REPLAY_RECEIPT_PROCESS_COUNTS")
    need(report["source_unchanged"] is True and report["retained_dependency_unchanged"] is True
         and report["upstream_scientific_producers_replayed"] is False
         and report["retained_dependency"] == dependency_record, "REPLAY_RECEIPT_PRESERVATION")
    children = report["children"]
    need(isinstance(children, list) and len(children) == 12, "REPLAY_TWELVE_ROWS")
    observed, outputs = set(), {}
    for row in children:
        role, optimized = row["role"], row["optimized"]
        need(role in ROLES and type(optimized) is bool and (role, optimized) not in observed, "REPLAY_ROLE_OR_MODE")
        observed.add((role, optimized))
        program, golden = ROLES[role]
        golden = golden.format(mode="OPTIMIZED" if optimized else "NORMAL")
        label = role+("_optimized" if optimized else "_normal")
        name = label+"/MANUSCRIPT_ARITHMETIC.json" if role == "manuscript_display" else label+".json"
        need(row["program"] == program and row["output"] == name
             and type(row["exit_code"]) is int and row["exit_code"] == 0
             and row["typed_and_raw_golden_equal"] is True, "REPLAY_CHILD_METADATA")
        raw = source[golden]
        need(row["bytes"] == len(raw) and row["sha256"] == sha(raw), "REPLAY_CHILD_GOLDEN_PIN")
        control_report(decode(raw), role)
        outputs[canonical(name)] = raw
        auxiliary = row["additional_raw_outputs"]
        if role == "manuscript_display":
            need(isinstance(auxiliary, list) and len(auxiliary) == 1, "MANUSCRIPT_AUXILIARY_COUNT")
            aux = auxiliary[0]
            tex_name, tex = label+"/CAUSAL_RESULT.tex", source["manuscript/CAUSAL_RESULT.tex"]
            need(aux["output"] == tex_name and aux["bytes"] == len(tex) and aux["sha256"] == sha(tex)
                 and aux["raw_golden_equal"] is True, "MANUSCRIPT_AUXILIARY_GOLDEN")
            outputs[tex_name] = tex
        else:
            need(auxiliary == [], "UNEXPECTED_CHILD_AUXILIARY")
    need(observed == {(role, mode) for role in ROLES for mode in (False, True)}, "REPLAY_COMPLETE_ROLE_MODE_CLOSURE")
    if replay_root is not None:
        for name, expected in outputs.items():
            raw, _ = read_regular(replay_root/name)
            need(raw == expected, "FRESH_CHILD_OUTPUT_BYTES:"+name)
    return outputs


def authenticate(repo, pins, pins_path, pins_sha):
    source = repo/"research"/CAPSULE
    payload, inventory = source_inventory(source, pins)
    archive_raw, archive_record = pinned_file(repo/"research"/(CAPSULE+".zip"), pins, "archive")
    payload = zip_payload(archive_raw, CAPSULE, len(payload), inventory["source_bytes"], payload)
    retained, retained_record = dependency_payload(payload["inputs/GEOMETRIC_SOURCE_DEPENDENCY.zip"])
    freeze_raw, freeze_record = pinned_file(repo/EVIDENCE/"SOURCE_FREEZE.json", pins, "source_freeze")
    freeze = decode(freeze_raw)
    need(freeze["capsule"] == CAPSULE and freeze["manifest_sha256"] == pins["manifest_sha256"]
         and freeze["source_files_including_manifest"] == len(payload)
         and freeze["source_bytes"] == inventory["source_bytes"]
         and freeze["archive"] == {"name": CAPSULE+".zip", "bytes": pins["archive_bytes"], "sha256": pins["archive_sha256"]}
         and freeze["source_files_read_only"] is True, "FREEZE_IDENTITY_OR_CLOSURE")
    official_raw, official_record = pinned_file(repo/EVIDENCE/"REPRODUCTION_RECEIPT.json", pins,
                                               "official_fresh_reproduction_receipt")
    validate_receipt(decode(official_raw), payload, pins, retained_record)
    _, verifier_record = pinned_file(Path(__file__).absolute(), pins, "verifier")
    pins_raw, pins_record = read_regular(pins_path)
    need(sha(pins_raw) == pins_sha and decode(pins_raw) == pins, "PINS_CHANGED")
    record = {"source": inventory, "archive": archive_record, "source_freeze": freeze_record,
              "official_receipt": official_record, "verifier": verifier_record, "pins": pins_record,
              "retained_dependency": retained_record}
    return payload, archive_raw, official_raw, retained, record


def output_location(path, repo):
    need(path.is_absolute() and str(path) == str(path.resolve()), "CANONICAL_ABSOLUTE_OUTPUT_REQUIRED")
    need(not path.exists() and not path.is_symlink() and not path.is_relative_to(repo), "FRESH_OUTPUT_OUTSIDE_REPOSITORY")
    for parent in path.parents:
        need(not parent.is_symlink(), "OUTPUT_SYMLINK_PARENT")
        need(not (parent/".git").exists() and not (parent/".git").is_symlink(), "OUTPUT_INSIDE_GIT_WORKTREE")
    return path


def extract_source(payload, destination):
    need(not destination.exists(), "EXTRACTION_DESTINATION_EXISTS")
    destination.mkdir(parents=True)
    for name, raw in sorted(payload.items()):
        target = destination/canonical(name)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as stream:
            stream.write(raw)
        target.chmod(0o444)
    for path in sorted((p for p in destination.rglob("*") if p.is_dir()), key=lambda p: len(p.parts), reverse=True):
        path.chmod(0o555)
    destination.chmod(0o555)
    destination.parent.chmod(0o555)


def run_replay(command, output, timeout):
    need(os.name == "posix", "FULL_REPLAY_REQUIRES_POSIX_PROCESS_GROUP_TIMEOUT")
    process = subprocess.Popen(command, cwd=output, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               start_new_session=True)
    timed_out = False
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
    (output/"reproduce.stdout").write_bytes(stdout)
    (output/"reproduce.stderr").write_bytes(stderr)
    need(not timed_out, "AGGREGATE_REPLAY_TIMEOUT")
    need(process.returncode == 0, "AUTHENTICATED_REPLAY_FAILED_SEE_EXTERNAL_LOGS")


def replay_output_closure(root, outputs, retained):
    root = real_path(root, directory=True)
    expected = dict(outputs)
    expected.update({"retained-source/"+DEPENDENCY+"/"+name: raw for name, raw in retained.items()})
    expected_names = set(expected) | {"REPRODUCTION_RECEIPT.json"}
    expected_dirs = {"."}
    for name in expected_names:
        expected_dirs.update(p.as_posix() for p in PurePosixPath(name).parents)
    actual, actual_dirs = set(), {"."}
    for path in root.rglob("*"):
        need(not path.is_symlink() and (path.is_file() or path.is_dir()), "REPLAY_OUTPUT_NONREGULAR")
        if path.is_dir():
            actual_dirs.add(canonical(path.relative_to(root).as_posix()))
        else:
            name = canonical(path.relative_to(root).as_posix())
            actual.add(name)
            if name != "REPRODUCTION_RECEIPT.json":
                need(name in expected and path.read_bytes() == expected[name], "REPLAY_OUTPUT_MEMBER_PIN")
                if name.startswith("retained-source/"):
                    need(not path.stat().st_mode & 0o222, "RETAINED_EXTRACTED_FILE_NOT_READONLY")
    need(actual == expected_names, "REPLAY_OUTPUT_EXACT_CLOSURE")
    need(actual_dirs == expected_dirs, "REPLAY_OUTPUT_EXACT_DIRECTORY_CLOSURE")
    return len(actual)


def main():
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode, "USE_PYTHON_I_S_B")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--pins", type=Path, required=True)
    parser.add_argument("--pins-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    need(concrete_hash(args.pins_sha256), "EXTERNAL_PINS_DIGEST_REQUIRED")
    repo = real_path(args.repo_root, directory=True)
    pins_path = real_path(args.pins)
    pins_raw, _ = read_regular(pins_path)
    need(sha(pins_raw) == args.pins_sha256, "PUBLICATION_PINS_EXTERNAL_DIGEST")
    pins = validate_pins(decode(pins_raw))
    output = output_location(args.output_dir, repo)
    payload, archive_raw, official_raw, retained, before = authenticate(repo, pins, pins_path, args.pins_sha256)
    output.mkdir(parents=True, exist_ok=False)
    replay_raw, extracted_before, replay_files = None, None, 0
    if not args.verify_only:
        extracted = output/"source"/CAPSULE
        extract_source(payload, extracted)
        _, extracted_before = source_inventory(extracted, pins, readonly=True)
        command = [sys.executable, "-I", "-S", "-B", str(extracted/"reproduce.py"),
                   "--manifest-sha256", pins["manifest_sha256"], "--output-dir", str(output/"replay")]
        run_replay(command, output, pins["replay_timeout_seconds"])
        replay_raw, _ = read_regular(output/"replay/REPRODUCTION_RECEIPT.json")
        need(replay_raw == official_raw, "FRESH_RECEIPT_RAW_OFFICIAL_AGREEMENT")
        outputs = validate_receipt(decode(replay_raw), payload, pins, before["retained_dependency"], output/"replay")
        replay_files = replay_output_closure(output/"replay", outputs, retained)
        _, extracted_after = source_inventory(extracted, pins, readonly=True)
        need(extracted_after == extracted_before, "EXTRACTED_SOURCE_CHANGED")
    _, _, _, _, after = authenticate(repo, pins, pins_path, args.pins_sha256)
    need(after == before, "PUBLICATION_SOURCE_ZIP_OR_PINNED_EVIDENCE_CHANGED")
    report = {"status": "PASS_PINNED_SOURCE_ZIP_AND_DEPENDENCY_AUTHENTICATION" if args.verify_only else
              "PASS_PINNED_SOURCE_ZIP_DEPENDENCY_AND_FRESH_REPLAY", "capsule": CAPSULE,
              "verify_only": args.verify_only, "publication_pins_sha256": args.pins_sha256,
              "verifier_sha256": pins["verifier_sha256"], "runner_sha256": pins["runner_sha256"],
              "manifest_sha256": pins["manifest_sha256"], "archive_sha256": pins["archive_sha256"],
              "source_freeze_sha256": pins["source_freeze_sha256"],
              "official_fresh_reproduction_receipt_sha256": pins["official_fresh_reproduction_receipt_sha256"],
              "source_files_including_manifest": len(payload), "matching_ZIP_members": len(payload),
              "retained_dependency": before["retained_dependency"],
              "scientific_child_processes": 0 if args.verify_only else 12,
              "targeted_in_process_controls_per_mode": None if args.verify_only else CONTROLS,
              "fresh_receipt_raw_and_typed_equal_official": None if args.verify_only else True,
              "fresh_reproduction_receipt_sha256": None if replay_raw is None else sha(replay_raw),
              "complete_replay_output_files": replay_files,
              "source_ZIP_evidence_PINS_and_wrapper_unchanged": True,
              "extracted_source_readonly_and_unchanged": None if args.verify_only else True,
              "retained_169_files_authenticated": True,
              "aggregate_replay_timeout_seconds": pins["replay_timeout_seconds"],
              "upstream_full_scientific_producers_replayed": False,
              "scope": "Externally pinned byte/dependency closure and the specified twelve-child replay; no full-domain waveform, full-domain residual, external peer review or additional theorem claim."}
    target = output/"VERIFICATION_RECORD.json"
    with target.open("x") as stream:
        stream.write(json.dumps(report, indent=2, sort_keys=True)+"\n")
    target.chmod(0o444)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError, zipfile.BadZipFile) as error:
        sys.exit("FORCED_RESPONSE_VERIFICATION_FAILED:"+str(error))
