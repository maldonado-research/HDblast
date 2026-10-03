"""Verify frozen payload, safe ZIP inventory, and registered/repair provenance.

Run this before generating new runtime outputs. It reads files without extraction.
The manifest is an integrity record; it is not a signature or independent proof.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import stat
import zipfile


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Duplicate JSON key: " + key)
        result[key] = value
    return result


def load_json(data):
    def reject_constant(value):
        raise ValueError("Nonfinite JSON constant: " + value)
    return json.loads(data, object_pairs_hook=unique_object, parse_constant=reject_constant)


def safe_relative(value):
    if not isinstance(value, str) or not value or "\\" in value or "\0" in value:
        raise ValueError("Invalid payload path")
    path = PurePosixPath(value)
    if path.is_absolute() or any(part in (".", "..") for part in path.parts) or path.as_posix() != value:
        raise ValueError("Unsafe/noncanonical payload path: " + value)
    return value


def digest_bytes(data):
    return hashlib.sha256(data).hexdigest()


def digest_file(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def check_digest(value):
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError("Invalid SHA-256 digest")
    return value


def verify_provenance(root, payload):
    sources = load_json((root / "PROTOCOL_SOURCES.json").read_bytes())
    repairs = load_json((root / "REPAIRS.json").read_bytes())
    original = sources["sha256"]
    if not isinstance(original, dict) or not original:
        raise ValueError("Empty prospective source manifest")
    if sources["stage"] != "before physical mode evolution":
        raise ValueError("Prospective source stage changed")
    current = dict(original)
    for relative, expected in current.items():
        safe_relative(relative)
        check_digest(expected)
    history = []
    for repair in repairs["repairs"]:
        relative = safe_relative(repair["file"])
        before = check_digest(repair["original_sha256"])
        after = check_digest(repair["repaired_sha256"])
        preserved = safe_relative(repair["original_preserved_at"])
        if relative not in current or current[relative] != before:
            raise ValueError("Repair chain does not match original source: " + relative)
        if preserved not in payload or digest_file(root / preserved) != before:
            raise ValueError("Original repair source not preserved: " + preserved)
        if repair["matrix_parameters_or_gates_changed"] is not False:
            raise ValueError("Repair changed registered settings/gates; needs a new protocol")
        current[relative] = after
        history.append(dict(file=relative, original_preserved_at=preserved, original_sha256=before, current_sha256=after))
    for relative, expected in current.items():
        if relative not in payload or digest_file(root / relative) != expected:
            raise ValueError("Prospective/repaired source digest mismatch: " + relative)
    protocol = digest_file(root / "REGISTRATION.md")
    if protocol != original["REGISTRATION.md"] or protocol != repairs["protocol_sha256"]:
        raise ValueError("Registered protocol differs from original or repair record")
    followup = load_json((root / "FOLLOWUP_SOURCES.json").read_bytes())
    if followup["stage"] != "before K384 physical mode evolution" or followup["original_matrix_status"] != "FAIL":
        raise ValueError("Follow-up chronology/original outcome changed")
    critical = followup["sha256"]
    if not isinstance(critical, dict) or not critical:
        raise ValueError("Empty follow-up source manifest")
    for relative, expected in critical.items():
        safe_relative(relative)
        check_digest(expected)
        if relative not in payload or digest_file(root / relative) != expected:
            raise ValueError("Prospective follow-up source mismatch: " + relative)
    constants = {}
    for statement in ast.parse((root / "code/smooth_frw_followup.py").read_text()).body:
        if isinstance(statement, ast.Assign):
            for target in statement.targets:
                if isinstance(target, ast.Name) and target.id in ("FROZEN_CORE_SHA256", "BASELINE_HASHES"):
                    constants[target.id] = ast.literal_eval(statement.value)
    if constants["FROZEN_CORE_SHA256"] != current["code/smooth_frw_control.py"]:
        raise ValueError("Follow-up wrapper core pin differs from repaired source")
    baseline_pins = constants["BASELINE_HASHES"]
    expected_baseline_paths = {"summary.json", "A0.2_tight_K192/arrays.npz", "A0.2_quadrature_K192/arrays.npz"}
    if set(baseline_pins) != expected_baseline_paths:
        raise ValueError("Follow-up baseline inventory changed")
    for relative, expected in baseline_pins.items():
        safe_relative(relative)
        check_digest(expected)
        mapped = "followup_baseline/" + relative
        if mapped not in critical or critical[mapped] != expected or digest_file(root / mapped) != expected:
            raise ValueError("Portable baseline differs from wrapper's frozen pin: " + relative)
    baseline_manifest = load_json((root / "followup_baseline/BASELINE_MANIFEST.json").read_bytes())
    if baseline_manifest["summary_sha256"] != baseline_pins["summary.json"]:
        raise ValueError("Baseline summary differs from curation record")
    keys = {"t", "K192_rho", "K192_p", "K192_Q"}
    if set(baseline_manifest["artifacts"]) != {"A0.2_tight_K192", "A0.2_quadrature_K192"}:
        raise ValueError("Portable baseline curation inventory changed")
    for label, artifact in baseline_manifest["artifacts"].items():
        relative = label + "/arrays.npz"
        path = root / "followup_baseline" / relative
        if artifact["minimal_sha256"] != baseline_pins[relative] or artifact["minimal_bytes"] != path.stat().st_size:
            raise ValueError("Portable baseline curation digest/size differs: " + label)
        if set(artifact["keys"]) != keys or set(artifact["selected_column_exact_equality"]) != keys or not all(value is True for value in artifact["selected_column_exact_equality"].values()) or artifact["finite"] is not True:
            raise ValueError("Portable baseline curation did not preserve all consumed columns")
        check_digest(artifact["original_raw_sha256"])
        with zipfile.ZipFile(path) as arrays:
            names = arrays.namelist()
            if len(names) != len(set(names)) or set(names) != {key + ".npy" for key in keys} or arrays.testzip():
                raise ValueError("Portable baseline numeric archive inventory/CRC differs")
    baseline_summary = load_json((root / "followup_baseline/summary.json").read_bytes())
    if baseline_summary["status"] != "FAIL" or baseline_summary["provenance"]["source_sha256"] != constants["FROZEN_CORE_SHA256"] or baseline_summary["provenance"]["protocol_sha256"] != protocol:
        raise ValueError("Frozen historical baseline no longer records original failure")
    return dict(original_sources_checked=len(original), protocol_sha256=protocol, repairs=history,
                original_registration_commit=repairs["original_registration_commit"],
                followup_sources_checked=len(critical),
                followup_protocol_sha256=critical["FOLLOWUP_REGISTRATION.md"],
                historical_baseline_status="FAIL", portable_baseline_pins=baseline_pins)


def verify(root, archive_path):
    root = root.resolve()
    archive_path = archive_path.resolve()
    manifest_path = root / "MANIFEST.sha256.json"
    manifest_bytes = manifest_path.read_bytes()
    manifest = load_json(manifest_bytes)
    payload = manifest["sha256"]
    if not isinstance(payload, dict) or not payload:
        raise ValueError("Empty scientific payload manifest")
    for relative, expected in payload.items():
        safe_relative(relative)
        check_digest(expected)
    excluded = {"MANIFEST.sha256.json", root.name + ".zip", root.name + ".zip.sha256"}
    if set(payload) & excluded:
        raise ValueError("Manifest includes itself or outer archive")
    files = set()
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError("Symlink in checked-out payload: " + str(path))
        if path.is_file():
            relative = path.relative_to(root).as_posix()
            if relative not in excluded:
                files.add(relative)
    if files != set(payload):
        raise ValueError("Checked-out inventory differs from manifest: " + repr(dict(missing=sorted(set(payload)-files), extra=sorted(files-set(payload)))))
    prefix = root.name + "/"
    manifest_name = prefix + "MANIFEST.sha256.json"
    with zipfile.ZipFile(archive_path) as archive:
        infos = archive.infolist()
        names = [info.filename for info in infos]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate ZIP paths")
        for info in infos:
            safe_relative(info.filename)
            mode = stat.S_IFMT(info.external_attr >> 16)
            if not info.filename.startswith(prefix) or info.is_dir() or mode not in (0, stat.S_IFREG):
                raise ValueError("Unexpected ZIP path/type: " + info.filename)
        expected_names = {prefix + path for path in payload} | {manifest_name}
        if set(names) != expected_names:
            raise ValueError("ZIP inventory differs from manifest")
        if archive.read(manifest_name) != manifest_bytes:
            raise ValueError("ZIP and checked-out manifests differ")
        damaged = archive.testzip()
        if damaged:
            raise ValueError("ZIP CRC failure: " + damaged)
        for relative, expected in payload.items():
            if digest_bytes(archive.read(prefix + relative)) != expected:
                raise ValueError("ZIP SHA-256 mismatch: " + relative)
            if digest_file(root / relative) != expected:
                raise ValueError("Checked-out SHA-256 mismatch: " + relative)
    sidecar = root / (root.name + ".zip.sha256")
    archive_digest = digest_file(archive_path)
    if sidecar.exists() and sidecar.read_text().split()[0] != archive_digest:
        raise ValueError("Outer ZIP checksum sidecar differs")
    provenance = verify_provenance(root, payload)
    return dict(status="PASS", files_checked=len(payload), zip=str(archive_path),
                archive_sha256=archive_digest, manifest_sha256=digest_bytes(manifest_bytes),
                algorithm="sha256", provenance=provenance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("zip", type=Path)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args()
    print(json.dumps(verify(args.root, args.zip), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
