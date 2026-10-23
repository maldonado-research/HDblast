#!/usr/bin/env python3
"""Authenticate and safely extract the exact published geometric source capsule."""
import argparse
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import stat
import sys
import zipfile

CAPSULE = "HDBLAST_GEOMETRIC_BULK_SOURCES_20261016"
ZIP_SHA256 = "a300d732fddc8c6e3d309a7ab955860e0835b3d3c0cafb9ca4621767fbeebe49"
MANIFEST_SHA256 = "3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2"
ZIP_BYTES = 2251635
FILES_INCLUDING_MANIFEST = 169
TOTAL_UNCOMPRESSED_BYTES = 15173440
ROOT = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def unique_object(pairs):
    out = {}
    for key, value in pairs:
        need(key not in out, "DEPENDENCY_DUPLICATE_JSON_KEY:"+key)
        out[key] = value
    return out


def canonical(name):
    pure = PurePosixPath(name)
    need(name == pure.as_posix() and not pure.is_absolute() and
         name not in (".", "") and ".." not in pure.parts and "\\" not in name,
         "DEPENDENCY_NONCANONICAL_PATH")
    return name


def authenticate_archive(path):
    path = Path(path)
    need(path.is_file() and not path.is_symlink(), "DEPENDENCY_ZIP_NOT_REGULAR")
    raw = path.read_bytes()
    need(len(raw) == ZIP_BYTES and sha(raw) == ZIP_SHA256, "DEPENDENCY_ZIP_EXTERNAL_PIN")
    payload = {}
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        infos = archive.infolist()
        need(len(infos) == FILES_INCLUDING_MANIFEST, "DEPENDENCY_ZIP_MEMBER_COUNT")
        need(sum(row.file_size for row in infos) == TOTAL_UNCOMPRESSED_BYTES,
             "DEPENDENCY_ZIP_TOTAL_SIZE")
        for row in infos:
            name = canonical(row.filename)
            need(name.startswith(CAPSULE+"/"), "DEPENDENCY_ZIP_PREFIX")
            relative = canonical(name[len(CAPSULE)+1:])
            need(relative not in payload, "DEPENDENCY_ZIP_DUPLICATE_MEMBER")
            need(not row.is_dir() and not row.flag_bits & 1 and
                 stat.S_IFMT(row.external_attr >> 16) == stat.S_IFREG,
                 "DEPENDENCY_ZIP_NONREGULAR_MEMBER")
            data = archive.read(row)
            need(len(data) == row.file_size, "DEPENDENCY_ZIP_MEMBER_SIZE")
            payload[relative] = data
    need("MANIFEST.json" in payload and sha(payload["MANIFEST.json"]) == MANIFEST_SHA256,
         "DEPENDENCY_MANIFEST_EXTERNAL_PIN")
    manifest = json.loads(payload["MANIFEST.json"], object_pairs_hook=unique_object)
    need(manifest["capsule"] == CAPSULE, "DEPENDENCY_MANIFEST_CAPSULE")
    expected = {}
    for row in manifest["files"]:
        name = canonical(row["path"])
        need(name not in expected and name != "MANIFEST.json", "DEPENDENCY_MANIFEST_DUPLICATE_MEMBER")
        expected[name] = row
    need(set(payload) == set(expected) | {"MANIFEST.json"}, "DEPENDENCY_ZIP_MANIFEST_CLOSURE")
    for name, row in expected.items():
        data = payload[name]
        need(len(data) == row["bytes"] and sha(data) == row["sha256"],
             "DEPENDENCY_MEMBER_PIN:"+name)
    record = {"capsule": CAPSULE, "archive_bytes": len(raw), "archive_sha256": sha(raw),
              "manifest_sha256": MANIFEST_SHA256, "files_including_manifest": len(payload),
              "uncompressed_bytes": sum(len(data) for data in payload.values()),
              "complete_archive_and_manifest_authenticated": True}
    return payload, record


def extracted_inventory(root):
    root = Path(root).resolve()
    manifest_path = root/"MANIFEST.json"
    need(manifest_path.is_file() and not manifest_path.is_symlink(), "EXTRACTED_MANIFEST_NOT_REGULAR")
    raw = manifest_path.read_bytes()
    need(sha(raw) == MANIFEST_SHA256, "EXTRACTED_MANIFEST_EXTERNAL_PIN")
    manifest = json.loads(raw, object_pairs_hook=unique_object)
    need(manifest["capsule"] == CAPSULE, "EXTRACTED_CAPSULE_NAME")
    expected = {canonical(row["path"]): row for row in manifest["files"]}
    need(len(expected) == len(manifest["files"]) and "MANIFEST.json" not in expected,
         "EXTRACTED_MANIFEST_DUPLICATE_MEMBER")
    actual = set()
    for path in root.rglob("*"):
        need(not path.is_symlink() and (path.is_file() or path.is_dir()), "EXTRACTED_NONREGULAR_MEMBER")
        if path.is_file():
            name = path.relative_to(root).as_posix()
            actual.add(name)
            if name != "MANIFEST.json":
                need(name in expected, "EXTRACTED_UNEXPECTED_MEMBER:"+name)
                data = path.read_bytes();row = expected[name]
                need(len(data) == row["bytes"] and sha(data) == row["sha256"], "EXTRACTED_MEMBER_PIN:"+name)
    need(actual == set(expected) | {"MANIFEST.json"} and len(actual) == FILES_INCLUDING_MANIFEST,
         "EXTRACTED_SOURCE_CLOSURE")
    return {"capsule": CAPSULE, "manifest_sha256": MANIFEST_SHA256,
            "files_including_manifest": len(actual), "complete_source_authenticated": True}


def extract_authenticated(archive_path, destination):
    payload, record = authenticate_archive(archive_path)
    destination = Path(destination)
    need(destination.is_absolute(), "DEPENDENCY_EXTRACTION_NEEDS_ABSOLUTE_DESTINATION")
    destination = destination.resolve()
    need(not destination.exists() and not destination.is_relative_to(ROOT),
         "DEPENDENCY_EXTRACTION_DESTINATION")
    destination.mkdir(parents=True)
    root = destination/CAPSULE
    root.mkdir()
    for name, data in sorted(payload.items()):
        output = root/name
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("xb") as handle:
            handle.write(data)
        output.chmod(0o444)
    inventory = extracted_inventory(root)
    need(inventory["files_including_manifest"] == record["files_including_manifest"], "EXTRACTION_COUNT")
    return root, record


def main():
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode,
         "USE_PYTHON_I_S_B")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ROOT/"inputs/GEOMETRIC_SOURCE_DEPENDENCY.zip")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root, record = extract_authenticated(args.archive, args.output_dir)
    print(json.dumps({"status": "PASS_FULL_RETAINED_DEPENDENCY_AUTHENTICATION_AND_EXTRACTION",
                      "source_root_name": root.name, "dependency": record}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError, zipfile.BadZipFile) as error:
        sys.exit("RETAINED_DEPENDENCY_FAILED:"+str(error))
