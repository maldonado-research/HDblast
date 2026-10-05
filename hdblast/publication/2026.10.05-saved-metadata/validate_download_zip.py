#!/usr/bin/env python3
"""Independent Python zipfile/CRC/MD5/SHA validation against the frozen candidate."""
import hashlib
import json
import struct
import sys
import zipfile
from pathlib import Path

BASE = Path(__file__).resolve().parent
PACKET = BASE.parent / "2026.10.03-verified-source-operator"

def main():
    filename = sys.argv[1] if len(sys.argv) > 1 else "HDBLAST_ZENODO_UPLOAD_ADDITIONS_23114217.zip"
    source = Path(filename).expanduser().resolve()
    receipt_path = Path(sys.argv[2] if len(sys.argv) > 2 else "PYTHON_ZIP_VALIDATION.json").expanduser().resolve()
    frozen = json.loads((PACKET / "FILE_MANIFEST.json").read_text())
    manifest = json.loads((BASE / "DOWNLOAD_MANIFEST.json").read_text())
    assert hashlib.sha256((PACKET / "FILE_MANIFEST.json").read_bytes()).hexdigest() == "00f0ed73c27e846603680fdc52667392935e59bbe18ac942f30e8e66248937bb"
    files = frozen["new_files"]
    assert len(files) == 13 and source.stat().st_size == manifest["zip"]["bytes"]
    digest = hashlib.sha256()
    with source.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
        stream.seek(-22, 2)
        end = struct.unpack("<4s4H2LH", stream.read(22))
    assert end[:5] == (b"PK\x05\x06", 0, 0, 13, 13) and end[-1] == 0
    assert end[5] + end[6] + 22 == source.stat().st_size
    assert digest.hexdigest() == manifest["zip"]["sha256"]
    payload_bytes = 0
    with zipfile.ZipFile(source) as archive:
        assert archive.comment == b""
        assert archive.namelist() == ["additions/" + item["filename"] for item in files]
        for item, info in zip(files, archive.infolist()):
            assert info.file_size == item["bytes"] and info.compress_size == item["bytes"]
            assert info.compress_type == zipfile.ZIP_STORED and info.flag_bits == 0x0800
            assert info.date_time == (1980, 1, 1, 0, 0, 0)
            assert info.extra == b"" and info.comment == b"" and info.create_system == 0
            assert info.extract_version == 20 and info.external_attr == 0
            with source.open("rb") as raw:
                raw.seek(info.header_offset)
                header = struct.unpack("<4s5H3L2H", raw.read(30))
                name = raw.read(header[-2])
            assert header[:6] == (b"PK\x03\x04", 20, 0x0800, 0, 0, 0x21)
            assert header[6:9] == (info.CRC, item["bytes"], item["bytes"]) and header[-1] == 0
            assert name == info.filename.encode("utf-8")
            sha, md5 = hashlib.sha256(), hashlib.md5()
            with archive.open(info) as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    sha.update(block); md5.update(block)
            assert sha.hexdigest() == item["sha256"] and md5.hexdigest() == item["md5"]
            payload_bytes += info.file_size
    receipt = {"status": "PASS_INDEPENDENT_PYTHON_ZIPFILE_VALIDATION", "members": 13,
               "payload_bytes": payload_bytes, "zip_bytes": source.stat().st_size,
               "zip_sha256": digest.hexdigest(), "all_frozen_sizes_md5_sha256_match": True,
               "zip_stored": True, "utf8_flags": True, "fixed_dos_epoch": True, "zip64": False,
               "network_requests": 0, "science_executed": False}
    with receipt_path.open("x") as target:
        target.write(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt))

if __name__ == "__main__":
    main()
