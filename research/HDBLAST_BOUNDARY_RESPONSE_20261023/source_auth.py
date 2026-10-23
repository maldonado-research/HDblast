#!/usr/bin/env python3
"""Strict source closure checked before importing scientific capsule code."""
import hashlib
import json
from pathlib import Path, PurePosixPath
import stat

ROOT = Path(__file__).resolve().parent
CAPSULE = 'HDBLAST_BOUNDARY_RESPONSE_20261023'


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def unique_object(pairs):
    result = {}
    for name, value in pairs:
        need(name not in result, 'DUPLICATE_JSON_KEY:' + name)
        result[name] = value
    return result


def decode(raw):
    return json.loads(raw, object_pairs_hook=unique_object)


def authenticate(expected_sha256, readonly=True):
    manifest_path = ROOT / 'MANIFEST.json'
    need(manifest_path.is_file() and not manifest_path.is_symlink(), 'SOURCE_MANIFEST_NOT_REGULAR')
    raw = manifest_path.read_bytes()
    need(sha(raw) == expected_sha256, 'SOURCE_MANIFEST_EXTERNAL_PIN')
    manifest = decode(raw)
    need(manifest['capsule'] == CAPSULE, 'SOURCE_CAPSULE_NAME')
    expected = {}
    for row in manifest['files']:
        name = row['path']; pure = PurePosixPath(name)
        need(name == pure.as_posix() and not pure.is_absolute() and name not in ('', '.', 'MANIFEST.json')
             and '..' not in pure.parts and '\\' not in name, 'NONCANONICAL_SOURCE_PATH')
        need(name not in expected and type(row['bytes']) is int and row['bytes'] >= 0,
             'DUPLICATE_OR_INVALID_SOURCE_ROW')
        expected[name] = row
    actual = set(); modes = {}; total = 0
    for path in [ROOT] + sorted(ROOT.rglob('*')):
        need(not path.is_symlink() and (path.is_file() or path.is_dir()), 'SOURCE_NONREGULAR_MEMBER')
        name = '.' if path == ROOT else path.relative_to(ROOT).as_posix()
        modes[name] = stat.S_IMODE(path.stat().st_mode)
        if readonly:
            need(not modes[name] & 0o222, 'SOURCE_WRITABLE_MEMBER:' + name)
        if path.is_file():
            actual.add(name); data = path.read_bytes(); total += len(data)
            if name != 'MANIFEST.json':
                need(name in expected, 'SOURCE_UNEXPECTED_MEMBER:' + name)
                row = expected[name]
                need(len(data) == row['bytes'] and sha(data) == row['sha256'], 'SOURCE_MEMBER_PIN:' + name)
    need(actual == set(expected) | {'MANIFEST.json'}, 'SOURCE_CLOSURE')
    return {'capsule': CAPSULE, 'manifest_sha256': expected_sha256,
            'files_including_manifest': len(actual), 'source_bytes': total,
            'inventory_and_modes_sha256': sha(json.dumps(modes, sort_keys=True).encode()),
            'source_read_only': readonly, 'complete_source_authenticated': True}
