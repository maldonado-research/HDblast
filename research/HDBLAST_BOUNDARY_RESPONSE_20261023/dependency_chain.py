#!/usr/bin/env python3
"""Authenticate the complete retained 49-file and nested 169-file capsules.

Only standard-library code runs here. Both archives are authenticated before
either is extracted or any retained numerical module can be imported.
"""
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import stat
import zipfile

ROOT = Path(__file__).resolve().parent
PREVIOUS = {
    'capsule': 'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023',
    'archive_sha256': 'f310241dd0fc7894fd737bbe4ba48e25f2efc45be64577a05433503015857f41',
    'archive_bytes': 3132931,
    'manifest_sha256': '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76',
    'files_including_manifest': 49, 'uncompressed_bytes': 5718317,
}
GEOMETRIC = {
    'capsule': 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
    'archive_sha256': 'a300d732fddc8c6e3d309a7ab955860e0835b3d3c0cafb9ca4621767fbeebe49',
    'archive_bytes': 2251635,
    'manifest_sha256': '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2',
    'files_including_manifest': 169, 'uncompressed_bytes': 15173440,
}


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


def canonical(name):
    pure = PurePosixPath(name)
    need(name == pure.as_posix() and not pure.is_absolute() and
         name not in ('', '.') and '..' not in pure.parts and '\\' not in name,
         'NONCANONICAL_PATH')
    return name


def manifest_rows(raw, spec):
    need(sha(raw) == spec['manifest_sha256'], 'MANIFEST_EXTERNAL_PIN')
    manifest = decode(raw)
    need(manifest['capsule'] == spec['capsule'], 'MANIFEST_CAPSULE')
    expected = {}
    for row in manifest['files']:
        name = canonical(row['path'])
        need(name not in expected and name != 'MANIFEST.json', 'MANIFEST_DUPLICATE_MEMBER')
        need(type(row['bytes']) is int and row['bytes'] >= 0, 'MANIFEST_SIZE')
        expected[name] = row
    need(len(expected) + 1 == spec['files_including_manifest'], 'MANIFEST_COUNT')
    return expected


def authenticate_bytes(raw, spec):
    need(len(raw) == spec['archive_bytes'] and sha(raw) == spec['archive_sha256'],
         'ARCHIVE_EXTERNAL_PIN:' + spec['capsule'])
    payload = {}
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        rows = archive.infolist()
        need(len(rows) == spec['files_including_manifest'], 'ARCHIVE_MEMBER_COUNT')
        need(sum(row.file_size for row in rows) == spec['uncompressed_bytes'], 'ARCHIVE_TOTAL_SIZE')
        for row in rows:
            name = canonical(row.filename)
            need(name.startswith(spec['capsule'] + '/'), 'ARCHIVE_PREFIX')
            relative = canonical(name[len(spec['capsule']) + 1:])
            need(relative not in payload, 'ARCHIVE_DUPLICATE_MEMBER')
            need(not row.is_dir() and not row.flag_bits & 1 and
                 stat.S_IFMT(row.external_attr >> 16) == stat.S_IFREG, 'ARCHIVE_NONREGULAR_MEMBER')
            data = archive.read(row)
            need(len(data) == row.file_size, 'ARCHIVE_MEMBER_SIZE')
            payload[relative] = data
    need('MANIFEST.json' in payload, 'ARCHIVE_MANIFEST_MISSING')
    expected = manifest_rows(payload['MANIFEST.json'], spec)
    need(set(payload) == set(expected) | {'MANIFEST.json'}, 'ARCHIVE_MANIFEST_CLOSURE')
    for name, row in expected.items():
        data = payload[name]
        need(len(data) == row['bytes'] and sha(data) == row['sha256'], 'ARCHIVE_MEMBER_PIN:' + name)
    return payload


def inventory(root, spec, readonly=True):
    root = Path(root)
    need(root.is_dir() and not root.is_symlink(), 'EXTRACTED_ROOT_NOT_DIRECTORY')
    path = root / 'MANIFEST.json'
    need(path.is_file() and not path.is_symlink(), 'EXTRACTED_MANIFEST_NOT_REGULAR')
    expected = manifest_rows(path.read_bytes(), spec)
    actual = set()
    modes = {}
    total = 0
    for path in [root] + sorted(root.rglob('*')):
        need(not path.is_symlink() and (path.is_file() or path.is_dir()), 'EXTRACTED_NONREGULAR_MEMBER')
        relative = '.' if path == root else path.relative_to(root).as_posix()
        mode = stat.S_IMODE(path.stat().st_mode)
        modes[relative] = mode
        if readonly:
            need(not mode & 0o222, 'EXTRACTED_MEMBER_WRITABLE:' + relative)
        if path.is_file():
            actual.add(relative)
            data = path.read_bytes()
            total += len(data)
            if relative != 'MANIFEST.json':
                need(relative in expected, 'EXTRACTED_UNEXPECTED_MEMBER:' + relative)
                row = expected[relative]
                need(len(data) == row['bytes'] and sha(data) == row['sha256'], 'EXTRACTED_MEMBER_PIN:' + relative)
    need(actual == set(expected) | {'MANIFEST.json'} and len(actual) == spec['files_including_manifest'],
         'EXTRACTED_SOURCE_CLOSURE')
    need(total == spec['uncompressed_bytes'], 'EXTRACTED_TOTAL_SIZE')
    return {'capsule': spec['capsule'], 'manifest_sha256': spec['manifest_sha256'],
            'files_including_manifest': len(actual), 'uncompressed_bytes': total,
            'inventory_and_modes_sha256': sha(json.dumps(modes, sort_keys=True).encode()),
            'complete_source_authenticated': True, 'source_read_only': readonly}


def extract_payload(payload, destination, spec):
    destination.mkdir()
    for name, raw in sorted(payload.items()):
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as handle:
            handle.write(raw)
        target.chmod(0o444)
    for path in sorted(destination.rglob('*'), reverse=True):
        if path.is_dir():
            path.chmod(0o555)
    destination.chmod(0o555)
    inventory(destination, spec)


def extract_chain(archive, destination):
    archive = Path(archive)
    need(archive.is_file() and not archive.is_symlink(), 'PREVIOUS_ARCHIVE_NOT_REGULAR')
    previous = authenticate_bytes(archive.read_bytes(), PREVIOUS)
    geometric = authenticate_bytes(previous['inputs/GEOMETRIC_SOURCE_DEPENDENCY.zip'], GEOMETRIC)
    destination = Path(destination)
    need(destination.is_absolute() and not destination.exists() and not destination.resolve().is_relative_to(ROOT),
         'FRESH_DEPENDENCY_DESTINATION_REQUIRED')
    destination.mkdir(parents=True)
    previous_root = destination / PREVIOUS['capsule']
    geometric_root = destination / GEOMETRIC['capsule']
    extract_payload(previous, previous_root, PREVIOUS)
    extract_payload(geometric, geometric_root, GEOMETRIC)
    destination.chmod(0o555)
    record = {'previous': dict(PREVIOUS), 'geometric': dict(GEOMETRIC),
              'both_complete_archives_authenticated_before_extraction': True,
              'upstream_scientific_producers_replayed': False}
    return previous_root, geometric_root, record
