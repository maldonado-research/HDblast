#!/usr/bin/env python3
"""Build or verify a deterministic, complete, standard-library capsule ZIP.

The capsule is authenticated before any archive is opened for writing.  No
project module is imported.  ZIP_STORED, a fixed timestamp, a sorted roster
and normalized regular-file permissions make the archive bytes repeatable.
The external manifest SHA-256 is mandatory in both build and verify modes.
"""
import sys
sys.dont_write_bytecode = True

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import zipfile

CAPSULE_NAME = 'HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026'
MANIFEST_SCHEMA = 'HDBLAST_SCREEN_MAGNETIC_CURVATURE_MANIFEST_V1'


def need(condition, reason):
    if not condition:
        raise ValueError(reason)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def valid_sha(value):
    return isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) is not None


def no_symlinks(path, reason):
    absolute = Path(os.path.abspath(path))
    need(all(not part.is_symlink() for part in (absolute,) + tuple(absolute.parents)), reason)
    return absolute


def member_path(value):
    need(isinstance(value, str) and value and '\\' not in value and '\x00' not in value,
         'MANIFEST_MEMBER_PATH')
    path = PurePosixPath(value)
    need(not path.is_absolute() and path.as_posix() == value and ':' not in value and value not in ('.', '..')
         and all(part not in ('.', '..') for part in path.parts), 'MANIFEST_MEMBER_PATH')


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        need(key not in result, 'JSON_DUPLICATE_KEY')
        result[key] = value
    return result


def authenticate(root, pin):
    need(valid_sha(pin), 'MANIFEST_PIN_FORMAT')
    root = no_symlinks(root, 'ROOT_SYMLINK_OR_TYPE')
    need(root.is_dir() and root.name == CAPSULE_NAME, 'ROOT_SYMLINK_OR_TYPE')
    manifest = root / 'MANIFEST.json'
    need(manifest.is_file() and not manifest.is_symlink() and sha(manifest) == pin, 'MANIFEST_PIN')
    data = json.loads(manifest.read_text(encoding='utf-8'), object_pairs_hook=unique_object)
    need(isinstance(data, dict) and data.get('schema') == MANIFEST_SCHEMA
         and isinstance(data.get('files'), list), 'MANIFEST_SCHEMA')
    expected = {'MANIFEST.json'}
    for row in data['files']:
        need(isinstance(row, dict), 'MANIFEST_ENTRY')
        relative = row.get('path')
        member_path(relative)
        need(relative not in expected, 'MANIFEST_MEMBER_PATH')
        expected.add(relative)
        need(type(row.get('bytes')) is int and row['bytes'] >= 0 and valid_sha(row.get('sha256')),
             'MANIFEST_ENTRY_FIELDS')
        path = root / relative
        need(path.is_file() and not path.is_symlink(), 'MANIFEST_MEMBER_TYPE')
        need(path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], 'MANIFEST_MEMBER_BYTES')
    actual = set()
    for path in root.rglob('*'):
        need(not path.is_symlink(), 'TREE_SYMLINK')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
        else:
            need(path.is_dir(), 'TREE_SPECIAL_MEMBER')
    need(actual == expected, 'MANIFEST_COMPLETE_CLOSURE')
    need(sha(manifest) == pin, 'MANIFEST_CHANGED_DURING_AUTHENTICATION')
    return sorted(expected)


def verify(root, archive, pin, roster):
    need(valid_sha(pin), 'ARCHIVE_PIN_FORMAT')
    archive = no_symlinks(archive, 'ARCHIVE_SYMLINK')
    need(archive.is_file() and sha(archive) == pin, 'ARCHIVE_PIN')
    expected = {CAPSULE_NAME + '/' + relative: root / relative for relative in roster}
    with zipfile.ZipFile(archive) as packed:
        names = packed.namelist()
        need(len(names) == len(set(names)) and set(names) == set(expected), 'ARCHIVE_EXACT_MEMBERS')
        for entry in packed.infolist():
            member_path(entry.filename)
            mode = entry.external_attr >> 16
            need(not entry.is_dir() and not stat.S_ISLNK(mode)
                 and stat.S_IFMT(mode) in (0, stat.S_IFREG), 'ARCHIVE_MEMBER_TYPE')
            need(not entry.flag_bits & 1, 'ARCHIVE_ENCRYPTED_MEMBER')
            source = expected[entry.filename]
            need(entry.file_size == source.stat().st_size
                 and packed.read(entry) == source.read_bytes(), 'ARCHIVE_MEMBER_BYTES')
        need(packed.testzip() is None, 'ARCHIVE_CRC')
    need(sha(archive) == pin, 'ARCHIVE_CHANGED_DURING_AUTHENTICATION')
    return {'sha256': pin, 'bytes': archive.stat().st_size, 'members': len(expected),
            'crc_pass': True, 'exact_member_set': True, 'all_member_bytes_equal': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--capsule-root', type=Path, required=True)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--verify', action='store_true', help='Verify an existing ZIP; create no files.')
    parser.add_argument('--archive-sha256', help='External ZIP pin, required with --verify.')
    args = parser.parse_args()
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode,
         'ENTRY_REQUIRES_ISOLATED_NO_SITE_NO_BYTECODE_FLAGS')
    entry = no_symlinks(__file__, 'ENTRY_SOURCE_SYMLINK')
    entry_pin = sha(entry)
    root = no_symlinks(args.capsule_root, 'ROOT_SYMLINK_OR_TYPE')
    roster = authenticate(root, args.manifest_sha256)
    archive = no_symlinks(args.archive, 'ARCHIVE_SYMLINK')
    if args.verify:
        need(args.archive_sha256 is not None, 'EXTERNAL_ARCHIVE_PIN_REQUIRED')
        result = verify(root, archive, args.archive_sha256, roster)
    else:
        need(args.archive_sha256 is None, 'BUILD_DOES_NOT_ACCEPT_ARCHIVE_PIN')
        need(not archive.exists(), 'FRESH_EXTERNAL_ARCHIVE_REQUIRED')
        need(archive.parent.is_dir(), 'ARCHIVE_PARENT_REQUIRED')
        need(archive != root and root not in archive.parents, 'ARCHIVE_INSIDE_CAPSULE')
        with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_STORED, allowZip64=True) as packed:
            for relative in roster:
                info = zipfile.ZipInfo(CAPSULE_NAME + '/' + relative, date_time=(1980, 1, 1, 0, 0, 0))
                info.create_system = 3
                info.compress_type = zipfile.ZIP_STORED
                info.external_attr = (stat.S_IFREG | 0o644) << 16
                info.extra = b''
                info.comment = b''
                packed.writestr(info, (root / relative).read_bytes())
        result = verify(root, archive, sha(archive), roster)
    need(authenticate(root, args.manifest_sha256) == roster
         and sha(entry) == entry_pin, 'SOURCE_CHANGED_DURING_PACKAGING')
    print(json.dumps({'status': 'PASS_DETERMINISTIC_CAPSULE_ZIP', 'mode': 'verify' if args.verify else 'build',
        'manifest_sha256': args.manifest_sha256, 'packager_sha256': entry_pin,
        'project_modules_imported': False, 'source_unchanged': True, 'archive': result},
        sort_keys=True, separators=(',', ':')), flush=True)


if __name__ == '__main__':
    main()
