"""Check every frozen capsule file and ZIP against the Git index or an immutable commit.

Read-only for the repository. Expected bytes come from Git objects, so files
present only in a working directory cannot satisfy this check.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import zipfile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(name):
    path = PurePosixPath(name)
    require(bool(name) and not path.is_absolute() and '..' not in path.parts
            and '\\' not in name and path.as_posix() == name, 'Noncanonical path: ' + name)
    return name


def command(repo, *args):
    run = subprocess.run(['git', '-C', str(repo), *args], capture_output=True)
    require(run.returncode == 0, 'Git object lookup failed: ' + run.stderr.decode(errors='replace').strip())
    return run.stdout


def inventory(repo, revision, capsule):
    entries = {}
    if revision is None:
        rows = command(repo, 'ls-files', '--stage', '-z', '--', capsule)
    else:
        rows = command(repo, 'ls-tree', '-r', '-z', revision, '--', capsule)
    for raw in rows.split(b'\0'):
        if not raw:
            continue
        header, encoded = raw.split(b'\t', 1)
        name = canonical(encoded.decode())
        parts = header.decode().split()
        if revision is None:
            mode, oid, stage = parts
            require(stage == '0', 'Unresolved index entry: ' + name)
        else:
            mode, kind, oid = parts
            require(kind == 'blob', 'Regular Git blob required: ' + name)
        require(mode in ('100644', '100755'), 'Symlink or special Git entry forbidden: ' + name)
        require(name not in entries, 'Duplicate Git path: ' + name)
        entries[name] = oid
    return entries


def verify(repo, capsule, manifest_sha256, archive_sha256, revision=None):
    canonical(capsule)
    for pin in (manifest_sha256, archive_sha256):
        require(re.fullmatch('[0-9a-f]{64}', pin) is not None, 'External complete SHA256 pins required')
    if revision is not None:
        require(re.fullmatch('[0-9a-f]{40}', revision) is not None, 'Immutable full commit SHA required')
        require(command(repo, 'rev-parse', revision + '^{commit}').decode().strip() == revision,
                'Revision must name the exact commit, not a tag')
    def blob(name):
        canonical(name)
        return command(repo, 'show', (revision + ':' if revision else ':') + name)
    manifest_name = capsule + '/MANIFEST.json'
    manifest_raw = blob(manifest_name)
    require(sha(manifest_raw) == manifest_sha256, 'Git manifest differs from the external pin')
    manifest = json.loads(manifest_raw)
    require(manifest['capsule'] == PurePosixPath(capsule).name, 'Manifest capsule identity differs')
    expected = {manifest_name: manifest_raw}
    for row in manifest['files']:
        relative = canonical(row['path'])
        name = capsule + '/' + relative
        require(name not in expected, 'Duplicate or reserved manifest source path')
        require(type(row['bytes']) is int and row['bytes'] >= 0
                and re.fullmatch('[0-9a-f]{64}', row['sha256']) is not None, 'Invalid manifest record')
        expected[name] = None
    tracked = inventory(repo, revision, capsule)
    missing, extra = sorted(set(expected) - set(tracked)), sorted(set(tracked) - set(expected))
    require(not missing and not extra,
            'Git/manifest inventory differs: ' + json.dumps({'missing': missing, 'extra': extra}))
    require(len(expected) > 1, 'Empty frozen capsule')
    for row in manifest['files']:
        name = capsule + '/' + row['path']
        raw = blob(name)
        require(len(raw) == row['bytes'] and sha(raw) == row['sha256'], 'Frozen Git bytes differ: ' + name)
        expected[name] = raw
    archive_name = capsule + '.zip'
    archive_entry = inventory(repo, revision, archive_name)
    require(set(archive_entry) == {archive_name}, 'Exactly one regular archive Git entry required')
    archive_raw = blob(archive_name)
    require(sha(archive_raw) == archive_sha256, 'Git archive differs from the external pin')
    prefix = str(PurePosixPath(capsule).parent) + '/'
    zipped_expected = {name.removeprefix(prefix): raw for name, raw in expected.items()}
    with zipfile.ZipFile(io.BytesIO(archive_raw)) as archive:
        names = archive.namelist()
        require(len(names) == len(set(names)) == len(zipped_expected)
                and set(names) == set(zipped_expected), 'ZIP/Git member inventory differs')
        for name, raw in zipped_expected.items():
            require(archive.read(name) == raw, 'ZIP/Git bytes differ: ' + name)
    # A second inventory observes an intervening index update instead of silently
    # certifying a mixture. Commit mode is immutable by construction.
    require(inventory(repo, revision, capsule) == tracked, 'Git source inventory changed during verification')
    require(inventory(repo, revision, archive_name) == archive_entry, 'Git archive entry changed during verification')
    require(blob(manifest_name) == manifest_raw and blob(archive_name) == archive_raw,
            'Manifest/archive Git objects changed during verification')
    return {'status': 'PASS', 'scope': 'Complete frozen-source and ZIP membership in Git; no scientific or publication claim.',
            'observed_utc': datetime.now(timezone.utc).isoformat(), 'repository': str(repo.resolve()),
            'basis': 'immutable_commit' if revision else 'staged_index', 'commit': revision,
            'tree': command(repo, 'rev-parse', revision + '^{tree}').decode().strip() if revision else None,
            'capsule_path': capsule, 'source_files': len(manifest['files']),
            'tracked_source_and_manifest': len(expected), 'archive_members': len(expected),
            'all_git_bytes_equal_frozen_manifest_and_zip': True,
            'manifest_sha256': manifest_sha256, 'archive_sha256': archive_sha256,
            'manifest': {'bytes': len(manifest_raw), 'sha256': sha(manifest_raw)},
            'archive': {'bytes': len(archive_raw), 'sha256': sha(archive_raw)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--capsule', required=True, help='Repository-relative capsule directory')
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--archive-sha256', required=True)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--index', action='store_true')
    choice.add_argument('--revision', help='Complete immutable commit SHA')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify(args.repo, args.capsule, args.manifest_sha256, args.archive_sha256, args.revision)
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
