"""Independently capture a complete externally pinned data-only witness closure."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import os
import stat

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
SOURCE = BASE / 'publication-delivery-candidate/hdblast/publication/2026.10.08-retained-endpoint-flow/proof'
DEST = HERE / 'captured-witness'
EXTERNAL = '534a07bcd2d7e6c5ce1d21642724cccbb51fcdd93f47ff55ba5f526217bae398'


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def read(path):
    require(all(not p.is_symlink() for p in (path, *path.parents)), 'no witness aliases')
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        before = os.fstat(fd)
        require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and before.st_size <= 20 * 1024 * 1024, 'bounded regular witness')
        with os.fdopen(os.dup(fd), 'rb') as stream:
            raw = stream.read(20 * 1024 * 1024 + 1)
        after = os.fstat(fd)
        require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns), 'witness changed')
        require(len(raw) == before.st_size, 'witness size changed')
        return raw
    finally:
        os.close(fd)


def pin(raw):
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def walk(root):
    files, directories = set(), set()
    stack = [root]
    while stack:
        parent = stack.pop()
        with os.scandir(parent) as entries:
            for entry in entries:
                path = Path(entry.path)
                metadata = entry.stat(follow_symlinks=False)
                name = str(path.relative_to(root))
                if stat.S_ISDIR(metadata.st_mode):
                    directories.add(name)
                    stack.append(path)
                else:
                    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1, 'unsafe witness closure entry')
                    files.add(name)
    return files, directories


def main():
    manifest_raw = read(SOURCE / 'WITNESS_MANIFEST.json')
    require(pin(manifest_raw)['sha256'] == EXTERNAL, 'external witness pin')
    manifest = json.loads(manifest_raw)
    names = set(manifest['files'])
    parents = set()
    for name in names:
        path = PurePosixPath(name)
        require(path.as_posix() == name and not path.is_absolute() and all(p not in ('', '.', '..') for p in name.split('/')), 'canonical relative witness path')
        parents.update(str(p) for p in path.parents if str(p) != '.')
    require(len(names) == 31 and walk(SOURCE) == (names | {'WITNESS_MANIFEST.json'}, parents), 'exact complete witness closure')
    captured = {}
    for name, expected in manifest['files'].items():
        raw = read(SOURCE / name)
        require(pin(raw) == expected, 'witness bytes differ:' + name)
        captured[name] = raw
    roles = manifest['witnesses']
    expected_roles = {
        'controller_journal': '20386fdd5734494bbb2dd92516af3099ab4c8eb8ade8d3a0da7d616e36cf6d1f',
        'verified_public': '71804ff951b4a0aa6643b939fc133afc3352f9545162a141443c1ab26d267e0b',
        'inventory': '2c6d0a3d4c4e227c3a8aa5f6e48dd3ff8494ac4e2f1b9045c603cf35617ced3c',
        'metadata': '43feec5ca6bd9ab8c5c6a2cd925f6f42ff13418bd90b334de7ab604367dc385c',
    }
    for role, expected in expected_roles.items():
        require(pin(captured[roles[role]])['sha256'] == expected, 'terminal external role pin:' + role)
    for rid, group in manifest['preserved_records'].items():
        for role, suffix in (('record', 'record'), ('files', 'files')):
            require(captured[group[role]] == read(HERE / f'observation-003/{rid}-{suffix}.json'), 'witness differs from own independent live observation')
    observed_record = json.loads(read(HERE / 'new-public-observation-001/23244754-record.json'))
    observed_files = json.loads(read(HERE / 'new-public-observation-001/23244754-files.json'))
    public_record = json.loads(captured[roles['new_public_record']])
    public_files = json.loads(captured[roles['new_public_files']])
    for key in ('id', 'metadata', 'custom_fields', 'access', 'pids', 'is_published', 'is_draft', 'files'):
        require(public_record[key] == observed_record[key], 'new witness/live record disagreement:' + key)
    require(public_files == observed_files, 'new witness/live complete file listing disagreement')
    DEST.mkdir()
    for name, raw in captured.items():
        path = DEST / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(raw)
    with (DEST / 'WITNESS_MANIFEST.json').open('xb') as stream:
        stream.write(manifest_raw)
    require(read(SOURCE / 'WITNESS_MANIFEST.json') == manifest_raw and walk(SOURCE) == (names | {'WITNESS_MANIFEST.json'}, parents), 'source witness changed during capture')
    for name, raw in captured.items():
        require(read(SOURCE / name) == raw == read(DEST / name), 'witness changed during capture')
    receipt = {'status': 'PASS_INDEPENDENT_EXACT_TERMINAL_WITNESS_AUTHENTICATION',
               'witness_manifest': pin(manifest_raw), 'witness_files': manifest['files'],
               'complete_witness_files_excluding_manifest': 31, 'witness_roles': 11,
               'complete_parent_directory_roster_checked': True,
               'root_terminal_source_role_sha256': expected_roles,
               'all_five_preserved_groups_use_exact_own_postpublication_GET_bytes': True,
               'new_record_and_complete_file_listing_match_own_live_GET': True,
               'captured_once_source_and_copy_rechecked': True,
               'network_requests': 0, 'remote_writes': 0, 'attachment_streams': 0}
    with (HERE / 'TERMINAL_WITNESS_AUTHENTICATION.json').open('x') as stream:
        json.dump(receipt, stream, sort_keys=True, indent=2)
        stream.write('\n')
    print(json.dumps({key: receipt[key] for key in ('status', 'complete_witness_files_excluding_manifest', 'witness_roles')}))


if __name__ == '__main__':
    main()
