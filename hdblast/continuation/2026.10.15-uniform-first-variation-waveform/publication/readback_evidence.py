"""Authenticate an additive evidence-only commit, then read its merged bytes.

Normal operation observes an already merged PR and reads HTTP bytes. It never
fetches or mutates Git. The result is local outside the evidence worktree to
avoid recursively including a receipt for its own publication.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import urllib.parse
import urllib.request
import uuid

REPOSITORY = 'maldonado-research/HDblast'
DELIVERY = 'hdblast/continuation/2026.10.15-uniform-first-variation-waveform'


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(name):
    require(isinstance(name, str), 'String path required')
    relative = PurePosixPath(name)
    require(name and not relative.is_absolute() and '..' not in relative.parts
            and '\\' not in name and relative.as_posix() == name, 'Canonical relative path required')
    return name


def git(root, *args):
    return subprocess.check_output(['git', '-C', str(root), *args])


def textgit(root, *args):
    return git(root, *args).decode().strip()


def check_commit(root, prepared):
    for key in ('base', 'head_sha', 'reviewed_tree'):
        require(re.fullmatch('[0-9a-f]{40}', prepared[key]) is not None, 'Full immutable evidence pins required')
    require(prepared['delivery'] == DELIVERY and type(prepared['added_files']) is int
            and prepared['added_files'] > 0 and type(prepared['total_evidence_files']) is int
            and re.fullmatch('[0-9a-f]{64}', prepared['manifest_sha256']) is not None,
            'Reviewed evidence metadata differs')
    head = prepared['head_sha']
    require(textgit(root, 'rev-parse', 'HEAD') == head and not textgit(root, 'status', '--porcelain'),
            'Clean evidence worktree must remain at reviewed head')
    require(textgit(root, 'rev-parse', head+'^{tree}') == prepared['reviewed_tree'], 'Evidence reviewed tree differs')
    changes = git(root, 'diff', '--no-renames', '--name-status', '-z', prepared['base'], head).split(b'\0')
    changes = [c for c in changes if c]
    require(len(changes) == 2*prepared['added_files'], 'Evidence addition count differs')
    additions = []
    for status, name in zip(changes[::2], changes[1::2]):
        name = canonical(name.decode())
        require(status == b'A' and name.startswith(DELIVERY+'/'), 'Evidence commit changes existing research/site bytes')
        additions.append(name)
    require(len(set(additions)) == len(additions), 'Duplicate evidence additions')
    def blob(name):
        raw = git(root, 'show', head+':'+name)
        path = root/name
        require(path.is_file() and not any(p.is_symlink() for p in [path, *path.parents])
                and path.read_bytes() == raw, 'Evidence worktree/committed bytes differ: '+name)
        return raw
    manifest_name = DELIVERY+'/PUBLICATION_EVIDENCE_MANIFEST.json'
    manifest_raw = blob(manifest_name)
    require(sha(manifest_raw) == prepared['manifest_sha256'], 'External evidence manifest pin differs')
    manifest = json.loads(manifest_raw)
    expected = {}
    for row in manifest['files']:
        name = canonical(row['path'])
        require(name not in expected and name != 'PUBLICATION_EVIDENCE_MANIFEST.json'
                and type(row['bytes']) is int and row['bytes'] >= 0
                and re.fullmatch('[0-9a-f]{64}', row['sha256']) is not None, 'Invalid/duplicate final evidence record')
        raw = blob(DELIVERY+'/'+name)
        require(len(raw) == row['bytes'] and sha(raw) == row['sha256'], 'Manifest-pinned evidence bytes differ')
        expected[name] = {'path': name, 'bytes': len(raw), 'sha256': sha(raw), 'raw': raw}
    require(bool(expected), 'Empty evidence manifest')
    expected['PUBLICATION_EVIDENCE_MANIFEST.json'] = {'path': 'PUBLICATION_EVIDENCE_MANIFEST.json',
        'bytes': len(manifest_raw), 'sha256': sha(manifest_raw), 'raw': manifest_raw}
    require(len(expected) == prepared['total_evidence_files'], 'Evidence file count differs')
    tracked = set()
    for line in git(root, 'ls-tree', '-r', '-z', head, '--', DELIVERY).split(b'\0'):
        if not line:
            continue
        header, name = line.split(b'\t', 1)
        mode, kind, oid = header.decode().split()
        name = canonical(name.decode())
        require(kind == 'blob' and mode in ('100644', '100755'), 'Regular committed evidence required')
        tracked.add(name)
    require(tracked == {DELIVERY+'/'+n for n in expected}, 'Complete committed evidence inventory differs')
    return list(expected.values()), additions


def main():
    task = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pr', type=int, help='Actual evidence-only PR number; required outside preflight')
    parser.add_argument('--prepared', type=Path, default=task/'PREPARED_EVIDENCE_COMMIT.json')
    parser.add_argument('--worktree', type=Path, default=task/'publish-evidence')
    parser.add_argument('--output', type=Path, default=task/'FINAL_EVIDENCE_READBACK.json')
    parser.add_argument('--workers', type=int, default=6, choices=range(1, 9))
    parser.add_argument('--preflight-only', action='store_true', help='Authenticate local evidence commit, without PR/API/HTTP observation')
    args = parser.parse_args()
    root = args.worktree.resolve()
    require(not args.output.resolve().is_relative_to(root), 'Final evidence readback receipt must stay outside published worktree')
    prepared_raw = args.prepared.read_bytes()
    prepared = json.loads(prepared_raw)
    expected, additions = check_commit(root, prepared)
    helper_pin = sha(Path(__file__).read_bytes())
    if args.preflight_only:
        receipt = {'status': 'PREPARED_NOT_EXECUTED', 'prepared_evidence_sha256': sha(prepared_raw),
                   'helper_sha256': helper_pin, 'target_count': len(expected), 'new_files': len(additions),
                   'manifest_sha256': prepared['manifest_sha256'],
                   'scope': 'Local immutable evidence commit and addition-only scope authentication; no PR or HTTP observation.'}
        args.output.write_text(json.dumps(receipt, indent=2)+'\n')
        print(json.dumps(receipt, indent=2))
        return
    require(type(args.pr) is int and args.pr > 0, 'Actual positive evidence PR number required')
    raw_pr = subprocess.check_output(['gh', 'api', f'repos/{REPOSITORY}/pulls/{args.pr}'], timeout=60)
    pr = json.loads(raw_pr)
    require(pr['number'] == args.pr and pr['merged'] is True and pr['head']['sha'] == prepared['head_sha']
            and pr['html_url'] == f'https://github.com/{REPOSITORY}/pull/{args.pr}', 'Expected reviewed evidence PR must be merged')
    merge = pr['merge_commit_sha']
    require(re.fullmatch('[0-9a-f]{40}', merge) is not None, 'Full evidence merge commit required')
    require(subprocess.run(['git', '-C', str(root), 'cat-file', '-e', merge+'^{commit}'],
                           capture_output=True).returncode == 0, 'Evidence merge object missing locally; obtain it before readback')
    require(textgit(root, 'rev-parse', merge+'^{tree}') == prepared['reviewed_tree'], 'Merged evidence tree differs')
    require(git(root, 'diff', '--no-renames', '--name-status', '-z', prepared['base'], merge)
            == git(root, 'diff', '--no-renames', '--name-status', '-z', prepared['base'], prepared['head_sha']),
            'Evidence merge changed the reviewed additive scope')
    nonce = uuid.uuid4().hex
    def verify(item):
        url = ('https://raw.githubusercontent.com/'+REPOSITORY+'/main/'+DELIVERY+'/'
               +urllib.parse.quote(item['path'], safe='/'))
        result = {'path': item['path'], 'url': url, 'expected_bytes': item['bytes'], 'expected_sha256': item['sha256']}
        try:
            request = urllib.request.Request(url+'?hdblast_evidence_readback='+nonce, headers={
                'Cache-Control': 'no-cache', 'User-Agent': 'HDBLAST-evidence-readback'})
            with urllib.request.urlopen(request, timeout=45) as response:
                raw = response.read()
                result.update(http_status=response.status, bytes=len(raw), sha256=sha(raw),
                              byte_equal=raw == item['raw'], final_url=response.geturl())
        except Exception as error:
            result.update(http_status=None, byte_equal=False, error=str(error))
        return result
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        targets = list(pool.map(verify, expected))
    success = all(t['http_status'] == 200 and t['byte_equal'] for t in targets)
    receipt = {'schema': 'hdblast-uniform-waveform-final-evidence-readback-v1',
               'observed_utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS' if success else 'FAIL',
               'target_count': len(targets), 'pull_request': {'repository': REPOSITORY, 'number': args.pr,
                   'url': pr['html_url'], 'merged': True, 'head_sha': pr['head']['sha'],
                   'merge_commit_sha': merge, 'merged_at': pr['merged_at']},
               'merged_tree_equals_reviewed_tree': True, 'all_existing_research_and_site_files_unchanged': True,
               'added_files': len(additions), 'manifest_sha256': prepared['manifest_sha256'],
               'prepared_evidence_sha256': sha(prepared_raw), 'observed_pr_response_sha256': sha(raw_pr),
               'helper_sha256': helper_pin, 'targets': targets,
               'scope': 'Final evidence publication byte readback retained locally, outside the published evidence tree, to avoid recursive receipts.'}
    args.output.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': receipt['status'], 'targets': len(targets), 'merge_commit': merge,
                      'failed_paths': [t['path'] for t in targets if not t['byte_equal']]}))
    raise SystemExit(0 if success else 1)


if __name__ == '__main__':
    main()
