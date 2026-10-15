"""Read back merged HDBLAST source, summaries and live pages after deployment."""
import argparse
import concurrent.futures
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import urllib.parse
import urllib.request
import uuid

CAPSULE = 'research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015'
DELIVERY = 'hdblast/continuation/2026.10.15-quasilinear-shell-evolution'
REVIEWED = {
    'maldonado-research/HDblast': '3292f5dd6eaafdf513d14486cded1dcbdaab5726',
    'maldonado-research/maldonado-research.github.io': '9225d99ff9231ded1538dd0dbc56fe1be5601f37',
}
MANIFEST_SHA256 = '1758d6479eaa4ed3637ef8d0aaf9f61fcbf3d8b0518e5f95fe5ecdd6fd8b979a'


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def git(root, *args):
    return subprocess.check_output(['git', *args], cwd=root).decode().strip()


def validate_merge_evidence(evidence, roots):
    """Require observed merges and successful deployment at both exact merges."""
    prs = {row['repository']: row for row in evidence['research_and_site_pull_requests']}
    require(set(REVIEWED) <= set(prs), 'Both merged repository observations required')
    snapshots = []
    for repository, reviewed in REVIEWED.items():
        row = prs[repository]
        merged = row['merge_commit_sha']
        require(row['merged'] is True and row['head_sha'] == reviewed,
                'Expected reviewed head must be observed merged: ' + repository)
        require(re.fullmatch('[0-9a-f]{40}', merged) is not None, 'Full merge SHA required')
        deployments = [run for run in evidence['deployments'] if run['head_sha'] == merged]
        require(any(run['status'] == 'completed' and run['conclusion'] == 'success'
                    for run in deployments), 'Successful exact-merge deployment required: ' + repository)
        root = roots[repository]
        require(git(root, 'rev-parse', 'HEAD') == reviewed,
                'Readback worktree must remain at reviewed source head: ' + repository)
        require(not git(root, 'diff', '--name-only', reviewed, '--'),
                'Reviewed tracked files changed locally: ' + repository)
        reviewed_tree = git(root, 'rev-parse', reviewed + '^{tree}')
        merged_tree = git(root, 'rev-parse', merged + '^{tree}')
        require(reviewed_tree == merged_tree, 'Reviewed and merged trees differ: ' + repository)
        snapshots.append({'repository': repository, 'reviewed_head': reviewed,
                          'merge_commit_sha': merged, 'reviewed_and_merged_tree': reviewed_tree})
    return snapshots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hdblast-root', required=True, type=Path)
    parser.add_argument('--shared-root', required=True, type=Path)
    parser.add_argument('--merge-evidence', required=True, type=Path,
                        help='Observed merged PRs and successful exact-merge deployments')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--workers', type=int, default=6, choices=range(1, 9))
    args = parser.parse_args()
    public, shared = args.hdblast_root.resolve(), args.shared_root.resolve()
    roots = dict(zip(REVIEWED, [public, shared]))
    merge_evidence_raw = args.merge_evidence.read_bytes()
    merge_evidence = json.loads(merge_evidence_raw)
    snapshots = validate_merge_evidence(merge_evidence, roots)
    manifest_raw = (public / CAPSULE / 'MANIFEST.json').read_bytes()
    require(hashlib.sha256(manifest_raw).hexdigest() == MANIFEST_SHA256,
            'Reviewed source manifest differs')
    manifest = json.loads(manifest_raw)
    require(len(manifest['files']) == 51, 'Expected 51 frozen source members')
    for row in manifest['files']:
        raw = (public / CAPSULE / row['path']).read_bytes()
        require(len(raw) == row['bytes'] and hashlib.sha256(raw).hexdigest() == row['sha256'],
                'Local frozen source differs: ' + row['path'])
    targets = []
    nonce = uuid.uuid4().hex
    def add(root, name, url, kind):
        raw = (root / name).read_bytes()
        targets.append({'path': name, 'url': url,
                        'request_url': url + '?hdblast_readback=' + nonce,
                        'kind': kind, 'expected_bytes': len(raw),
                        'expected_sha256': hashlib.sha256(raw).hexdigest(),
                        'expected_raw': raw})
    public_names = [CAPSULE + '/' + row['path'] for row in manifest['files']]
    public_names += [CAPSULE + '/MANIFEST.json', CAPSULE + '.zip',
                     '.github/workflows/hdblast_quasilinear_shell_evolution.yml',
                     DELIVERY + '/FROZEN_SOURCE_RECEIPT.json',
                     DELIVERY + '/NEXT_SHELL_RESPONSE_TARGET.md',
                     'README.md', 'hdblast/README.md', 'hdblast/CURRENT_STATUS.md',
                     'hdblast/PROGRAM_MAP.md', 'site/index.html', 'site/es/index.html']
    for name in public_names:
        add(public, name, 'https://raw.githubusercontent.com/maldonado-research/HDblast/main/'
            + urllib.parse.quote(name, safe='/'), 'GitHub public source')
    for name in ['README.md', 'index.html', 'projects/hdblast/index.html']:
        add(shared, name, 'https://raw.githubusercontent.com/maldonado-research/maldonado-research.github.io/main/'
            + urllib.parse.quote(name, safe='/'), 'GitHub shared-site source')
    for root, name, url in [
        (public, 'site/index.html', 'https://maldonado-research.github.io/HDblast/'),
        (public, 'site/es/index.html', 'https://maldonado-research.github.io/HDblast/es/'),
        (shared, 'index.html', 'https://maldonado-research.github.io/'),
        (shared, 'projects/hdblast/index.html', 'https://maldonado-research.github.io/projects/hdblast/'),
    ]:
        add(root, name, url, 'Live Pages')
    require(len(targets) == 69, 'Readback target inventory differs')
    require(len({t['url'] for t in targets}) == len(targets), 'Duplicate readback URL')
    def verify(target):
        row = {k: v for k, v in target.items() if k != 'expected_raw'}
        request = urllib.request.Request(target['request_url'], headers={
            'User-Agent': 'HDBLAST-scientific-publication-readback', 'Cache-Control': 'no-cache'})
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                actual = response.read()
                row.update(http_status=response.status, bytes=len(actual),
                           sha256=hashlib.sha256(actual).hexdigest(),
                           byte_equal=actual == target['expected_raw'],
                           final_url=response.geturl())
        except Exception as error:
            row.update(http_status=None, byte_equal=False, error=str(error))
        return row
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(verify, targets))
    success = all(row['http_status'] == 200 and row['byte_equal'] for row in rows)
    receipt = {'schema': 'hdblast-quasilinear-shell-publication-readback-v1',
               'checked_utc': datetime.now(timezone.utc).isoformat(),
               'status': 'PASS' if success else 'FAIL', 'target_count': len(rows),
               'merge_and_deployment_prerequisites': snapshots,
               'merge_evidence_sha256': hashlib.sha256(merge_evidence_raw).hexdigest(),
               'manifest_sha256': MANIFEST_SHA256, 'concurrent_workers': args.workers,
               'targets': rows,
               'scope': 'Dated byte-equality witness for reviewed merged source and deployed HTML; not analytic proof certification.'}
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': receipt['status'], 'targets': len(rows),
                      'failed_urls': [row['url'] for row in rows if not row['byte_equal']],
                      'output': str(args.output)}))
    raise SystemExit(0 if success else 1)


if __name__ == '__main__':
    main()
