"""Read reviewed source, initial evidence, HTML and figure bytes after five exact-merge workflows pass.

The preflight mode is local only. Normal mode gates every HTTP request on
passing reviewed-head checks, exact merge trees and five successful workflows.
This helper never fetches, checks out, merges or otherwise mutates Git state.
"""
import argparse
import concurrent.futures
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import urllib.parse
import urllib.request
import uuid
import zipfile

PUBLIC = 'maldonado-research/HDblast'
SHARED = 'maldonado-research/maldonado-research.github.io'
CAPSULE = 'research/HDBLAST_QUADRATIC_SHELL_FORCING_20261015'
DELIVERY = 'hdblast/continuation/2026.10.15-quadratic-shell-forcing'
WORKFLOW = '.github/workflows/hdblast_quadratic_shell_forcing.yml'
SUMMARIES = {
    PUBLIC: ['README.md', 'hdblast/README.md', 'hdblast/CURRENT_STATUS.md', 'hdblast/PROGRAM_MAP.md',
             'site/index.html', 'site/es/index.html'],
    SHARED: ['README.md', 'index.html', 'projects/hdblast/index.html'],
}
PREMERGE_WORKFLOWS = {'Registered stored-state comparison exact replay',
                      'Registered retained-endpoint source replay',
                      'Quadratic material shell forcing replay'}
POSTMERGE_WORKFLOWS = {
    PUBLIC: PREMERGE_WORKFLOWS | {'Deploy website'}, SHARED: {'pages build and deployment'},
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def load_input(path, pins):
    raw = path.read_bytes()
    require(path.name not in pins, 'Input receipt basename collision')
    pins[path.name] = {'path': str(path), 'bytes': len(raw), 'sha256': sha(raw)}
    return json.loads(raw)


def canonical(name):
    require(isinstance(name, str), 'String path required')
    path = PurePosixPath(name)
    require(bool(name) and not path.is_absolute() and '..' not in path.parts
            and '\\' not in name and path.as_posix() == name, 'Noncanonical source path: '+name)
    return name


def git(root, *args):
    return subprocess.check_output(['git', '-C', str(root), *args])


def textgit(root, *args):
    return git(root, *args).decode().strip()


def tree_files(root, revision, prefix):
    result = {}
    for line in git(root, 'ls-tree', '-r', '-z', revision, '--', prefix).split(b'\0'):
        if not line:
            continue
        entry, path = line.split(b'\t', 1)
        mode, kind, oid = entry.decode().split()
        name = canonical(path.decode())
        require(kind == 'blob' and mode in ('100644', '100755'), 'Regular committed file required: '+name)
        require(name not in result, 'Duplicate committed path')
        result[name] = oid
    return result


def prepared_pins(prepared, require_pr=False):
    require(isinstance(prepared, list) and len(prepared) == 2
            and {r['repository'] for r in prepared} == {PUBLIC, SHARED}, 'Exactly two prepared heads required')
    pins = {r['repository']: r for r in prepared}
    for row in prepared:
        require(all(re.fullmatch('[0-9a-f]{40}', row[k]) is not None
                    for k in ('base', 'head_sha', 'reviewed_tree')), 'Full commit/tree pins required')
        if require_pr:
            require(type(row.get('pr_number')) is int and row['pr_number'] > 0, 'Dynamic prepared pr_number required')
    return pins


def build_targets(roots, prepared, frozen, figures):
    """Expected bytes come only from reviewed immutable Git objects, not untracked files."""
    reviewed = prepared_pins(prepared)
    for repo, pin in reviewed.items():
        root = roots[repo]
        require(textgit(root, 'rev-parse', 'HEAD') == pin['head_sha'], 'Worktree must remain at reviewed head')
        require(textgit(root, 'rev-parse', pin['head_sha']+'^{tree}') == pin['reviewed_tree'], 'Reviewed tree differs')
        require(not textgit(root, 'diff', '--name-only', pin['head_sha'], '--'), 'Tracked worktree differs from reviewed commit')
    cache = {}
    def blob(repo, name):
        key = (repo, canonical(name))
        if key not in cache:
            require(set(tree_files(roots[repo], reviewed[repo]['head_sha'], name)) == {name},
                    'Exactly one regular committed target required: '+name)
            raw = git(roots[repo], 'show', reviewed[repo]['head_sha']+':'+name)
            local = roots[repo]/name
            require(local.is_file() and not local.is_symlink() and local.read_bytes() == raw,
                    'Local/committed target bytes differ: '+name)
            cache[key] = raw
        return cache[key]
    require(frozen['status'] == 'PASS' and frozen['capsule'] == PurePosixPath(CAPSULE).name,
            'Passing matching external source freeze required')
    manifest_raw = blob(PUBLIC, CAPSULE+'/MANIFEST.json')
    require(len(manifest_raw) == frozen['manifest']['bytes']
            and sha(manifest_raw) == frozen['manifest']['sha256'], 'External source-manifest pin differs')
    manifest = json.loads(manifest_raw)
    require(manifest['capsule'] == PurePosixPath(CAPSULE).name, 'Manifest capsule differs')
    names = set()
    for row in manifest['files']:
        source_parts = PurePosixPath(row['path']).parts
        require(not ({p.lower() for p in source_parts} & {'internal', 'private', '_internal', '_private'})
                and (not row['path'].startswith('millennium/')
                     or row['path'] == 'millennium/HDBLAST_METHODS_TRANSFER_PUBLIC_SAFE.md'),
                'Private audit or paper-body path in publication source')
        name = CAPSULE+'/'+canonical(row['path'])
        require(name not in names and row['path'] != 'MANIFEST.json', 'Duplicate/reserved manifest path')
        require(type(row['bytes']) is int and row['bytes'] >= 0
                and re.fullmatch('[0-9a-f]{64}', row['sha256']) is not None, 'Invalid source manifest record')
        names.add(name)
        raw = blob(PUBLIC, name)
        require(len(raw) == row['bytes'] and sha(raw) == row['sha256'], 'Frozen source bytes differ: '+name)
    require(bool(names), 'Empty source capsule')
    source_names = names | {CAPSULE+'/MANIFEST.json'}
    require(set(tree_files(roots[PUBLIC], reviewed[PUBLIC]['head_sha'], CAPSULE)) == source_names,
            'Complete committed source inventory differs (including ignored files)')
    require(frozen['source_files'] == len(names) and frozen['archive_members'] == len(source_names),
            'External frozen inventory count differs')
    archive_raw = blob(PUBLIC, CAPSULE+'.zip')
    require(len(archive_raw) == frozen['archive']['bytes'] and sha(archive_raw) == frozen['archive']['sha256'],
            'External archive pin differs')
    zip_prefix = PurePosixPath(CAPSULE).name+'/'
    expected_zip = {zip_prefix+n.removeprefix(CAPSULE+'/'): blob(PUBLIC, n) for n in source_names}
    with zipfile.ZipFile(io.BytesIO(archive_raw)) as archive:
        actual = archive.namelist()
        require(len(actual) == len(set(actual)) == len(expected_zip) and set(actual) == set(expected_zip),
                'ZIP/committed-source member inventory differs')
        for name, raw in expected_zip.items():
            require(archive.read(name) == raw, 'ZIP/committed source bytes differ: '+name)
    workflow_raw = blob(PUBLIC, WORKFLOW)
    initial_manifest_name = DELIVERY+'/INITIAL_REPLAY_MANIFEST.json'
    initial_raw = blob(PUBLIC, initial_manifest_name)
    initial = json.loads(initial_raw)
    delivery_names = {initial_manifest_name}
    for row in initial['files']:
        name = DELIVERY+'/'+canonical(row['path'])
        require(name not in delivery_names, 'Duplicate/reserved initial evidence path')
        require(type(row['bytes']) is int and row['bytes'] >= 0
                and re.fullmatch('[0-9a-f]{64}', row['sha256']) is not None, 'Invalid initial evidence record')
        delivery_names.add(name)
        raw = blob(PUBLIC, name)
        require(len(raw) == row['bytes'] and sha(raw) == row['sha256'], 'Initial evidence bytes differ: '+name)
    require(set(tree_files(roots[PUBLIC], reviewed[PUBLIC]['head_sha'], DELIVERY)) == delivery_names,
            'Complete initial evidence Git inventory differs')
    require(figures['status'] == 'PASS' and figures['all_asset_bytes_equal_to_source'] is True
            and len(figures['site_assets']) in (0, 2), 'Zero or two explicitly declared figure assets required')
    if figures['site_assets']:
        source_figure = blob(PUBLIC, CAPSULE+'/'+canonical(figures['source']['path']))
        require(len(source_figure) == figures['source']['bytes'] and sha(source_figure) == figures['source']['sha256'],
                'Figure source differs from integration receipt')
    else:
        require(figures.get('mode') == 'no_figure_assets' and figures['source'] is None,
                'Explicit no-figure receipt required when no assets are selected')
        source_figure = None
    asset_keys = set()
    for asset in figures['site_assets']:
        repo, path, url = asset['repository'], canonical(asset['path']), urllib.parse.urlsplit(asset['url'])
        require(repo in reviewed and (repo, path) not in asset_keys, 'Duplicate or unknown figure repository/path')
        asset_keys.add((repo, path))
        require(url.scheme == 'https' and url.netloc == 'maldonado-research.github.io'
                and not url.query and not url.fragment
                and (repo != PUBLIC or url.path.startswith('/HDblast/')), 'Canonical live Pages asset URL required')
        raw = blob(repo, path)
        require(raw == source_figure and len(raw) == asset['bytes'] and sha(raw) == asset['sha256'],
                'Committed figure asset differs from pinned source')
    targets, nonce = [], uuid.uuid4().hex
    def add(repo, name, url, kind):
        raw = blob(repo, name)
        targets.append({'repository': repo, 'path': name, 'url': url,
                        'request_url': url+'?hdblast_readback='+nonce, 'kind': kind,
                        'expected_bytes': len(raw), 'expected_sha256': sha(raw), 'expected_raw': raw})
    public_names = sorted(source_names)+[CAPSULE+'.zip', WORKFLOW]+sorted(delivery_names)+SUMMARIES[PUBLIC]
    for repo, paths in ((PUBLIC, public_names), (SHARED, SUMMARIES[SHARED])):
        for name in paths:
            add(repo, name, 'https://raw.githubusercontent.com/'+repo+'/main/'+urllib.parse.quote(name, safe='/'), 'GitHub source')
    for repo, name, url in (
        (PUBLIC, 'site/index.html', 'https://maldonado-research.github.io/HDblast/'),
        (PUBLIC, 'site/es/index.html', 'https://maldonado-research.github.io/HDblast/es/'),
        (SHARED, 'index.html', 'https://maldonado-research.github.io/'),
        (SHARED, 'projects/hdblast/index.html', 'https://maldonado-research.github.io/projects/hdblast/')):
        add(repo, name, url, 'Live Pages HTML')
    for asset in figures['site_assets']:
        repo, name = asset['repository'], asset['path']
        add(repo, name, 'https://raw.githubusercontent.com/'+repo+'/main/'+urllib.parse.quote(name, safe='/'), 'GitHub figure asset')
        add(repo, name, asset['url'], 'Live Pages figure asset')
    counts = {'source_files': len(names), 'source_with_manifest': len(source_names),
              'zip_and_workflow': 2, 'initial_delivery_files': len(delivery_names),
              'raw_summaries': 9, 'rendered_html_pages': 4,
              'raw_figure_assets': len(figures['site_assets']), 'live_figure_assets': len(figures['site_assets'])}
    total = sum(v for k, v in counts.items() if k != 'source_files')
    require(len(targets) == total and len({r['url'] for r in targets}) == total, 'Readback target inventory differs')
    pins = {'manifest': sha(manifest_raw), 'archive': sha(archive_raw), 'workflow': sha(workflow_raw),
            'runner': sha(blob(PUBLIC, CAPSULE+'/reproduce.py')),
            'verifier': sha(blob(PUBLIC, CAPSULE+'/verify_manifest.py')),
            'initial_delivery_manifest': sha(initial_raw), 'figure': sha(source_figure) if source_figure is not None else None}
    return targets, {**counts, 'target_count': total}, pins


def validate_observations(pre, post, merges, prepared, roots, input_pins):
    reviewed = prepared_pins(prepared, require_pr=True)
    require(pre['status'] == 'PASS' and len(pre['pull_requests']) == 2
            and pre['prepared_commits_sha256'] == input_pins['PREPARED_COMMITS.json']['sha256'],
            'Passing premerge observation of current prepared heads required')
    pre_rows = {r['repository']: r for r in pre['pull_requests']}
    require(set(pre_rows) == set(reviewed), 'Both premerge repositories required')
    for repo, pin in reviewed.items():
        row = pre_rows[repo]
        require(row['status'] == 'PASS' and row['number'] == pin['pr_number']
                and row['head_sha'] == pin['head_sha'] and row['base_sha'] == pin['base'], 'Reviewed premerge PR differs')
        expected = PREMERGE_WORKFLOWS if repo == PUBLIC else set()
        require(len(row['checks']) == len(row['workflow_runs']) == len(expected)
                and len({c['id'] for c in row['checks']}) == len(expected)
                and len({r['id'] for r in row['workflow_runs']}) == len(expected)
                and {r['name'] for r in row['workflow_runs']} == expected, 'Distinct required premerge checks/workflows differ')
        for run in row['checks']+row['workflow_runs']:
            require(run['head_sha'] == pin['head_sha'] and run['status'] == 'completed'
                    and run['conclusion'] == 'success', 'Successful reviewed-head premerge check required')
        require(all(s['state'] == 'success' for s in row['legacy_statuses']), 'Legacy premerge status failed')
    require(isinstance(merges, list) and len(merges) == 2
            and {r['repository'] for r in merges} == set(reviewed), 'Two merged PR observations required')
    require(post['status'] == 'PASS' and len(post['runs']) == 5
            and len({(r['repository'], r['id']) for r in post['runs']}) == 5
            and post['prepared_commits_sha256'] == input_pins['PREPARED_COMMITS.json']['sha256']
            and post['merge_results_sha256'] == input_pins['MERGE_RESULTS.json']['sha256'],
            'Five successful exact-merge workflows with current receipt pins required')
    snapshots = []
    for row in merges:
        repo, merged = row['repository'], row['merge_commit_sha']
        pin = reviewed[repo]
        require(row['merged'] is True and row['head_sha'] == pin['head_sha'] and row['number'] == pin['pr_number']
                and row['url'] == f"https://github.com/{repo}/pull/{pin['pr_number']}"
                and re.fullmatch('[0-9a-f]{40}', merged) is not None, 'Reviewed merged PR differs')
        require(row['reviewed_tree'] == row['merged_tree'] == pin['reviewed_tree']
                and row['merged_tree_equals_reviewed_tree'] is True, 'Observed merge tree differs')
        matching = [r for r in post['runs'] if r['repository'] == repo]
        require({r['name'] for r in matching} == POSTMERGE_WORKFLOWS[repo], 'Required exact-merge workflow names differ')
        for run in matching:
            require(run['head_sha'] == merged and run['status'] == 'completed' and run['conclusion'] == 'success'
                    and run['event'] == ('push' if repo == PUBLIC else 'dynamic')
                    and run['html_url'] == f"https://github.com/{repo}/actions/runs/{run['id']}",
                    'Successful exact-merge workflow/deployment required')
        # Root may read/fetch these immutable merge objects beforehand. This helper never mutates Git.
        require(subprocess.run(['git', '-C', str(roots[repo]), 'cat-file', '-e', merged+'^{commit}'],
                               capture_output=True).returncode == 0,
                'Exact merge object missing locally; obtain it before readback')
        require(textgit(roots[repo], 'rev-parse', merged+'^{tree}') == pin['reviewed_tree'],
                'Actual Git merge tree differs from reviewed tree')
        snapshots.append({'repository': repo, 'reviewed_head': pin['head_sha'],
                          'merge_commit_sha': merged, 'reviewed_and_merged_tree': pin['reviewed_tree']})
    return snapshots


def main():
    task = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared-commits', type=Path, default=task/'PREPARED_COMMITS.json')
    parser.add_argument('--source-freeze-receipt', type=Path, default=task/'frozen-source/SOURCE_FREEZE_RECEIPT.json')
    parser.add_argument('--figure-integration', type=Path, default=task/'publication/FIGURE_INTEGRATION_RECEIPT.json')
    parser.add_argument('--premerge-checks', type=Path, default=task/'PRE_MERGE_CHECKS.json')
    parser.add_argument('--merge-results', type=Path, default=task/'MERGE_RESULTS.json')
    parser.add_argument('--postmerge-workflows', type=Path, default=task/'POST_MERGE_WORKFLOW_RUNS.json')
    parser.add_argument('--hdblast-root', type=Path)
    parser.add_argument('--shared-root', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=6, choices=range(1, 9))
    parser.add_argument('--preflight-only', action='store_true', help='Authenticate local committed inventory; no HTTP or merge observation')
    args = parser.parse_args()
    inputs = {}
    prepared = load_input(args.prepared_commits, inputs)
    reviewed = prepared_pins(prepared)
    roots = {PUBLIC: (args.hdblast_root or Path(reviewed[PUBLIC]['worktree'])).resolve(),
             SHARED: (args.shared_root or Path(reviewed[SHARED]['worktree'])).resolve()}
    frozen, figures = load_input(args.source_freeze_receipt, inputs), load_input(args.figure_integration, inputs)
    targets, inventory, pins = build_targets(roots, prepared, frozen, figures)
    if args.preflight_only:
        receipt = {'status': 'PREPARED_NOT_EXECUTED', 'reviewed': reviewed, 'inventory': inventory,
                   'pins': pins, 'input_receipts': inputs, 'helper_sha256': sha(Path(__file__).read_bytes()),
                   'target_inventory': [{k: v for k, v in r.items() if k not in ('expected_raw', 'request_url')} for r in targets],
                   'scope': 'Reviewed committed-object and figure target authentication only; no HTTP, Git mutation or deployment observation.'}
        args.output.write_text(json.dumps(receipt, indent=2)+'\n')
        print(json.dumps({'status': receipt['status'], 'inventory': inventory}, indent=2))
        return
    pre = load_input(args.premerge_checks, inputs)
    merges = load_input(args.merge_results, inputs)
    post = load_input(args.postmerge_workflows, inputs)
    snapshots = validate_observations(pre, post, merges, prepared, roots, inputs)
    def fetch(target):
        row = {k: v for k, v in target.items() if k != 'expected_raw'}
        request = urllib.request.Request(target['request_url'], headers={
            'User-Agent': 'HDBLAST-scientific-publication-readback', 'Cache-Control': 'no-cache'})
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                raw = response.read()
                row.update(http_status=response.status, bytes=len(raw), sha256=sha(raw),
                           byte_equal=raw == target['expected_raw'], final_url=response.geturl())
        except Exception as error:
            row.update(http_status=None, byte_equal=False, error=str(error))
        return row
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(fetch, targets))
    success = all(r['http_status'] == 200 and r['byte_equal'] for r in rows)
    receipt = {'schema': 'hdblast-quadratic-shell-forcing-publication-readback-v1',
               'checked_utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS' if success else 'FAIL',
               'target_count': len(rows), 'inventory': inventory, 'reviewed': reviewed,
               'merge_and_deployment_prerequisites': snapshots, 'pins': pins, 'input_receipts': inputs,
               'helper_sha256': sha(Path(__file__).read_bytes()), 'concurrent_workers': args.workers, 'targets': rows,
               'scope': 'Dated byte-equality witness for reviewed merged source, initial evidence, deployed HTML and figures; not analytic proof certification.'}
    args.output.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': receipt['status'], 'targets': len(rows),
                      'failed_urls': [r['url'] for r in rows if not r['byte_equal']], 'output': str(args.output)}))
    raise SystemExit(0 if success else 1)


if __name__ == '__main__':
    main()
