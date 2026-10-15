"""Synthetic, network-free gate checks for newly prepared publication observers."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parent
PUBLIC, SHARED = 'maldonado-research/HDblast', 'maldonado-research/maldonado-research.github.io'


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def require(ok, message):
    if not ok:
        raise ValueError(message)


def rejects(action):
    try:
        action()
    except RuntimeError:
        return True
    raise ValueError('Synthetic invalid observation was accepted')


def main():
    pre, post, rb = (load(n) for n in ('monitor_premerge', 'monitor_postmerge', 'readback_published'))
    prepared = [{'repository': repo, 'base': 'b'*40, 'head_sha': 'a'*40,
                 'reviewed_tree': 'c'*40, 'pr_number': number, 'worktree': '/synthetic/'+str(number)}
                for repo, number in ((PUBLIC, 94), (SHARED, 81))]
    states = {PUBLIC: [], SHARED: []}
    workflows = {}
    for k, name in enumerate(sorted(pre.WORKFLOWS)):
        run_id, check_id = 1000+k, 2000+k
        states[PUBLIC].append({'id': check_id, 'head_sha': 'a'*40, 'status': 'completed', 'conclusion': 'success',
                               'name': 'synthetic-job', 'details_url': f'https://github.com/{PUBLIC}/actions/runs/{run_id}/job/{check_id}'})
        workflows[run_id] = {'id': run_id, 'name': name, 'head_sha': 'a'*40, 'event': 'pull_request',
                              'status': 'completed', 'conclusion': 'success',
                              'html_url': f'https://github.com/{PUBLIC}/actions/runs/{run_id}'}
    def pre_api(path):
        repo = PUBLIC if '/HDblast/' in path else SHARED
        item = next(p for p in prepared if p['repository'] == repo)
        if '/pulls/' in path:
            return {'number': item['pr_number'], 'head': {'sha': item['head_sha']}, 'base': {'sha': item['base']},
                    'mergeable': True, 'mergeable_state': 'clean', 'draft': False, 'state': 'open'}
        if '/check-runs?' in path:
            return {'check_runs': states[repo], 'total_count': len(states[repo])}
        if '/status?' in path:
            return {'statuses': [], 'total_count': 0}
        return workflows[int(path.rsplit('/', 1)[1])]
    pre_rows = [pre.observe(p, pre_api) for p in prepared]
    require(all(r['status'] == 'PASS' for r in pre_rows), 'Synthetic positive premerge failed')
    controls = {}
    def pre_modified(modify):
        def get(path):
            value = copy.deepcopy(pre_api(path))
            modify(path, value)
            return value
        return get
    controls['moved_pr_head'] = rejects(lambda: pre.observe(prepared[0], pre_modified(
        lambda path, value: value['head'].update(sha='d'*40) if '/pulls/' in path else None)))
    controls['failed_premerge_check'] = rejects(lambda: pre.observe(prepared[0], pre_modified(
        lambda path, value: value['check_runs'][0].update(conclusion='failure') if '/check-runs?' in path else None)))
    controls['wrong_workflow_name'] = rejects(lambda: pre.observe(prepared[0], pre_modified(
        lambda path, value: value.update(name='unrelated workflow') if '/actions/runs/' in path else None)))
    pending = pre.observe(prepared[0], pre_modified(
        lambda path, value: value['check_runs'][0].update(status='in_progress', conclusion=None) if '/check-runs?' in path else None))
    require(pending['status'] == 'PENDING', 'Incomplete premerge promoted to PASS')
    merges = [{'repository': p['repository'], 'number': p['pr_number'], 'url': f"https://github.com/{p['repository']}/pull/{p['pr_number']}",
               'merged': True, 'head_sha': p['head_sha'], 'merge_commit_sha': ('e' if p['repository'] == PUBLIC else 'f')*40,
               'reviewed_tree': p['reviewed_tree'], 'merged_tree': p['reviewed_tree'], 'merged_tree_equals_reviewed_tree': True}
              for p in prepared]
    post.validate_merges(merges, prepared)
    post_runs = {}
    for merge in merges:
        repo = merge['repository']
        post_runs[repo] = [{'id': 3000+i+(10 if repo == SHARED else 0), 'name': name,
                            'head_sha': merge['merge_commit_sha'], 'event': 'push' if repo == PUBLIC else 'dynamic',
                            'status': 'completed', 'conclusion': 'success',
                            'html_url': f'https://github.com/{repo}/actions/runs/{3000+i+(10 if repo == SHARED else 0)}'}
                           for i, name in enumerate(sorted(post.REQUIRED[repo]))]
    def post_api(path):
        repo = PUBLIC if '/HDblast/' in path else SHARED
        return {'total_count': len(post_runs[repo]), 'workflow_runs': post_runs[repo]}
    rows = [r for m in merges for r in post.observe(m, post_api)]
    require(len(rows) == 5, 'Synthetic exact-merge observation failed')
    corrupt_merge = copy.deepcopy(merges)
    corrupt_merge[0]['head_sha'] = 'd'*40
    controls['unreviewed_merge_head'] = rejects(lambda: post.validate_merges(corrupt_merge, prepared))
    def failed_post(path):
        value = copy.deepcopy(post_api(path))
        value['workflow_runs'][0]['conclusion'] = 'failure'
        return value
    controls['failed_exact_merge_workflow'] = rejects(lambda: post.observe(merges[0], failed_post))
    pre_receipt = {'status': 'PASS', 'pull_requests': pre_rows, 'prepared_commits_sha256': '1'*64}
    post_receipt = {'status': 'PASS', 'runs': [{**r, 'repository': m['repository']} for m in merges for r in post_runs[m['repository']]],
                    'prepared_commits_sha256': '1'*64, 'merge_results_sha256': '2'*64}
    inputs = {'PREPARED_COMMITS.json': {'sha256': '1'*64}, 'MERGE_RESULTS.json': {'sha256': '2'*64}}
    roots = {r: Path('/synthetic/'+str(i)) for i, r in enumerate((PUBLIC, SHARED))}
    # Stub Git-read observations only; no Git command or HTTP request is executed.
    rb.subprocess = SimpleNamespace(run=lambda *a, **k: SimpleNamespace(returncode=0))
    rb.textgit = lambda root, *args: 'c'*40
    require(len(rb.validate_observations(pre_receipt, post_receipt, merges, prepared, roots, inputs)) == 2,
            'Synthetic complete prerequisite gate failed')
    bad = copy.deepcopy(post_receipt)
    bad['status'] = 'PENDING'
    controls['pending_postmerge_blocks_readback'] = rejects(lambda: rb.validate_observations(pre_receipt, bad, merges, prepared, roots, inputs))
    bad = copy.deepcopy(post_receipt)
    bad['runs'].pop()
    controls['missing_fifth_workflow_blocks_readback'] = rejects(lambda: rb.validate_observations(pre_receipt, bad, merges, prepared, roots, inputs))
    bad = copy.deepcopy(post_receipt)
    bad['runs'][1]['head_sha'] = 'd'*40
    controls['wrong_merge_workflow_head_blocks_readback'] = rejects(lambda: rb.validate_observations(pre_receipt, bad, merges, prepared, roots, inputs))
    bad = copy.deepcopy(post_receipt)
    bad['runs'][1]['id'] = bad['runs'][0]['id']
    controls['duplicate_run_blocks_readback'] = rejects(lambda: rb.validate_observations(pre_receipt, bad, merges, prepared, roots, inputs))
    source_pins = {}
    for name in ('monitor_premerge.py', 'monitor_postmerge.py', 'readback_published.py'):
        raw = (ROOT/name).read_bytes()
        compile(raw, str(ROOT/name), 'exec')
        run = subprocess.run(['python', '-I', '-B', str(ROOT/name), '--help'], capture_output=True)
        require(run.returncode == 0, 'CLI help smoke failed: '+name)
        source_pins[name] = hashlib.sha256(raw).hexdigest()
    receipt = {'status': 'PASS_SYNTHETIC_NETWORK_FREE_TOOL_GATES', 'positive_schema_cases': 6,
               'pending_observation_remains_pending': True, 'invalid_observations_rejected': controls,
               'source_sha256': source_pins, 'http_requests': 0, 'git_commands': 0,
               'scope': 'Preparation smoke tests only. No actual PR, merge, CI, inventory or live readback observed.'}
    (ROOT/'MONITOR_READBACK_TOOL_VALIDATION.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
