"""Require five scientific/regression/Pages runs at exact reviewed merge commits."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

PUBLIC = 'maldonado-research/HDblast'
SHARED = 'maldonado-research/maldonado-research.github.io'
REQUIRED = {
    PUBLIC: {'Registered stored-state comparison exact replay', 'Registered retained-endpoint source replay',
             'Coordinate shell source coefficients replay', 'Deploy website'},
    SHARED: {'pages build and deployment'},
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def api(path):
    return json.loads(subprocess.check_output(['gh', 'api', path], text=True, timeout=60))


def validate_merges(merges, prepared):
    require(isinstance(merges, list) and len(merges) == 2
            and {r['repository'] for r in merges} == set(REQUIRED), 'Exactly two merged PR records required')
    pins = {r['repository']: r for r in prepared}
    require(len(prepared) == len(pins) == 2 and set(pins) == set(REQUIRED), 'Exactly two prepared heads required')
    for row in merges:
        p = pins[row['repository']]
        require(type(p['pr_number']) is int and row['number'] == p['pr_number']
                and row['merged'] is True and row['head_sha'] == p['head_sha']
                and row['reviewed_tree'] == row['merged_tree'] == p['reviewed_tree']
                and row['merged_tree_equals_reviewed_tree'] is True
                and re.fullmatch('[0-9a-f]{40}', row['merge_commit_sha']) is not None,
                'Merged PR does not match reviewed dynamic PR/head/tree')


def observe(item, get=api):
    repo, head = item['repository'], item['merge_commit_sha']
    payload = get(f'repos/{repo}/actions/runs?head_sha={head}&per_page=100')
    require(payload['total_count'] == len(payload['workflow_runs']), 'Workflow pagination requires explicit inspection')
    selected = [r for r in payload['workflow_runs'] if r['name'] in REQUIRED[repo]]
    require(len({r['name'] for r in selected}) == len(selected), 'Multiple required workflow runs; inspect reruns explicitly')
    rows = []
    for run in selected:
        require(run['head_sha'] == head and run['event'] == ('push' if repo == PUBLIC else 'dynamic'),
                'Unexpected postmerge workflow head or trigger')
        require(not (run['status'] == 'completed' and run['conclusion'] != 'success'),
                'Required postmerge workflow failed: '+run['html_url'])
        require(run['html_url'] == f"https://github.com/{repo}/actions/runs/{run['id']}", 'Workflow URL differs')
        rows.append({k: run[k] for k in ('id', 'name', 'head_sha', 'event', 'status', 'conclusion', 'html_url')}
                    | {'repository': repo})
    return rows


def main():
    task = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared-commits', type=Path, default=task/'PREPARED_COMMITS.json')
    parser.add_argument('--merge-results', type=Path, default=task/'MERGE_RESULTS.json')
    parser.add_argument('--output', type=Path, default=task/'POST_MERGE_WORKFLOW_RUNS.json')
    parser.add_argument('--once', action='store_true', help='One observation; pending exits with code 2')
    parser.add_argument('--max-observations', type=int, default=90)
    parser.add_argument('--interval-seconds', type=int, default=20, choices=range(1, 61))
    args = parser.parse_args()
    require(args.max_observations > 0, 'Positive observation budget required')
    raw, prepared_raw = args.merge_results.read_bytes(), args.prepared_commits.read_bytes()
    merges, prepared = json.loads(raw), json.loads(prepared_raw)
    validate_merges(merges, prepared)
    previous = None
    for attempt in range(1 if args.once else args.max_observations):
        rows = [r for item in merges for r in observe(item)]
        complete = (len(rows) == 5 and len({(r['repository'], r['id']) for r in rows}) == 5
                    and all(r['status'] == 'completed' and r['conclusion'] == 'success' for r in rows))
        receipt = {'status': 'PASS' if complete else 'PENDING',
                   'observed_utc': datetime.now(timezone.utc).isoformat(), 'runs': rows,
                   'merge_results_sha256': hashlib.sha256(raw).hexdigest(),
                   'prepared_commits_sha256': hashlib.sha256(prepared_raw).hexdigest()}
        args.output.write_text(json.dumps(receipt, indent=2)+'\n')
        summary = [(r['repository'], r['name'], r['status'], r['conclusion']) for r in rows]
        if summary != previous:
            print(json.dumps(summary), flush=True)
            previous = summary
        if complete:
            print('ALL_EXACT_MERGE_WORKFLOWS_AND_DEPLOYMENTS_PASSED', flush=True)
            return
        if args.once:
            raise SystemExit(2)
        if attempt + 1 < args.max_observations:
            time.sleep(args.interval_seconds)
    raise RuntimeError('Postmerge observation window exceeded')


if __name__ == '__main__':
    main()
