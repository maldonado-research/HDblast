"""Observe exact reviewed PR heads and their three named workflows, without mutations."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

PUBLIC = 'maldonado-research/HDblast'
SHARED = 'maldonado-research/maldonado-research.github.io'
WORKFLOWS = {
    'Registered stored-state comparison exact replay',
    'Registered retained-endpoint source replay',
    'Uniform first-variation waveform enclosure replay',
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def api(path):
    return json.loads(subprocess.check_output(['gh', 'api', path], text=True, timeout=60))


def load_prepared(path):
    raw = path.read_bytes()
    rows = json.loads(raw)
    require(isinstance(rows, list) and len(rows) == 2
            and {r['repository'] for r in rows} == {PUBLIC, SHARED}, 'Exactly two prepared repositories required')
    for row in rows:
        require(all(re.fullmatch('[0-9a-f]{40}', row[k]) for k in ('base', 'head_sha', 'reviewed_tree')),
                'Complete reviewed commit/tree pins required')
        require(type(row['pr_number']) is int and row['pr_number'] > 0, 'Dynamic positive pr_number required')
    return rows, hashlib.sha256(raw).hexdigest()


def observe(item, get=api):
    repo, head, number = item['repository'], item['head_sha'], item['pr_number']
    queries = {'pr': f'repos/{repo}/pulls/{number}',
               'checks': f'repos/{repo}/commits/{head}/check-runs?per_page=100',
               'legacy': f'repos/{repo}/commits/{head}/status?per_page=100'}
    with ThreadPoolExecutor(max_workers=3) as pool:
        result = dict(pool.map(lambda pair: (pair[0], get(pair[1])), queries.items()))
    pr, payload = result['pr'], result['checks']
    checks, statuses = payload['check_runs'], result['legacy']['statuses']
    require(payload['total_count'] == len(checks), 'Check-run pagination requires explicit inspection')
    require(result['legacy']['total_count'] == len(statuses), 'Legacy-status pagination requires explicit inspection')
    require(pr['number'] == number and pr['head']['sha'] == head and pr['base']['sha'] == item['base'],
            'Reviewed PR head/base moved; preserve and reassess before merging')
    require(all(c['head_sha'] == head for c in checks)
            and len({c['id'] for c in checks}) == len(checks), 'Distinct exact-head check runs required')
    bad = [c for c in checks if c['status'] == 'completed' and c['conclusion'] not in ('success', 'neutral', 'skipped')]
    require(not bad and not any(s['state'] in ('failure', 'error') for s in statuses),
            'A premerge check or legacy status failed')
    workflows = []
    for check in checks:
        url = check.get('details_url') or check.get('html_url') or ''
        found = re.fullmatch(r'https://github\.com/' + re.escape(repo) + r'/actions/runs/(\d+)(?:/job/\d+)?', url)
        require(found is not None, 'Expected GitHub Actions check provenance required')
        run = get(f'repos/{repo}/actions/runs/{found.group(1)}')
        require(run['head_sha'] == head and run['event'] == 'pull_request', 'Workflow head or event differs')
        workflows.append({k: run[k] for k in ('id', 'name', 'head_sha', 'event', 'status', 'conclusion', 'html_url')})
    expected = WORKFLOWS if repo == PUBLIC else set()
    require(len({r['id'] for r in workflows}) == len(workflows), 'Duplicate premerge workflow run')
    require({r['name'] for r in workflows} <= expected, 'Unexpected premerge workflow; inspect explicitly')
    passed = (len(checks) == len(expected) and {r['name'] for r in workflows} == expected
              and all(c['status'] == 'completed' and c['conclusion'] == 'success' for c in checks)
              and all(r['status'] == 'completed' and r['conclusion'] == 'success' for r in workflows)
              and all(s['state'] == 'success' for s in statuses)
              and pr['mergeable'] is True and pr['mergeable_state'] == 'clean'
              and pr['draft'] is False and pr['state'] == 'open')
    return {'repository': repo, 'number': number, 'head_sha': head, 'base_sha': pr['base']['sha'],
            'mergeable': pr['mergeable'], 'mergeable_state': pr['mergeable_state'],
            'checks': checks, 'legacy_statuses': statuses, 'workflow_runs': workflows,
            'status': 'PASS' if passed else 'PENDING'}


def main():
    task = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared-commits', type=Path, default=task/'PREPARED_COMMITS.json')
    parser.add_argument('--output', type=Path, default=task/'PRE_MERGE_CHECKS.json')
    parser.add_argument('--once', action='store_true', help='One observation; pending exits with code 2')
    parser.add_argument('--max-observations', type=int, default=90)
    parser.add_argument('--interval-seconds', type=int, default=20, choices=range(1, 61))
    args = parser.parse_args()
    require(args.max_observations > 0, 'Positive observation budget required')
    prepared, pin = load_prepared(args.prepared_commits)
    previous = None
    for attempt in range(1 if args.once else args.max_observations):
        rows = [observe(item) for item in prepared]
        complete = all(r['status'] == 'PASS' for r in rows)
        receipt = {'observed_utc': datetime.now(timezone.utc).isoformat(),
                   'status': 'PASS' if complete else 'PENDING', 'pull_requests': rows,
                   'prepared_commits_sha256': pin}
        args.output.write_text(json.dumps(receipt, indent=2)+'\n')
        summary = [(r['repository'], r['number'], r['status'],
                    [(c['name'], c['status'], c['conclusion']) for c in r['checks']]) for r in rows]
        if summary != previous:
            print(json.dumps(summary), flush=True)
            previous = summary
        if complete:
            print('ALL_REVIEWED_PREMERGE_WORKFLOWS_PASSED', flush=True)
            return
        if args.once:
            raise SystemExit(2)
        if attempt + 1 < args.max_observations:
            time.sleep(args.interval_seconds)
    raise RuntimeError('Premerge observation window exceeded')


if __name__ == '__main__':
    main()
