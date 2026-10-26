#!/usr/bin/env python3
"""Retain complete observed GitHub CI and PR state for an exact commit."""
import argparse
import json
from pathlib import Path
import subprocess


def api(path):
    return json.loads(subprocess.check_output(['gh', 'api', path], text=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', required=True)
    parser.add_argument('--sha', required=True)
    parser.add_argument('--phase', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--pr', type=int)
    args = parser.parse_args()
    repo = args.repository
    name = repo.split('/')[-1]
    prefix = name+'-'+args.phase
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    queries = {
        'RUNS': f'repos/{repo}/actions/runs?head_sha={args.sha}&per_page=100',
        'CHECKS': f'repos/{repo}/commits/{args.sha}/check-runs?per_page=100',
        'STATUS': f'repos/{repo}/commits/{args.sha}/status?per_page=100',
    }
    if args.pr:
        queries['PR'] = f'repos/{repo}/pulls/{args.pr}'
    results = {key: api(path) for key, path in queries.items()}
    for key, value in results.items():
        if key in ('RUNS', 'CHECKS'):
            field = 'workflow_runs' if key == 'RUNS' else 'check_runs'
            if value['total_count'] != len(value[field]):
                raise ValueError('INCOMPLETE_API_PAGE:'+key)
        (output/(prefix+'_'+key+'.json')).write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    runs, checks = results['RUNS']['workflow_runs'], results['CHECKS']['check_runs']
    rows = runs+checks
    complete = bool(runs) and bool(checks) and all(r['head_sha'] == args.sha and r['status'] == 'completed'
                 and r['conclusion'] == 'success' for r in rows)
    legacy = results['STATUS']
    complete = complete and (legacy['total_count'] == 0 or legacy['state'] == 'success')
    print(json.dumps({'repository': repo, 'head': args.sha, 'all_observed_success': complete,
        'runs': [{'id': r['id'], 'name': r['name'], 'status': r['status'], 'conclusion': r['conclusion']} for r in runs],
        'checks': len(checks)}))


if __name__ == '__main__':
    main()
