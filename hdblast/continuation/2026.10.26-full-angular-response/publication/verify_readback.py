#!/usr/bin/env python3
"""Read the exact published changed files and complete live HTML pages."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen


def fetch(row):
    try:
        request = Request(row['url'], headers={'User-Agent': 'HDBLAST-publication-readback/1.0',
                                              'Cache-Control': 'no-cache'})
        with urlopen(request, timeout=35) as response:
            data, status, final = response.read(), response.status, response.url
        digest = hashlib.sha256(data).hexdigest()
        blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        ok = status == 200 and digest == row['sha256'] and len(data) == row['bytes']
        if 'git_blob_sha' in row:
            ok = ok and blob == row['git_blob_sha']
        return {**row, 'http_status': status, 'final_url': final,
                'actual_sha256': digest, 'actual_git_blob_sha': blob,
                'actual_bytes': len(data), 'matches': ok}
    except Exception as error:
        return {**row, 'matches': False, 'error': str(error)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--research-sha', required=True)
    parser.add_argument('--shared-sha', required=True)
    parser.add_argument('--expected', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--phase', choices=['source', 'live', 'all'], default='all')
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('FRESH_READBACK_OUTPUT_REQUIRED')
    expected = json.loads(args.expected.read_text())
    targets = []
    for row in expected:
        commit = args.research_sha if row['repository'].endswith('/HDblast') else args.shared_sha
        targets.append({**row, 'kind': 'immutable_source', 'commit': commit,
                        'url': 'https://raw.githubusercontent.com/'+row['repository']+'/'+commit+'/'+row['path']})
    live = [('maldonado-research/HDblast', 'site/index.html', 'https://maldonado-research.github.io/HDblast/'),
            ('maldonado-research/HDblast', 'site/es/index.html', 'https://maldonado-research.github.io/HDblast/es/'),
            ('maldonado-research/maldonado-research.github.io', 'projects/hdblast/index.html',
             'https://maldonado-research.github.io/projects/hdblast/')]
    for repo, path, url in live:
        row = next(r for r in expected if r['repository'] == repo and r['path'] == path)
        targets.append({**row, 'kind': 'complete_live_HTML',
                        'url': url+'?checkpoint=full-angular-response-20261026-'+args.research_sha[:12]})
    if args.phase != 'all':
        kind = 'immutable_source' if args.phase == 'source' else 'complete_live_HTML'
        targets = [row for row in targets if row['kind'] == kind]
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(fetch, targets))
    passed = all(row['matches'] for row in results)
    receipt = {'schema': 'HDBLAST_FULL_ANGULAR_PUBLIC_READBACK_V1',
               'status': 'PASS' if passed else 'INCOMPLETE',
               'research_commit': args.research_sha, 'shared_commit': args.shared_sha,
               'expected_list_sha256': hashlib.sha256(args.expected.read_bytes()).hexdigest(),
               'phase': args.phase,
               'immutable_source_targets': sum(r['kind'] == 'immutable_source' for r in results),
               'complete_live_HTML_targets': sum(r['kind'] == 'complete_live_HTML' for r in results),
               'matched': sum(row['matches'] for row in results), 'targets': results}
    args.output.write_text(json.dumps(receipt, sort_keys=True, indent=2)+'\n')
    print(json.dumps({'status': receipt['status'], 'matched': receipt['matched'],
                      'targets': len(results), 'failed': [r['url'] for r in results if not r['matches']]}))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
