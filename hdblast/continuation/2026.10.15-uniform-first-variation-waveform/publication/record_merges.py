"""Record already completed PR merges against fetched immutable Git objects."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent.parent
prepared = json.loads((ROOT / 'PREPARED_COMMITS.json').read_text())


def observe(row):
    repo, number = row['repository'], row['pr_number']
    pr = json.loads(subprocess.check_output(['gh', 'api', f'repos/{repo}/pulls/{number}'], text=True))
    if not (pr['merged'] is True and pr['head']['sha'] == row['head_sha']):
        raise RuntimeError('Reviewed PR must already be merged')
    merge = pr['merge_commit_sha']
    def git(*args):
        return subprocess.check_output(['git', '-C', row['worktree'], *args], text=True).strip()
    tree = git('rev-parse', merge + '^{tree}')
    if tree != row['reviewed_tree']:
        raise RuntimeError('Actual merged tree differs from the reviewed tree')
    if git('rev-parse', 'HEAD') != row['head_sha']:
        raise RuntimeError('Reviewed worktree must stay at its original head')
    parents = git('rev-list', '--parents', '-n', '1', merge).split()
    if parents != [merge, row['base'], row['head_sha']]:
        raise RuntimeError('Reviewed merge parents differ')
    return {'repository': repo, 'number': number, 'url': pr['html_url'], 'merged': True,
            'head_sha': row['head_sha'], 'merge_commit_sha': merge, 'merged_at': pr['merged_at'],
            'reviewed_tree': row['reviewed_tree'], 'merged_tree': tree,
            'merged_tree_equals_reviewed_tree': True}


with ThreadPoolExecutor(max_workers=2) as pool:
    results = list(pool.map(observe, prepared))
(ROOT / 'MERGE_RESULTS.json').write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps(results, indent=2))
