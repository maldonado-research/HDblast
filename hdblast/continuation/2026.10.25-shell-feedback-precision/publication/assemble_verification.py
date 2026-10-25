#!/usr/bin/env python3
"""Assemble a publication receipt from retained observed API/readback data."""
import datetime
import hashlib
import json
from pathlib import Path


def need(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    root = Path(__file__).resolve().parent
    def read(name):
        return json.loads((root/name).read_text())
    specs = [
        ('HDblast', '9a45ee0deaf3496fee93491c2803cb168bb36c92',
         'a97052a03a188e63d77cbe468b50f56e7a2e411d',
         '5388dd48cca21bde0388fa752b180992826878c8',
         'acb59120624b2a92939c52a4f00801dbdd005726', 'HDblast-MERGE_COMMIT.json',
         'HDblast-MERGED_PR.json', 85, 10),
        ('maldonado-research.github.io', 'eb97aef8ee7f72795dfffdb1668c942f8e09fcc6',
         'da1e5312d993c72eef5fab5b4895f73531d6a73a',
         '2c642cd8cbf63ef7479eaa6d0f0618ba8cf19b45',
         'fff3ff5f998d04d6f395149c671d05886c377ab9', 'shared-MERGE_COMMIT.json',
         'shared-MERGED_PR.json', 76, 1),
    ]
    final = []
    for name, merge, base, head, tree, commit_file, pr_file, number, minimum_runs in specs:
        commit, pr = read(commit_file), read(pr_file)
        need(commit['sha'] == merge and commit['tree']['sha'] == tree, 'CANONICAL_TREE')
        need([p['sha'] for p in commit['parents']] == [base, head], 'ORDERED_PARENTS')
        need(pr['merged'] and pr['merge_commit_sha'] == merge, 'MERGED_PR')
        pre_runs, pre_checks = read(name+'-PREMERGE_RUNS.json'), read(name+'-PREMERGE_CHECKS.json')
        post_runs, post_checks = read(name+'-POSTMERGE_RUNS.json'), read(name+'-POSTMERGE_CHECKS.json')
        for obj, collection, sha in [(pre_runs, 'workflow_runs', head), (pre_checks, 'check_runs', head),
                                      (post_runs, 'workflow_runs', merge), (post_checks, 'check_runs', merge)]:
            rows = obj[collection]
            need(obj['total_count'] == len(rows), 'COMPLETE_API_PAGE')
            need(all(r['head_sha'] == sha and r['status'] == 'completed' and r['conclusion'] == 'success'
                     for r in rows), 'EXACT_HEAD_CI_SUCCESS')
        need(len(post_runs['workflow_runs']) >= minimum_runs, 'POSTMERGE_RUN_COVERAGE')
        need(name != 'HDblast' or len(pre_runs['workflow_runs']) == len(pre_checks['check_runs']) == 9,
             'PREMERGE_RESEARCH_COVERAGE')
        for phase in ('PREMERGE', 'POSTMERGE'):
            status = read(name+'-'+phase+'_STATUS.json')
            need(status['total_count'] == 0 or status['state'] == 'success', 'LEGACY_STATUS')
        final.append({'repository': 'maldonado-research/'+name, 'pull_request': number,
                      'merge_commit': merge, 'ordered_parents': [base, head], 'tree': tree,
                      'premerge_successful_runs': len(pre_runs['workflow_runs']),
                      'premerge_successful_checks': len(pre_checks['check_runs']),
                      'postmerge_successful_runs': len(post_runs['workflow_runs']),
                      'postmerge_successful_checks': len(post_checks['check_runs'])})
    source, live = read('FINAL_SOURCE_READBACK.json'), read('FINAL_LIVE_READBACK.json')
    need(source['status'] == live['status'] == 'PASS', 'PUBLIC_READBACK')
    for observed in (source, live):
        need(observed['research_commit'] == specs[0][1] and observed['shared_commit'] == specs[1][1],
             'READBACK_MERGE_IDENTITIES')
    need(source['matched'] == source['immutable_source_targets'] == 46 and live['matched'] == 3,
         'PUBLIC_READBACK_COVERAGE')
    need(all(r['matches'] for r in source['targets']+live['targets']), 'ALL_READBACK_TARGETS')
    tree_review = read('STAGED_TREE_REVIEW.json')
    need(all(not r['deleted'] for r in tree_review), 'PRESERVATION_NO_DELETION')
    need([r['exact_old_entries_preserved'] for r in tree_review] == [11733, 18], 'OLD_TREE_PRESERVATION')
    native = read('PREMERGE_NATIVE_LOG_REPLAY.json')
    native_post = read('POSTMERGE_NATIVE_LOG_REPLAY.json')
    initial = read('INITIAL_PACKAGED_REPRODUCTION.json')
    need(native['golden_runs'] == initial['runs'] and len(initial['runs']) == 18, 'RUNNER_LOCAL_PARITY')
    need(native_post['golden_runs'] == initial['runs'] and native_post['head_sha'] == specs[0][1],
         'POSTMERGE_RUNNER_LOCAL_PARITY')
    need(native['completion']['corruption_controls'] == len(initial['corruption_controls_rejected']) == 7,
         'INTEGRITY_CONTROL_COVERAGE')
    result = {'schema': 'HDBLAST_SHELL_PRECISION_FINAL_PUBLICATION_VERIFICATION_V1',
              'status': 'PASS', 'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'publications': final, 'immutable_source_targets_matched': 46,
              'complete_live_HTML_pages_matched': 3, 'source_files_in_capsule': 35,
              'local_and_GitHub_golden_runs': 18, 'source_corruption_controls_rejected': 7,
              'old_tree_entries_preserved': [11733, 18],
              'claims': ['Conditional total N0 and J0 sample signs at x=9/8 and x=5/4.',
                         'Conditional anchored proper-time volume second parameter derivative negative at x=9/8.',
                         'All retained parent premises and smooth kinematical completion assumptions remain explicit.'],
              'limits': ['No full angular second-order physical solution, finite-amplitude error range, nonlinear family or cosmological evidence.',
                         'Internal review and computational error bounds are distinct from external peer review and formalization.',
                         'CI artifact storage was inaccessible to the runtime; exact result hashes and completion were read from the GitHub job log.',
                         'No new Zenodo version or journal submission.']}
    destination = root/'FINAL_PUBLICATION_VERIFICATION.json'
    need(not destination.exists(), 'FRESH_VERIFICATION_OUTPUT')
    destination.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
                      'source_targets': 46, 'live_pages': 3}))


if __name__ == '__main__':
    main()
