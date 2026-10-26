#!/usr/bin/env python3
"""Validate observed publication, replay, preservation and complete-byte readback."""
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
    identities = read('PUBLICATION_IDENTITIES.json')
    need(identities['phase'] == 'MERGED', 'MERGED_IDENTITIES')
    need(read('FINAL_POSTMERGE_REVIEW.json')['status'] == 'PASS', 'FINAL_INDEPENDENT_REVIEW')
    initial = read('INITIAL_PACKAGED_REPRODUCTION.json')
    expected_goldens = [{k:r[k] for k in ('actor','golden','sha256','byte_identical')}
                        for r in initial['golden_matches']]
    need(initial['status'] == 'PASS' and len(initial['direct_runs']) == 24
         and len(expected_goldens) == 35 and len(initial['source_corruption_controls']) == 7,
         'FRESH_PACKAGED_REPLAY')
    need(initial['complete_eight_upstream_roots_before'] == initial['complete_eight_upstream_roots_after']
         and len(initial['complete_eight_upstream_roots_after']) == 8, 'EIGHT_UPSTREAM_CLOSURES')
    final = []
    for item in identities['publications']:
        name, merge, head = item['name'], item['merge'], item['head']
        commit, pr = read(name+'-MERGE_COMMIT.json'), read(name+'-MERGED_PR.json')
        need(commit['sha'] == merge and commit['tree']['sha'] == item['tree'], 'CANONICAL_MERGE_TREE')
        need([p['sha'] for p in commit['parents']] == [item['base'],head], 'ORDERED_MERGE_PARENTS')
        need(pr['merged'] and pr['merge_commit_sha'] == merge and pr['head']['sha'] == head,
             'MERGED_PR_IDENTITY')
        coverage = {}
        for phase, sha in [('PREMERGE',head),('POSTMERGE',merge)]:
            runs, checks = read(name+'-'+phase+'_RUNS.json'), read(name+'-'+phase+'_CHECKS.json')
            for obj, field in [(runs,'workflow_runs'),(checks,'check_runs')]:
                need(obj['total_count'] == len(obj[field]), 'COMPLETE_API_PAGE')
                need(all(r['head_sha'] == sha and r['status'] == 'completed' and r['conclusion'] == 'success'
                         for r in obj[field]), 'EXACT_COMMIT_CI_SUCCESS')
            status = read(name+'-'+phase+'_STATUS.json')
            need(status['total_count'] == 0 or status['state'] == 'success', 'LEGACY_STATUS')
            if phase == 'PREMERGE':
                need(len(runs['workflow_runs']) == item['expected_premerge_runs'], 'PREMERGE_RUN_COVERAGE')
            else:
                need(len(runs['workflow_runs']) >= item['minimum_postmerge_runs'], 'POSTMERGE_RUN_COVERAGE')
            coverage[phase.lower()] = {'successful_runs':len(runs['workflow_runs']),
                                       'successful_checks':len(checks['check_runs'])}
        final.append({**item, **coverage})
    research_merge, shared_merge = [r['merge'] for r in identities['publications']]
    for phase, sha in [('PREMERGE',identities['publications'][0]['head']),('POSTMERGE',research_merge)]:
        native = read(phase+'_NATIVE_LOG_REPLAY.json')
        need(native['status'] == 'PASS' and native['head_sha'] == sha and native['goldens'] == expected_goldens,
             'RUNNER_LOCAL_GOLDEN_PARITY')
        need(set(native['passes']) == set(r['id'] for r in initial['direct_runs']), 'RUNNER_DIRECT_RUN_PARITY')
        need(len(native['summaries']) == 1 and native['summaries'][0]['status'] == 'PASS', 'RUNNER_PASS_SUMMARY')
    source, live = read('FINAL_SOURCE_READBACK.json'), read('FINAL_LIVE_READBACK.json')
    for observed in (source,live):
        need(observed['status'] == 'PASS' and observed['research_commit'] == research_merge
             and observed['shared_commit'] == shared_merge and all(r['matches'] for r in observed['targets']),
             'COMPLETE_BYTE_READBACK')
    need(source['matched'] == source['immutable_source_targets'] == 87 and live['matched'] == 3,
         'READBACK_COVERAGE')
    preservation = read('STAGED_TREE_REVIEW.json')
    need([r['exact_old_entries_preserved'] for r in preservation] == [11775,18]
         and all(not r['deleted'] for r in preservation), 'COMPLETE_OLD_TREE_PRESERVATION')
    need([r['staged_tree'] for r in preservation] == [r['tree'] for r in final], 'REVIEWED_TREES_PUBLISHED')
    refs = read('FINAL_MAIN_REFS.json')
    need(refs['research_main'] == research_merge and refs['shared_main'] == shared_merge, 'FINAL_MAIN_REFS')
    result = {'schema':'HDBLAST_FULL_ANGULAR_FINAL_PUBLICATION_VERIFICATION_V1','status':'PASS',
        'observed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'publications':final,'immutable_source_files_matched':87,'complete_live_HTML_pages_matched':3,
        'scientific_capsule_files':77,'direct_replay_runs':24,'scientific_golden_comparisons':35,
        'source_integrity_controls_rejected':7,'additional_nested_negative_control_runs':40,
        'upstream_complete_capsules_authenticated':8,'old_tree_entries_preserved_exactly':[11775,18],
        'local_premerge_postmerge_scientific_hashes_equal':True,
        'claims':['Continuous full23 formal second response approximation meets its stated C1 comparison-norm target.',
                  'Complete integral residual and residual/resolvent error bound are supplied.',
                  'Actual total equal TT coefficients are strictly negative at one exterior point before shell contact.'],
        'limits':['Retained conditional background, master, gauge, coefficient-existence and subsidiary premises apply.',
                  'No pointwise C2 Einstein residual, finite-amplitude nonlinear realization or observational origin evidence.',
                  'Registered TT coefficient is not a gravitational-wave observable, energy or detector signal.',
                  'Internal AI-assisted review is distinct from external specialist review and proof-assistant formalization.',
                  'CI scientific hashes and completion were checked in decoded GitHub job logs; full artifact storage was not downloaded.',
                  'No new Zenodo version or journal submission.']}
    destination = root/'FINAL_PUBLICATION_VERIFICATION.json'
    need(not destination.exists(), 'FRESH_VERIFICATION_OUTPUT')
    destination.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','source_files':87,'live_pages':3,
                      'sha256':hashlib.sha256(destination.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    main()
