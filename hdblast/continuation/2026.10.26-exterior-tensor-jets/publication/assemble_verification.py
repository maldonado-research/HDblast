#!/usr/bin/env python3
"""Validate retained exact-commit publication, replay and byte-readback evidence."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path


def need(ok, label):
    if not ok:
        raise ValueError(label)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-existing', action='store_true',
                        help='Validate an existing receipt without changing any file.')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    def read(name):
        return json.loads((root/name).read_text())
    identities = read('PUBLICATION_IDENTITIES.json')
    need(identities['phase'] == 'MERGED', 'MERGED_IDENTITIES')
    need(read('FINAL_POSTMERGE_REVIEW.json')['status'] == 'PASS', 'FINAL_INDEPENDENT_REVIEW')
    initial = read('INITIAL_PACKAGED_REPRODUCTION.json')
    positives = [r for r in initial['direct_runs'] if r['exit_code'] == 0]
    negatives = [r for r in initial['direct_runs'] if r['exit_code'] != 0]
    need(initial['status'] == 'PASS' and len(positives) == 16 and len(negatives) == 4
         and all(r['reason_verified'] and r['output_absent'] for r in negatives), 'EXPECTED_REPLAY_RESULTS')
    expected_goldens = [{k:r[k] for k in ('actor','golden','sha256','byte_identical')}
                        for r in initial['golden_matches']]
    need(len(expected_goldens) == 17 and len(initial['source_corruption_controls']) == 7,
         'REFERENCE_AND_SOURCE_CONTROL_COVERAGE')
    need(initial['complete_nine_upstream_roots_before'] == initial['complete_nine_upstream_roots_after']
         and len(initial['complete_nine_upstream_roots_after']) == 9, 'NINE_COMPLETE_UPSTREAM_CLOSURES')
    need(initial['complete_new_capsule_before'] == initial['complete_new_capsule_after']
         and initial['complete_new_capsule_after']['files_including_manifest'] == 46,
         'FROZEN_NEW_SOURCE_CLOSURE')
    publications = []
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
                need(obj['total_count'] == len(obj[field]), 'COMPLETE_CI_API_PAGE')
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
        publications.append({**item,**coverage})
    research_merge, shared_merge = [r['merge'] for r in publications]
    expected_passes = [{'id':r['id'],'expected_rejection':r.get('expected_rejection')}
                       for r in initial['direct_runs']]
    for phase, sha in [('PREMERGE',publications[0]['head']),('POSTMERGE',research_merge)]:
        native = read(phase+'_NATIVE_LOG_REPLAY.json')
        need(native['status'] == 'PASS' and native['head_sha'] == sha and native['goldens'] == expected_goldens,
             'RUNNER_LOCAL_REFERENCE_PARITY')
        need(sorted(native['passes'],key=lambda r:r['id']) == expected_passes,
             'RUNNER_POSITIVE_AND_EXPECTED_NEGATIVE_PARITY')
        need(len(native['summaries']) == 1 and native['summaries'][0]['status'] == 'PASS', 'NATIVE_RECEIPT_SUMMARY')
    source, live = read('FINAL_SOURCE_READBACK.json'), read('FINAL_LIVE_READBACK.json')
    for observed in (source,live):
        need(observed['status'] == 'PASS' and observed['research_commit'] == research_merge
             and observed['shared_commit'] == shared_merge and all(r['matches'] for r in observed['targets']),
             'COMPLETE_BYTE_READBACK')
    need(source['matched'] == source['immutable_source_targets'] == 56 and live['matched'] == 3,
         'SOURCE_AND_COMPLETE_LIVE_READBACK_COVERAGE')
    preservation = read('STAGED_TREE_REVIEW.json')
    need([r['exact_old_entries_preserved'] for r in preservation] == [11858,18]
         and all(not r['deleted'] for r in preservation), 'FULL_OLD_TREE_PRESERVATION')
    need([r['staged_tree'] for r in preservation] == [r['tree'] for r in publications], 'REVIEWED_TREES_PUBLISHED')
    refs = read('FINAL_MAIN_REFS.json')
    need(refs['research_main'] == research_merge and refs['shared_main'] == shared_merge, 'FINAL_MAIN_REFS')
    result = {'schema':'HDBLAST_EXTERIOR_TENSOR_JETS_FINAL_PUBLICATION_VERIFICATION_V1','status':'PASS',
        'observed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'publications':publications,
        'immutable_source_files_matched':56,'complete_live_HTML_pages_matched':3,'scientific_capsule_files':46,
        'positive_replay_runs':16,'expected_reason_verified_nonzero_runs':4,'reference_byte_comparisons':17,
        'source_integrity_controls_rejected':7,'upstream_complete_capsules_authenticated':9,
        'old_tree_entries_preserved_exactly':[11858,18],'local_premerge_postmerge_scientific_hashes_equal':True,
        'claims':['Actual total equal tensor coefficients and all six coordinate jet entries enclosed at one exterior point.',
                  'Coefficient and scaled time derivative strictly negative; value enclosure is a subset of the old interval.',
                  'Complete source radial derivative needs only the authenticated C3 inputs in this exact tensor row.'],
        'limits':['All three second-derivative signs remain inconclusive at the declared maximum budget.',
                  'Retained smooth background/master/gauge, differentiated-error, formal coefficient and subsidiary premises apply.',
                  'Ordinary coefficient derivatives are not covariant curvature or a gauge-invariant gravitational-wave observable.',
                  'No full-domain C2 Einstein residual, finite-amplitude nonlinear realization or Big Bang-origin evidence.',
                  'Internal AI-assisted review differs from external specialist review and proof-assistant formalization.',
                  'CI reference hashes and completion were checked in decoded GitHub job logs; full artifacts were not downloaded.',
                  'No new Zenodo version or journal submission.']}
    destination = root/'FINAL_PUBLICATION_VERIFICATION.json'
    if args.check_existing:
        retained = json.loads(destination.read_text())
        need(isinstance(retained.get('observed_UTC'), str), 'RETAINED_OBSERVATION_TIME')
        retained.pop('observed_UTC')
        comparison = dict(result)
        comparison.pop('observed_UTC')
        need(retained == comparison, 'EXISTING_RECEIPT_MATCHES_RECOMPUTED_EVIDENCE')
        print(json.dumps({'status':'PASS','mode':'read-only receipt verification',
                          'sha256':hashlib.sha256(destination.read_bytes()).hexdigest()}))
        return
    need(not destination.exists(), 'FRESH_PUBLICATION_RECEIPT')
    destination.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','source_files':56,'complete_live_pages':3,
                      'sha256':hashlib.sha256(destination.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    main()
