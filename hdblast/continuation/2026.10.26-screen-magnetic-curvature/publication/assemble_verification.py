#!/usr/bin/env python3
"""Reconcile exact-commit CI, frozen science, preservation and public readbacks."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path


def need(value, label):
    if not value:
        raise ValueError(label)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-existing', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    def read(name):
        return json.loads((root/name).read_text())
    identities = read('PUBLICATION_IDENTITIES.json')
    need(identities['phase'] == 'MERGED', 'CANONICAL_MERGED_IDENTITIES')
    need(read('FINAL_PREMERGE_REVIEW.json')['status'] == 'PASS'
         and read('FINAL_POSTMERGE_REVIEW.json')['status'] == 'PASS', 'INDEPENDENT_FINAL_PUBLICATION_REVIEWS')
    local,pins = read('INITIAL_PACKAGED_REPRODUCTION.json'),read('PACKAGE_PINS.json')
    need(local['status'] == 'PASS' and local['counts'] == {'positive_runs':26,'negative_runs':6,'golden_matches':28},
         'FRESH_RELOCATED_REPLAY')
    need(local['complete_capsule_before'] == local['complete_capsule_after']
         and local['complete_capsule_after']['files_including_manifest'] == 60
         and local['complete_capsule_after']['manifest_sha256'] == pins['manifest_sha256']
         and local['source_unchanged_after_replay'] is True, 'FROZEN_60_FILE_CAPSULE')
    need(local['complete_ten_upstream_roots_before'] == local['complete_ten_upstream_roots_after']
         and len(local['complete_ten_upstream_roots_after']) == 10
         and local['all_authentication_before_any_project_execution'] is True, 'TEN_COMPLETE_UPSTREAM_ROOTS')
    need(local['archive_before'] == local['archive_after'] and local['archive_after']['sha256'] == pins['archive_sha256']
         and local['archive_after']['members'] == 60 and local['entry_sha256'] == pins['reproduce_sha256'],
         'EXTERNAL_ARCHIVE_AND_ENTRY_PINS')
    need(len(local['source_integrity_controls']) == 16 and len(local['archive_integrity_controls']) == 3
         and all(r['reason_verified'] for r in local['source_integrity_controls']+local['archive_integrity_controls']),
         'REASON_SPECIFIC_INTEGRITY_CONTROLS')
    runs = local['direct_runs']
    need(len(runs) == len({r['id'] for r in runs}) == 32
         and sum(r['exit_code'] == 0 for r in runs) == 26
         and all(r['reason_verified'] and r['output_absent'] for r in runs if r['exit_code'] != 0)
         and len(local['golden_matches']) == 28 and len(local['program_index']) == 18,
         'EXACT_PROGRAM_RUN_REFERENCE_COVERAGE')
    publications = []
    for item in identities['publications']:
        name,merge,head = item['name'],item['merge'],item['head']
        commit,pr = read(name+'-MERGE_COMMIT.json'),read(name+'-MERGED_PR.json')
        need(commit['sha'] == merge and commit['tree']['sha'] == item['tree']
             and [p['sha'] for p in commit['parents']] == [item['base'],head], 'MERGE_TREE_AND_ORDERED_PARENTS')
        need(pr['merged'] is True and pr['merge_commit_sha'] == merge and pr['head']['sha'] == head,
             'MERGED_PR_IDENTITY')
        coverage = {}
        for phase,sha in [('PREMERGE',head),('POSTMERGE',merge)]:
            observed,checks = read(name+'-'+phase+'_RUNS.json'),read(name+'-'+phase+'_CHECKS.json')
            for obj,field in [(observed,'workflow_runs'),(checks,'check_runs')]:
                need(obj['total_count'] == len(obj[field]), 'COMPLETE_CI_API_PAGE')
                need(all(r['head_sha'] == sha and r['status'] == 'completed' and r['conclusion'] == 'success'
                         for r in obj[field]), 'EXACT_COMMIT_ALL_CI_SUCCESS')
            status = read(name+'-'+phase+'_STATUS.json')
            need(status['total_count'] == 0 or status['state'] == 'success', 'LEGACY_STATUS_SUCCESS')
            if phase == 'PREMERGE':
                need(len(observed['workflow_runs']) == item['expected_premerge_runs'], 'PREMERGE_TRIGGER_COVERAGE')
            else:
                need(len(observed['workflow_runs']) >= item['minimum_postmerge_runs'], 'POSTMERGE_TRIGGER_COVERAGE')
            coverage[phase.lower()] = {'successful_runs':len(observed['workflow_runs']),
                                       'successful_checks':len(checks['check_runs'])}
        publications.append({**item,**coverage})
    expected_passes = [{'id':r['id'],'expected_rejection':r.get('expected_rejection')} for r in runs]
    for phase,sha in [('PREMERGE',publications[0]['head']),('POSTMERGE',publications[0]['merge'])]:
        native = read(phase+'_NATIVE_LOG_REPLAY.json')
        need(native['status'] == 'PASS' and native['head_sha'] == sha
             and native['goldens'] == local['golden_matches'] and native['passes'] == expected_passes
             and len(native['integrity_controls']) == 19, 'LOCAL_AND_NATIVE_SCIENTIFIC_PARITY')
        need(native['summaries'][0]['status'] == 'PASS'
             and native['summaries'][0]['manifest_sha256'] == pins['manifest_sha256']
             and native['summaries'][0]['archive_sha256'] == pins['archive_sha256']
             and native['summaries'][0]['reproduce_sha256'] == pins['reproduce_sha256'], 'NATIVE_PINNED_RECEIPT')
    research_merge,shared_merge = [r['merge'] for r in publications]
    source,live = read('FINAL_SOURCE_READBACK.json'),read('FINAL_LIVE_READBACK.json')
    for result in (source,live):
        need(result['status'] == 'PASS' and result['research_commit'] == research_merge
             and result['shared_commit'] == shared_merge and all(r['matches'] for r in result['targets']),
             'PUBLISHED_COMPLETE_BYTE_READBACK')
    need(source['matched'] == source['immutable_source_targets'] == 69 and live['matched'] == 3,
         'ALL_69_CHANGED_FILES_AND_3_COMPLETE_HTML_PAGES')
    preservation = read('STAGED_TREE_REVIEW.json')
    need([r['exact_old_entries_preserved'] for r in preservation] == [11910,18]
         and all(not r['deleted'] for r in preservation)
         and [r['staged_tree'] for r in preservation] == [r['tree'] for r in publications],
         'FULL_OLD_TREE_PRESERVATION_AND_IDENTICAL_MERGE_TREES')
    refs = read('FINAL_MAIN_REFS.json')
    need([refs['research_main'],refs['shared_main']] == [research_merge,shared_merge], 'FINAL_MAIN_IDENTITIES')
    result = {'schema':'HDBLAST_SCREEN_MAGNETIC_FINAL_PUBLICATION_VERIFICATION_V1','status':'PASS',
      'observed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(), 'publications':publications,
      'immutable_source_files_matched':69,'complete_live_HTML_pages_matched':3,
      'scientific_capsule_files':60,'python_programs':18,'positive_replay_runs':26,
      'reason_verified_expected_nonzero_runs':6,'reference_byte_matches':28,
      'source_integrity_rejections':16,'archive_integrity_rejections':3,'complete_upstream_capsules':10,
      'unchanged_prior_leaf_entries':[11910,18], 'local_premerge_postmerge_scientific_hashes_equal':True,
      'package_pins':pins,
      'claims':['Complete registered-observer selected magnetic Weyl second-variation coefficient excludes zero on one exterior angular sphere.',
                'Full scalar quadratic, observer, radial-normal, connection and epsilon corrections retained.',
                'Exact angular closure and independent generic five-dimensional curvature fixtures verify the compact coefficient.'],
      'limits':['Registered t-normal/spatial-y-normal observers; no geometrically anchored clock, point or second-order gauge-invariance certificate.',
                'Selected three-screen projection, not full five-dimensional magnetic Weyl reconstruction; two great circles are angular zeros.',
                'Retained smooth actual first-field and formal-response/subsidiary premises remain conditional.',
                'Full compatible second-order Einstein–matter boundary solution and finite-amplitude nonlinear realization remain unfinished.',
                'No radiation, detector force, energy or cosmological-origin conclusion.',
                'Internal AI-assisted review differs from independent specialist review and proof-assistant formalization.',
                'Native decoded job evidence was checked; raw logs and full artifacts are not part of published proof.',
                'No new Zenodo version or journal submission.']}
    destination = root/'FINAL_PUBLICATION_VERIFICATION.json'
    if args.check_existing:
        observed = json.loads(destination.read_text());need(isinstance(observed.pop('observed_UTC',None),str),'RECEIPT_TIME')
        comparison = dict(result);comparison.pop('observed_UTC')
        need(observed == comparison, 'RECOMPUTED_PUBLICATION_RECEIPT')
    else:
        need(not destination.exists(), 'FRESH_PUBLICATION_RECEIPT')
        with destination.open('x') as handle:
            json.dump(result,handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps({'status':'PASS','source_files':69,'complete_live_pages':3,
                      'sha256':hashlib.sha256(destination.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    main()
