#!/usr/bin/env python3
"""Compare completed native job-log science with the fresh local ZIP replay."""
import argparse
import hashlib
import json
from pathlib import Path


def need(value, label):
    if not value:
        raise ValueError(label)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log', type=Path, required=True)
    parser.add_argument('--head', required=True)
    parser.add_argument('--run', type=int, required=True)
    parser.add_argument('--job', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(), 'FRESH_NATIVE_PROOF')
    root = Path(__file__).resolve().parent
    local = json.loads((root/'INITIAL_PACKAGED_REPRODUCTION.json').read_text())
    pins = json.loads((root/'PACKAGE_PINS.json').read_text())
    goldens,passes,integrity,completions,summaries = [],[],[],[],[]
    for line in args.log.read_text().splitlines():
        start = line.find('{')
        if start < 0:
            continue
        try:
            row = json.loads(line[start:])
        except ValueError:
            continue
        if not isinstance(row,dict):
            continue
        if row.get('phase') == 'golden':
            goldens.append({key:row[key] for key in ('run','output','reference','sha256','bytes','byte_identical')})
        if row.get('phase') == 'pass':
            passes.append({'id':row['run'],'expected_rejection':row.get('expected_rejection')})
        if row.get('phase') == 'integrity-control':
            need(row.get('pass') is True, 'NATIVE_INTEGRITY_REJECTION')
            integrity.append({'id':row.get('run',row.get('id')),
                              'reason':row.get('reason',row.get('expected_rejection'))})
        if row.get('status') == 'PASS' and 'positive_runs' in row:
            completions.append(row)
        if row.get('phase') == 'ci-receipt':
            summaries.append(row)
    need(local['status'] == 'PASS' and len(goldens) == 28 and goldens == local['golden_matches'],
         'ALL_28_NATIVE_REFERENCE_HASHES_AND_BYTES')
    expected_passes = [{'id':r['id'],'expected_rejection':r.get('expected_rejection')}
                       for r in local['direct_runs']]
    need(len(passes) == 32 and passes == expected_passes, 'ALL_32_RUN_IDS_AND_6_REASONS')
    expected_integrity = [{'id':r['id'],'reason':r['expected_rejection']}
                          for r in local['source_integrity_controls']+local['archive_integrity_controls']]
    need(len(integrity) == 19 and integrity == expected_integrity, 'ALL_19_NATIVE_INTEGRITY_CONTROLS')
    need(len(completions) == len(summaries) == 1, 'ONE_NATIVE_COMPLETION_AND_CI_RECEIPT')
    end,ci = completions[0],summaries[0]
    need({k:end[k] for k in ('positive_runs','negative_runs','golden_matches')} == local['counts']
         == {'positive_runs':26,'negative_runs':6,'golden_matches':28}, 'NATIVE_COMPLETION_COUNTS')
    need(end['source_integrity_controls'] == ci['source_integrity_controls'] == 16
         and end['archive_integrity_controls'] == ci['archive_integrity_controls'] == 3
         and end['upstream_capsules'] == ci['upstream_closures'] == 10, 'NATIVE_CLOSURE_AND_CONTROL_COUNTS')
    need(ci['status'] == 'PASS' and ci['direct_runs'] == 32 and ci['positive_direct_runs'] == 26
         and ci['expected_nonzero_direct_runs'] == 6 and ci['scientific_golden_matches'] == 28
         and ci['capsule_source_files'] == 60 and ci['indexed_programs'] == 18
         and ci['output_files'] == 88, 'NATIVE_INDEPENDENT_RECEIPT_VERIFICATION')
    for key in ('manifest_sha256','archive_sha256'):
        need(end[key] == ci[key] == pins[key], 'NATIVE_EXTERNAL_PIN:'+key)
    need(ci['reproduce_sha256'] == pins['reproduce_sha256']
         and ci['replay_plan_sha256'] == local['replay_plan_sha256']
         and end['receipt_sha256'] == ci['receipt_sha256'], 'NATIVE_ENTRY_PLAN_AND_RECEIPT_PINS')
    report = {'schema':'HDBLAST_SCREEN_MAGNETIC_NATIVE_LOG_PARITY_V1','status':'PASS',
      'repository':'maldonado-research/HDblast','head_sha':args.head,'run_id':args.run,'job_id':args.job,
      'source':'GitHub connector decoded completed job log; job/commit metadata checked separately.',
      'native_log_sha256':hashlib.sha256(args.log.read_bytes()).hexdigest(),
      'goldens':goldens,'passes':passes,'integrity_controls':integrity,
      'completion':end,'summaries':summaries,
      'all_28_reference_hashes_equal_local':True,'all_32_run_ids_and_6_expected_reasons_equal_local':True,
      'all_19_integrity_rejections_equal_local':True,
      'scope':'Scientific byte/hash parity and independent completion summary; raw job response and logs excluded from public evidence. Full workflow artifact not downloaded.'}
    with args.output.open('x') as handle:
        json.dump(report,handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps({'status':'PASS','goldens':28,'run_ids':32,'expected_nonzero':6,'integrity_controls':19}))


if __name__ == '__main__':
    main()
