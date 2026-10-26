#!/usr/bin/env python3
"""Compare decoded exact-commit CI job evidence with the fresh local replay."""
import argparse
import hashlib
import json
from pathlib import Path


def need(ok, label):
    if not ok:
        raise ValueError(label)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--log', type=Path, required=True)
    p.add_argument('--head', required=True)
    p.add_argument('--run', type=int, required=True)
    p.add_argument('--job', type=int, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    root = Path(__file__).resolve().parent
    local = json.loads((root/'INITIAL_PACKAGED_REPRODUCTION.json').read_text())
    pins = json.loads((root/'PACKAGE_PINS.json').read_text())
    goldens, passes, completions, summaries = [], [], [], []
    for line in args.log.read_text().splitlines():
        at = line.find('{')
        if at < 0:
            continue
        try:
            row = json.loads(line[at:])
        except ValueError:
            continue
        if row.get('phase') == 'golden':
            goldens.append({key:row[key] for key in ('actor','golden','sha256','byte_identical')})
        if row.get('phase') == 'pass':
            passes.append({'id':row['run'],'expected_rejection':row.get('expected_rejection')})
        if row.get('status') == 'PASS' and 'final_pins' in row:
            completions.append(row)
        if row.get('phase') == 'ci-receipt':
            summaries.append(row)
    goldens.sort(key=lambda x:(x['actor'],x['golden']))
    expected = [{key:row[key] for key in ('actor','golden','sha256','byte_identical')}
                for row in local['golden_matches']]
    need(len(goldens) == 17 and goldens == expected, 'ALL_SCIENTIFIC_GOLDEN_HASHES')
    need(len(passes) == 20 and len({row['id'] for row in passes}) == 20, 'EXACT_RUN_COUNT')
    need(sorted(passes,key=lambda r:r['id']) ==
         [{'id':r['id'],'expected_rejection':r.get('expected_rejection')} for r in local['direct_runs']],
         'ALL_POSITIVE_AND_EXPECTED_NEGATIVE_RUNS')
    need(len(completions) == len(summaries) == 1, 'ONE_NATIVE_COMPLETION_AND_RECEIPT')
    for row in completions+summaries:
        need(row['status'] == 'PASS' and row['direct_runs'] == 20
             and row['positive_direct_runs'] == 16 and row['expected_nonzero_direct_runs'] == 4
             and row['scientific_golden_matches'] == 17 and row['source_corruption_controls'] == 7
             and row['independent_math_negative_controls'] == 4, 'NATIVE_SUMMARY_COUNTS')
    completion, summary = completions[0], summaries[0]
    need(completion['complete_upstream_closures'] == summary['upstream_closures'] == 9
         and summary['capsule_source_files'] == 46 and summary['indexed_programs'] == 11,
         'CLOSURE_AND_PROGRAM_COVERAGE')
    for field in ('manifest_sha256','archive_sha256','reproduce_sha256'):
        need(summary[field] == pins[field], 'EXTERNAL_PIN:'+field)
    need(completion['final_pins']['capsule_manifest_sha256'] == pins['manifest_sha256'], 'COMPLETION_MANIFEST')
    result = {'schema':'HDBLAST_EXTERIOR_TENSOR_JETS_NATIVE_LOG_PARITY_V1','status':'PASS',
        'repository':'maldonado-research/HDblast','head_sha':args.head,'run_id':args.run,'job_id':args.job,
        'source':'GitHub connector decoded completed workflow job log',
        'native_log_sha256':hashlib.sha256(args.log.read_bytes()).hexdigest(),
        'goldens':goldens,'passes':passes,'completion':completion,'summaries':summaries,
        'all_golden_hashes_match_local':True,'all_direct_run_ids_and_expected_rejections_match_local':True,
        'scope':'Scientific byte hashes and checked completion read from native job log; full artifact not downloaded. Host paths/log metadata make execution receipt hashes host-dependent.'}
    need(not args.output.exists(), 'FRESH_NATIVE_PARITY_RECEIPT')
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','goldens':17,'positive_runs':16,'expected_negative_runs':4}))


if __name__ == '__main__':
    main()
