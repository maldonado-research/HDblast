#!/usr/bin/env python3
"""Pure byte/hash and synthetic receipt review; never decodes saved arrays."""
import hashlib
import json
from pathlib import Path

BASE = Path('/workspace/hdblast-research-work')
CPP = Path('/workspace/HDblast/research/HDBLAST_CHECKPOINT_20261002_SOURCE_FREE_LEDGER')
OUT = BASE / 'ledger-review-20261002'

def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            digest.update(block)
    return digest.hexdigest()

def read(path):
    return json.loads(path.read_text())

checks = []
def check(name, condition):
    checks.append({'name': name, 'pass': bool(condition)})

ind = BASE / 'ledger-independent-20261002'
prep = read(ind / 'PREPARATION_CHECKS.json')
for name, pin in prep['source_pins'].items():
    check('independent_source/' + name, sha(ind / name) == pin)
for name, pin in prep['evidence'].items():
    check('independent_receipt/' + name, sha(ind / name) == pin)
manifest = read(ind / 'MANIFEST.json')
for name, pin in manifest['files'].items():
    check('independent_manifest/' + name, sha(ind / name) == pin)
for name in ('diagnostic_independent.py', 'ledger_core.py', 'profile_kernel.py', 'exact_binary.py',
             'CONFIGURATION.json', 'MANIFEST.json', 'DESIGN.md', 'OUTPUT_SCHEMA.json'):
    check('canonical_independent/' + name, sha(CPP / 'independent' / name) == sha(ind / name))
check('independent_preparation_zero_retained_evaluations', prep['retained_physical_evaluations'] == 0
      and prep['retained_physical_payloads_decoded'] is False)
timing = read(ind / 'SYNTHETIC_PIPELINE_TIMING_FINAL.json')
check('independent_resource_projection_accurately_labeled',
      timing['projection_is_not_full_retained_data_resource_proof'] is True
      and prep['resource_projection_is_full_retained_data_proof'] is False)

primary = BASE / 'ledger-primary-20261002'
full = read(primary / 'proof/FINAL_FULL_SHAPE_NORMAL.json')
optimized = read(primary / 'proof/FINAL_PURE_OPTIMIZED.json')
for mode, receipt in (('normal', full), ('optimized', optimized)):
    check('primary_source_bound_' + mode, receipt['source_sha256'] == sha(primary / 'ledger_primary.py')
          and receipt['test_sha256'] == sha(primary / 'synthetic_preflight.py'))
    pure = receipt['pure_checks']
    check('primary_pure_checks_' + mode, pure['checks_count'] == 83 and pure['rejected_mutations'] == 6
          and pure['physical_input_files_opened'] == 0 and pure['physical_source_evaluations'] == 0)
check('canonical_primary', sha(CPP / 'code/diagnostic_primary.py') == sha(primary / 'ledger_primary.py'))
check('primary_full_shape_synthetic_only', full['physical_saved_data_loaded'] is False
      and full['full_shape_run']['fixed_cases'] == 12)
check('primary_complete_two_precision_synthetic_budget',
      full['full_shape_run']['resources']['seconds'] <= 900 and full['peak_rss_kib'] <= 262144)

tooling = BASE / 'ledger-replay-20261002'
ready = read(tooling / 'evidence/TOOLING_READINESS.json')
for name, pin in ready['files'].items():
    check('tooling_source/' + name, sha(tooling / name) == pin)
    if name.startswith('code/') or name in ('REPRODUCTION.md', 'requirements-replay.txt'):
        check('canonical_tooling/' + name, sha(CPP / name) == pin)
check('tooling_guards_and_zero_physical_evaluation', ready['normal_guard_tests'] == 27
      and ready['optimized_guard_tests'] == 27 and ready['physical_evaluations'] == 0
      and ready['kernel_process_creation_limit'] == 0)
import_receipt = read(tooling / 'evidence/import-resource-check/pinned_imports.exit.json')
check('actual_pinned_imports_under_single_process_limits', import_receipt['status'] == 'PASS_EXECUTION'
      and import_receipt['process_creation_limit'] == 0 and import_receipt['peak_rss_kib'] <= 262144
      and import_receipt['elapsed_seconds'] <= 900 and import_receipt['physical_route'] is False)
check('shared_integrity_helper_binding', ready['helper_sha256'] == prep['shared_integrity_helper_sha256'])
report = {
    'schema_version': 1,
    'status': 'PASS_SOURCE_AND_SYNTHETIC_EVIDENCE' if all(x['pass'] for x in checks) else 'FAIL_SOURCE_AND_SYNTHETIC_EVIDENCE',
    'checks_count': len(checks), 'checks': checks,
    'physical_arrays_decoded': False, 'physical_calculations_executed': 0,
    'validator_review': 'Separate pure synthetic mutation suite; pending final source-bound receipt.',
    'source_pins': {'primary': full['source_sha256'],
                    'independent_entry': prep['source_pins']['diagnostic_independent.py'],
                    'independent_core': prep['source_pins']['ledger_core.py'],
                    'shared_integrity': ready['helper_sha256'],
                    'replay_driver': ready['files']['code/replay_ledger.py']},
    'resources': {'primary_full_shape_synthetic_seconds': full['full_shape_run']['resources']['seconds'],
                  'primary_full_shape_synthetic_peak_rss_kib': full['peak_rss_kib'],
                  'independent_projected_full_two_precision_seconds': timing['projected_both_precision_seconds_including_node_comparison'],
                  'independent_synthetic_peak_rss_kib': timing['peak_rss_kib'],
                  'independent_projection_is_actual_resource_proof': False},
}
(OUT / 'SOURCE_EVIDENCE_REVIEW.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
print(json.dumps({k: report[k] for k in ('status', 'checks_count', 'source_pins', 'resources')}, sort_keys=True))
if report['status'].startswith('FAIL'):
    print(json.dumps([x for x in checks if not x['pass']]))
    raise SystemExit(1)
