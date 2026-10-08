#!/usr/bin/env python3
"""Independent complete opaque-byte review immediately before the root seal."""
import ast
import importlib.util
import json
from pathlib import Path
from datetime import datetime, timezone

W = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
OUT = W / 'replay-delivery/package-stage-independent-review'
spec = importlib.util.spec_from_file_location('original_static_byte_reviewer', OUT / 'review_unsealed_stage.py')
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)  # Original reviewer utility only; no package import/main().
P = check.P
STAGE_SHA = '4d24e9484f85fab6c68def9e982a650748e950f01524f7b6409a43f1c3865178'
RESULTS_PIN = {'bytes': 22147, 'sha256': '02ad104b9cb6e69881ce4ea3040fcf674465f64688e1da1dd1e8df5e904109e4'}


def main():
    before = check.opaque_pin(check.STAGE)
    check.require(before['sha256'] == STAGE_SHA, 'final unsealed staging snapshot differs')
    stage = check.load_json(check.STAGE)
    check.require(stage['sealed'] is False and stage['status'] == 'UNSEALED_PORTABLE_PACKAGE_STAGED', 'not expected final unsealed stage')
    previous = check.load_json(OUT / 'OBSERVED_REFRESHED_SNAPSHOT.json')['files']
    files, dirs = check.tree(P)
    check.require(files == stage['files'], 'final stage differs from complete external leaf map')
    check.require(check.summary(files) == {'file_count': 239, 'total_leaf_bytes': 258002288}, 'final unsealed leaf counts differ')
    added = sorted(set(files) - set(previous))
    expected_additions = sorted(check.C_REL + '/' + n for n in ['RESULTS.md', 'NEXT_CALCULATION_CONTRACT.md', 'evidence/RESULTS_PRESENTATION_REVIEW.json'])
    check.require(added == expected_additions and not set(previous) - set(files), 'unexpected final supplemental additions/removals')
    check.require(all(files[n] == pin for n, pin in previous.items()), 'previous approved leaf bytes changed')
    check.require(files[check.C_REL + '/RESULTS.md'] == RESULTS_PIN, 'externally approved RESULTS identity differs')
    check.require(files['replay_checkpoint.py']['sha256'] == check.HELPER_SHA and files['PUBLIC_GO.json']['sha256'] == check.GO_SHA, 'helper/GO identity differs')
    require_absent = [P / 'PAYLOAD_MANIFEST.json', P / 'REPLAY_PINS.json']
    check.require(all(not q.exists() for q in require_absent), 'root final controls already present')
    registration_path = P / check.C_REL / 'FULL_REGISTRATION.json'
    registration_pin = check.opaque_pin(registration_path)
    check.require(registration_pin['sha256'] == check.REG_SHA, 'final source registration differs')
    registration = check.load_json(registration_path)
    check.require(set(registration) == {'schema_version', 'files'} and type(registration['schema_version']) is int and registration['schema_version'] == 1 and len(registration['files']) == 117, '117-file registration schema differs')
    for n, pin in registration['files'].items():
        check.rel_safe(n)
        check.validate_pin(pin)
        check.require(files[check.C_REL + '/' + n] == pin, 'final registered leaf differs: ' + n)
    scientific_copy = check.subset(files, check.C_REL)
    scientific_source, _ = check.tree(check.SOURCE_C)
    check.require(scientific_copy == scientific_source, 'final complete C does not match current repository source')
    check.require(len(scientific_copy) == 161, 'final complete C file count differs')

    inputspec = check.load_json(P / check.C_REL / 'provenance/INPUT_SPEC.json')
    dependencies = {}
    for capsule in inputspec['capsules']:
        dependencies[check.rel_safe(capsule['repository_path'])] = {'bytes': capsule['reader_inputspec']['file_bytes'], 'sha256': capsule['reader_inputspec']['file_sha256']}
        original = capsule['original']
        dependencies[check.rel_safe(original['repository_path'])] = {'bytes': original['bytes'], 'sha256': original['sha256']}
    npz_bytes = sum(pin['bytes'] for pin in dependencies.values())
    check.require(len(dependencies) == 8 and npz_bytes == 212975098, 'eight NPZ dependency pins differ')
    for path_key, sha_key in [('manifest_repository_path', 'manifest_sha256'), ('producer_repository_path', 'producer_sha256')]:
        n = check.rel_safe(inputspec[path_key])
        dependencies[n] = files[n]
        check.require(dependencies[n]['sha256'] == inputspec[sha_key], 'non-array dependency SHA differs')
    audit = inputspec['prior_static_audit']
    dependencies[check.rel_safe(audit['repository_path'])] = {'bytes': audit['bytes'], 'sha256': audit['sha256']}
    check.require(len(dependencies) == 11 and dependencies == stage['input_dependencies_verified'], 'final 11-dependency map differs')
    for n, pin in dependencies.items():
        check.validate_pin(pin)
        check.require(files[n] == pin == check.opaque_pin(check.REPO / n), 'final full opaque dependency/source bytes differ: ' + n)
    modes = {}
    for mode in ['normal', 'optimized']:
        saved = check.subset(files, 'evidence/actual/' + mode)
        source, _ = check.tree(W / ('root/actual-' + mode + '-001'))
        check.require(saved == source, 'final saved mode differs from complete actual root copy: ' + mode)
        modes[mode] = check.summary(saved)
    presentation_path = check.C_REL + '/evidence/RESULTS_PRESENTATION_REVIEW.json'
    presentation = check.load_json(P / presentation_path)
    python_files = sorted(n for n in files if n.endswith('.py'))
    for n in python_files:
        ast.parse((P / n).read_text(encoding='utf-8'), filename=n)
    check.require(check.opaque_pin(check.STAGE) == before, 'final inventory changed during review')
    end, end_dirs = check.tree(P)
    check.require(end == files and end_dirs == dirs, 'final unsealed P changed during review')
    snapshot = {'schema_version': 1, 'status': 'INDEPENDENT_FINAL_UNSEALED_PRE_SEAL_SNAPSHOT', 'sealed': False, 'package_path': str(P), 'files': files, 'directories': dirs, **check.summary(files)}
    check.json_write(OUT / 'OBSERVED_FINAL_PRE_SEAL_SNAPSHOT.json', snapshot)
    report = {
        'schema_version': 1,
        'status': 'PASS_FINAL_UNSEALED_PRE_SEAL_PACKAGE_BYTE_AND_STATIC_REVIEW',
        'sealed': False,
        'utc': datetime.now(timezone.utc).isoformat(),
        'package_path': str(P),
        'staging_inventory': before,
        'prior_stage_receipts': {n: check.opaque_pin(OUT / n) for n in ['INITIAL_STAGE_REVIEW_RECEIPT.json', 'UPDATED_README_STATIC_REVIEW_RECEIPT.json', 'REFRESHED_STAGE_REVIEW_RECEIPT.json']},
        'observed_snapshot': check.opaque_pin(OUT / 'OBSERVED_FINAL_PRE_SEAL_SNAPSHOT.json'),
        **check.summary(files),
        'all_previous_approved_leaf_hashes_reverified': True,
        'added_supplemental_files': {n: files[n] for n in added},
        'previous_leaf_modifications': [],
        'previous_leaf_removals': [],
        'complete_C_current_repository_identity': check.summary(scientific_copy),
        'registered_files': 117,
        'source_registration': registration_pin,
        'input_dependencies': dependencies,
        'input_dependency_count': 11,
        'eight_NPZ_total_bytes': npz_bytes,
        'actual_modes_exact_complete_root_copy': modes,
        'helper': files['replay_checkpoint.py'],
        'public_go': files['PUBLIC_GO.json'],
        'RESULTS_external_approved_pin': RESULTS_PIN,
        'results_presentation_receipt': files[presentation_path],
        'results_presentation_receipt_status_observed': presentation.get('status'),
        'RESULTS_numerical_recalculation_by_this_reviewer': 'NOT_PERFORMED',
        'python_AST_count': len(python_files),
        'python_AST': 'PASS_NO_PACKAGE_IMPORTS_OR_EXECUTIONS',
        'safe_regular_single_link_leaf_files': len(files),
        'unsafe_links_specials_caches_or_empty_directories': 0,
        'root_additions_pending': stage['root_additions_pending'],
        'final_payload_manifest_exists': False,
        'final_replay_pins_exists': False,
        'final_seal': 'ROOT_PENDING_AFTER_THIS_BYTE_REVIEW',
        'fresh_portable_replay': 'NOT_RUN_BY_THIS_REVIEWER',
        'retained_array_decodes': 0,
        'scientific_output_numeric_reads': 0,
        'study_imports': 0,
        'source_callbacks': 0,
        'target_evaluations': 0,
        'generator_or_archived_script_executions': 0,
        'network_calls': 0,
        'repository_writes': 0,
        'package_writes': 0,
        'review_script': check.opaque_pin(Path(__file__)),
        'shared_original_reviewer_script': check.opaque_pin(OUT / 'review_unsealed_stage.py'),
        'limitations': ['This receipt authenticates the final unsealed 239-leaf snapshot before the root adds externally pinned final controls.', 'No final payload seal, archive identity, registered scientific readback rerun, fresh portable replay or publication outcome is established by this review.'],
    }
    check.json_write(OUT / 'FINAL_STAGE_REVIEW_RECEIPT.json', report)
    print(json.dumps({'status': report['status'], 'files': len(files), 'bytes': report['total_leaf_bytes'], 'C': check.summary(scientific_copy), 'registered': 117, 'dependencies': len(dependencies), 'NPZ_bytes': npz_bytes, 'receipt': check.opaque_pin(OUT / 'FINAL_STAGE_REVIEW_RECEIPT.json')}, sort_keys=True))


if __name__ == '__main__':
    main()
