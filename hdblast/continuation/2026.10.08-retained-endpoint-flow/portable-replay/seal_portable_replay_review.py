"""Seal already completed portable replay evidence; never invokes a worker."""
from pathlib import Path
import hashlib
import json
import stat

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
FRESH = HERE / 'fresh-pair-001'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    return json.loads(path.read_bytes())


def pin(path):
    metadata = path.lstat()
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1, 'unaliased regular evidence')
    raw = path.read_bytes()
    require(len(raw) == metadata.st_size, 'evidence size changed')
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def write_new(name, value):
    with (HERE / name).open('x') as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write('\n')


def main():
    package = load(HERE / 'PACKAGE_AND_OPAQUE_DEPENDENCIES.json')
    post = load(HERE / 'POST_EXECUTION_PACKAGE_AND_DEPENDENCIES.json')
    require(package == post, 'package/original dependency reauthentication differed')
    verify = load(HERE / 'VERIFY_ONLY_ACCEPTANCE.json')
    require(verify['status'] == 'PASS_INDEPENDENT_CAPTURED_BOOTSTRAP_VERIFY_ONLY', 'verify-only acceptance')
    audit = load(HERE / 'FRESH_PAIR_INDEPENDENT_AUDIT.json')
    require(audit['status'] == 'PASS_INDEPENDENT_FRESH_ACTUAL_PORTABLE_PAIR_CUSTODY_AND_IDENTITY', 'independent pair audit')
    replay = load(FRESH / 'PORTABLE_REPLAY_RECEIPT.json')
    require(load(HERE / 'FRESH_PAIR_STDOUT.json') == replay, 'helper final stdout differs from durable receipt')
    require((HERE / 'FRESH_PAIR_STDERR.log').read_bytes() == b'', 'helper stderr not empty')
    for mode in ('normal', 'optimized'):
        require((FRESH / (mode + '-custodian.stderr.log')).read_bytes() == b'', 'custodian stderr not empty')
        require(load(FRESH / (mode + '-custodian.stdout.log')) == load(FRESH / mode / 'EXECUTION.json'), 'custodian stdout differs from receipt')
    prior = BASE / 'actual-independent-review/FINAL_ACTUAL_REVIEW_RECEIPT.json'
    require(pin(prior)['sha256'] == package['actual_review_sha256'], 'prior actual review changed')
    science = load(prior)
    require(science['scientific_file_pins'] == audit['scientific_file_pins'], 'scientific authority differs')
    prior_checks = load(HERE / 'PRIOR_REVIEWS_UNCHANGED.json')
    require(prior_checks['status'] == 'PASS_ALL_PRIOR_SEALED_REVIEW_BYTES_UNCHANGED', 'prior evidence changed')
    process = {
        'status': 'PASS_CAPTURED_BOOTSTRAP_PROCESS_AND_DURABLE_RECEIPT',
        'observed_exec_session_exit_code': 0,
        'observed_exec_session': 31146,
        'invocation': 'One captured-bootstrap actual replay invocation after independent verify-only acceptance.',
        'fresh_workers': 2,
        'worker_modes': ['normal', 'optimized'],
        'helper_stdout_matches_durable_receipt': True,
        'helper_stderr_bytes': 0,
        'custodian_stdout_matches_execution_each_mode': True,
        'custodian_stderr_bytes_each_mode': 0,
        'pins': {name: pin(HERE / name) for name in ('invoke_captured_bootstrap.py', 'trusted_bootstrap_replay.py', 'FRESH_PAIR_STDOUT.json', 'FRESH_PAIR_STDERR.log')},
        'portable_replay_receipt': pin(FRESH / 'PORTABLE_REPLAY_RECEIPT.json'),
    }
    write_new('FRESH_PAIR_PROCESS_REVIEW.json', process)
    receipt = {
        'status': 'PASS_INDEPENDENT_ACTUAL_PORTABLE_PACKAGE_AND_SINGLE_FRESH_REPLAY_PAIR',
        'scope': 'Reproduction of the accepted finite retained-endpoint certificate using unchanged source005 and the existing public GO.',
        'authorizes_additional_physical_or_retained_operations': False,
        'public_freeze_commit': package['freeze_commit'],
        'source_registration_sha256': package['source_registration_sha256'],
        'public_GO_sha256': package['public_GO_sha256'],
        'trusted_bootstrap_sha256': package['bootstrap_sha256'],
        'external_package_control_pins': package['external_controls'],
        'package': {
            'payload_files': package['payload_files'], 'total_files': package['total_files'],
            'payload_bytes': package['payload_bytes'], 'frozen_core_files': package['source_files'],
            'exact_parent_directory_roster_checked': True,
            'all_bytes_unchanged_after_both_fresh_workers': True,
            'original_input_payloads_embedded': False,
        },
        'original_dependencies': {
            'files': package['opaque_dependency_files'], 'bytes': package['opaque_dependency_bytes'],
            'file_pins': package['opaque_dependency_pins'],
            'hashes_checked_before_and_after_fresh_workers': True,
            'opaque_audit_archive_containers_parsed': 0,
            'opaque_audit_scalars_decoded': 0,
            'origin_commit': '13a30ef46f5c90c1b01ae83b609db65fd0f8a709',
            'current_repository_HEAD_equality_required': False,
        },
        'independent_verify_only': {
            'receipt': pin(HERE / 'VERIFY_ONLY_ACCEPTANCE.json'),
            'saved_modes': 2, 'nodes_per_mode': 49152, 'endpoint_comparisons_per_mode': 98304,
            'finite_aggregates_per_mode': 24, 'new_source_callbacks': 0, 'new_array_decodes': 0,
            'fresh_execution_condition_met_before_launch': True,
        },
        'execution_authorization_and_accounting': {
            'authorization': 'Explicit parent/root instruction for one actual normal/optimized portable pair after independent package/dependency and verify-only acceptance.',
            'portable_pair_invocations': 1,
            'actual_worker_executions': 2,
            'genuine_modes': ['normal', 'optimized'],
            'physical_source_callbacks': 256,
            'selected_array_decodes': 72,
            'extra_source_free_audit_callbacks_or_decodes': 0,
            'retries': 0, 'adapted_parameters': False, 'new_GO_issued': False,
            'repository_or_remote_mutations_by_this_review': 0,
            'additional_actual_executions_performed': 0,
        },
        'unchanged_per_worker_limits': {
            'wall_seconds': 900, 'cpu_seconds': 900, 'address_space_bytes': 536870912,
            'aggregate_output_bytes': 134217728, 'file_size_bytes': 134217728,
            'custodian_receipt_reserve_bytes': 65536, 'process_limit': 0, 'core_bytes': 0,
        },
        'custody': {
            'both_workers_nonroot': True, 'kernel_offline_no_children_controls_bound': True,
            'wait4_exit_and_wrapper_statuses_zero': True,
            'directory_identity_logs_and_complete_output_hashes_checked': True,
            'aggregate_output_enforcement': 'POLL_AND_FINAL_ACCEPTANCE_CHECK_WITH_RECEIPT_RESERVE',
            'portable_helper_receipt': audit['portable_replay_receipt'],
            'independent_audit_receipt': pin(HERE / 'FRESH_PAIR_INDEPENDENT_AUDIT.json'),
            'process_receipt': pin(HERE / 'FRESH_PAIR_PROCESS_REVIEW.json'),
        },
        'fresh_runs': audit['runs'],
        'coverage_per_run': science['coverage_per_run'],
        'scientific_file_pins': audit['scientific_file_pins'],
        'nine_files_identical_across_fresh_modes_packaged_saved_modes_and_root_accepted_modes': True,
        'prior_independent_actual_scientific_review': pin(prior),
        'prior_allnode_arithmetic_applies_by_exact_export_byte_identity': True,
        'prior_all_48_signed_finite_R_P_intervals_still_exclude_zero': True,
        'prior_exact_interval_authority': science['compact_review_evidence_pins']['SIGNED_INTERVAL_ZERO_EXCLUSIONS_48.json'],
        'target_export_gate': science['target_export_gate'],
        'exact_global_finite_error_upper_bounds': science['exact_global_finite_error_upper_bounds'],
        'prior_review_preservation': prior_checks['checks'],
        'prior_helper_review_sha256': '2ba3deb951a1297c983e94d9263edeb04fc76c3018e3e6b2344c84fe7f9f4076',
        'prior_helper_independent_controls': 92,
        'preserved_scientific_scope': {
            'incoming_state': 'Exact represented t1=-4.5 saved state.',
            'target_endpoints': ['-4', '-7/2'],
            'error_scope': 'Total subsequent saved-endpoint error; incoming preparation error remains separate.',
            'metric_calibration': 'FAIL_UNCHANGED',
            'complete_continuous_momentum_time_contact_UV_certificate': 'UNRESOLVED',
            'unsaved_historical_per_step_trajectory_recovered': False,
            'new_independent_analytic_proof': False,
            'independent_source_free_audit_redecoded_original_values': False,
            'original_value_identity_basis': 'Pinned authenticated decoder executed by the two authorized workers; no additional independent decoder invocation.',
        },
        'failures_in_this_portable_acceptance': [],
        'prior_failed_controls_and_rejected_source_versions_preserved': True,
        'compact_manifest_policy': 'Relative selected evidence paths; large fresh source/node artifacts excluded from compact payload but completely hash-bound by fresh_runs.',
    }
    write_new('FINAL_PORTABLE_REPLAY_REVIEW_RECEIPT.json', receipt)
    selected = {p.name for p in HERE.iterdir() if p.is_file()}
    selected.discard('FINAL_PORTABLE_REPLAY_REVIEW_MANIFEST.json')
    selected.add('fresh-pair-001/PORTABLE_REPLAY_RECEIPT.json')
    for mode in ('normal', 'optimized'):
        for stream in ('stdout', 'stderr'):
            selected.add(f'fresh-pair-001/{mode}-custodian.{stream}.log')
        for name in ('ENTRY_RECEIPT.json', 'EXECUTION.json', 'RESOURCE.json', 'child.log', 'INPUT_AUTHENTICATION.json', 'OUTPUT_READBACK.json', 'SCIENCE_SUMMARY.json', 'SOURCE_ATTEMPTS.jsonl', 'DECODE_ATTEMPTS.jsonl', 'TARGET_BUDGETS.json'):
            selected.add(f'fresh-pair-001/{mode}/{name}')
    files = {name: pin(HERE / name) for name in sorted(selected)}
    manifest = {
        'status': 'SEALED_COMPACT_INDEPENDENT_ACTUAL_PORTABLE_REPLAY_REVIEW',
        'path_base': '.', 'files': files, 'file_count': len(files),
        'total_bytes': sum(item['bytes'] for item in files.values()),
        'excluded_large_fresh_outputs_bound_in_receipt': ['DATA.json', 'NODE_CERTIFICATES.jsonl.gz', 'SOURCE_CERTIFICATE.json'],
        'final_receipt_sha256': files['FINAL_PORTABLE_REPLAY_REVIEW_RECEIPT.json']['sha256'],
    }
    write_new('FINAL_PORTABLE_REPLAY_REVIEW_MANIFEST.json', manifest)
    for name, expected in files.items():
        require(pin(HERE / name) == expected, 'compact evidence changed during seal')
    print(json.dumps({'receipt': pin(HERE / 'FINAL_PORTABLE_REPLAY_REVIEW_RECEIPT.json'), 'manifest': pin(HERE / 'FINAL_PORTABLE_REPLAY_REVIEW_MANIFEST.json'), 'files': len(files), 'bytes': manifest['total_bytes']}, sort_keys=True))


if __name__ == '__main__':
    main()
