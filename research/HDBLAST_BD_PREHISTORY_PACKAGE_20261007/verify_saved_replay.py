"""Read-only independent review of saved portable replay outputs.

This reviewer imports only authenticated JSON validators, never numerical
routes, physical source builders, or saved quantum arrays.
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
PINS = 'a839dc8ccd4efb2be41711093987eabd2e7ee39bdf4fc2fd9be646db82c72243'
MANIFEST = '1bb65ecadc6ea11b9ffc58c4e88021b82728bb3696f50a7941468b4428a1b83a'
HELPER = '37df86c5850c71b8c8cc4cb1926dbf21a6e841ff42358a272e5eb2f6bf293c22'
REGISTRATION = '1d55b6310631fcf600a8165ec1d326b9880a0de346b531c9d9659336778f765a'
GO = '5e28dbeb90b3eef8d3cf429aed69b6791e4c4529b1852928dfb72868642001d6'
COMMIT = 'a253f50158051881f18bbdec81c12dbbb4f49eba'

def require(ok, label):
    if not ok:
        raise ValueError(label)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def identity(path):
    return {'bytes': path.stat().st_size, 'sha256': digest(path)}

def load(path):
    def pairs(items):
        value = {}
        for key, item in items:
            require(key not in value, 'Duplicate JSON key')
            value[key] = item
        return value
    def nonfinite(value):
        raise ValueError('Nonfinite JSON ' + value)
    return json.loads(path.read_bytes(), object_pairs_hook=pairs, parse_constant=nonfinite)

def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()

def safe(root, relative):
    require(type(relative) is str and relative and '\\' not in relative, 'Filename required')
    name = PurePosixPath(relative)
    require(not name.is_absolute() and name.as_posix() == relative and
            all(part not in ('', '.', '..') for part in relative.split('/')), 'Unsafe relative name')
    path = root / relative
    require(path.is_file() and all(not p.is_symlink() for p in (path, *path.parents)), 'Real regular file required')
    return path

def verify_manifest(root, manifest):
    require(set(manifest) == {'schema_version', 'files'} and type(manifest['schema_version']) is int and
            manifest['schema_version'] == 1, 'Manifest schema')
    files = manifest['files']
    for name, item in files.items():
        require(item == identity(safe(root, name)), 'Manifest leaf bytes differ ' + name)
    leaves = set()
    for path in root.rglob('*'):
        require(not path.is_symlink(), 'Manifest tree symlink')
        if path.is_file():
            leaves.add(path.relative_to(root).as_posix())
    require(leaves == set(files), 'Manifest leaf inventory differs')
    return len(leaves)

def pinned_module(root, path, label):
    # The caller authenticated the whole 101-leaf checkpoint first.
    spec = importlib.util.spec_from_file_location(label, safe(root, path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def budget_core(budget, route):
    core = dict(budget)
    if route == 'primary':
        removed = {key: core.pop(key) for key in ('wall_seconds', 'peak_rss_kib')}
    else:
        require(type(core.get('resource')) is dict and set(core['resource']) == {'wall_seconds', 'peak_rss_kib'},
                'Only declared resource fields may differ')
        removed = core['resource']
        core['resource'] = {}
    require(type(removed['wall_seconds']) in (int, float) and removed['wall_seconds'] >= 0 and
            type(removed['peak_rss_kib']) is int and removed['peak_rss_kib'] >= 0, 'Resource measurements')
    return core

def review(extracted, replay):
    package = extracted / 'HDBLAST_BD_PREHISTORY_PACKAGE_20261007'
    root = extracted / 'HDBLAST_CHECKPOINT_20261007_BD_PREHISTORY'
    for name, expected in [('REPLAY_PINS.json', PINS), ('PAYLOAD_MANIFEST.json', MANIFEST),
                           ('PUBLIC_GO.json', GO), ('replay_checkpoint.py', HELPER)]:
        require(digest(safe(package, name)) == expected, 'External package byte pin differs ' + name)
    pins = load(package / 'REPLAY_PINS.json')
    require(pins['payload_manifest_sha256'] == MANIFEST and pins['registration_sha256'] == REGISTRATION and
            pins['public_go_sha256'] == GO, 'Release pins differ')
    manifest = load(package / 'PAYLOAD_MANIFEST.json')
    count = verify_manifest(root, manifest)
    require(count == 101 and digest(root / 'FULL_REGISTRATION.json') == REGISTRATION, 'Full payload/registration differs')
    full = load(root / 'FULL_REGISTRATION.json')
    for name, item in full['files'].items():
        require(item == identity(safe(root, name)), 'Registered leaf bytes differ ' + name)
    go = load(package / 'PUBLIC_GO.json')
    require(go['freeze_commit'] == COMMIT and go['registration_sha256'] == REGISTRATION and
            go['all_registered_public_blobs_independently_downloaded'] is True and
            go['new_target_evaluations_before_public_verification'] ==
            {'likelihood_evaluations': 0, 'physical_trajectories': 0, 'registered_source_callbacks': 0,
             'retained_array_decodes': 0, 'stored_state_comparisons': 0}, 'Public source freeze metadata differs')
    receipt = load(replay / 'REPLAY_RECEIPT.json')
    require(receipt['status'] == 'PASS_COMPLETE_PORTABLE_BD_PREHISTORY_REPLAY' and
            receipt['helper_sha256'] == HELPER and receipt['registration_sha256'] == REGISTRATION and
            receipt['public_go_sha256'] == GO and receipt['replay_pins_sha256'] == PINS and
            receipt['payload_manifest_sha256'] == MANIFEST and receipt['payload_files'] == 101 and
            receipt['source_callbacks_in_four_registered_replays'] == 88 and receipt['saved_array_decodes'] == 0,
            'Final replay receipt differs')
    require(receipt['original_binary80_state_error'] == 'NOT_ENCLOSED' and
            receipt['full_twelve_case_certificate'] == 'UNRESOLVED' and receipt['metric_calibration'] == 'FAIL' and
            receipt['higher_dimensional_Big_Bang_origin'] == 'NOT_ESTABLISHED', 'Scientific scope widened')
    require(len(receipt['steps']) == 18, 'All four routes and retained proof/control steps required')
    for step in receipt['steps']:
        require(step['exit_code'] == 0 and step['log_sha256'] == digest(safe(replay, step['step'] + '.log')),
                'Step status/log pin differs ' + step['step'])
    sys.path.insert(0, str(root / 'execution'))
    validate = pinned_module(root, 'execution/validate_outputs.py', 'replay_review_saved_validator')
    post = pinned_module(root, 'analysis/certify_continuous_model.py', 'replay_review_saved_postprocessor')
    data_pins = {}
    core_pins = {}
    resources = {}
    callback_total = 0
    mode_rows = []
    original_certificate = load(safe(root, pins['reference_continuous_certificate']))
    original_certificate_core = dict(original_certificate)
    original_certificate_core.pop('provenance')
    for mode, optimized in [('normal', False), ('optimized', True)]:
        collected = {}
        for route in ('primary', 'independent'):
            relative = route + '-' + mode
            path = replay / relative
            data = load(path / 'DATA.json')
            budget = load(path / 'BUDGET.json')
            entry = load(path / 'ENTRY_RECEIPT.json')
            execution = load(path / 'EXECUTION.json')
            require(execution['status'] == 'PASS_BOUNDED_REGISTERED_ENTRY' and execution['exit_code'] == 0 and
                    execution['fabricated_only'] is False and execution['optimized'] is optimized and
                    execution['route'] == route and execution['registration_sha256'] == REGISTRATION and
                    execution['remote_receipt_sha256'] == GO and execution['uid'] != 0 and execution['stop_reason'] is None,
                    'Bounded execution metadata differs ' + relative)
            require(execution['outer_limits'] == {'wall_seconds': 900, 'peak_rss_kib': 524288, 'output_bytes': 20971520} and
                    execution['wall_seconds'] <= 900 and execution['peak_rss_kib_wait4'] <= 524288 and
                    execution['worker_output_bytes'] <= 20971520, 'Outer limits exceeded ' + relative)
            require(execution['runtime']['package_versions'] == {'python-flint': '0.9.0', 'sympy': '1.14.0', 'mpmath': '1.3.0'} and
                    execution['runtime']['sys_version'].startswith('3.12.14 ') and execution['runtime']['platform'] == 'linux',
                    'Pinned runtime differs ' + relative)
            for name, item in execution['files'].items():
                require(item == identity(safe(path, name)), 'Bounded result file pin differs ' + relative + '/' + name)
            require(entry['status'] == 'PASS_AUTHENTICATED_OFFLINE_REGISTERED_ENTRY' and entry['route'] == route and
                    entry['fabricated_only'] is False and entry['registration_sha256'] == REGISTRATION and
                    entry['receipt_sha256'] == GO and entry['freeze_commit'] == COMMIT, 'Entry authentication differs')
            for name in ('DATA.json', 'BUDGET.json'):
                require(entry['files'][name] == identity(path / name), 'Entry scientific-byte binding differs')
            events = entry['source_callbacks']['events']
            journal = [json.loads(line) for line in (path / 'SOURCE_ATTEMPTS.jsonl').read_bytes().splitlines()]
            require(entry['source_callbacks']['physical_source_callbacks'] == 22 and
                    entry['source_callbacks']['archive_arrays_decoded'] == 0 and len(events) == 22 and
                    canonical(events) == canonical(journal), 'Callback inventory/journal differs')
            callback_total += 22
            original = safe(root, pins['reference_data'][route])
            require((path / 'DATA.json').read_bytes() == original.read_bytes(), 'Original exact DATA differs ' + relative)
            core = budget_core(budget, route)
            original_core = budget_core(load(original.with_name('BUDGET.json')), route)
            require(canonical(core) == canonical(original_core), 'Original canonical scientific budget differs ' + relative)
            core_sha = hashlib.sha256(canonical(core)).hexdigest()
            require(receipt['exact_scientific_budget_core_sha256'][relative] == core_sha, 'Budget core receipt pin differs')
            data_pins[relative] = digest(path / 'DATA.json')
            core_pins[relative] = core_sha
            resources[relative] = {key: execution[key] for key in ('wall_seconds', 'peak_rss_kib_wait4', 'worker_output_bytes')}
            collected[route] = (data, budget)
        intersections = validate.compare_outputs(*collected['primary'], *collected['independent'], False, True)
        require(intersections['overlap_rectangles'] == 36 and intersections['component_intersections'] == 72,
                'Saved independent intersections incomplete ' + mode)
        retained_intersections = load(replay / ('INTERSECTIONS_' + mode + '.json'))
        require(canonical(intersections) == canonical(retained_intersections), 'Intersection receipt differs')
        certificate = load(replay / ('CONTINUOUS_MODEL_' + mode + '.json'))
        provenance = certificate.pop('provenance')
        independent = replay / ('independent-' + mode)
        require(provenance['independent_budget_sha256'] == digest(independent / 'BUDGET.json') and
                provenance['independent_entry_receipt_sha256'] == digest(independent / 'ENTRY_RECEIPT.json') and
                provenance['registration_sha256'] == REGISTRATION and provenance['freeze_commit'] == COMMIT and
                provenance['postprocessor_sha256'] == digest(root / 'analysis/certify_continuous_model.py'),
                'Fresh certificate execution provenance differs ' + mode)
        derived = post.build_certificate(collected['independent'][1])
        require(canonical(derived) == canonical(certificate) == canonical(original_certificate_core),
                'Continuous exact mathematical core differs ' + mode)
        mode_rows.append({'mode': mode, 'rectangle_intersections': 36, 'component_intersections': 72,
                          'continuous_certificate_core_sha256': hashlib.sha256(canonical(derived)).hexdigest()})
    for route in ('primary', 'independent'):
        require(data_pins[route + '-normal'] == data_pins[route + '-optimized'] and
                core_pins[route + '-normal'] == core_pins[route + '-optimized'] and
                receipt['exact_original_and_normal_optimized_data_sha256'][route] == data_pins[route + '-normal'],
                'Normal/-O exact data/budget scientific core differs ' + route)
    proof_tasks = [('representation', 'evidence/preparation/PRIMARY_EXACT_'),
                   ('independent-cap', 'evidence/actual/INDEPENDENT_CAP_CURRENT_'),
                   ('continuous-selftest', 'evidence/preparation/CONTINUOUS_MODEL_SYNTHETIC_'),
                   ('primary-fabricated', 'evidence/preparation/PRIMARY_FABRICATED_'),
                   ('outer-controls', 'evidence/preparation/OUTER_CONTROLS_')]
    proof_pins = {}
    for mode in ('normal', 'optimized'):
        for label, prefix in proof_tasks:
            name = label + '-' + mode + '.json'
            reference = safe(root, prefix + mode.upper() + '.json')
            require(safe(replay, name).read_bytes() == reference.read_bytes(), 'Retained proof/control bytes differ ' + name)
            proof_pins[name] = digest(replay / name)
    require(callback_total == 88 and verify_manifest(root, manifest) == 101, 'Final callback/payload reconciliation differs')
    return {'status': 'PASS_AUTHOR_INDEPENDENT_CLEAN_REPLAY_SAVED_EVIDENCE_REVIEW',
            'external_pins': {'pins': PINS, 'manifest': MANIFEST, 'helper': HELPER, 'registration': REGISTRATION, 'public_GO': GO},
            'payload_leaves_verified': 101, 'registered_files_verified': len(full['files']),
            'input_final_REPLAY_RECEIPT_sha256': digest(replay / 'REPLAY_RECEIPT.json'),
            'data_sha256': data_pins, 'canonical_scientific_budget_core_sha256': core_pins,
            'retained_proof_control_sha256': proof_pins, 'modes': mode_rows, 'bounded_execution_resources': resources,
            'source_callbacks_in_review': 0, 'recorded_source_callbacks_in_four_replays': callback_total,
            'numerical_route_imports_in_review': 0, 'saved_quantum_arrays_decoded': 0,
            'scope': 'Saved payload/result authentication and exact rational JSON validators only; original scientific limits unchanged.'}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extracted', type=Path, required=True)
    parser.add_argument('--replay', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = review(args.extracted, args.replay)
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': result['status'], 'payload_leaves_verified': 101,
                      'physical_source_callbacks_in_review': 0, 'receipt_sha256': digest(args.output)}))

if __name__ == '__main__':
    main()
