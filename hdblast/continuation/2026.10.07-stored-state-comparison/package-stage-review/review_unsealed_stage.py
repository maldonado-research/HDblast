#!/usr/bin/env python3
"""Original stdlib-only byte/static review; imports no package/scientific code."""
import ast
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
from datetime import datetime, timezone

W = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
P = W / 'package-candidate/HDBLAST_STORED_STATE_COMPARISON_PACKAGE_20261007'
C_REL = 'research/HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON'
C = P / C_REL
REPO = Path('/workspace/HDblast')
SOURCE_C = REPO / C_REL
OUT = W / 'replay-delivery/package-stage-independent-review'
STAGE = W / 'package-candidate/STAGING_INVENTORY.json'
REG_SHA = '05e9943f0b4c8134252a2fecef7631ddba8bd398d18553d6e126fbf50fc3aacc'
GO_SHA = '27784552d9d7e10167066927d54f1757351b0f925b499e951b79849a44b62769'
HELPER_SHA = 'c4d04d22aff10d04dd93ee70ca94e45fcd04d9bbafd69fbcff3a919bafc10377'
FREEZE = 'bedcca7e86995da1230c12a20fae3755b31f4a94'
HEX = re.compile(r'[0-9a-f]{64}\Z')
FORBIDDEN_DIRS = {'.git', '.aws', '.codex', '.agents', '__pycache__', '.pytest_cache', '.mypy_cache', '.cache', 'node_modules'}
FORBIDDEN_SUFFIXES = {'.pyc', '.pyo', '.so', '.dll', '.dylib', '.o', '.a'}


def fail(message):
    raise RuntimeError(message)


def require(condition, message):
    if not condition:
        fail(message)


def rel_safe(name):
    require(isinstance(name, str), 'non-string relative path')
    pp = PurePosixPath(name)
    require(name and not pp.is_absolute() and str(pp) == name and all(x not in ('', '.', '..') for x in pp.parts), 'unsafe relative path: ' + str(name))
    require(not any(ord(x) < 32 or ord(x) == 127 for x in name) and '\\' not in name, 'unsafe path characters: ' + name)
    return name


def ancestor_check(path):
    for q in [path, *path.parents]:
        st = q.lstat()
        require(not stat.S_ISLNK(st.st_mode), 'symlink ancestor: ' + str(q))


def identity(st):
    return (st.st_dev, st.st_ino, st.st_mode, st.st_nlink, st.st_size, st.st_mtime_ns, st.st_ctime_ns)


def opaque_pin(path):
    ancestor_check(path)
    before = path.lstat()
    require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1, 'unsafe leaf: ' + str(path))
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        opened = os.fstat(fd)
        require(identity(before) == identity(opened), 'identity changed at open: ' + str(path))
        h = hashlib.sha256()
        n = 0
        while True:
            chunk = os.read(fd, 1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
            n += len(chunk)
        after = os.fstat(fd)
        require(identity(opened) == identity(after) == identity(path.lstat()), 'identity changed during stream: ' + str(path))
        require(n == before.st_size, 'short stream: ' + str(path))
        return {'bytes': n, 'sha256': h.hexdigest()}
    finally:
        os.close(fd)


def tree(root):
    ancestor_check(root)
    require(stat.S_ISDIR(root.lstat().st_mode), 'not directory: ' + str(root))
    files = {}
    dirs = []
    def descend(parent):
        before = identity(parent.lstat())
        children = sorted(parent.iterdir(), key=lambda x: x.name)
        require(children, 'empty directory: ' + str(parent))
        for child in children:
            name = rel_safe(child.relative_to(root).as_posix())
            st = child.lstat()
            require(not stat.S_ISLNK(st.st_mode), 'symlink entry: ' + name)
            if stat.S_ISDIR(st.st_mode):
                require(child.name not in FORBIDDEN_DIRS, 'cache/private directory: ' + name)
                dirs.append(name)
                descend(child)
            else:
                require(stat.S_ISREG(st.st_mode) and st.st_nlink == 1, 'hardlink or special entry: ' + name)
                require(child.suffix.lower() not in FORBIDDEN_SUFFIXES, 'cached/native artifact: ' + name)
                files[name] = opaque_pin(child)
        require(before == identity(parent.lstat()), 'directory changed during scan: ' + str(parent))
    descend(root)
    return dict(sorted(files.items())), sorted(dirs)


def duplicate_guard(pairs):
    result = {}
    for k, v in pairs:
        require(k not in result, 'duplicate JSON key: ' + k)
        result[k] = v
    return result


def load_json(path):
    raw = path.read_bytes()
    require(raw.decode('utf-8').encode('utf-8') == raw, 'not exact UTF8 JSON')
    return json.loads(raw, object_pairs_hook=duplicate_guard, parse_constant=lambda x: fail('nonfinite JSON constant: ' + x))


def validate_pin(value):
    require(type(value) is dict and set(value) == {'bytes', 'sha256'}, 'invalid pin schema')
    require(type(value['bytes']) is int and value['bytes'] >= 0 and type(value['sha256']) is str and HEX.fullmatch(value['sha256']), 'invalid typed byte pin')


def json_write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def summary(files):
    return {'file_count': len(files), 'total_leaf_bytes': sum(row['bytes'] for row in files.values())}


def subset(files, prefix):
    prefix += '/'
    return {name[len(prefix):]: row for name, row in files.items() if name.startswith(prefix)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    stage_pin_before = opaque_pin(STAGE)
    stage = load_json(STAGE)
    require(stage['sealed'] is False and stage['status'] == 'UNSEALED_PORTABLE_PACKAGE_STAGED', 'not expected unsealed stage')
    expected = stage['files']
    require(type(expected) is dict, 'stage files map required')
    for name, pin in expected.items():
        rel_safe(name)
        validate_pin(pin)
    require(len(expected) == 227 and stage['file_count'] == 227 and stage['total_leaf_bytes'] == 257464451, 'initial snapshot changed')
    observed, dirs = tree(P)
    require(observed == expected, 'P does not equal initial staging inventory')
    require(summary(observed) == {'file_count': 227, 'total_leaf_bytes': 257464451}, 'initial counts differ')
    require(not (P / 'PAYLOAD_MANIFEST.json').exists() and not (P / 'REPLAY_PINS.json').exists(), 'root controls prematurely present')
    require(observed['replay_checkpoint.py']['sha256'] == HELPER_SHA, 'helper identity differs')
    require(observed['PUBLIC_GO.json']['sha256'] == GO_SHA, 'public GO identity differs')
    require(opaque_pin(W / 'root/public-readback/PUBLIC_GO.json') == observed['PUBLIC_GO.json'], 'GO copy differs from actual source')
    reg_path = C / 'FULL_REGISTRATION.json'
    reg_pin = opaque_pin(reg_path)
    require(reg_pin['sha256'] == REG_SHA, 'FULL_REGISTRATION identity differs')
    require(opaque_pin(SOURCE_C / 'FULL_REGISTRATION.json') == reg_pin, 'repository registration differs')
    registration = load_json(reg_path)
    require(set(registration) == {'files', 'schema_version'} and type(registration['schema_version']) is int and registration['schema_version'] == 1, 'registration schema differs')
    require(type(registration['files']) is dict and len(registration['files']) == 117, '117 registered files required')
    for name, pin in registration['files'].items():
        rel_safe(name)
        validate_pin(pin)
        require(observed[C_REL + '/' + name] == pin, 'registered P leaf differs: ' + name)
        require(opaque_pin(SOURCE_C / name) == pin, 'registered source leaf differs: ' + name)

    spec = load_json(C / 'provenance/INPUT_SPEC.json')
    inputs = {}
    for capsule in spec['capsules']:
        inputs[rel_safe(capsule['repository_path'])] = {'bytes': capsule['reader_inputspec']['file_bytes'], 'sha256': capsule['reader_inputspec']['file_sha256']}
        original = capsule['original']
        inputs[rel_safe(original['repository_path'])] = {'bytes': original['bytes'], 'sha256': original['sha256']}
    npz_bytes = sum(v['bytes'] for v in inputs.values())
    require(len(inputs) == 8 and npz_bytes == 212975098, 'eight NPZ pins/count differ')
    manifest_path = rel_safe(spec['manifest_repository_path'])
    producer_path = rel_safe(spec['producer_repository_path'])
    inputs[manifest_path] = opaque_pin(P / manifest_path)
    require(inputs[manifest_path]['sha256'] == spec['manifest_sha256'], 'prior manifest SHA differs')
    inputs[producer_path] = opaque_pin(P / producer_path)
    require(inputs[producer_path]['sha256'] == spec['producer_sha256'], 'producer SHA differs')
    audit = spec['prior_static_audit']
    inputs[rel_safe(audit['repository_path'])] = {'bytes': audit['bytes'], 'sha256': audit['sha256']}
    require(len(inputs) == 11 and inputs == stage['input_dependencies_verified'], '11 input dependency pins differ')
    for name, pin in inputs.items():
        validate_pin(pin)
        require(observed[name] == pin, 'package dependency mismatch: ' + name)
        require(opaque_pin(REPO / name) == pin, 'dependency source mismatch: ' + name)

    modes = {}
    stable_names = ['DATA.json', 'SCIENCE_SUMMARY.json', 'SOURCE_CERTIFICATE.json', 'DIAGNOSTIC_CROSSCHECK.json', 'NODE_TARGETS.jsonl.gz', 'SOURCE_ATTEMPTS.jsonl', 'DECODE_ATTEMPTS.jsonl', 'OUTPUT_READBACK.json']
    limits = {'address_space_bytes': 536870912, 'aggregate_output_bytes': 134217728, 'core_bytes': 0, 'cpu_seconds': 900, 'custodian_receipt_reserve_bytes': 65536, 'file_size_bytes': 134217728, 'process_limit': 0, 'wall_seconds': 900}
    for mode in ['normal', 'optimized']:
        dest = P / 'evidence/actual' / mode
        copied = subset(observed, 'evidence/actual/' + mode)
        source, _ = tree(W / ('root/actual-' + mode + '-001'))
        require(copied == source, 'actual mode is not a complete exact source copy: ' + mode)
        require(len(copied) == 12 and all(name in copied for name in ['ENTRY_RECEIPT.json', 'EXECUTION.json', 'RESOURCE.json', 'child.log', *stable_names]), 'actual mode leaves incomplete: ' + mode)
        execution = load_json(dest / 'EXECUTION.json')
        entry = load_json(dest / 'ENTRY_RECEIPT.json')
        require(execution['status'] == 'PASS_BOUNDED_REGISTERED_EXECUTION' and execution['fabricated_only'] is False and execution['optimized'] is (mode == 'optimized'), 'actual custodian status/mode differs')
        require(all(type(execution[k]) is int and execution[k] == 0 for k in ['exit_code', 'wrapper_exit_code', 'wait_status']), 'actual outcomes are not integer zero')
        require(execution['stop_reason'] is None and execution['monitor_error'] is None, 'actual monitor did not accept')
        require(execution['limits'] == limits and all(type(v) is int for v in execution['limits'].values()), 'fixed custody limits differ')
        require(execution['entry_receipt'] == copied['ENTRY_RECEIPT.json'], 'custodian entry link differs')
        require(execution['captured_child_log_sha256'] == execution['observed_child_log_sha256'] == copied['child.log']['sha256'] and execution['preserved_child_log'] is None, 'custodian log links differ')
        require(type(execution['uid']) is int and execution['uid'] > 0, 'actual nonroot UID missing')
        require(type(execution['peak_rss_kib_wait4']) is int and 0 < execution['peak_rss_kib_wait4'] <= 524288, 'actual RSS outside bound')
        for key in ['wall_seconds', 'cpu_seconds_wait4']:
            val = execution[key]
            require(type(val) in (int, float) and math.isfinite(val) and 0 <= val <= 900, 'actual timing outside bound')
        worker = execution['worker_output']
        require(type(worker['bytes']) is int and 0 < worker['bytes'] <= 134217728 - 65536 and type(worker['file_count']) is int and worker['file_count'] == 11, 'worker custody size/count differs')
        require(entry['status'] == 'PASS_AUTHENTICATED_STORED_COMPARISON' and entry['fabricated_only'] is False and entry['registration_sha256'] == REG_SHA and entry['public_go_sha256'] == GO_SHA and entry['freeze_commit'] == FREEZE, 'entry identity differs')
        for name, pin in entry['files'].items():
            rel_safe(name)
            validate_pin(pin)
            require(copied[name] == pin, 'entry opaque file link differs: ' + name)
        modes[mode] = {'file_count': len(copied), 'total_leaf_bytes': sum(v['bytes'] for v in copied.values()), 'entry_sha256': copied['ENTRY_RECEIPT.json']['sha256'], 'execution_sha256': copied['EXECUTION.json']['sha256'], 'child_log_sha256': copied['child.log']['sha256'], 'stable_files': {n: copied[n] for n in stable_names}, 'static_custody_links': 'PASS'}
    require(modes['normal']['stable_files'] == modes['optimized']['stable_files'], 'eight stable opaque output byte identities differ')

    whitelist_root = P / 'preparation/root-independent-review'
    whitelist_source = W / 'actual-independent-review'
    artifact_pin = opaque_pin(whitelist_root / 'ARTIFACT_PINS.json')
    require(opaque_pin(whitelist_source / 'ARTIFACT_PINS.json') == artifact_pin, 'review whitelist manifest differs')
    artifacts = load_json(whitelist_root / 'ARTIFACT_PINS.json')
    require(len(artifacts['files']) == 9 and 'review_serialized_norms.py' in artifacts['files'], 'nine-artifact whitelist differs')
    require(artifacts['registration_sha256'] == REG_SHA and artifacts['public_go_sha256'] == GO_SHA, 'review provenance differs')
    require(set(subset(observed, 'preparation/root-independent-review')) == set(artifacts['files']) | {'ARTIFACT_PINS.json'}, 'review folder exceeds public-safe whitelist')
    for name, pin in artifacts['files'].items():
        rel_safe(name)
        validate_pin(pin)
        require(observed['preparation/root-independent-review/' + name] == pin and opaque_pin(whitelist_source / name) == pin, 'review whitelist bytes differ: ' + name)

    python_files = [name for name in observed if name.endswith('.py')]
    for name in python_files:
        ast.parse((P / name).read_text(encoding='utf-8'), filename=name)
    require(not any(name.lower().endswith('.pdf') for name in observed), 'initial snapshot unexpectedly contains PDF body')
    readme = (P / 'README.md').read_text(encoding='utf-8')
    require(not re.search(r'\b[0-9a-fA-F]{64}\b', readme), 'README contains concrete control/self SHA')
    for phrase in ['49,152', '12 finite', '20 selected', 'U_saved = u_1 / epsilon_represented', 'W_saved = w_1 / epsilon_represented', '1e-18', 'UNRESOLVED', '**FAIL**', 'NOT_ESTABLISHED', 'NOT_ASSESSED', '--verify-only', '--output-root']:
        require(phrase in readme, 'README essential scope/unit/command missing: ' + phrase)

    source_now, _ = tree(SOURCE_C)
    c_observed = subset(observed, C_REL)
    extra_now = sorted(set(source_now) - set(c_observed))
    missing_now = sorted(set(c_observed) - set(source_now))
    changed_now = sorted(name for name in set(source_now) & set(c_observed) if source_now[name] != c_observed[name])
    require(not missing_now and not changed_now, 'copied C source leaf changed or disappeared')
    require(opaque_pin(STAGE) == stage_pin_before, 'STAGING_INVENTORY changed during review')
    end_observed, end_dirs = tree(P)
    require(end_observed == observed and end_dirs == dirs, 'P changed during review')
    snapshot = {'schema_version': 1, 'status': 'INDEPENDENT_INITIAL_UNSEALED_SNAPSHOT', 'sealed': False, 'package_path': str(P), 'files': observed, 'directories': dirs, **summary(observed)}
    json_write(OUT / 'OBSERVED_INITIAL_SNAPSHOT.json', snapshot)
    report = {
        'schema_version': 1,
        'status': 'PASS_INITIAL_UNSEALED_STAGE_BYTE_AND_STATIC_REVIEW',
        'sealed': False,
        'utc': datetime.now(timezone.utc).isoformat(),
        'package_path': str(P),
        'staging_inventory': stage_pin_before,
        'observed_snapshot': opaque_pin(OUT / 'OBSERVED_INITIAL_SNAPSHOT.json'),
        **summary(observed),
        'safe_regular_single_link_leaf_files': len(observed),
        'directory_count': len(dirs),
        'empty_or_unsafe_directories': 0,
        'cached_native_or_private_cache_files': 0,
        'registered_files': 117,
        'registration': reg_pin,
        'public_go': observed['PUBLIC_GO.json'],
        'helper': observed['replay_checkpoint.py'],
        'input_dependencies': inputs,
        'input_dependency_count': 11,
        'eight_npz_bytes': npz_bytes,
        'actual_modes': modes,
        'actual_eight_stable_byte_identity': 'PASS_OPAQUE_HASHES_ONLY',
        'independent_review_whitelist_count': 9,
        'independent_review_artifact_manifest': artifact_pin,
        'python_AST_syntax_count': len(python_files),
        'python_AST_syntax': 'PASS_NO_IMPORT_OR_EXECUTION',
        'README': observed['README.md'],
        'README_scope_and_generic_control_commands': 'PASS_INITIAL_TEXT',
        'README_R_P_normalized_units_and_time_interval': 'PENDING_PARENT_SUPPLEMENTAL_TEXT',
        'repository_C_snapshot_comparison': {'source_current': summary(source_now), 'package_copy': summary(c_observed), 'source_extra_pending_additions': extra_now, 'source_missing': missing_now, 'source_changed': changed_now},
        'final_payload_manifest_exists': False,
        'final_replay_pins_exists': False,
        'fresh_portable_replay': 'NOT_RUN_IN_THIS_REVIEW',
        'publication': 'NOT_PERFORMED_IN_THIS_REVIEW',
        'retained_array_decodes': 0,
        'scientific_output_numeric_reads': 0,
        'study_imports': 0,
        'source_callbacks': 0,
        'target_evaluations': 0,
        'network_calls': 0,
        'repository_writes': 0,
        'package_writes': 0,
        'archive_script_executions': 0,
        'limitations': ['Initial stage snapshot only; expected later RESULTS/figures/approved generator additions are outside this receipt.', 'README normalized R/P units and time interval are pending parent update.', 'Opaque saved-output identity and static custody links do not replace registered scientific readback or a fresh portable replay.', 'No final control seal, external publication identity, remote custody or publication conclusion is established.'],
        'review_script': opaque_pin(Path(__file__)),
    }
    json_write(OUT / 'INITIAL_STAGE_REVIEW_RECEIPT.json', report)
    print(json.dumps({'status': report['status'], 'file_count': report['file_count'], 'bytes': report['total_leaf_bytes'], 'registration_files': 117, 'dependencies': 11, 'NPZ_bytes': npz_bytes, 'python_AST': len(python_files), 'source_current_extra': extra_now, 'receipt': opaque_pin(OUT / 'INITIAL_STAGE_REVIEW_RECEIPT.json')}, sort_keys=True))


if __name__ == '__main__':
    main()
