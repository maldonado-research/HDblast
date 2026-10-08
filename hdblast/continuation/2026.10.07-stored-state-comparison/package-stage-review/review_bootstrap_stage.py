#!/usr/bin/env python3
"""Final one-leaf bootstrap correction; static AST and complete opaque bytes."""
import ast
import importlib.util
import json
import re
from pathlib import Path
from datetime import datetime, timezone

W = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
OUT = W / 'replay-delivery/package-stage-independent-review'
spec = importlib.util.spec_from_file_location('original_static_byte_reviewer', OUT / 'review_unsealed_stage.py')
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)  # Only original stdlib reviewer, outside P.
P = check.P
STAGE_SHA = 'b38d88a36b3ed7cf269f8e6d01cded24c675a41b364a1405f55c19d075b5026e'
README_PIN = {'bytes': 12883, 'sha256': '468a1585aea4e7f93cd39fd9a8929cb482d90fb484b4f86a4739e33274138b2b'}


def main():
    before = check.opaque_pin(check.STAGE)
    check.require(before['sha256'] == STAGE_SHA, 'bootstrap-final staging pin differs')
    stage = check.load_json(check.STAGE)
    check.require(stage['sealed'] is False and stage['status'] == 'UNSEALED_PORTABLE_PACKAGE_STAGED', 'bootstrap stage is not unsealed')
    previous = check.load_json(OUT / 'OBSERVED_FINAL_PRE_SEAL_SNAPSHOT.json')['files']
    files, dirs = check.tree(P)
    check.require(files == stage['files'], 'bootstrap-final complete leaf map differs')
    check.require(check.summary(files) == {'file_count': 239, 'total_leaf_bytes': 258004208}, 'bootstrap-final counts differ')
    check.require(set(files) == set(previous), 'bootstrap update added or removed a leaf')
    changed = sorted(n for n in files if files[n] != previous[n])
    check.require(changed == ['README.md'], 'bootstrap update changed a leaf besides README')
    check.require(files['README.md'] == README_PIN == check.opaque_pin(W / 'replay-delivery/PACKAGE_README_UNSEALED.md'), 'bootstrap README/source pin differs')
    text = (P / 'README.md').read_text(encoding='utf-8')
    check.require(re.search(r'\b[0-9a-fA-F]{64}\b', text) is None, 'README embeds a concrete control/self SHA256')
    shell = re.findall(r'```sh\n(.*?)\n```', text, flags=re.S)
    check.require(len(shell) == 1, 'unexpected additional shell invocation')
    lines = shell[0].splitlines()
    check.require(lines[0] == 'python -I -B - "$PWD" "$PAYLOAD_SHA256" "$HELPER_SHA256" "$REPLAY_PINS_SHA256" --verify-only <<\'PY\'', 'bootstrap entry/pin argv differs')
    check.require(lines[-1] == 'PY', 'bootstrap here-document terminator differs')
    code = '\n'.join(lines[1:-1]) + '\n'
    parsed = ast.parse(code, filename='README captured-byte stdlib bootstrap')
    imports = []
    for node in ast.walk(parsed):
        if isinstance(node, ast.Import):
            imports.extend(a.name for a in node.names)
        if isinstance(node, ast.ImportFrom):
            imports.append(node.module)
    check.require(set(imports) == {'hashlib', 'os', 'pathlib', 're', 'stat', 'sys'}, 'bootstrap imports non-stdlib or unexpected modules')
    exec_calls = [n for n in ast.walk(parsed) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == 'exec']
    check.require(len(exec_calls) == 1, 'bootstrap must execute helper exactly once')
    last = parsed.body[-1]
    check.require(isinstance(last, ast.Expr) and last.value is exec_calls[0], 'helper execution is not final after all checks')
    call = exec_calls[0]
    check.require(len(call.args) == 2 and isinstance(call.args[0], ast.Call) and isinstance(call.args[0].func, ast.Name) and call.args[0].func.id == 'compile', 'bootstrap does not compile captured code')
    compile_call = call.args[0]
    check.require(ast.unparse(compile_call.args[0]) == "captured['replay_checkpoint.py']", 'bootstrap execution reloads helper instead of captured bytes')
    check.require(ast.literal_eval(compile_call.args[2]) == 'exec', 'bootstrap compile mode differs')
    check.require(ast.unparse(call.args[1]) == "{'__name__': '__main__', '__file__': entry}", 'helper runtime entry context differs')
    for required in [
        'controls = dict(zip(("PAYLOAD_MANIFEST.json", "replay_checkpoint.py", "REPLAY_PINS.json"), sys.argv[2:5]))',
        'for name, expected in controls.items():',
        're.fullmatch(r"[0-9a-f]{64}", expected)',
        'os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK',
        'not stat.S_ISREG(before.st_mode) or before.st_nlink != 1',
        'before.st_size > 8 * 1024 * 1024',
        'stream.read(8 * 1024 * 1024 + 1)',
        'before.st_size, before.st_mtime_ns, before.st_ctime_ns',
        'after.st_size, after.st_mtime_ns, after.st_ctime_ns',
        'len(raw) != before.st_size or hashlib.sha256(raw).hexdigest() != expected',
        'captured[name] = raw',
        '*sys.argv[5:]',
    ]:
        check.require(required in code, 'required bootstrap capture/check/argv structure missing: ' + required)
    check.require(code.index('captured[name] = raw') < code.index('exec(compile('), 'helper executes before complete capture loop')
    check.require('replacing its trailing `--verify-only` argument before the here-document with `--output-root /absolute/external/new-replay-directory`' in text, 'fresh replay does not reuse captured-byte bootstrap')
    check.require('Perform the bootstrap checks before every invocation' in text, 'bootstrap is not required for each invocation')
    for phrase in ['a0(eta)^4 / epsilon_represented', 'a0(eta)=L(eta)=-1/eta', 'mu_j=dk_j*k_j^2/(2*Pi_represented^2)', 'eta in [-4.5,-3.5]', 'UNRESOLVED', '**FAIL**', 'NOT_ESTABLISHED', 'NOT_ASSESSED']:
        check.require(phrase in text, 'accepted normalization/scope text lost: ' + phrase)
    scientific = check.subset(files, check.C_REL)
    source, _ = check.tree(check.SOURCE_C)
    check.require(scientific == source and len(scientific) == 161, 'unchanged final C differs from current source')
    check.require(not (P / 'PAYLOAD_MANIFEST.json').exists() and not (P / 'REPLAY_PINS.json').exists(), 'final controls present before root seal')
    check.require(check.opaque_pin(check.STAGE) == before, 'bootstrap inventory changed during review')
    again, again_dirs = check.tree(P)
    check.require(again == files and again_dirs == dirs, 'bootstrap package changed during review')
    snapshot = {'schema_version': 1, 'status': 'INDEPENDENT_BOOTSTRAP_FINAL_UNSEALED_SNAPSHOT', 'sealed': False, 'package_path': str(P), 'files': files, 'directories': dirs, **check.summary(files)}
    check.json_write(OUT / 'OBSERVED_BOOTSTRAP_FINAL_SNAPSHOT.json', snapshot)
    report = {
        'schema_version': 1,
        'status': 'PASS_BOOTSTRAP_FINAL_UNSEALED_PACKAGE_ROOT_SEAL_READY',
        'sealed': False,
        'utc': datetime.now(timezone.utc).isoformat(),
        'package_path': str(P),
        'staging_inventory': before,
        'prior_final_stage_receipt': check.opaque_pin(OUT / 'FINAL_STAGE_REVIEW_RECEIPT.json'),
        'observed_snapshot': check.opaque_pin(OUT / 'OBSERVED_BOOTSTRAP_FINAL_SNAPSHOT.json'),
        **check.summary(files),
        'modified_files': changed,
        'added_files': [],
        'removed_files': [],
        'all_other_238_leaf_hashes_reverified': True,
        'README': README_PIN,
        'README_bootstrap_AST': 'PASS_STATIC_ONLY',
        'capture_semantics': 'ALL_THREE_STANDARD_LIBRARY_CAPTURE_AND_EXTERNAL_HASH_CHECKS_COMPLETE_BEFORE_SINGLE_EXEC_OF_CAPTURED_AUTHENTICATED_HELPER_BYTES',
        'no_helper_code_second_load': True,
        'generic_pin_argument_and_trailing_option_forwarding': 'PASS_STATIC',
        'fresh_mode_requires_same_bootstrap': True,
        'no_concrete_control_or_self_SHA256_in_README': True,
        'complete_C_current_source_identity': check.summary(scientific),
        'previous_117_frozen_and_11_dependency_pins': 'UNCHANGED_REAUTHENTICATED_OPAQUE_PACKAGE_BYTES_FROM_PRIOR_FINAL_REVIEW',
        'safe_regular_single_link_leaf_files': len(files),
        'unsafe_links_specials_caches_or_empty_dirs': 0,
        'final_payload_manifest_exists': False,
        'final_replay_pins_exists': False,
        'root_seal_readiness': 'READY_FOR_ROOT_TO_ADD_AND_EXTERNALLY_PIN_FINAL_CONTROLS',
        'bootstrap_or_helper_execution_by_this_reviewer': 0,
        'retained_array_decodes': 0,
        'scientific_output_numeric_reads': 0,
        'study_imports': 0,
        'source_callbacks': 0,
        'target_evaluations': 0,
        'network_calls': 0,
        'repository_writes': 0,
        'package_writes': 0,
        'fresh_portable_replay': 'NOT_RUN_BY_THIS_REVIEWER',
        'publication': 'NOT_PERFORMED_BY_THIS_REVIEWER',
        'review_script': check.opaque_pin(Path(__file__)),
        'shared_original_reviewer_script': check.opaque_pin(OUT / 'review_unsealed_stage.py'),
        'limitations': ['Static bootstrap inspection and opaque byte authentication do not claim actual bootstrap execution, a completed seal, a fresh portable replay or publication.', 'Final manifest/replay controls will increase the package leaf count and require root authentication under external pins.'],
    }
    check.json_write(OUT / 'BOOTSTRAP_FINAL_STAGE_REVIEW_RECEIPT.json', report)
    print(json.dumps({'status': report['status'], 'files': len(files), 'bytes': report['total_leaf_bytes'], 'changes': changed, 'receipt': check.opaque_pin(OUT / 'BOOTSTRAP_FINAL_STAGE_REVIEW_RECEIPT.json')}, sort_keys=True))


if __name__ == '__main__':
    main()
