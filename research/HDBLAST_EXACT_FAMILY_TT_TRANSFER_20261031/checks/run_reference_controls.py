#!/usr/bin/env python3
"""Replay real source/receipt corruptions for the bounded reference checker."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def need(value, witness):
    if not value:
        raise RuntimeError(witness)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cap = Path(__file__).resolve().parent.parent
    repo = args.repo_root.resolve(strict=True)
    output = args.output.absolute()
    need(not output.exists() and not output.is_symlink(), 'FRESH_CONTROL_RECEIPT')
    need(output.parent.is_dir() and all(not p.is_symlink() for p in output.parents), 'REGULAR_OUTPUT_PARENT')
    need(repo != output and repo not in output.parents and cap not in output.parents, 'OUTPUT_OUTSIDE_SOURCE')
    script = cap / 'check_reference_identities.py'
    receipt = cap / 'PUBLIC_SOURCE_RECEIPT.json'
    code, receipt_bytes = script.read_bytes(), receipt.read_bytes()
    pins = json.loads(receipt_bytes)
    original = {}
    for row in pins['sources']:
        relative = Path(row['path'])
        need(not relative.is_absolute() and '..' not in relative.parts, 'SAFE_SOURCE_RELATIVE_PATH')
        source = repo / relative
        need(source.is_file() and not source.is_symlink(), 'REGULAR_SOURCE')
        digest = sha(source)
        need(digest == row['sha256'], 'ORIGINAL_SOURCE_PIN')
        original[row['path']] = digest
    cases = []
    with tempfile.TemporaryDirectory(prefix='reference-source-controls-', dir=output.parent) as temp_name:
        temp = Path(temp_name)
        for optimized in (False, True):
            mode = 'optimized' if optimized else 'normal'
            for name in ('actual_source_bytes', 'receipt_expected_hash'):
                case = temp / (mode + '-' + name)
                case.mkdir()
                copied_script = case / script.name
                copied_script.write_bytes(code)
                copied_receipt = case / receipt.name
                copied_receipt.write_bytes(receipt_bytes)
                source_root = repo
                if name == 'actual_source_bytes':
                    source_root = case / 'repo'
                    for row in pins['sources']:
                        relative = Path(row['path'])
                        destination = source_root / relative
                        destination.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copyfile(repo / relative, destination)
                    changed = source_root / pins['sources'][0]['path']
                    data = changed.read_bytes()
                    need(bool(data), 'NONEMPTY_SOURCE_MUTATION')
                    changed.write_bytes(bytes([data[0] ^ 1]) + data[1:])
                else:
                    changed_pins = json.loads(receipt_bytes)
                    digest = changed_pins['sources'][0]['sha256']
                    changed_pins['sources'][0]['sha256'] = ('0' if digest[0] != '0' else '1') + digest[1:]
                    copied_receipt.write_text(json.dumps(changed_pins, sort_keys=True, indent=2) + '\n')
                # Negative runs emit their explicit failure JSON on stdout. A
                # successful report would be status=pass; none may be emitted.
                command = [sys.executable, '-I', '-B', *(['-O'] if optimized else []),
                           str(copied_script), '--repo-root', str(source_root)]
                result = subprocess.run(command, capture_output=True, text=True, timeout=60)
                need(result.returncode == 2, 'EXPECTED_SOURCE_REJECTION_EXIT:' + name)
                document = json.loads(result.stdout)
                need(document['status'] == 'fail' and document['stage'] == 'source_authentication',
                     'EXPLICIT_SOURCE_FAILURE:' + name)
                need(document['sympy_imported'] is False, 'FAILURE_PRECEDES_SYMPY_IMPORT:' + name)
                need(document['error'].startswith('SHA256 mismatch: '), 'ACTUAL_HASH_FAILURE_WITNESS:' + name)
                need('identities' not in document and 'mutant_controls' not in document,
                     'NO_SUCCESSFUL_IDENTITY_REPORT:' + name)
                cases.append({'case': name, 'mode': mode, 'actual_returncode': result.returncode,
                              'actual_witness': document['error'], 'failure_status': document['status'],
                              'failure_stage': document['stage'], 'sympy_imported': document['sympy_imported'],
                              'successful_stdout_report_emitted': document.get('status') == 'pass', 'executed_code_sha256': sha(copied_script)})
    need(script.read_bytes() == code and receipt.read_bytes() == receipt_bytes, 'ORIGINAL_REFERENCE_CODE_UNCHANGED')
    need(all(sha(repo / relative) == digest for relative, digest in original.items()), 'ORIGINAL_REFERENCE_SOURCES_UNCHANGED')
    report = {'schema': 'HDBLAST_REFERENCE_SOURCE_CONTROLS_V1',
              'status': 'PASS_REFERENCE_SOURCE_AUTHENTICATION_CONTROLS',
              'source_authentication_subprocess_rejections': len(cases),
              'original_sources_and_code_unchanged': True, 'controls': cases,
              'scope': 'Real copied-source and copied-receipt corruptions before symbolic imports. The five in-process formula fixtures in each canonical reference report are separate.'}
    with output.open('x') as stream:
        stream.write(json.dumps(report, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'source_authentication_subprocess_rejections': len(cases)}, sort_keys=True))


if __name__ == '__main__':
    main()
