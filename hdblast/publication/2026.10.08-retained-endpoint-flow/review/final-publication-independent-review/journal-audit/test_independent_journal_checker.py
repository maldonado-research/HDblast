#!/usr/bin/env python3
"""Synthetic corruption controls over the already preserved partial journal."""
import contextlib
import copy
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys

OWN = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('independent_local_audit', OWN / 'audit_publication_journal.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
original = (OWN / 'partial-observation-002/JOURNAL.jsonl').read_bytes()
events = [json.loads(x) for x in original.splitlines()]
first_stream = next(i for i, x in enumerate(events) if x['kind'] == 'FULL_CONTENT_STREAM_VERIFIED')
first_upload = next(i for i, x in enumerate(events) if x['kind'] == 'WRITE_INTENT' and x['operation'] == 'put_content')
publish = next(i for i, x in enumerate(events) if x['kind'] == 'WRITE_INTENT' and x['operation'] == 'publish')
cases = []

def case(name, edit):
    changed = copy.deepcopy(events)
    edit(changed)
    cases.append((name, b''.join(json.dumps(x).encode() + b'\n' for x in changed)))

case('wrong-md5', lambda e: e[first_stream]['observed'].__setitem__('md5', '0' * 32))
case('wrong-sha256', lambda e: e[first_stream]['observed'].__setitem__('sha256', '0' * 64))
case('wrong-size', lambda e: e[first_stream]['observed'].__setitem__('bytes', 1))
case('bool-size', lambda e: e[first_stream]['observed'].__setitem__('bytes', True))
case('unknown-stream-name', lambda e: e[first_stream].__setitem__('filename', 'unknown.zip'))
case('duplicate-stream-name', lambda e: e.__setitem__(first_stream + 1, copy.deepcopy(e[first_stream])))
case('wrong-draft-list-endpoint', lambda e: e[first_stream - 1].__setitem__('url', audit.PUBLIC + '/files'))
case('wrong-draft-list-authentication', lambda e: e[first_stream - 1].__setitem__('authenticated_request', False))
case('unscoped-stream', lambda e: e[first_stream - 1].__setitem__('url', audit.ORIGIN + audit.PRIOR))
case('wrong-upload-sha256', lambda e: e[first_upload].__setitem__('body_sha256', '0' * 64))
case('extra-publish-intent', lambda e: e.insert(publish + 1, copy.deepcopy(e[publish])))
case('unknown-write-outcome', lambda e: e[first_upload + 1].__setitem__('kind', 'WRITE_OUTCOME_UNKNOWN'))
case('truncated-write-response', lambda e: e[first_upload + 1].__setitem__('body_truncated', True))
case('utc-order-regression', lambda e: e[first_stream].__setitem__('utc', '2026-10-08T00:00:00+00:00'))
cases.append(('truncated-jsonl-line', original[:-1]))
cases.append(('duplicate-json-key', original.replace(b'"kind": "FULL_CONTENT_STREAM_VERIFIED"',
                                                   b'"kind": "GET", "kind": "FULL_CONTENT_STREAM_VERIFIED"', 1)))
root = OWN / ('controls-optimized' if sys.flags.optimize else 'controls-normal')
root.mkdir()
results = []
for name, raw, expected_pass in [('baseline-partial', original, True)] + [(n, r, False) for n, r in cases]:
    folder = root / name
    folder.mkdir()
    evidence = folder / 'evidence'
    evidence.mkdir()
    (evidence / 'JOURNAL.jsonl').write_bytes(raw)
    (folder / 'JOURNAL_PARTIAL_OBSERVATION_001.jsonl').write_bytes(raw.splitlines(keepends=True)[0])
    audit.OWN, audit.EVIDENCE = folder, evidence
    sys.argv = ['audit_publication_journal.py', '--snapshot', 'observed', '--partial']
    output = io.StringIO()
    error = None
    with contextlib.redirect_stdout(output):
        try:
            audit.main()
        except Exception as exc:
            error = type(exc).__name__ + ': ' + str(exc)
    (folder / 'RESULT.json').write_text(json.dumps({'expected_pass': expected_pass, 'error': error,
                                                'stdout': output.getvalue()}, indent=2) + '\n')
    passed = (error is None) is expected_pass
    results.append({'case': name, 'control_passed': passed, 'observed_error': error})
receipt = {'status': 'PASS' if all(x['control_passed'] for x in results) else 'FAIL',
           'cases': len(results), 'optimized': bool(sys.flags.optimize),
           'checker_sha256': hashlib.sha256((OWN / 'audit_publication_journal.py').read_bytes()).hexdigest(),
           'all_inputs_manufactured_from_saved_partial_journal': True,
           'network_calls': 0, 'remote_mutations': 0, 'attachment_streams': 0, 'results': results}
(root / 'CONTROL_RECEIPT.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
print(json.dumps({'status': receipt['status'], 'cases': receipt['cases'], 'output': str(root / 'CONTROL_RECEIPT.json')}))
raise SystemExit(0 if receipt['status'] == 'PASS' else 1)
