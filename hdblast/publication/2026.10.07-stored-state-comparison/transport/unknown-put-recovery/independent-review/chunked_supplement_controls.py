#!/usr/bin/env python3
"""Thin manufactured supplement; accepted eighty-control source stays unchanged."""
import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import socket
import tempfile
from unittest.mock import patch


def check(ok, reason):
    if not ok:
        raise AssertionError(reason)


def deny(*args, **kwargs):
    raise AssertionError('NETWORK_FORBIDDEN_MANUFACTURED_CONTROLS')


def load(path):
    spec = importlib.util.spec_from_file_location('_manufactured_chunked_recovery', path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def reject(m, fn):
    try:
        fn()
    except (m.Stop, OSError, ValueError):
        return
    raise AssertionError('invalid manufactured prior evidence accepted')


def encoded(obj):
    return (json.dumps(obj, sort_keys=True) + '\n').encode()


def fixture(m, base):
    prior = base / 'nonexecutable-manufactured-prior-source.txt'
    prior.write_bytes(b'MANUFACTURED_NOT_A_REAL_RECOVERY_SOURCE_OR_AUTHORIZATION\n')
    folder = base / 'manufactured-failed-attempt'
    folder.mkdir()
    intent = {'operation': 'ONE_IDENTICAL_PUT_RECOVERY', 'url': m.URL,
              'asset': copy.deepcopy(m.ASSET_PIN), 'controller_sha256': m.CONTROLLER_SHA}
    result = {'original_state_journal_unchanged': True, 'pending_cleared': False,
              'commits': 0, 'publishes': 0,
              'result': {'metrics': {'http_status': 401, 'uploaded_bytes': 24051375,
                                    'ssl_verify_result': 0, 'redirects': 0}}}
    raws = {'ONE_PUT_ATTEMPT_LATCH.json': encoded(intent), 'RECOVERY_RESULT.json': encoded(result),
            'curl_body.raw': b'Unauthorized: authentication failed for integration codex-secret-zenodo.org\n',
            'curl_headers.raw': b'HTTP/1.1 401 Unauthorized\r\nx-at-upstream-error: false\r\n\r\n',
            'curl_stderr.raw': b'MANUFACTURED_STDERR_NOT_ACTUAL_EVIDENCE\n'}
    for name, raw in raws.items():
        (folder/name).write_bytes(raw)
    return prior, folder, raws


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    rows = []
    def test(name, call):
        try:
            call()
        except Exception as error:
            rows.append({'name': name, 'status': 'FAIL', 'exception': type(error).__name__, 'reason': str(error)})
        else:
            rows.append({'name': name, 'status': 'PASS'})
    with patch.object(socket, 'socket', deny), patch.object(socket, 'create_connection', deny):
        m = load(args.source)
        def transport_headers():
            argv = m.curl_command(51, 52)
            headers = [argv[i+1] for i, value in enumerate(argv) if value == '--header']
            check(headers.count('Content-Length:') == 1, 'Content-Length not suppressed exactly once')
            check(headers.count('Transfer-Encoding: chunked') == 1, 'chunked header not exact/unique')
            check(not any(x.lower().startswith('content-length:') and x != 'Content-Length:' for x in headers), 'fixed length emitted')
            check(headers.count('Expect:') == 1 and not any('100-continue' in x.lower() for x in headers), 'Expect replay route')
            check(argv.count(m.URL) == 1 and argv[-1] == m.URL, 'multiple or changed URLs')
            check('/proc/self/fd/51' in argv and '--upload-file' in argv, 'whole captured descriptor not used')
            check(not set(argv) & {'--range', '--continue-at', '-C', '--data', '--data-binary'}, 'partial/resume or alternate body enabled')
        test('exact_empty_length_chunked_expect_whole_asset_headers', transport_headers)
        test('production_asset_identity_unchanged', lambda: check(m.ASSET_PIN == {'filename': 'HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON.zip', 'bytes': 240135519, 'sha256': 'ead5ffa3987c760d99da7a79189ff5ad766fb9bef6e032c489f4efb2316cd189', 'md5': 'e424502f12dbaacd36abaaabdd3fb74e'}, 'original sealed asset pin changed'))
        test('original_state_and_controller_pins_unchanged', lambda: check(m.ORIGINAL_STATE_SHA == '9927d046277e4c808e43fff49d3c6add2495df8795cff7b166941694830132ad' and m.CONTROLLER_SHA == 'aa7e6cfe271559731ed3b45c5339f8eb612ebe6f017ddafa5fdadad9579e084b', 'original guard pins changed'))
        test('updated_journal_and_fresh_fixed_side', lambda: check(m.ORIGINAL_JOURNAL_SHA == 'ad2b052940b711225dd7af7ef475b362f4cd6fdf866d429e89c5434ebc2b677a' and m.SIDE.name == 'live-chunked-attempt-001', 'journal or isolated side differs'))
        test('accepted_prior_source_pin_and_exact_five_evidence_names', lambda: check(m.PRIOR_RECOVERY_SHA == '39ec701ae3133b699bbad1f74031c549a0aedbc2116109cf5be199d3c7b0916d' and set(m.PRIOR_FAILURE_PINS) == {'ONE_PUT_ATTEMPT_LATCH.json', 'RECOVERY_RESULT.json', 'curl_body.raw', 'curl_headers.raw', 'curl_stderr.raw'}, 'prior source or closure differs'))
        def prior_case(label, mutate=None, repin=False):
            with tempfile.TemporaryDirectory(prefix='manufactured-chunked-prior-') as temporary:
                prior, folder, raws = fixture(m, Path(temporary))
                pins = {name: hashlib.sha256(raw).hexdigest() for name, raw in raws.items()}
                prior_sha = hashlib.sha256(prior.read_bytes()).hexdigest()
                if mutate:
                    mutate(prior, folder, raws)
                    if repin:
                        pins = {name: hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in raws}
                with patch.object(m, 'PRIOR_RECOVERY', prior), patch.object(m, 'PRIOR_RECOVERY_SHA', prior_sha), patch.object(m, 'PRIOR_FAILURE', folder), patch.object(m, 'PRIOR_FAILURE_PINS', pins):
                    if label == 'valid':
                        proof = m.verify_prior_failure()
                        check(proof['prior_http_status'] == 401 and proof['prior_uploaded_bytes'] == 24051375 and proof['prior_failure_preserved'] is True, 'prior proof differs')
                        check(proof['failure_pins'] == pins and proof['prior_source_sha256'] == prior_sha, 'prior proof pins not linked')
                    else:
                        reject(m, m.verify_prior_failure)
        test('manufactured_valid_prior_failure_proof', lambda: prior_case('valid'))
        test('reject_stale_prior_source_pin', lambda: prior_case('bad', lambda prior, folder, raws: prior.write_bytes(b'changed source')))
        for name in ['ONE_PUT_ATTEMPT_LATCH.json', 'RECOVERY_RESULT.json', 'curl_body.raw', 'curl_headers.raw', 'curl_stderr.raw']:
            test('reject_changed_raw_pin_' + name, lambda name=name: prior_case('bad', lambda prior, folder, raws: (folder/name).write_bytes(raws[name]+b'changed')))
            test('reject_missing_prior_leaf_' + name, lambda name=name: prior_case('bad', lambda prior, folder, raws: (folder/name).unlink()))
        def alter_json(name, field, value):
            def mutate(prior, folder, raws):
                obj = json.loads(raws[name]); obj[field] = value
                (folder/name).write_bytes(encoded(obj))
            return mutate
        for field, value in [('operation', 'another-operation'), ('url', m.URL+'?wrong=1'), ('asset', {**m.ASSET_PIN, 'bytes': m.ASSET_PIN['bytes']-1}), ('controller_sha256', 'wrong')]:
            test('prior_intent_reject_' + field, lambda field=field, value=value: prior_case('bad', alter_json('ONE_PUT_ATTEMPT_LATCH.json', field, value), True))
        for field, value in [('original_state_journal_unchanged', False), ('pending_cleared', True), ('commits', 1), ('commits', False), ('publishes', 1), ('publishes', False)]:
            test('prior_scope_reject_' + field + '_' + type(value).__name__, lambda field=field, value=value: prior_case('bad', alter_json('RECOVERY_RESULT.json', field, value), True))
        def metric_change(field, value):
            def mutate(prior, folder, raws):
                obj=json.loads(raws['RECOVERY_RESULT.json']); obj['result']['metrics'][field]=value
                (folder/'RECOVERY_RESULT.json').write_bytes(encoded(obj))
            return mutate
        for field, value in [('http_status', 200), ('http_status', 401.0), ('uploaded_bytes', 24051374), ('uploaded_bytes', 24051375.0), ('ssl_verify_result', 1), ('redirects', 1)]:
            test('prior_metrics_reject_' + field + '_' + type(value).__name__, lambda field=field, value=value: prior_case('bad', metric_change(field,value), True))
        test('prior_body_reject_near_match_even_resealed', lambda: prior_case('bad', lambda prior, folder, raws: (folder/'curl_body.raw').write_bytes(raws['curl_body.raw']+b' '), True))
        test('prior_headers_reject_missing_upstream_guard_even_resealed', lambda: prior_case('bad', lambda prior, folder, raws: (folder/'curl_headers.raw').write_bytes(b'HTTP/1.1 401 Unauthorized\r\n\r\n'), True))
        test('prior_headers_reject_wrong_status_even_resealed', lambda: prior_case('bad', lambda prior, folder, raws: (folder/'curl_headers.raw').write_bytes(b'HTTP/1.1 200 OK\r\nx-at-upstream-error: false\r\n\r\n'), True))
    r={'status':'PASS' if all(x['status']=='PASS' for x in rows) else 'FAIL', 'optimized':not __debug__, 'source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(), 'manufactured_only':True, 'authorization_receipt':False, 'main_invoked':False, 'live_transport_calls':0, 'real_prior_failure_files_read':False, 'pass':sum(x['status']=='PASS' for x in rows), 'fail':sum(x['status']=='FAIL' for x in rows), 'controls':rows}
    args.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:r[k] for k in ['status','source_sha256','optimized','pass','fail']}))


if __name__ == '__main__':
    main()
