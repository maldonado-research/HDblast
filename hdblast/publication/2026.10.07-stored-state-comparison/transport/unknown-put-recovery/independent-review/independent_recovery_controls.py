#!/usr/bin/env python3
"""Manufactured controls only. Never invokes helper main or a live transport."""
import argparse
import copy
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import socket
import tempfile
import threading
import types
from unittest.mock import patch


def check(ok, reason):
    if not ok:
        raise AssertionError(reason)


def deny_network(*args, **kwargs):
    raise AssertionError('NETWORK_FORBIDDEN_MANUFACTURED_CONTROLS')


def load(source):
    spec = importlib.util.spec_from_file_location('_manufactured_put_recovery', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def state(m):
    return {'draft_id': m.DRAFT, 'published': False,
            'continuation_binding': {'controller_sha256': m.CONTROLLER_SHA},
            'inventory_sha256': m.INVENTORY_SHA, 'metadata_sha256': m.METADATA_SHA,
            'create_post_attempted': True, 'publish_post_attempted': False,
            'content_verified_files': [],
            'pending': {'kind': 'put_content', 'url': m.URL,
                        'body_sha256': m.ASSET_PIN['sha256'],
                        'body_bytes': m.ASSET_PIN['bytes'],
                        'context': {'filename': m.NAME}, 'utc': 'manufactured-time'}}


def stopped(m, call):
    try:
        call()
    except (m.Stop, OSError, ValueError):
        return
    raise AssertionError('invalid manufactured input was accepted')


class FakeModule:
    @staticmethod
    def immutable_file_identity(f):
        return (f.get('file_id'), f.get('version_id'))

    @staticmethod
    def assert_pins(files, pins, exact=False):
        for p in pins:
            check(files[p['filename']]['size'] == p['bytes'], 'inherited size')
            check(files[p['filename']]['checksum'] == 'md5:' + p['md5'], 'inherited checksum')


def remote(m):
    old = [{'filename': 'manufactured-old-a', 'bytes': 3, 'md5': 'a'*32},
           {'filename': 'manufactured-old-b', 'bytes': 4, 'md5': 'b'*32}]
    files = {p['filename']: {'key': p['filename'], 'size': p['bytes'],
             'checksum': 'md5:' + p['md5'], 'bucket_id': m.BUCKET_ID,
             'file_id': 'fake-file-' + str(i), 'version_id': 'fake-version-' + str(i)}
             for i, p in enumerate(old)}
    files[m.NAME] = {'key': m.NAME, 'status': 'pending', 'size': None, 'checksum': None,
                     'file_id': m.FILE_ID, 'version_id': m.VERSION_ID, 'bucket_id': m.BUCKET_ID,
                     'links': {'content': m.URL, 'commit': m.URL.removesuffix('/content') + '/commit'}}
    return types.SimpleNamespace(module=FakeModule, old=old), files, copy.deepcopy(files)


class Response:
    def __init__(self, m, status=400, payload=None, truncated=False, failed=False):
        self.status = status
        self.body_truncated = truncated
        self.body_read_failed = failed
        self.payload = ({'status': 400, 'message': 'File with key "' + m.NAME + '" is not available.'}
                        if payload is None else payload)
    def json(self):
        return self.payload


class TokenInput(io.BytesIO):
    def __init__(self, process):
        super().__init__()
        self.process = process
    def close(self):
        if not self.closed:
            self.process.config = self.getvalue()
            self.process.thread.start()
        super().close()


class FakeProcess:
    """Actual anonymous pipes, manufactured bytes; never launches a program."""
    def __init__(self, argv, kwargs, body, stderr, headers, calls):
        calls.append({'argv': list(argv), 'kwargs': kwargs, 'process': self})
        self.config = None
        self.returncode = None
        self.kills = 0
        self.write_fds = []
        out_r, out_w = os.pipe()
        err_r, err_w = os.pipe()
        self.stdout = os.fdopen(out_r, 'rb', buffering=0)
        self.stderr = os.fdopen(err_r, 'rb', buffering=0)
        header_w = os.dup(kwargs['pass_fds'][1])
        self.write_fds = [out_w, err_w, header_w]
        def produce():
            try:
                for fd, raw in zip(self.write_fds, [body, stderr, headers]):
                    view = memoryview(raw)
                    while view:
                        try:
                            n = os.write(fd, view[:65536])
                        except OSError:
                            break
                        view = view[n:]
            finally:
                for fd in self.write_fds:
                    try:
                        os.close(fd)
                    except OSError:
                        pass
                if self.returncode is None:
                    self.returncode = 0
        self.thread = threading.Thread(target=produce, daemon=True)
        self.stdin = TokenInput(self)
    def kill(self):
        self.kills += 1
        self.returncode = -9
        for fd in self.write_fds:
            try:
                os.close(fd)
            except OSError:
                pass
    def poll(self):
        return self.returncode
    def wait(self, timeout=None):
        self.thread.join(timeout)
        check(not self.thread.is_alive(), 'fake collector did not terminate')
        return self.returncode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = []
    def test(name, call):
        try:
            call()
        except Exception as error:
            rows.append({'name': name, 'status': 'FAIL', 'exception': type(error).__name__, 'reason': str(error)})
        else:
            rows.append({'name': name, 'status': 'PASS'})
    with patch.object(socket, 'socket', deny_network), patch.object(socket, 'create_connection', deny_network):
        m = load(args.source)
        test('valid_exact_state', lambda: m.validate_state(state(m)))
        for field, value in [('draft_id', 'wrong'), ('published', True), ('inventory_sha256', 'wrong'),
                             ('metadata_sha256', 'wrong'), ('create_post_attempted', False),
                             ('publish_post_attempted', True), ('content_verified_files', [m.NAME])]:
            def changed(field=field, value=value):
                s = state(m); s[field] = value
                stopped(m, lambda: m.validate_state(s))
            test('state_reject_' + field, changed)
        test('state_reject_controller_pin', lambda: stopped(m, lambda: m.validate_state({**state(m), 'continuation_binding': {'controller_sha256': 'wrong'}})))
        for field, value in [('kind', 'post'), ('url', m.URL + '?wrong=1'), ('body_sha256', 'wrong'),
                             ('body_bytes', m.ASSET_PIN['bytes'] + 1), ('body_bytes', True),
                             ('context', {'filename': 'wrong'})]:
            def pending_changed(field=field, value=value):
                s = state(m); s['pending'][field] = value
                stopped(m, lambda: m.validate_state(s))
            test('pending_reject_' + field + '_' + type(value).__name__, pending_changed)
        test('pending_reject_extra_field', lambda: stopped(m, lambda: m.validate_state({**state(m), 'pending': {**state(m)['pending'], 'extra': 1}})))
        for value in [None, 0, '']:
            test('state_strict_publish_false_' + repr(value), lambda value=value: stopped(m, lambda: m.validate_state({**state(m), 'publish_post_attempted': value})))
        def missing_publish():
            s = state(m); del s['publish_post_attempted']; m.validate_state(s)
        test('state_accepts_original_absent_publish_flag', missing_publish)
        test('valid_remote_inventory', lambda: m.validate_remote(*remote(m)))
        for key, value in [('file_id', 'wrong'), ('version_id', 'wrong'), ('bucket_id', 'wrong'),
                            ('status', 'completed'), ('size', 0), ('checksum', 'md5:' + 'c'*32), ('key', 'wrong')]:
            def changed_remote(key=key, value=value):
                c, files, initialized = remote(m); files[m.NAME][key] = value
                stopped(m, lambda: m.validate_remote(c, files, initialized))
            test('remote_new_reject_' + key, changed_remote)
        def bad_links():
            c, files, initialized = remote(m); files[m.NAME]['links']['content'] += '?wrong=1'
            stopped(m, lambda: m.validate_remote(c, files, initialized))
        test('remote_reject_destination_override', bad_links)
        for key in ['file_id', 'version_id', 'bucket_id']:
            def old_changed(key=key):
                c, files, initialized = remote(m); files['manufactured-old-a'][key] = 'wrong'
                stopped(m, lambda: m.validate_remote(c, files, initialized))
            test('remote_inherited_reject_' + key, old_changed)
        def extra_remote():
            c, files, initialized = remote(m); files['extra'] = {}
            stopped(m, lambda: m.validate_remote(c, files, initialized))
        test('remote_reject_extra_file', extra_remote)
        def changed_initialization():
            c, files, initialized = remote(m); initialized[m.NAME]['file_id'] = 'different-initialized-file'
            stopped(m, lambda: m.validate_remote(c, files, initialized))
        test('remote_bind_new_initialized_identity', changed_initialization)
        test('valid_exact_missing_body', lambda: m.validate_missing(Response(m)))
        for label, kw in [('status', {'status': 404}), ('truncated', {'truncated': True}), ('read_failed', {'failed': True}),
                          ('extra_body_field', {'payload': {**Response(m).payload, 'extra': True}}),
                          ('near_match', {'payload': {'status': 400, 'message': Response(m).payload['message'] + ' '}})]:
            test('missing_reject_' + label, lambda kw=kw: stopped(m, lambda: m.validate_missing(Response(m, **kw))))
        test('missing_strict_integer_status', lambda: stopped(m, lambda: m.validate_missing(Response(m, status=400.0))))
        with tempfile.TemporaryDirectory(prefix='put-recovery-fabricated-') as temporary:
            base = Path(temporary)
            evidence = base / 'manufactured-pinned.json'
            evidence.write_bytes(b'{"manufactured":true}\n')
            evidence_sha = hashlib.sha256(evidence.read_bytes()).hexdigest()
            test('read_pinned_exact_complete_hash', lambda: check(m.read_pinned(evidence, evidence_sha) == evidence.read_bytes(), 'pinned bytes differ'))
            test('read_pinned_reject_hash_mismatch', lambda: stopped(m, lambda: m.read_pinned(evidence, '0'*64)))
            evidence_link = base / 'pinned-symlink'; evidence_link.symlink_to(evidence)
            test('read_pinned_reject_symlink', lambda: stopped(m, lambda: m.read_pinned(evidence_link, evidence_sha)))
            evidence_hard = base / 'pinned-hardlink'; os.link(evidence, evidence_hard)
            test('read_pinned_reject_hardlink', lambda: stopped(m, lambda: m.read_pinned(evidence, evidence_sha)))
            evidence_hard.unlink()
            def escaped_json_redaction():
                token = 'MANUFACTURED_"\\_CREDENTIAL'
                path = base / 'reflected.json'
                escaped = json.dumps(token)[1:-1]
                m.write_json(path, {'echo': token, 'nested': [token], 'escaped_echo': escaped, token: 'key'}, token)
                parsed = json.loads(path.read_bytes())
                check(token not in repr(parsed) and escaped not in repr(parsed), 'credential survived JSON redaction')
                check(parsed['echo'] == '[REDACTED]' and parsed['nested'] == ['[REDACTED]'], 'recursive redaction failed')
            test('write_json_redacts_escaped_recursive_token', escaped_json_redaction)
            raw = b'fabricated immutable asset\x00\xff' * 30
            asset = base / 'asset.bin'; asset.write_bytes(raw)
            pin = {'filename': 'asset.bin', 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest(), 'md5': hashlib.md5(raw).hexdigest()}
            def capture_valid():
                fd, receipt = m.capture_asset(asset, pin)
                try:
                    check(os.read(fd, len(raw)+1) == raw, 'captured bytes differ')
                    check(receipt['sha256'] == pin['sha256'] and receipt['md5'] == pin['md5'], 'capture digests differ')
                    stopped(m, lambda: os.write(fd, b'x'))
                    stopped(m, lambda: os.ftruncate(fd, 0))
                    check(receipt['stat_identity'] == list(m.stat_identity(asset.stat())), 'capture stat identity differs')
                finally:
                    os.close(fd)
            test('capture_exact_hashes_stat_and_seals', capture_valid)
            for field, value in [('bytes', len(raw)+1), ('sha256', '0'*64), ('md5', '0'*32)]:
                test('capture_reject_' + field, lambda field=field, value=value: stopped(m, lambda: m.capture_asset(asset, {**pin, field: value})))
            symlink = base / 'asset-link'; symlink.symlink_to(asset)
            test('capture_reject_symlink', lambda: stopped(m, lambda: m.capture_asset(symlink, pin)))
            hardlink = base / 'asset-hard'; os.link(asset, hardlink)
            test('capture_reject_hardlink', lambda: stopped(m, lambda: m.capture_asset(asset, pin)))
            hardlink.unlink()
            def snapshot_survives_path_change():
                fd, receipt = m.capture_asset(asset, pin)
                try:
                    asset.write_bytes(b'replaced after valid capture')
                    check(os.read(fd, len(raw)+1) == raw, 'later path changed captured upload')
                finally:
                    os.close(fd); asset.write_bytes(raw)
            test('capture_preserves_original_after_path_change', snapshot_survives_path_change)
            def partial_capture_write():
                real_write = os.write
                with patch.object(m.os, 'write', lambda fd, view: real_write(fd, view[:3])):
                    fd, receipt = m.capture_asset(asset, pin)
                try:
                    check(os.read(fd, len(raw)+1) == raw, 'partial writes lost bytes')
                finally:
                    os.close(fd)
            test('capture_handles_partial_memfd_writes', partial_capture_write)
            def during_change():
                real_read = os.read
                changed = False
                def read(fd, size):
                    nonlocal changed
                    chunk = real_read(fd, size)
                    if chunk and not changed:
                        changed = True
                        s = asset.stat(); os.utime(asset, ns=(s.st_atime_ns, s.st_mtime_ns + 1))
                    return chunk
                with patch.object(m.os, 'read', read):
                    stopped(m, lambda: m.capture_asset(asset, pin))
            test('capture_rejects_changed_stat_during_read', during_change)
            def during_path_change():
                real_read = os.read
                changed = False
                def read(fd, size):
                    nonlocal changed
                    chunk = real_read(fd, size)
                    if chunk and not changed:
                        changed = True
                        replacement = base / 'replacement.bin'; replacement.write_bytes(raw)
                        os.replace(replacement, asset)
                    return chunk
                with patch.object(m.os, 'read', read):
                    stopped(m, lambda: m.capture_asset(asset, pin))
            test('capture_rejects_path_replaced_during_read', during_path_change)
            def latch_repeat():
                side = base / 'latch-test'; side.mkdir()
                m.durable_latch(side, {'manufactured': True})
                first = (side / 'ONE_PUT_ATTEMPT_LATCH.json').read_bytes()
                stopped(m, lambda: m.durable_latch(side, {'manufactured': False}))
                check((side / 'ONE_PUT_ATTEMPT_LATCH.json').read_bytes() == first, 'persistent latch overwritten')
            test('durable_exclusive_latch_refuses_repeat', latch_repeat)
            def concurrent_latch():
                side = base / 'concurrent-latch'; side.mkdir()
                statuses = []
                def reserve():
                    try:
                        m.durable_latch(side, {'manufactured': True})
                    except OSError:
                        statuses.append('refused')
                    else:
                        statuses.append('reserved')
                threads = [threading.Thread(target=reserve) for _ in range(2)]
                for thread in threads: thread.start()
                for thread in threads: thread.join(2)
                check(sorted(statuses) == ['refused', 'reserved'], 'concurrent latch allowed multiple attempts')
            test('durable_latch_concurrent_single_reservation', concurrent_latch)
            for kind in ['symlink', 'directory']:
                def bad_latch(kind=kind):
                    side = base / ('latch-' + kind); side.mkdir()
                    latch = side / 'ONE_PUT_ATTEMPT_LATCH.json'
                    if kind == 'symlink': latch.symlink_to(asset)
                    else: latch.mkdir()
                    stopped(m, lambda: m.durable_latch(side, {'manufactured': True}))
                test('durable_latch_reject_' + kind, bad_latch)
        dummy = 'MANUFACTURED_TOKEN_NEVER_A_REAL_CREDENTIAL'
        def transport_command():
            command = m.curl_command(91, 92)
            check(command[1] == '-q', 'curl default config enabled')
            check(command[-1] == m.URL and command.count(m.URL) == 1, 'curl destination differs')
            check(dummy not in repr(command), 'token in argv')
            for option, value in [('--request', 'PUT'), ('--config', '-'), ('--retry', '0'), ('--max-redirs', '0'), ('--proto', '=https'), ('--proto-redir', '=https')]:
                check(command[command.index(option)+1] == value, 'transport option differs: ' + option)
            check('--http1.1' in command and 'Expect:' in command and not any('100-continue' in x for x in command), 'Expect fallback enabled')
            check(not set(command) & {'-L', '--location', '-k', '--insecure', '--anyauth', '--retry-all-errors'}, 'unsafe transport option')
            check('/proc/self/fd/91' in command, 'captured descriptor not upload source')
        test('curl_command_fixed_no_repeat_no_redirect_verified_tls', transport_command)
        test('token_only_stdin_config', lambda: check(dummy.encode() in m.token_config(dummy), 'config token missing'))
        for bad in ['', 'line\nbreak', 'line\rbreak', 'NUL\0', 'non-ascii-é']:
            test('token_reject_' + repr(bad), lambda bad=bad: stopped(m, lambda: m.token_config(bad)))
        metric = b'\nHDBLAST_CURL_METRICS|400|900|37|0.1|0|0\n'
        test('metrics_valid_status', lambda: check(m.parse_metrics(metric)['http_status'] == 400, 'status parse'))
        test('metrics_reject_absent', lambda: stopped(m, lambda: m.parse_metrics(b'400 in ordinary error text')))
        test('metrics_reject_trailing_untrusted_bytes', lambda: stopped(m, lambda: m.parse_metrics(metric + b'stale')))
        test('metrics_reject_conflicting_lines', lambda: stopped(m, lambda: m.parse_metrics(metric + metric.replace(b'|400|', b'|200|'))))
        test('metrics_reject_malformed_duration', lambda: stopped(m, lambda: m.parse_metrics(metric.replace(b'|0.1|', b'|1..2|'))))
        test('metrics_zero_records_no_http_outcome', lambda: check(m.parse_metrics(metric.replace(b'|400|', b'|000|'))['http_status'] == 0, 'curl no-HTTP outcome not preserved'))
        test('metrics_reject_out_of_range_status', lambda: stopped(m, lambda: m.parse_metrics(metric.replace(b'|400|', b'|999|'))))
        def collector(body, stderr, headers, overflow=False):
            calls = []
            def popen(command, **kwargs):
                return FakeProcess(command, kwargs, body, stderr, headers, calls)
            with patch.object(m.subprocess, 'Popen', popen):
                result, outputs = m.run_curl(99, dummy)
            check(len(calls) == 1, 'collector launched repeated process')
            check(calls[0]['process'].config == m.token_config(dummy), 'token not sent only via config')
            check(dummy not in repr(calls[0]['argv']), 'token in captured argv')
            check(all(dummy.encode() not in value for value in outputs.values()), 'token reflected in output')
            check(all(len(value) <= m.MAX_JSON for value in outputs.values()), 'response buffer unbounded')
            check(result['response_bound_exceeded'] is overflow, 'collector bound outcome differs')
            if overflow: check(calls[0]['process'].kills > 0, 'overflow did not kill process')
            else: check(result['metrics']['http_status'] == 400, 'collector status differs')
            for stream in [calls[0]['process'].stdout, calls[0]['process'].stderr]: stream.close()
        test('collector_one_process_redacts_three_channels', lambda: collector(dummy.encode(), dummy.encode() + metric, b'HTTP/1.1 400 Bad Request\r\n' + dummy.encode()))
        test('collector_body_bound_kills_without_repeat', lambda: collector(b'x'*(m.MAX_JSON+1), metric, b'HTTP/1.1 400 Bad Request\r\n', True))
        test('collector_header_bound_kills_without_repeat', lambda: collector(b'{}', metric, b'x'*(m.MAX_JSON+1), True))
        test('collector_stderr_bound_kills_without_repeat', lambda: collector(b'{}', metric + b'x'*(m.MAX_JSON+1), b'HTTP/1.1 400 Bad Request\r\n', True))
        def manufactured_child_environment():
            calls = []
            env = {'ZENODO_ACCESS_TOKEN': dummy, 'SSLKEYLOGFILE': '/manufactured/keylog',
                   'HTTPS_PROXY': 'http://manufactured-proxy.invalid',
                   'CURL_CA_BUNDLE': '/manufactured/ca-bundle'}
            def popen(command, **kwargs):
                return FakeProcess(command, kwargs, b'{}', metric, b'HTTP/1.1 400 Bad Request\r\n', calls)
            with patch.object(m.os, 'environ', env), patch.object(m.subprocess, 'Popen', popen):
                result, outputs = m.run_curl(99, dummy)
            child = calls[0]['kwargs']['env']
            check('ZENODO_ACCESS_TOKEN' not in child and 'SSLKEYLOGFILE' not in child, 'credential or TLS keylog environment forwarded')
            check(child == {k: value for k, value in env.items() if k not in {'ZENODO_ACCESS_TOKEN', 'SSLKEYLOGFILE'}}, 'proxy/CA environment unexpectedly changed')
            check(dummy not in repr(child) and dummy not in repr(calls[0]['argv']), 'credential outside stdin config')
            calls[0]['process'].stdout.close(); calls[0]['process'].stderr.close()
        test('collector_child_environment_token_and_keylog_removed', manufactured_child_environment)
        def bounded_kill_drain():
            calls = []
            def popen(command, **kwargs):
                process = FakeProcess(command, kwargs, b'', b'', b'', calls)
                process.thread = threading.Thread(target=lambda: None, daemon=True)
                def kill_without_pipe_eof():
                    process.kills += 1
                    process.returncode = -9
                process.kill = kill_without_pipe_eof
                return process
            times = iter([0.0, 191.0, 197.0])
            try:
                with patch.object(m.subprocess, 'Popen', popen), patch.object(m.time, 'monotonic', lambda: next(times)):
                    result, outputs = m.run_curl(99, dummy)
                check(result['wall_deadline_exceeded'] is True, 'wall deadline not recorded')
                check(result.get('pipe_drain_deadline_exceeded') is True, 'held pipe drain unbounded')
                check(len(calls) == 1 and calls[0]['process'].kills == 1, 'kill or process repeated')
            finally:
                for call in calls:
                    process = call['process']
                    for fd in process.write_fds:
                        try: os.close(fd)
                        except OSError: pass
                    process.stdout.close(); process.stderr.close()
        test('collector_wall_and_held_pipe_drain_are_bounded', bounded_kill_drain)
    result = {'status': 'PASS' if all(r['status'] == 'PASS' for r in rows) else 'FAIL',
              'manufactured_only': True, 'authorization_receipt': False, 'main_invoked': False,
              'live_transport_calls': 0, 'source_path': str(args.source),
              'source_sha256': hashlib.sha256(args.source.read_bytes()).hexdigest(),
              'optimized': not __debug__, 'controls': rows,
              'pass': sum(r['status'] == 'PASS' for r in rows), 'fail': sum(r['status'] == 'FAIL' for r in rows)}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'source_sha256', 'optimized', 'pass', 'fail']}))


if __name__ == '__main__':
    main()
