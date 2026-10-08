#!/usr/bin/env python3
"""Socket-free manufactured early-response probe; never imports controller."""
from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
import hashlib
import http.client
import io
import json
import mmap
from pathlib import Path
import re
import socket
import ssl
import sys
import urllib.error as error
import urllib.parse as parse
import urllib.request as request
from unittest.mock import patch

CONTROLLER = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007/publication-planning/new_edition_controller.py')
CONTROLLER_SHA256 = 'aa7e6cfe271559731ed3b45c5339f8eb612ebe6f017ddafa5fdadad9579e084b'
BODY_BYTES = 240_135_519
URL = 'https://zenodo.org/api/offline-manufactured-test/content'
ERROR_BODY = b'{"error":"manufactured rejection"}'
REJECTION = (b'HTTP/1.1 413 Payload Too Large\r\nContent-Type: application/json\r\n'
             + b'Content-Length: ' + str(len(ERROR_BODY)).encode()
             + b'\r\nConnection: close\r\n\r\n' + ERROR_BODY)
REDIRECTION = (b'HTTP/1.1 302 Found\r\nLocation: https://example.invalid/no-network\r\n'
               b'Content-Length: 0\r\nConnection: close\r\n\r\n')


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def deny_network(*args, **kwargs):
    raise AssertionError('NETWORK_CALL_FORBIDDEN_IN_OFFLINE_PROBE')


class SyntheticSocket:
    def __init__(self, response, body_error):
        self.response = response
        self.body_error = body_error
        self.send_lengths = []
        self.makefile_calls = 0
        self.close_calls = 0

    def sendall(self, data):
        self.send_lengths.append(len(data))
        if len(self.send_lengths) == 2 and self.body_error is not None:
            raise self.body_error

    def makefile(self, mode):
        check(mode == 'rb', 'Unexpected response stream mode')
        self.makefile_calls += 1
        return io.BytesIO(self.response)

    def close(self):
        self.close_calls += 1


class OfflineHTTPSConnection(http.client.HTTPSConnection):
    def __init__(self, host, timeout, response, body_error, context=None):
        super().__init__(host, timeout=timeout, context=context or ssl.create_default_context())
        self.synthetic_socket = SyntheticSocket(response, body_error)
        self.sock = self.synthetic_socket
        self.request_calls = 0
        self.response_calls = 0

    def connect(self):
        raise AssertionError('Synthetic connection must not connect')

    def request(self, *args, **kwargs):
        self.request_calls += 1
        return super().request(*args, **kwargs)

    def getresponse(self):
        self.response_calls += 1
        return super().getresponse()


class OfflineHTTPSHandler(request.HTTPSHandler):
    def __init__(self, response, body_error):
        self.context = ssl.create_default_context()
        super().__init__(context=self.context)
        self.response = response
        self.body_error = body_error
        self.connections = []

    def https_open(self, req):
        def factory(host, timeout, **kwargs):
            connection = OfflineHTTPSConnection(host, timeout, self.response,
                                                self.body_error, **kwargs)
            self.connections.append(connection)
            return connection
        return self.do_open(factory, req, context=self.context)


def source_request_method():
    raw = CONTROLLER.read_bytes()
    check(hashlib.sha256(raw).hexdigest() == CONTROLLER_SHA256, 'Controller pin differs')
    tree = ast.parse(raw.decode(), filename=str(CONTROLLER))
    selected = [n for n in tree.body if getattr(n, 'name', None)
                in {'Stop', 'require', 'safe_url', 'Response', 'NoRedirect'}]
    transport = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'Transport')
    method = next(n for n in transport.body if isinstance(n, ast.FunctionDef) and n.name == 'request')
    selected.append(ast.ClassDef(name='SourceRequestOnly', bases=[], keywords=[],
                                 body=[method], decorator_list=[], type_params=[]))
    namespace = {'__name__': __name__, 'request': request, 'error': error,
                 'parse': parse, 're': re, 'dataclass': dataclass, 'MAX_JSON': 8 * 1024 * 1024}
    isolated = ast.fix_missing_locations(ast.Module(body=selected, type_ignores=[]))
    exec(compile(isolated, str(CONTROLLER), 'exec'), namespace)
    return namespace, namespace['SourceRequestOnly']


def summarize(connection):
    sock = connection.synthetic_socket
    return {'request_calls': connection.request_calls,
            'getresponse_calls': connection.response_calls,
            'makefile_calls': sock.makefile_calls,
            'send_lengths': sock.send_lengths,
            'socket_close_calls': sock.close_calls,
            'manufactured_queued_response_bytes': len(sock.response)}


def run_probe():
    namespace, source_transport_class = source_request_method()
    cases = []
    errors = [
        ('broken_pipe', BrokenPipeError('MANUFACTURED_ONLY')),
        ('connection_reset', ConnectionResetError('MANUFACTURED_ONLY')),
        ('tls_eof', ssl.SSLEOFError('MANUFACTURED_ONLY')),
        ('timeout', TimeoutError('MANUFACTURED_ONLY')),
    ]
    # A lazy anonymous mapping provides the exact bytes-like body length. Fake
    # sendall only measures its memoryview; no upload bytes or real ZIP are read.
    with mmap.mmap(-1, BODY_BYTES) as mapping:
        body = memoryview(mapping)
        try:
            for label, manufactured_error in errors:
                handler = OfflineHTTPSHandler(REJECTION, manufactured_error)
                opener = request.build_opener(request.ProxyHandler({}), namespace['NoRedirect'](), handler)
                req = request.Request(URL, data=body, method='PUT',
                                      headers={'Content-Type': 'application/octet-stream',
                                               'Content-Length': str(BODY_BYTES)})
                try:
                    opener.open(req, timeout=60)
                except error.URLError as exc:
                    check(type(exc.reason) is type(manufactured_error), 'URLError reason differs')
                else:
                    raise AssertionError('Failed send unexpectedly returned HTTP status')
                check(len(handler.connections) == 1, 'urllib attempted more than one request')
                conn = handler.connections[0]
                check(conn.response_calls == 0 and conn.synthetic_socket.makefile_calls == 0,
                      'urllib unexpectedly parsed queued response')
                check(conn.synthetic_socket.send_lengths[1] == BODY_BYTES, 'Body send length differs')
                cases.append({'name': 'urllib_queued_413_send_' + label, 'status': 'PASS',
                              'observed_exception': 'URLError', 'reason_type': type(manufactured_error).__name__,
                              'http_status_exposed': False, **summarize(conn)})

                handler = OfflineHTTPSHandler(REJECTION, manufactured_error)
                transport = source_transport_class()
                transport.token = 'SYNTHETIC_OFFLINE_CREDENTIAL_ONLY'
                transport.opener = request.build_opener(request.ProxyHandler({}), namespace['NoRedirect'](), handler)
                try:
                    transport.request('PUT', URL, data=body, headers={'Content-Type': 'application/octet-stream'})
                except namespace['Stop'] as exc:
                    check(str(exc) == 'TRANSPORT_OUTCOME_UNKNOWN', 'Controller exception differs')
                else:
                    raise AssertionError('Source transport exposed status on failed send')
                conn = handler.connections[0]
                check(len(handler.connections) == 1 and conn.response_calls == 0,
                      'Source transport retried or read response')
                cases.append({'name': 'source_request_queued_413_send_' + label, 'status': 'PASS',
                              'observed_exception': 'Stop', 'reason': 'TRANSPORT_OUTCOME_UNKNOWN',
                              'http_status_exposed': False, **summarize(conn)})

                conn = OfflineHTTPSConnection('zenodo.org', 60, REJECTION, manufactured_error)
                try:
                    conn.request('PUT', '/api/offline-manufactured-test/content', body,
                                 {'Content-Length': str(BODY_BYTES), 'Content-Type': 'application/octet-stream'})
                except OSError:
                    pass
                else:
                    raise AssertionError('Direct send unexpectedly succeeded')
                response = conn.getresponse()
                check(response.status == 413 and response.read() == ERROR_BODY,
                      'Same connection did not expose manufactured final response')
                check(conn.request_calls == 1 and conn.response_calls == 1, 'Direct rescue repeated PUT')
                response.close()
                conn.close()
                cases.append({'name': 'direct_same_connection_rescue_' + label, 'status': 'PASS',
                              'http_status_exposed': 413, **summarize(conn)})

            handler = OfflineHTTPSHandler(REJECTION, None)
            transport = source_transport_class()
            transport.token = 'SYNTHETIC_OFFLINE_CREDENTIAL_ONLY'
            transport.opener = request.build_opener(request.ProxyHandler({}), namespace['NoRedirect'](), handler)
            response = transport.request('PUT', URL, body, {'Content-Type': 'application/octet-stream'})
            check(response.status == 413 and response.body == ERROR_BODY,
                  'Completed body did not return definitive 413')
            check(len(handler.connections) == 1, 'Completed body retried')
            cases.append({'name': 'source_request_complete_send_definitive_413', 'status': 'PASS',
                          'http_status_exposed': 413, **summarize(handler.connections[0])})

            handler = OfflineHTTPSHandler(REDIRECTION, None)
            transport.opener = request.build_opener(request.ProxyHandler({}), namespace['NoRedirect'](), handler)
            response = transport.request('PUT', URL, body, {'Content-Type': 'application/octet-stream'})
            check(response.status == 302 and len(handler.connections) == 1,
                  'Redirect was followed or replayed')
            cases.append({'name': 'source_request_redirect_not_followed', 'status': 'PASS',
                          'http_status_exposed': 302, **summarize(handler.connections[0])})

            conn = OfflineHTTPSConnection('zenodo.org', 60, b'', BrokenPipeError('MANUFACTURED_ONLY'))
            try:
                conn.request('PUT', '/api/offline-manufactured-test/content', body,
                             {'Content-Length': str(BODY_BYTES)})
            except OSError:
                pass
            try:
                conn.getresponse()
            except http.client.RemoteDisconnected:
                pass
            else:
                raise AssertionError('Missing response must remain unknown')
            check(conn.request_calls == 1, 'Missing response retried')
            conn.close()
            cases.append({'name': 'direct_same_connection_no_response_remains_unknown', 'status': 'PASS',
                          'http_status_exposed': False, **summarize(conn)})
        finally:
            body.release()
    return {'schema_version': 1, 'status': 'PASS_SOCKET_FREE_MANUFACTURED_MECHANISM',
            'python': sys.version, 'optimized': bool(sys.flags.optimize),
            'controller_sha256': CONTROLLER_SHA256, 'synthetic_body_bytes': BODY_BYTES,
            'case_count': len(cases), 'cases': cases,
            'scope': 'No real sockets, TLS handshakes, uploads, controller imports/main, saved assets, or credentials. '
                     'Manufactured queued 413 responses do not identify the real failure cause.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    with patch.object(socket, 'socket', deny_network), patch.object(socket, 'create_connection', deny_network):
        result = run_probe()
    output = Path(args.output)
    with output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({'status': result['status'], 'case_count': result['case_count'],
                      'optimized': result['optimized'], 'output': str(output)}))


if __name__ == '__main__':
    main()
