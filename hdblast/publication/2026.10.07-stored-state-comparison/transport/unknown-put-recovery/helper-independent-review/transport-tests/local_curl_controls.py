#!/usr/bin/env python3
"""Exercise the installed curl collector against loopback manufactured HTTP only.

Never call helper.main or any real API; production HTTPS/proxy guards are reviewed
separately. Tests derive the accepted command and replace only fake endpoint,
protocol and bounded timeouts. Child environment contains no real credentials.
"""
from pathlib import Path
import argparse
import hashlib
import http.server
import importlib.util
import json
import os
import socket
import subprocess
import sys
import tempfile
import threading
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent.parent / 'recover_put_once.py'
PAYLOAD = b'independent fabricated opaque upload bytes\n' * 2048
TOKEN = 'OFFLINE_SYNTHETIC_TOKEN_9ec8'

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

class ManufacturedServer(http.server.ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = False
    def __init__(self):
        super().__init__(('127.0.0.1', 0), ManufacturedHandler)
        self.observations = []
        self.scenario = 'success'

class ManufacturedHandler(http.server.BaseHTTPRequestHandler):
    protocol_version = 'HTTP/1.1'
    def log_message(self, *args):
        pass
    def handle_expect_100(self):
        self.server.observations.append({'unexpected_expect_handshake': True})
        return False
    def do_PUT(self):
        row = {'method': 'PUT', 'path': self.path, 'expect_present': 'Expect' in self.headers,
               'authorization_correct': self.headers.get('Authorization') == 'Bearer ' + TOKEN,
               'content_length': self.headers.get('Content-Length'), 'received_bytes': 0}
        self.server.observations.append(row)
        scenario = self.server.scenario
        if scenario in ('417', '500', 'redirect'):
            status = {'417': 417, '500': 500, 'redirect': 307}[scenario]
            body = json.dumps({'status': status, 'manufactured': True}).encode()
            self.send_response(status)
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Connection', 'close')
            if scenario == 'redirect':
                self.send_header('Location', '/must-never-follow')
            self.end_headers()
            self.wfile.write(body)
            self.wfile.flush()
            self.close_connection = True
            return
        body = self.rfile.read(len(PAYLOAD))
        row['received_bytes'] = len(body)
        row['received_sha256'] = sha(body)
        if scenario == 'disconnect':
            self.connection.shutdown(socket.SHUT_RDWR)
            self.connection.close()
            self.close_connection = True
            return
        if scenario == 'timeout':
            time.sleep(3)
            self.close_connection = True
            return
        if scenario == 'malformed':
            self.wfile.write(b'NOT_AN_HTTP_STATUS\r\n\r\n')
            self.wfile.flush()
            self.close_connection = True
            return
        if scenario == 'partial':
            self.wfile.write(b'HTTP/1.1 200 OK\r\nContent-Length: 1000\r\nConnection: close\r\n\r\n{}')
            self.wfile.flush()
            self.close_connection = True
            return
        if scenario == 'oversize':
            result = b'x' * 16384
        elif scenario == 'wrong-content-receipt':
            result = json.dumps({'sha256': '0' * 64, 'bytes': len(PAYLOAD)}).encode()
        elif scenario == 'reflected-token':
            result = json.dumps({'reflection': TOKEN}).encode()
        else:
            result = json.dumps({'sha256': sha(body), 'bytes': len(body)}).encode()
        self.send_response(200)
        self.send_header('Content-Length', str(len(result)))
        self.send_header('Connection', 'close')
        self.end_headers()
        try:
            self.wfile.write(result)
            self.wfile.flush()
        except (BrokenPipeError, ConnectionResetError):
            pass
        self.close_connection = True

def load(pin):
    raw = SOURCE.read_bytes()
    require(sha(raw) == pin, 'helper source pin differs at load')
    snapshot = HERE / ('helper_transport_snapshot_' + pin + '.py')
    if not snapshot.exists():
        snapshot.write_bytes(raw)
    require(snapshot.read_bytes() == raw, 'existing snapshot differs')
    spec = importlib.util.spec_from_file_location('manufactured_curl_helper', snapshot)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--helper-sha256', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    helper = load(args.helper_sha256)
    helper.ASSET_PIN = {'filename': 'manufactured.bin', 'bytes': len(PAYLOAD),
                        'sha256': sha(PAYLOAD), 'md5': hashlib.md5(PAYLOAD).hexdigest()}
    helper.MAX_JSON = 8192
    server = ManufacturedServer()
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    helper.URL = 'http://127.0.0.1:' + str(server.server_port) + '/manufactured-upload'
    original_command = helper.curl_command
    command_records = []
    def local_command(asset_fd, header_fd):
        command = original_command(asset_fd, header_fd)
        require(command[0:2] == ['/usr/bin/curl','-q'], 'implicit curlrc not disabled')
        require(command.count('PUT') == 1 and command[command.index('--request')+1] == 'PUT', 'not one PUT command')
        require('--http1.1' in command and '--location' not in command and '-L' not in command, 'HTTP version/redirect boundary differs')
        require(command[command.index('--retry')+1] == '0' and command[command.index('--max-redirs')+1] == '0', 'retry or redirect boundary differs')
        require(command.count('Expect:') == 1 and '100-continue' not in ' '.join(command), 'Expect header not suppressed')
        require(command[command.index('--config')+1] == '-', 'config not via stdin')
        require(command[-1] == helper.URL and TOKEN not in ' '.join(command), 'test endpoint or credential argument boundary differs')
        command[command.index('--proto')+1] = '=http'
        command[command.index('--proto-redir')+1] = '=http'
        command[command.index('--connect-timeout')+1] = '2'
        command[command.index('--max-time')+1] = '2'
        command_records.append({'curlrc_disabled':True,'http1_1':True,'empty_expect':True,
                                'retry_0':True,'no_redirects':True,'config_stdin':True,
                                'fabricated_loopback_url':True,'production_https_override_for_test_only':True})
        return command
    helper.curl_command = local_command
    original_popen = subprocess.Popen
    def sanitized_popen(command, **kwargs):
        kwargs['env'] = {'PATH':'/usr/bin:/bin','LANG':'C.UTF-8','LC_ALL':'C.UTF-8'}
        return original_popen(command, **kwargs)
    helper.subprocess.Popen = sanitized_popen
    rows = []
    try:
        with tempfile.TemporaryDirectory(prefix='hdblast-local-curl-') as temp:
            asset = Path(temp) / 'manufactured.bin'
            asset.write_bytes(PAYLOAD)
            for scenario in ['success','417','500','partial','wrong-content-receipt','redirect','disconnect','malformed','timeout','oversize','reflected-token']:
                server.scenario = scenario
                start = len(server.observations)
                captured = None
                try:
                    captured, capture_receipt = helper.capture_asset(asset, helper.ASSET_PIN)
                    outcome, outputs = helper.run_curl(captured, TOKEN)
                    time.sleep(0.1)
                    observations = server.observations[start:]
                    require(len(observations) == 1, 'unexpected repeated wire request: '+str(observations))
                    request = observations[0]
                    require(request['method'] == 'PUT' and request['path'] == '/manufactured-upload', 'unexpected method/path')
                    require(not request['expect_present'] and request['authorization_correct'], 'Expect or fabricated auth mismatch')
                    require(request['content_length'] == str(len(PAYLOAD)), 'Content-Length mismatch')
                    require(all(len(x) <= helper.MAX_JSON for x in outputs.values()), 'unbounded collector output')
                    require(all(TOKEN.encode() not in x for x in outputs.values()), 'synthetic token not redacted')
                    status = outcome.get('metrics',{}).get('http_status')
                    expected_status = {'success':200,'417':417,'500':500,'partial':200,'wrong-content-receipt':200,
                                       'redirect':307,'disconnect':0,'malformed':0,'timeout':0,'oversize':200,'reflected-token':200}[scenario]
                    require(status == expected_status, 'final status mismatch: '+str(outcome))
                    require(outcome['metrics']['redirects'] == 0, 'redirect was followed')
                    if scenario in ['success','wrong-content-receipt','reflected-token']:
                        require(outcome['curl_exit'] == 0 and request['received_sha256'] == sha(PAYLOAD), 'successful upload byte custody mismatch')
                        require(outcome['metrics']['uploaded_bytes'] == len(PAYLOAD), 'actual size_upload mismatch')
                    if scenario in ['partial','disconnect','malformed','timeout','oversize']:
                        require(outcome['curl_exit'] != 0, 'incomplete/failed result reported zero curl exit')
                    if scenario == 'wrong-content-receipt':
                        require(json.loads(outputs['body'])['sha256'] != sha(PAYLOAD), 'mismatch fixture did not mutate')
                    rows.append({'scenario':scenario,'status':'PASS','outcome':outcome,'request':request,
                                 'outputs':{k:{'bytes':len(v),'sha256':sha(v)} for k,v in outputs.items()},
                                 'capture':capture_receipt,'wire_puts':1,'helper_main_invoked':False})
                except Exception as exc:
                    rows.append({'scenario':scenario,'status':'FAIL','exception':type(exc).__name__,'message':str(exc),
                                 'observations':server.observations[start:]})
                finally:
                    if captured is not None:
                        os.close(captured)
    finally:
        helper.subprocess.Popen = original_popen
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
    result={'schema_version':1,'status':'PASS' if all(x['status']=='PASS' for x in rows) else 'FAIL',
            'optimized':bool(sys.flags.optimize),'helper_sha256':args.helper_sha256,
            'suite_sha256':sha(Path(__file__).read_bytes()),'case_count':len(rows),
            'passed':sum(x['status']=='PASS' for x in rows),'cases':rows,
            'derived_command_guards':command_records,'child_environment':'MINIMAL_NO_REAL_CREDENTIALS_NO_PROXY',
            'network':'LOOPBACK_MANUFACTURED_HTTP_ONLY','real_api_requests':0,
            'live_helper_main_invocations':0,'live_state_edits':0,'scientific_callbacks':0,
            'limits':'Tests observe selected actual curl outcomes; source audit separately verifies production HTTPS/proxy/default TLS and persistent latch. No server CAS claim.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','optimized','helper_sha256','case_count','passed']},indent=2))
    for row in rows:
        if row['status']!='PASS':
            print(json.dumps(row,indent=2))
    return 0 if result['status']=='PASS' else 1

if __name__=='__main__':
    raise SystemExit(main())
