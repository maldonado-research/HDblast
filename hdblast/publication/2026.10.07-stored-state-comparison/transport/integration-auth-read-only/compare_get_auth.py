"""Read-only, same-endpoint urllib/curl Authorization representation comparison."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import ssl
import subprocess
import time
from urllib import request, error

HERE = Path(__file__).resolve().parent
URL = 'https://zenodo.org/api/records/23228395/draft'
VENDOR = 'application/vnd.inveniordm.v1+json'
CAP = 8 * 1024 * 1024

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def redact(raw, token):
    escaped = json.dumps(token, ensure_ascii=False)[1:-1].encode()
    return raw.replace(token.encode(), b'[REDACTED]').replace(escaped, b'[REDACTED]')

def capture(label, row, raw, headers, token):
    raw = redact(raw, token)
    headers = redact(headers, token)
    (HERE / (label + '.body')).write_bytes(raw)
    (HERE / (label + '.headers')).write_bytes(headers)
    row.update(body_bytes=len(raw), body_sha256=digest(raw), headers_bytes=len(headers), headers_sha256=digest(headers))
    try:
        obj = json.loads(raw)
    except Exception:
        row['body_excerpt'] = raw[:256].decode(errors='replace')
    else:
        if row.get('status') == 200:
            row['draft'] = {k: obj.get(k) for k in ('id', 'is_draft', 'is_published', 'revision_id')}
            row['draft']['parent'] = obj.get('parent', {}).get('id')
            row['draft']['owner'] = obj.get('parent', {}).get('access', {}).get('owned_by', {}).get('user')
        else:
            row['error'] = obj
    (HERE / (label + '.receipt.json')).write_text(json.dumps(row, sort_keys=True, indent=2) + '\n')
    print(json.dumps(row, sort_keys=True), flush=True)

def main():
    token = os.environ['ZENODO_ACCESS_TOKEN']
    assert token and all(32 <= ord(x) < 127 for x in token)
    class NoRedirect(request.HTTPRedirectHandler):
        def redirect_request(self, *args, **kwargs):
            return None
    opener = request.build_opener(NoRedirect(), request.HTTPSHandler(context=ssl.create_default_context()))
    start = time.monotonic()
    req = request.Request(URL, method='GET', headers={'Authorization': 'Bearer ' + token,
        'Accept': VENDOR, 'Accept-Encoding': 'identity', 'User-Agent': 'HDBLAST-GET-only-integration-auth-diagnostic/20261008'})
    try:
        response = opener.open(req, timeout=30)
    except error.HTTPError as exception:
        response = exception
    with response:
        raw = response.read(CAP + 1)
        assert len(raw) <= CAP
        headers = json.dumps({k: v for k, v in response.headers.items()
            if k.lower() not in ('set-cookie', 'authorization', 'proxy-authorization')}, sort_keys=True).encode()
        row = {'label': 'urllib_vendor', 'method': 'GET', 'status': response.code,
            'final_url_unchanged': response.geturl() == URL, 'elapsed_seconds': round(time.monotonic() - start, 3)}
        capture('urllib_vendor', row, raw, headers, token)
    config = ('header = "Authorization: Bearer ' + token.replace('\\', '\\\\').replace('"', '\\"') + '"\n').encode()
    for retained in (False, True):
        label = 'curl_vendor_z_retained' if retained else 'curl_vendor_z_removed'
        environment = os.environ.copy()
        environment.pop('SSLKEYLOGFILE', None)
        if not retained:
            environment.pop('ZENODO_ACCESS_TOKEN', None)
        headers_path = HERE / (label + '.headers.tmp')
        command = ['/usr/bin/curl', '-q', '--config', '-', '--http1.1', '--request', 'GET',
            '--header', 'Accept: ' + VENDOR, '--header', 'Accept-Encoding: identity',
            '--proto', '=https', '--proto-redir', '=https', '--max-redirs', '0', '--retry', '0',
            '--connect-timeout', '30', '--max-time', '45', '--max-filesize', str(CAP),
            '--dump-header', str(headers_path), '--output', '-', '--silent', '--show-error',
            '--write-out', '%{stderr}\nAUTH_GET_METRICS|%{http_code}|%{ssl_verify_result}|%{num_redirects}\n', URL]
        start = time.monotonic()
        result = subprocess.run(command, input=config, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env=environment, timeout=50, check=False)
        assert len(result.stdout) <= CAP and len(result.stderr) <= 65536
        headers = headers_path.read_bytes()
        assert len(headers) <= 65536
        headers_path.unlink()
        stderr = redact(result.stderr, token)
        (HERE / (label + '.stderr')).write_bytes(stderr)
        parts = stderr.rsplit(b'AUTH_GET_METRICS|', 1)[-1].strip().split(b'|')
        assert len(parts) == 3
        row = {'label': label, 'method': 'GET', 'status': int(parts[0]), 'curl_exit': result.returncode,
            'ssl_verify_result': int(parts[1]), 'redirects': int(parts[2]),
            'zenodo_variable_retained_in_child': retained, 'ssl_keylog_removed': True,
            'proxy_trust_environment_preserved': True, 'token_delivery': 'STDIN_CONFIG_AUTHORIZATION_HEADER',
            'elapsed_seconds': round(time.monotonic() - start, 3)}
        capture(label, row, result.stdout, headers, token)
    print(json.dumps({'status': 'GET_ONLY_COMPARISON_COMPLETED', 'live_mutations': 0}), flush=True)

if __name__ == '__main__':
    main()
