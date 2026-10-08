"""GET-only diagnostic; never imports or changes the publication controller/state."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import ssl
import time
from urllib import request, error

BASE = Path(__file__).resolve().parent
EXECUTION = BASE.parent / 'new-edition-execution'
ORIGIN = 'https://zenodo.org'
DRAFT = '23228395'
NAME = 'HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON.zip'
EXPECTED_BYTES = 240135519
EXPECTED_SHA = 'ead5ffa3987c760d99da7a79189ff5ad766fb9bef6e032c489f4efb2316cd189'
CAP = 8 * 1024 * 1024

class NoRedirect(request.HTTPRedirectHandler):
    def redirect_request(self, *args, **kwargs):
        return None

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    token = os.environ['ZENODO_ACCESS_TOKEN']
    assert token and '\n' not in token and '\r' not in token
    preserved = {x: digest(EXECUTION / x) for x in ('STATE.json', 'JOURNAL.jsonl')}
    opener = request.build_opener(NoRedirect(), request.HTTPSHandler(context=ssl.create_default_context()))
    results = []
    targets = [('draft', '/api/records/' + DRAFT + '/draft'),
               ('files', '/api/records/' + DRAFT + '/draft/files'),
               ('file', '/api/records/' + DRAFT + '/draft/files/' + NAME),
               ('content', '/api/records/' + DRAFT + '/draft/files/' + NAME + '/content')]
    for label, path in targets:
        url = ORIGIN + path
        row = {'label': label, 'method': 'GET', 'path': path,
               'utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
        started = time.monotonic()
        req = request.Request(url, method='GET', headers={
            'Authorization': 'Bearer ' + token, 'Accept-Encoding': 'identity',
            'Accept': 'application/json' if label != 'content' else '*/*',
            'User-Agent': 'HDBLAST-GET-only-unknown-upload-diagnostic/20261008'})
        try:
            try:
                response = opener.open(req, timeout=45)
            except error.HTTPError as http_error:
                response = http_error
            with response:
                row['status'] = response.code
                row['final_url_unchanged'] = response.geturl() == url
                row['headers'] = {k.lower(): v for k, v in response.headers.items()
                    if k.lower() in ('content-type', 'content-length', 'content-range',
                                     'content-encoding', 'etag', 'retry-after', 'server',
                                     'via', 'x-request-id', 'x-ratelimit-limit',
                                     'x-ratelimit-remaining', 'x-ratelimit-reset')}
                if label == 'content' and response.code == 200:
                    sha = hashlib.sha256()
                    md5 = hashlib.md5()
                    count = 0
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        count += len(chunk)
                        if count > EXPECTED_BYTES:
                            raise ValueError('CONTENT_EXCEEDS_EXPECTED_BYTES')
                        sha.update(chunk)
                        md5.update(chunk)
                    row.update(bytes=count, sha256=sha.hexdigest(), md5=md5.hexdigest(),
                               matches_expected=count == EXPECTED_BYTES and sha.hexdigest() == EXPECTED_SHA)
                else:
                    limit = CAP if label in ('draft', 'files', 'file') else 65536
                    raw = response.read(limit + 1)
                    row['body_truncated'] = len(raw) > limit
                    raw = raw[:limit]
                    # Even an unexpected credential echo must not enter saved evidence.
                    raw = raw.replace(token.encode(), b'[REDACTED]')
                    filename = label + '.body'
                    (BASE / filename).write_bytes(raw)
                    row.update(body_file=filename, body_bytes=len(raw), body_sha256=hashlib.sha256(raw).hexdigest())
                    try:
                        obj = json.loads(raw)
                    except Exception:
                        row['body_excerpt'] = raw[:512].decode('utf-8', errors='replace')
                    else:
                        if label == 'draft':
                            row['draft'] = {k: obj.get(k) for k in ('id', 'is_draft', 'is_published', 'revision_id')}
                            row['draft']['parent_id'] = obj.get('parent', {}).get('id')
                            row['draft']['owner'] = obj.get('parent', {}).get('access', {}).get('owned_by', {}).get('user')
                        elif label == 'files':
                            entries = obj.get('entries', [])
                            row['entry_count'] = len(entries)
                            row['completed_count'] = sum(x.get('status') == 'completed' for x in entries)
                            selected = next((x for x in entries if x.get('key') == NAME), None)
                            if selected is not None:
                                row['archive'] = {k: selected.get(k) for k in ('key', 'status', 'size', 'checksum', 'file_id', 'version_id', 'bucket_id', 'updated')}
                        elif label == 'file':
                            row['file'] = {k: obj.get(k) for k in ('key', 'status', 'size', 'checksum', 'file_id', 'version_id', 'bucket_id', 'updated')}
                        else:
                            row['error'] = obj
        except Exception as exception:
            row['exception_type'] = type(exception).__name__
            row['exception'] = str(exception).replace(token, '[REDACTED]')[:1024]
            row['errno'] = getattr(exception, 'errno', None)
            reason = getattr(exception, 'reason', None)
            if reason is not None:
                row['reason_type'] = type(reason).__name__
                row['reason_errno'] = getattr(reason, 'errno', None)
        row['elapsed_seconds'] = round(time.monotonic() - started, 3)
        results.append(row)
        (BASE / (label + '.receipt.json')).write_text(json.dumps(row, indent=2, sort_keys=True) + '\n')
        print(json.dumps(row, sort_keys=True), flush=True)
    after = {x: digest(EXECUTION / x) for x in preserved}
    assert preserved == after, 'CONTROLLER_STATE_OR_JOURNAL_CHANGED_DURING_DIAGNOSTIC'
    receipt = {'status': 'GET_ONLY_DIAGNOSTIC_COMPLETED', 'controller_state_journal_unchanged': True,
               'preserved_hashes': preserved, 'results': results, 'live_mutations': 0}
    (BASE / 'GET_ONLY_DIAGNOSIS.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')

if __name__ == '__main__':
    main()
