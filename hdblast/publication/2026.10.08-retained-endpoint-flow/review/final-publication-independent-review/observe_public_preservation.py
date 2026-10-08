"""Anonymous GET-only historical Zenodo record/list preservation observation."""
from pathlib import Path
import concurrent.futures
import datetime
import hashlib
import json
import ssl
import urllib.error
import urllib.request
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
SOURCE = BASE / 'release-assets/publication-verifier'
MANIFEST_SHA = 'f0f094a5e68f2ed02b7d5810cb928ed178783327eadaaf4fac167976fbb34ca1'
IDS = {'23228395': '17088132', '23225288': '17088132', '23114217': '17088132', '22347452': '17088132', '23111008': '22922927'}
CAP = 8 * 1024 * 1024


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def unique(items):
    result = {}
    for key, value in items:
        require(key not in result, 'duplicate JSON key')
        result[key] = value
    return result


def loads(raw):
    return json.loads(raw, object_pairs_hook=unique)


def write(path, raw):
    with path.open('xb') as stream:
        stream.write(raw)


def dump(path, obj):
    write(path, (json.dumps(obj, sort_keys=True, indent=2) + '\n').encode())


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *args, **kwargs):
        return None


def get(output, rid, kind):
    url = 'https://zenodo.org/api/records/' + rid + ('/files' if kind == 'files' else '')
    start = now()
    context = ssl.create_default_context()
    opener = urllib.request.build_opener(NoRedirect(), urllib.request.HTTPSHandler(context=context))
    accept = 'application/json' if kind == 'files' else 'application/vnd.inveniordm.v1+json'
    request = urllib.request.Request(url, method='GET', headers={
        'Accept': accept, 'Accept-Encoding': 'identity',
        'User-Agent': 'HDBLAST-independent-read-only-preservation/20261008'})
    try:
        with opener.open(request, timeout=45) as response:
            raw = response.read(CAP + 1)
            status, final, headers = response.status, response.url, dict(response.headers.items())
    except urllib.error.HTTPError as error:
        raw = error.read(CAP + 1)
        status, final, headers = error.code, error.url, dict(error.headers.items())
    require(len(raw) <= CAP, 'JSON response too large')
    name = rid + '-' + kind
    write(output / (name + '.json'), raw)
    receipt = {'url': url, 'final_url': final, 'method': 'GET', 'http_status': status,
               'request_started_utc': start, 'response_completed_utc': now(),
               'bytes': len(raw), 'sha256': digest(raw), 'authenticated_request': False,
               'accept': accept,
               'TLS_certificate_and_hostname_verification': True, 'redirects_followed': 0,
               'response_date': headers.get('Date', headers.get('date')),
               'etag': headers.get('ETag', headers.get('etag'))}
    dump(output / (name + '.GET.json'), receipt)
    require(status == 200 and final == url, 'public JSON GET failed:' + name)
    return loads(raw), receipt


def files(obj):
    entries = obj['files']['entries'] if 'files' in obj else obj['entries']
    if isinstance(entries, dict):
        for name, row in entries.items():
            require(row.get('key', row.get('filename')) == name, 'file dictionary key mismatch')
        values = list(entries.values())
    else:
        require(type(entries) is list, 'complete file list')
        values = entries
    result = {}
    for row in values:
        name = row.get('key', row.get('filename'))
        require(type(name) is str and name and name not in result, 'unique file name')
        result[name] = row
    return result


def main():
    require(len(sys.argv) == 2 and sys.argv[1] in ('observation-002', 'observation-003'), 'fixed observation output name')
    out = HERE / sys.argv[1]
    out.mkdir()
    raw = (SOURCE / 'MANIFEST.json').read_bytes()
    require(digest(raw) == MANIFEST_SHA, 'baseline manifest external pin changed')
    manifest = loads(raw)
    baseline = {}
    for rid in IDS:
        for kind in ('record', 'files'):
            name = f'manufactured-seeds/{rid}-{kind}.json'
            raw = (SOURCE / name).read_bytes()
            require({'bytes': len(raw), 'sha256': digest(raw)} == manifest['files'][name], 'baseline pin changed')
            baseline[(rid, kind)] = loads(raw)
    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        futures = {(rid, kind): pool.submit(get, out, rid, kind) for rid in IDS for kind in ('record', 'files')}
        for key, future in futures.items():
            results[key] = future.result()
    checks = {}
    for rid, concept in IDS.items():
        record, rr = results[(rid, 'record')]
        listing, lr = results[(rid, 'files')]
        old = baseline[(rid, 'record')]
        old_files = files(baseline[(rid, 'files')])
        current_files, embedded = files(listing), files(record)
        require(record['id'] == rid and record['parent']['id'] == concept, 'record/concept identity')
        require(record['is_published'] is True and record['is_draft'] is False, 'public state')
        require(record['pids']['doi']['identifier'] == '10.5281/zenodo.' + rid, 'version DOI')
        require(record['parent']['pids']['doi']['identifier'] == '10.5281/zenodo.' + concept, 'concept DOI')
        require(str(record['parent']['access']['owned_by']['user']) == '1386319', 'owner identity')
        require(all(record[key] == old[key] for key in ('metadata', 'custom_fields', 'access')), 'full exact editable metadata changed')
        require(set(current_files) == set(old_files) == set(embedded) == set(files(old)), 'file membership changed')
        identities = {}
        for name, previous in old_files.items():
            current = current_files[name]
            require(type(current['size']) is int and current['size'] == previous['size'], 'file size changed')
            require(current['checksum'] == previous['checksum'] and current['status'] == 'completed', 'file checksum/status changed')
            require(current['file_id'] == previous['file_id'] and current['version_id'] == previous['version_id'], 'immutable file/version identity changed')
            e = embedded[name]
            require(e.get('id', e.get('file_id')) == current['file_id'], 'embedded file identity changed')
            require(e.get('size', e.get('filesize')) == current['size'] and e['checksum'] == current['checksum'], 'embedded file pins changed')
            identities[name] = {'bytes': current['size'], 'checksum': current['checksum'], 'file_id': current['file_id'], 'version_id': current['version_id']}
        checks[rid] = {'concept_id': concept, 'file_count': len(current_files), 'total_file_bytes': sum(row['size'] for row in current_files.values()),
                       'full_metadata_custom_fields_access_exact_equal': True, 'complete_file_listing_exact_equal': current_files == old_files,
                       'immutable_file_version_ids_equal': True, 'file_pins': identities,
                       'record_snapshot': rid + '-record.json', 'files_snapshot': rid + '-files.json',
                       'record_started_utc': rr['request_started_utc'], 'files_completed_utc': lr['response_completed_utc']}
    receipt = {'status': 'PASS_INDEPENDENT_PUBLIC_HISTORICAL_PRESERVATION_GETS', 'observation': sys.argv[1],
               'completed_utc': now(), 'baseline_manifest_sha256': MANIFEST_SHA,
               'timing_scope': 'Dated independent observations; before/after publication classification must be determined from the publisher journal, not assumed from this label.',
               'records': checks, 'GET_requests': 10, 'remote_writes': 0, 'attachment_streams': 0,
               'new_attachment_SHA256_verification': False, 'internal_independent_review_not_external_peer_review': True}
    dump(out / 'PUBLIC_PRESERVATION_REVIEW.json', receipt)
    print(json.dumps({'status': receipt['status'], 'observation': sys.argv[1], 'completed_utc': receipt['completed_utc'], 'records': list(checks)}))


if __name__ == '__main__':
    main()
