"""Read back the reviewed research packet, summaries and four live pages."""
import argparse
import concurrent.futures
import datetime
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hdblast-root', type=Path, required=True)
    parser.add_argument('--shared-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    public = args.hdblast_root.resolve()
    shared = args.shared_root.resolve()
    capsule = 'research/HDBLAST_ANGULAR_UNIFORM_CONTACT_20261015'
    delivery = 'hdblast/continuation/2026.10.15-angular-uniform-contact'
    manifest = json.loads((public / capsule / 'MANIFEST.json').read_text())
    targets = []
    def add(root, name, url, kind):
        raw = (root / name).read_bytes()
        targets.append({'path': name, 'url': url, 'kind': kind,
                        'expected_bytes': len(raw),
                        'expected_sha256': hashlib.sha256(raw).hexdigest(),
                        'expected_raw': raw})
    prefix = 'https://raw.githubusercontent.com/maldonado-research/HDblast/main/'
    public_paths = [capsule + '/' + row['path'] for row in manifest['files']]
    public_paths += [capsule + '/MANIFEST.json', capsule + '.zip']
    public_paths += [delivery + '/' + name for name in [
        'FROZEN_SOURCE_RECEIPT.json', 'SITE_INTEGRATION_PREFLIGHT.json',
        'NEXT_QUASILINEAR_CONTACT_TARGET.md']]
    public_paths += ['README.md', 'hdblast/README.md', 'hdblast/CURRENT_STATUS.md',
                     'hdblast/PROGRAM_MAP.md', 'site/index.html', 'site/es/index.html']
    for name in public_paths:
        add(public, name, prefix + urllib.parse.quote(name, safe='/'), 'GitHub public source')
    prefix_shared = 'https://raw.githubusercontent.com/maldonado-research/maldonado-research.github.io/main/'
    for name in ['README.md', 'index.html', 'projects/hdblast/index.html']:
        add(shared, name, prefix_shared + urllib.parse.quote(name, safe='/'), 'GitHub shared-site source')
    for root, name, url in [
        (public, 'site/index.html', 'https://maldonado-research.github.io/HDblast/'),
        (public, 'site/es/index.html', 'https://maldonado-research.github.io/HDblast/es/'),
        (shared, 'index.html', 'https://maldonado-research.github.io/'),
        (shared, 'projects/hdblast/index.html', 'https://maldonado-research.github.io/projects/hdblast/'),
    ]:
        add(root, name, url, 'Live Pages')
    def verify(target):
        row = {key: value for key, value in target.items() if key != 'expected_raw'}
        request = urllib.request.Request(target['url'], headers={
            'User-Agent': 'HDBLAST-scientific-publication-readback', 'Cache-Control': 'no-cache'})
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                actual = response.read()
                row.update(http_status=response.status, bytes=len(actual),
                           sha256=hashlib.sha256(actual).hexdigest(),
                           byte_equal=actual == target['expected_raw'])
        except Exception as error:
            row.update(http_status=None, byte_equal=False, error=str(error))
        return row
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        rows = list(pool.map(verify, targets))
    success = all(row['http_status'] == 200 and row['byte_equal'] for row in rows)
    result = {'schema': 'hdblast-angular-uniform-publication-readback-v1',
              'checked_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'status': 'PASS' if success else 'FAIL', 'target_count': len(rows),
              'targets': rows,
              'scope': 'Dated byte-equality witness for reviewed source and deployed HTML; not analytic proof verification.'}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'targets': len(rows), 'output': str(args.output)}))
    if not success:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
