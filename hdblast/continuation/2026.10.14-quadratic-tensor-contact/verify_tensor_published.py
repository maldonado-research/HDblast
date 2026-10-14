"""Read back every frozen tensor-packet file and the four live static pages."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
from pathlib import Path, PurePosixPath
import urllib.request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hdblast-root', required=True)
    parser.add_argument('--shared-root', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    public = Path(args.hdblast_root).resolve()
    shared = Path(args.shared_root).resolve()
    prefix = 'research/HDBLAST_QUADRATIC_TENSOR_CONTACT_20261014'
    external = 'hdblast/continuation/2026.10.14-quadratic-tensor-contact'
    frozen = json.loads((public / external / 'FROZEN_SOURCE_RECEIPT.json').read_text())
    manifest_path = public / prefix / 'MANIFEST.json'
    if hashlib.sha256(manifest_path.read_bytes()).hexdigest() != frozen['manifest']['sha256']:
        raise RuntimeError('Local external manifest pin differs')
    manifest = json.loads(manifest_path.read_text())
    base = 'https://raw.githubusercontent.com/maldonado-research/HDblast/main/'
    pairs = [
        ('https://maldonado-research.github.io/HDblast/', public / 'site/index.html'),
        ('https://maldonado-research.github.io/HDblast/es/', public / 'site/es/index.html'),
        ('https://maldonado-research.github.io/projects/hdblast/', shared / 'projects/hdblast/index.html'),
        ('https://maldonado-research.github.io/', shared / 'index.html'),
    ]
    seen = set()
    for entry in manifest['files']:
        name = entry['path']
        relative = PurePosixPath(name)
        if relative.is_absolute() or '..' in relative.parts or relative.as_posix() != name or name in seen:
            raise RuntimeError('Noncanonical manifest path')
        seen.add(name)
        local = public / prefix / name
        body = local.read_bytes()
        if len(body) != entry['bytes'] or hashlib.sha256(body).hexdigest() != entry['sha256']:
            raise RuntimeError('Local source differs: ' + name)
        pairs.append((base + prefix + '/' + name, local))
    for path in [prefix + '/MANIFEST.json', prefix + '.zip',
                 external + '/FROZEN_SOURCE_RECEIPT.json',
                 external + '/SITE_INTEGRATION_PREFLIGHT_20261014.json']:
        pairs.append((base + path, public / path))
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')

    def read(pair):
        url, local = pair
        request = urllib.request.Request(url + '?hdblast_readback=' + stamp,
                                         headers={'Cache-Control': 'no-cache', 'User-Agent': 'HDBLAST-verification'})
        try:
            with urllib.request.urlopen(request, timeout=35) as response:
                body = response.read()
                expected = local.read_bytes()
                return {'url': url, 'http_status': response.status,
                        'bytes': len(body), 'sha256': hashlib.sha256(body).hexdigest(),
                        'expected_sha256': hashlib.sha256(expected).hexdigest(),
                        'byte_equal_to_reviewed_source': body == expected,
                        'last_modified': response.headers.get('Last-Modified')}
        except Exception as error:
            return {'url': url, 'error': type(error).__name__ + ': ' + str(error),
                    'byte_equal_to_reviewed_source': False}

    with ThreadPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(read, pairs))
    status = 'PASS' if all(row.get('http_status') == 200 and row['byte_equal_to_reviewed_source']
                           for row in results) else 'FAIL'
    receipt = {'checked_utc': stamp, 'status': status, 'results': results,
               'scope': 'Public file and website byte equality; not a mathematical or cosmological proof.',
               'zenodo_version_created': False, 'aps_submission_created': False}
    Path(args.output).write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': status, 'targets_verified': len(results), 'receipt': args.output}, indent=2))
    if status != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
