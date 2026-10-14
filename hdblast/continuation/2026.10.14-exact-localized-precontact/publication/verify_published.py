"""Read back the frozen research files and live static pages after deployment."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
from pathlib import Path
import urllib.request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hdblast-root', default=str(Path(__file__).resolve().parents[4]),
                        help='HDblast checkout to compare with published main')
    parser.add_argument('--shared-root', required=True,
                        help='Shared website checkout to compare with the live pages')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    public = Path(args.hdblast_root).resolve()
    localized = public
    material = public
    shared = Path(args.shared_root).resolve()
    base = 'https://raw.githubusercontent.com/maldonado-research/HDblast/main/'
    pairs = [
        ('https://maldonado-research.github.io/HDblast/', public / 'site/index.html'),
        ('https://maldonado-research.github.io/HDblast/es/', public / 'site/es/index.html'),
        ('https://maldonado-research.github.io/projects/hdblast/', shared / 'projects/hdblast/index.html'),
        ('https://maldonado-research.github.io/', shared / 'index.html'),
        (base + 'research/HDBLAST_SECOND_ORDER_LOCALIZED_20261014/MANIFEST.json',
         localized / 'research/HDBLAST_SECOND_ORDER_LOCALIZED_20261014/MANIFEST.json'),
        (base + 'research/HDBLAST_SECOND_ORDER_LOCALIZED_20261014/continuation/REFLECTION_FORMAL_INITIAL_DATA.md',
         localized / 'research/HDBLAST_SECOND_ORDER_LOCALIZED_20261014/continuation/REFLECTION_FORMAL_INITIAL_DATA.md'),
        (base + 'research/HDBLAST_SECOND_ORDER_LOCALIZED_20261014/continuation/PRECONTACT_SECOND_ORDER_DYNAMICS.md',
         localized / 'research/HDBLAST_SECOND_ORDER_LOCALIZED_20261014/continuation/PRECONTACT_SECOND_ORDER_DYNAMICS.md'),
        (base + 'research/HDBLAST_SECOND_ORDER_LOCALIZED_20261014.zip',
         localized / 'research/HDBLAST_SECOND_ORDER_LOCALIZED_20261014.zip'),
        (base + 'research/HDBLAST_MATERIAL_JUNCTIONS_20261014/MANIFEST.json',
         material / 'research/HDBLAST_MATERIAL_JUNCTIONS_20261014/MANIFEST.json'),
        (base + 'research/HDBLAST_MATERIAL_JUNCTIONS_20261014/MATERIAL_SECOND_ORDER_JUNCTIONS.md',
         material / 'research/HDBLAST_MATERIAL_JUNCTIONS_20261014/MATERIAL_SECOND_ORDER_JUNCTIONS.md'),
        (base + 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/MANIFEST.json',
         public / 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/MANIFEST.json'),
        (base + 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md',
         public / 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md'),
        (base + 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014.zip',
         public / 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014.zip'),
        (base + 'hdblast/continuation/2026.10.14-second-order-formal/fresh-full-replay/AUDIT_COMPLETION_RECEIPT.json',
         public / 'hdblast/continuation/2026.10.14-second-order-formal/fresh-full-replay/AUDIT_COMPLETION_RECEIPT.json'),
        (base + 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/README.md',
         public / 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/README.md'),
        (base + 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/CAUCHY_SOURCE_LEDGER.md',
         public / 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/CAUCHY_SOURCE_LEDGER.md'),
        (base + 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/VERIFICATION_AND_SCOPE.json',
         public / 'research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/VERIFICATION_AND_SCOPE.json'),
        (base + 'hdblast/continuation/2026.10.14-exact-localized-precontact/publication/ASTRA_FRESH_REFEREE.md',
         public / 'hdblast/continuation/2026.10.14-exact-localized-precontact/publication/ASTRA_FRESH_REFEREE.md'),
        (base + 'hdblast/continuation/2026.10.14-exact-localized-precontact/publication/SOL_FRESH_CAUCHY_REFEREE.md',
         public / 'hdblast/continuation/2026.10.14-exact-localized-precontact/publication/SOL_FRESH_CAUCHY_REFEREE.md'),
        (base + 'hdblast/continuation/2026.10.14-exact-localized-precontact/publication/FRESH_REFEREE_REVIEWS.md',
         public / 'hdblast/continuation/2026.10.14-exact-localized-precontact/publication/FRESH_REFEREE_REVIEWS.md'),
        (base + 'hdblast/continuation/2026.10.14-exact-localized-precontact/publication/PUBLICATION_RECEIPT.json',
         public / 'hdblast/continuation/2026.10.14-exact-localized-precontact/publication/PUBLICATION_RECEIPT.json'),
    ]
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    def read(pair):
        url, local = pair
        request = urllib.request.Request(url + '?hdblast_readback=' + stamp,
                                         headers={'Cache-Control': 'no-cache', 'User-Agent': 'HDBLAST-verification'})
        with urllib.request.urlopen(request, timeout=35) as response:
            body = response.read()
            expected = local.read_bytes()
            return {'url': url, 'http_status': response.status,
                    'bytes': len(body), 'sha256': hashlib.sha256(body).hexdigest(),
                    'expected_sha256': hashlib.sha256(expected).hexdigest(),
                    'byte_equal_to_reviewed_source': body == expected,
                    'last_modified': response.headers.get('Last-Modified')}
    with ThreadPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(read, pairs))
    receipt = {'checked_utc': stamp, 'status': 'PASS' if all(
        item['http_status'] == 200 and item['byte_equal_to_reviewed_source'] for item in results) else 'FAIL',
        'results': results, 'zenodo_version_created': False, 'aps_submission_created': False}
    Path(args.output).write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': receipt['status'], 'targets_verified': len(results),
                      'receipt': str(Path(args.output).resolve())}, indent=2))
    if receipt['status'] != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
