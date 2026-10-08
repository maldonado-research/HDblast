#!/usr/bin/env python3
"""Prepare scoped publication-hold candidates outside all repositories.

Pure standard-library rendering of already completed results and dated publication
facts. This performs no numerical producer import, input decode or remote request.
"""
from pathlib import Path
from urllib.parse import quote
import hashlib
import html
import json
import re

W = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
D = W / 'sites-publication-hold-draft'
ORIGINAL = W / 'sites-draft'
OUTPUT = D / 'candidates'
GUIDE = 'https://github.com/maldonado-research/HDblast/tree/main/hdblast/publication/2026.10.07-stored-state-comparison'
PACKAGE = 'https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_STORED_STATE_COMPARISON_PACKAGE_20261007'
OLD_PROOF = 'https://github.com/maldonado-research/HDblast/tree/main/hdblast/publication/2026.10.07-bd-prehistory'
ARCHIVE = 'HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON.zip'
ARCHIVE_SIZE = 240135519
ARCHIVE_SHA = 'ead5ffa3987c760d99da7a79189ff5ad766fb9bef6e032c489f4efb2316cd189'
ARCHIVE_MD5 = 'e424502f12dbaacd36abaaabdd3fb74e'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def one(text, old, new):
    if text.count(old) != 1:
        raise ValueError(('Replacement occurrence count', old[:100], text.count(old)))
    return text.replace(old, new, 1)


def sub_one(text, expression, new):
    result, count = re.subn(expression, lambda _: new, text, flags=re.S)
    if count != 1:
        raise ValueError(('Replacement occurrence count', expression, count))
    return result


def hold_md():
    return ('The completed comparison and standalone replay are available on GitHub. '
            'The sealed Zenodo comparison archive is **240,135,519 bytes**, but its whole-file upload '
            'is still **UNVERIFIED** in owned same-family draft **23228395**. That draft remains '
            '**unsubmitted**; no new comparison edition or new version DOI is claimed. '
            'The latest independently verified published edition remains '
            '[23225288](https://zenodo.org/records/23225288), DOI **10.5281/zenodo.23225288**, '
            'with **27 files / 448,387,915 bytes**. '
            '[Prepared comparison edition and publication hold](hdblast/publication/2026.10.07-stored-state-comparison/README.md).')


def hold_html(spanish=False):
    if spanish:
        return ('<p data-hdb-stored-edition="publication-hold">La comparación completa y su reproducción '
                f'independiente están disponibles en <a href="{PACKAGE}">GitHub</a>. El archivo '
                'sellado para Zenodo tiene <strong>240.135.519 bytes</strong>, pero su carga íntegra '
                'sigue <strong>UNVERIFIED</strong> en el borrador propio de la misma familia '
                '<strong>23228395</strong>. El borrador sigue sin enviar; no se afirma una nueva '
                'edición publicada ni un DOI de versión nuevo. La última edición publicada y '
                'verificada sigue siendo <a href="https://zenodo.org/records/23225288">23225288</a>, '
                'con <strong>27 archivos / 448.387.915 bytes</strong>. '
                f'<a href="{GUIDE}">Preparación y publicación pendiente (inglés)</a>.</p>')
    return ('<p data-hdb-stored-edition="publication-hold">The completed comparison and independent '
            f'portable replay are available on <a href="{PACKAGE}">GitHub</a>. The sealed Zenodo archive '
            'is <strong>240,135,519 bytes</strong>, but its whole-file upload remains '
            '<strong>UNVERIFIED</strong> in owned same-family draft <strong>23228395</strong>. '
            'That draft is unsubmitted; no new comparison publication or new version DOI is claimed. '
            'The latest independently verified published edition remains '
            '<a href="https://zenodo.org/records/23225288">23225288</a>, '
            'with <strong>27 files / 448,387,915 bytes</strong>. '
            f'<a href="{GUIDE}">Prepared edition and publication hold</a>.</p>')


def published_file_page(inventory):
    rows = []
    for row in sorted(inventory['files'], key=lambda r: r['filename']):
        url = 'https://zenodo.org/records/23225288/files/' + quote(row['filename'], safe='') + '?download=1'
        rows.append('      <tr><th scope="row"><a href="' + html.escape(url, quote=True) + '">' +
                    html.escape(row['filename']) + '</a></th><td>' + format(row['bytes'], ',') +
                    '</td><td><code>SHA256 ' + row['sha256'] + '</code><code>MD5 ' + row['md5'] +
                    '</code></td></tr>')
    return '''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'self'; base-uri 'none'; object-src 'none'; form-action 'none'">
  <meta name="referrer" content="no-referrer">
  <meta name="description" content="Dated verified HDBLAST published inventory: 27 files in record 23225288; the later comparison edition remains unpublished.">
  <title>HDBLAST verified publication files</title>
  <link rel="canonical" href="https://maldonado-research.github.io/HDblast/zenodo-files/">
  <link rel="stylesheet" href="./style.css">
</head>
<body>
  <main>
    <h1>HDBLAST verified publication files</h1>
    <p>The latest independently verified published same-family edition is <a href="https://zenodo.org/records/23225288">23225288</a>, DOI <a href="https://doi.org/10.5281/zenodo.23225288">10.5281/zenodo.23225288</a>, semantic version <strong>2026.10.07-bd-prehistory-target</strong>.</p>
    <p>Its complete inventory contains <strong>27 files / 448,387,915 bytes</strong>. Every file was fully streamed and matched its size, MD5 and SHA256 before and after that publication. The dated public verification completed on <strong>8 October 2026 at 00:14:31 UTC</strong>, or 7 October Pacific. This page reproduces that preserved verified inventory; it does not claim another fresh full readback.</p>
    HOLD_HTML
    <p><a href="OLD_PROOF">Published edition 23225288: frozen witnesses and offline verifier</a> · <a href="https://doi.org/10.5281/zenodo.17088132">Preserved concept DOI</a> · <a href="../">Project website</a></p>
    <h2>Complete published attachment inventory</h2>
    <p>Checksums identify the downloaded bytes. Each linked file is hosted by Zenodo in the preserved published record. The later all-node comparison archive is not included in this 27-file edition.</p>
    <div class="inventory"><table><thead><tr><th scope="col">File</th><th scope="col">Bytes</th><th scope="col">SHA256 / MD5</th></tr></thead><tbody>
ROWS
    </tbody></table></div>
    <h2>Earlier delivery records</h2>
    <p>The <a href="https://github.com/maldonado-research/HDblast/tree/main/hdblast/publication/2026.10.05-saved-metadata">5 October browser-delivery packet</a> preserves the historical thirteen-file upload procedure. Edition <a href="https://zenodo.org/records/23114217">23114217</a> is published with its verified 23-file inventory; those old draft instructions are archival history. Conditional calculations and successful reproduction do not establish a higher-dimensional origin of the Big Bang.</p>
  </main>
</body>
</html>
'''.replace('HOLD_HTML', hold_html()).replace('OLD_PROOF', OLD_PROOF).replace('ROWS', '\n'.join(rows))


def render(entry, data, inventory):
    repo, rel = entry['repository'], entry['target_relative_path']
    if Path(rel).suffix not in {'.md', '.html', '.json', '.css'}:
        return data
    text = data.decode()
    if rel == 'README.md' and repo in {'HDblast', 'HDblast-archive'}:
        text = sub_one(text, r'The verified same-family edition \[NEW_RECORD_ID\].*?\[Publication proof\]\(hdblast/publication/2026\.10\.07-stored-state-comparison/README\.md\)\.', hold_md())
        if repo == 'HDblast-archive':
            text = sub_one(text, r'\*\*Verified main publication:\*\*.*?(?=\n)',
                           '**Latest verified published main edition:** [Zenodo 23225288](https://zenodo.org/records/23225288), '
                           'semantic version **2026.10.07-bd-prehistory-target**, DOI **10.5281/zenodo.23225288**, '
                           '**27 files / 448,387,915 bytes** in concept family 17088132. '
                           'The later comparison edition remains unpublished in owned draft 23228395; its archive upload is unverified.')
    elif rel == 'hdblast/CURRENT_STATUS.md':
        start = text.index('## Latest verified Zenodo publication\n')
        end = text.index('## Preserved verified edition 23225288: BD prehistory\n', start)
        replacement = '''## Comparison edition: prepared, with Zenodo publication on hold

The new all-node comparison and standalone package are complete. The sealed archive **HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON.zip** is **240,135,519 bytes**, SHA256 **ead5ffa3987c760d99da7a79189ff5ad766fb9bef6e032c489f4efb2316cd189**, MD5 **e424502f12dbaacd36abaaabdd3fb74e**. Fresh extraction and bounded normal/optimized portable replay reproduce all eight scientific files exactly. This is local package and replay evidence, not a completed Zenodo upload.

Owned same-family draft **23228395**, based on published 23225288 and concept **10.5281/zenodo.17088132**, remains **unsubmitted**. The whole comparison archive upload is **UNVERIFIED**. The approved candidate edition has **29 files / 688,531,967 bytes**, but those are prepared local inventory values, not a verified published inventory. No new version DOI or new comparison publication is claimed. Full saved metadata, file names, sizes, MD5 and SHA256 still require independent verification before any publication. Preserve the existing draft and journals; do not rerun completed historical publishers or publish an incomplete inherited-only draft.

The [prepared comparison edition and publication-hold guide](publication/2026.10.07-stored-state-comparison/README.md) records the sealed local assets and unresolved upload. The recovery outcome is not inferred from a timeout or the passage of time. Earlier publications remain preserved. Automatic GitHub-to-Zenodo archiving was not toggled; no fresh account-side integration-setting inspection is claimed. Client-side state checks do not establish an atomic server publication compare-and-swap guarantee.

'''
        text = text[:start] + replacement + text[end:]
        text = one(text, '## Preserved verified edition 23225288: BD prehistory',
                   '## Latest independently verified published edition 23225288: BD prehistory')
    elif rel in {'site/index.html', 'site/es/index.html'} and repo == 'HDblast':
        spanish = rel == 'site/es/index.html'
        text = sub_one(text, r'<p data-hdb-stored-edition="verified">.*?</p>', hold_html(spanish))
        text = one(text, '"version": "NEW_SEMANTIC_VERSION"', '"version": "2026.10.07-stored-state-comparison"')
        text = one(text, '"url": "https://zenodo.org/records/NEW_RECORD_ID"',
                   '"url": "https://zenodo.org/records/23225288", "name": "Latest verified published edition: BD prehistory; later comparison upload unresolved"')
        text = text.replace('The preserved previous main-family edition <a href="https://zenodo.org/records/23225288">',
                            'The latest independently verified published main-family edition <a href="https://zenodo.org/records/23225288">')
        text = text.replace('La edición principal anterior conservada <a href="https://zenodo.org/records/23225288">',
                            'La última edición principal publicada y verificada <a href="https://zenodo.org/records/23225288">')
    elif rel == 'site/zenodo-files/index.html':
        text = published_file_page(inventory)
    elif rel == 'README.md' and repo == 'maldonado-research.github.io':
        text = sub_one(text, r'The new verified edition NEW_RECORD_ID, DOI NEW_VERSION_DOI,.*?All contents were independently streamed before and after publication\.',
                       'The completed comparison and portable replay are available on GitHub. The sealed comparison archive is 240,135,519 bytes, '
                       'but its whole-file Zenodo upload remains UNVERIFIED in owned same-family draft 23228395. The draft remains unsubmitted; '
                       'no new comparison publication or new version DOI is claimed. Latest independently verified published record 23225288, '
                       'DOI 10.5281/zenodo.23225288, retains its 27 files / 448,387,915 bytes in concept 17088132. '
                       'The dated full public readback of that published edition completed on 8 October 2026 at 00:14:31 UTC.')
    elif rel == 'projects/hdblast/index.html':
        text = one(text, '"url": "https://zenodo.org/records/NEW_RECORD_ID", "name": "Verified stored-state comparison edition NEW_SEMANTIC_VERSION"',
                   '"url": "https://zenodo.org/records/23225288", "name": "Latest verified published edition: BD prehistory; later comparison upload unresolved"')
        text = sub_one(text, r'<p data-hdb-stored-edition="verified">.*?</p>', hold_html())
        text = text.replace('The preserved previous main-family edition <a href="https://zenodo.org/records/23225288">',
                            'The latest independently verified published main-family edition <a href="https://zenodo.org/records/23225288">')
    elif rel.endswith('HDBLAST_STORED_STATE_COMPARISON_CONTINUATION_20261007.md'):
        text = sub_one(text, r'The new verified edition \*\*NEW_RECORD_ID\*\*.*?(?=\n\n)',
                       'The comparison edition remains **unpublished** in owned same-family draft **23228395**, '
                       'preserving concept **10.5281/zenodo.17088132**. Its locally sealed archive is **240,135,519 bytes**, '
                       'SHA256 **ead5ffa3987c760d99da7a79189ff5ad766fb9bef6e032c489f4efb2316cd189**, '
                       'MD5 **e424502f12dbaacd36abaaabdd3fb74e**. The prepared 29-file candidate totals **688,531,967 bytes**; '
                       'these are local prepared values, not a verified remote or published inventory. Whole-file upload remains **UNVERIFIED** '
                       'and the draft is **unsubmitted**. No new version DOI or new published edition is claimed. '
                       'Latest independently verified published record **23225288**, DOI **10.5281/zenodo.23225288**, '
                       'retains **27 files / 448,387,915 bytes** with its dated full public readback at **2026-10-08 00:14:31 UTC**. '
                       'Preserve initialized draft file state and every controller/recovery journal; an uncertain upload outcome requires read-only '
                       'reconciliation and independently reviewed recovery, not a repeated publication attempt or removal of the draft. '
                       'Root owns any later guarded recovery result; no such result or failure timestamp is inferred here. '
                       'The controller does not rerun older completed publishers or change automatic GitHub-to-Zenodo settings. '
                       'Atomic server compare-and-swap remains NOT_ESTABLISHED. Earlier main records, static companion and separate DOI family remain unchanged.')
    elif rel.endswith('HDBLAST_CURRENT_STATE_STORED_STATE_COMPARISON_20261007.json'):
        state = json.loads(text)
        state['publication'] = {
            'status': 'LATEST_INDEPENDENTLY_VERIFIED_PUBLISHED_EDITION',
            'record_id': '23225288',
            'version_doi': '10.5281/zenodo.23225288',
            'concept_doi': '10.5281/zenodo.17088132',
            'semantic_version': '2026.10.07-bd-prehistory-target',
            'file_count': 27,
            'total_bytes': 448387915,
            'public_full_stream_verification_utc': '2026-10-08 00:14:31 UTC',
            'fresh_full_stream_readback_claimed_by_hold_candidates': False,
        }
        state['comparison_edition'] = {
            'status': 'UNSUBMITTED_ARCHIVE_UPLOAD_UNVERIFIED',
            'draft_id': '23228395',
            'parent_record_id': '23225288',
            'concept_doi': '10.5281/zenodo.17088132',
            'semantic_version': '2026.10.07-stored-state-comparison',
            'prepared_local_file_count': 29,
            'prepared_local_total_bytes': 688531967,
            'new_publication_claimed': False,
            'new_version_doi_claimed': False,
            'all_29_remote_streams_verified_before_publication': False,
            'archive': {'filename': ARCHIVE, 'bytes': ARCHIVE_SIZE, 'sha256': ARCHIVE_SHA, 'md5': ARCHIVE_MD5},
            'whole_archive_upload_verified': False,
            'guarded_recovery_result_claimed': False,
        }
        state['portable_fresh_extraction_replay'] = {
            'status': 'PASS_COMPLETE_PORTABLE_STORED_STATE_COMPARISON_REPLAY',
            'receipt_sha256': '0abe89a95cac5af9a0ebdf3afe736b8acc97c8bd5c4ff5729ef4a1b1aafd532a',
            'normal_optimized_scientific_files_exact': 8,
            'publication_upload_evidence': False,
        }
        text = json.dumps(state, indent=2) + '\n'
    if 'NEW_' in text or 'VERIFIED_INVENTORY_ROWS' in text:
        raise ValueError(('Unresolved placeholder', repo, rel))
    return text.encode()


def main():
    if OUTPUT.exists():
        raise ValueError('Fresh outside-repository output required')
    for repo in ('HDblast', 'HDblast-archive', 'maldonado-research.github.io'):
        if OUTPUT.resolve().is_relative_to((Path('/workspace') / repo).resolve()):
            raise ValueError('Output must be outside repositories')
    original_raw = (ORIGINAL / 'CANDIDATE_MAPPING.json').read_bytes()
    mapping = json.loads(original_raw)
    inventory_path = W / 'publication-planning/read-only-observation/CURRENT_PUBLISHED_INVENTORY.json'
    inventory_raw = inventory_path.read_bytes()
    inventory = json.loads(inventory_raw)
    if (inventory['current_record'], inventory['expected_files'], inventory['expected_bytes']) != (23225288, 27, 448387915):
        raise ValueError('Dated published inventory identity')
    if len(inventory['files']) != 27 or sum(r['bytes'] for r in inventory['files']) != 448387915:
        raise ValueError('Dated published inventory count or bytes')
    archive = W / 'publication-final/assets' / ARCHIVE
    archive_data = archive.read_bytes()
    if (len(archive_data), sha(archive_data), hashlib.md5(archive_data).hexdigest()) != (ARCHIVE_SIZE, ARCHIVE_SHA, ARCHIVE_MD5):
        raise ValueError('Sealed local archive pin')
    del archive_data
    results = []
    for row in mapping['mapping']:
        original = Path(row['candidate_path'])
        data = original.read_bytes()
        if (len(data), sha(data)) != (row['candidate_bytes'], row['candidate_sha256']):
            raise ValueError('Original candidate bytes changed')
        target = OUTPUT / row['repository'] / row['target_relative_path']
        target.parent.mkdir(parents=True, exist_ok=True)
        transformed = render(row, data, inventory)
        target.write_bytes(transformed)
        results.append({**row, 'original_candidate_path': str(original),
                        'original_candidate_sha256': row['candidate_sha256'],
                        'candidate_path': str(target), 'candidate_bytes': len(transformed),
                        'candidate_sha256': sha(transformed)})
    receipt = {
        'schema_version': 1,
        'status': 'PREPARED_PUBLICATION_HOLD_CANDIDATES_PENDING_ROOT_REVIEW',
        'use_only_if_root_confirms_current_publication_hold': True,
        'original_mapping_sha256': sha(original_raw),
        'preserved_published_inventory_sha256': sha(inventory_raw),
        'latest_independently_verified_published_record': '23225288',
        'latest_verified_published_files': 27,
        'latest_verified_published_bytes': 448387915,
        'comparison_draft_id': '23228395',
        'comparison_archive_upload_verified': False,
        'new_publication_claimed': False,
        'new_version_doi_claimed': False,
        'recovery_result_or_timestamp_invented': False,
        'network_requests': 0,
        'repository_mutations': 0,
        'remote_mutations': 0,
        'new_retained_array_decodes': 0,
        'new_physical_source_constructions': 0,
        'mapping': results,
    }
    (D / 'CANDIDATE_MAPPING.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ('status', 'comparison_draft_id', 'new_publication_claimed', 'repository_mutations')}))


if __name__ == '__main__':
    main()
