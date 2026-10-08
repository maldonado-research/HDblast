#!/usr/bin/env python3
"""Read-only local checks of fallback rendering and its preserved scope."""
from pathlib import Path
from fractions import Fraction
from html.parser import HTMLParser
from urllib.parse import urlparse, unquote
import hashlib
import json
import re
import subprocess

W = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
D = W / 'sites-publication-hold-draft'
ORIGINAL = W / 'sites-draft'


def sha(data):
    return hashlib.sha256(data).hexdigest()


class Page(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.links = []
        self.images = []
        self.ids = []
        self.canonicals = []
        self.table_rows = 0

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs:
            self.ids.append(attrs['id'])
        if tag == 'a' and 'href' in attrs:
            self.links.append(attrs['href'])
        if tag == 'img' and 'src' in attrs:
            self.images.append(attrs['src'])
        if tag == 'link' and attrs.get('rel') == 'canonical':
            self.canonicals.append(attrs['href'])
        if tag == 'tr':
            self.table_rows += 1


def require(condition, label):
    if not condition:
        raise ValueError(label)
    checks.append(label)


checks = []
mapping_raw = (D / 'CANDIDATE_MAPPING.json').read_bytes()
mapping = json.loads(mapping_raw)
old_mapping_raw = (ORIGINAL / 'CANDIDATE_MAPPING.json').read_bytes()
require(sha(old_mapping_raw) == mapping['original_mapping_sha256'], 'original candidate mapping preserved')
require(len(mapping['mapping']) == 15, 'exactly fifteen intended target candidates')
require(mapping['new_publication_claimed'] is False and mapping['new_version_doi_claimed'] is False,
        'no new publication or new version DOI claim')
require(mapping['comparison_archive_upload_verified'] is False, 'archive upload remains unverified')
require(mapping['repository_mutations'] == mapping['remote_mutations'] == mapping['network_requests'] == 0,
        'candidate preparation performs no repository or remote mutations')

per_file = []
for row in mapping['mapping']:
    candidate = Path(row['candidate_path'])
    data = candidate.read_bytes()
    require((len(data), sha(data)) == (row['candidate_bytes'], row['candidate_sha256']),
            'candidate bytes pinned: ' + row['repository'] + '/' + row['target_relative_path'])
    original = Path(row['original_candidate_path']).read_bytes()
    require(sha(original) == row['original_candidate_sha256'],
            'original candidate preserved: ' + row['repository'] + '/' + row['target_relative_path'])
    source = Path(row['source_path']) if row['source_path'] is not None else None
    target = Path('/workspace') / row['repository'] / row['target_relative_path']
    if source:
        source_data = source.read_bytes()
        require((len(source_data), sha(source_data)) == (row['source_bytes'], row['source_sha256']),
                'source still matches: ' + row['repository'] + '/' + row['target_relative_path'])
    else:
        require(not target.exists(), 'new target remains absent: ' + row['repository'] + '/' + row['target_relative_path'])
    if candidate.suffix in {'.md', '.html', '.json', '.css'}:
        text = data.decode()
        require('NEW_' not in text and 'VERIFIED_INVENTORY_ROWS' not in text,
                'no unresolved publication placeholder: ' + row['repository'] + '/' + row['target_relative_path'])
        require('10.5281/zenodo.23228395' not in text and 'https://zenodo.org/records/23228395' not in text,
                'unsubmitted draft is not given a published DOI link: ' + row['repository'] + '/' + row['target_relative_path'])
    if candidate.suffix == '.html':
        parser = Page()
        parser.feed(data.decode())
        require(len(parser.ids) == len(set(parser.ids)), 'unique HTML anchors: ' + row['target_relative_path'])
        for href in parser.links:
            require(not href.startswith(('javascript:', 'data:')), 'ordinary navigation link: ' + row['target_relative_path'])
            if href.startswith('#'):
                require(href[1:] in parser.ids, 'local fragment resolves: ' + row['target_relative_path'])
        old_parser = Page()
        old_parser.feed(original.decode())
        require(parser.canonicals == old_parser.canonicals, 'canonical unchanged: ' + row['target_relative_path'])
        for image in parser.images:
            if not urlparse(image).scheme:
                relative = (candidate.parent / image).resolve()
                source_relative = (target.parent / image).resolve()
                require(relative.is_file() or source_relative.is_file(), 'local figure resolves: ' + row['target_relative_path'])
        for ld in re.findall(r'<script type="application/ld\+json">(.*?)</script>', data.decode(), re.S):
            parsed = json.loads(ld)
            require(isinstance(parsed, dict), 'structured data parses: ' + row['target_relative_path'])
    if candidate.suffix == '.json':
        state = json.loads(data)
        require(state['publication']['record_id'] == '23225288' and state['publication']['file_count'] == 27 and
                state['publication']['total_bytes'] == 448387915,
                'current-state latest published inventory remains previous twenty-seven files')
        require(state['comparison_edition']['draft_id'] == '23228395' and
                state['comparison_edition']['status'] == 'UNSUBMITTED_ARCHIVE_UPLOAD_UNVERIFIED' and
                state['comparison_edition']['whole_archive_upload_verified'] is False,
                'current-state draft/upload flags remain unresolved')
        require(state['saved_incoming_state_error'] == 'ENCLOSED_AT_49152_CAPSULE_NODE_OCCURRENCES' and
                state['full_twelve_case_pressure_contacts'] == 'UNRESOLVED' and
                state['higher_dimensional_Big_Bang_origin'] == 'NOT_ESTABLISHED',
                'current-state mathematical and physical claims remain scoped')
    per_file.append({'target': row['repository'] + '/' + row['target_relative_path'],
                     'candidate_sha256': sha(data), 'source_still_matches': True})

for repo in ('HDblast', 'HDblast-archive'):
    original = (ORIGINAL / 'candidates' / repo / 'hdblast/CURRENT_STATUS.md').read_text()
    candidate = (D / 'candidates' / repo / 'hdblast/CURRENT_STATUS.md').read_text()
    require(original.split('## Latest verified Zenodo publication')[0] ==
            candidate.split('## Comparison edition: prepared, with Zenodo publication on hold')[0],
            'all current science and literature status text exactly preserved: ' + repo)
    old_tail = original.split('## Preserved verified edition 23225288: BD prehistory\n', 1)[1]
    new_tail = candidate.split('## Latest independently verified published edition 23225288: BD prehistory\n', 1)[1]
    require(old_tail == new_tail, 'all historical publication and next-calculation text preserved: ' + repo)

for rel in ('site/index.html', 'site/es/index.html'):
    old = (ORIGINAL / 'candidates/HDblast' / rel).read_text()
    new = (D / 'candidates/HDblast' / rel).read_text()
    strip_ld = lambda s: re.sub(r'<script type="application/ld\+json">.*?</script>', '', s, flags=re.S)
    require(strip_ld(old).split('<section id="cite">')[0] == strip_ld(new).split('<section id="cite">')[0],
            'all scientific page content exactly preserved: ' + rel)

shared = Path('/workspace/maldonado-research.github.io')
shared_candidate = D / 'candidates/maldonado-research.github.io'
old_index = (shared / 'index.html').read_text()
new_index = (shared_candidate / 'index.html').read_text()
cards = re.findall(r'<article class="project-card">.*?</article>', old_index, re.S)
new_cards = re.findall(r'<article class="project-card">.*?</article>', new_index, re.S)
require(len(cards) == len(new_cards) == 7, 'shared site retains seven cards')
unrelated = []
for old, new in zip(cards, new_cards):
    if './projects/hdblast/' not in old:
        require(old == new, 'unrelated shared project card exactly preserved')
        unrelated.append(sha(old.encode()))
require(len(unrelated) == 6, 'all six unrelated shared project cards preserved')
require((shared_candidate / 'README.md').read_bytes().startswith((shared / 'README.md').read_bytes()),
        'entire shared original README preserved verbatim as prefix')
require(not any(r['repository'] == 'maldonado-research.github.io' and r['target_relative_path'].startswith('projects/tmd/')
                for r in mapping['mapping']), 'TMD project page excluded from candidate mutation scope')

inventory_path = W / 'publication-planning/read-only-observation/CURRENT_PUBLISHED_INVENTORY.json'
inventory = json.loads(inventory_path.read_bytes())
require(sha(inventory_path.read_bytes()) == mapping['preserved_published_inventory_sha256'], 'dated published inventory input preserved')
page = (D / 'candidates/HDblast/site/zenodo-files/index.html').read_text()
parser = Page()
parser.feed(page)
require(parser.table_rows == 28, 'published inventory table has one header and twenty-seven attachment rows')
attachment_links = [link for link in parser.links if link.startswith('https://zenodo.org/records/23225288/files/')]
require(len(attachment_links) == len(set(attachment_links)) == 27, 'exactly twenty-seven unique published attachment links')
linked_names = {unquote(urlparse(link).path.rsplit('/', 1)[-1]) for link in attachment_links}
require(linked_names == {r['filename'] for r in inventory['files']}, 'exact dated published filenames linked')
for row in inventory['files']:
    require(row['sha256'] in page and row['md5'] in page, 'published SHA256 and MD5 present for ' + row['filename'])
require(sum(r['bytes'] for r in inventory['files']) == 448387915, 'published inventory exact byte sum')
require('2026 at 00:14:31 UTC' in page and 'does not claim another fresh full readback' in page,
        'published full-stream evidence explicitly dated')

data = json.loads((W / 'root/actual-normal-001/DATA.json').read_bytes())
fraction = lambda r: Fraction(int(r['numerator']), int(r['denominator']))
maximum_u = max(fraction(r['suprema']['delta_U_L1_upper']) for r in data['rows'])
coarse_w = max(fraction(r['suprema']['delta_W_L1_upper']) for r in data['rows'] if r['capsule'].endswith('/coarse'))
fine_w = max(fraction(r['suprema']['delta_W_L1_upper']) for r in data['rows'] if r['capsule'].endswith('/fine'))
require(maximum_u < Fraction('1.571e-16') and coarse_w < Fraction('1.914e-14') and fine_w < Fraction('5.776e-16'),
        'scientific display bounds round outward under exact rational arithmetic')
require(data['distinct_capsule_nodes'] == 49152 and data['row_count'] == 12, 'retained-node and prefix counts match actual outputs')

heads = {repo: subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=Path('/workspace') / repo, text=True).strip()
         for repo in ('HDblast', 'HDblast-archive', 'maldonado-research.github.io')}
receipt = {
    'schema_version': 1,
    'status': 'PASS_LOCAL_PUBLICATION_HOLD_CANDIDATE_REPRESENTATION_AND_SCOPE',
    'checks_passed': len(checks),
    'checks': checks,
    'candidate_mapping_sha256': sha(mapping_raw),
    'source_heads_observed': heads,
    'sources_and_candidates': per_file,
    'non_HDBLAST_project_cards_byte_preserved': 6,
    'non_HDBLAST_project_card_sha256': unrelated,
    'shared_README_original_bytes_preserved_as_prefix': True,
    'TMD_page_observed_sha256': sha((shared / 'projects/tmd/index.html').read_bytes()),
    'scientific_data_sha256': sha((W / 'root/actual-normal-001/DATA.json').read_bytes()),
    'new_publication_claimed': False,
    'new_version_doi_claimed': False,
    'archive_upload_verified': False,
    'remote_state_verified_by_this_review': False,
    'fresh_full_stream_readback_claimed': False,
    'required_root_delivery_dependency': 'Create and review hdblast/publication/2026.10.07-stored-state-comparison/README.md as an UNPUBLISHED publication-hold guide before application.',
    'recovery_outcome_or_timestamp_invented': False,
    'repository_mutations': 0,
    'remote_mutations': 0,
    'new_retained_array_decodes': 0,
    'new_physical_source_constructions': 0,
}
(D / 'LOCAL_CANDIDATE_REVIEW.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({k: receipt[k] for k in ('status', 'checks_passed', 'candidate_mapping_sha256', 'new_publication_claimed')}))
