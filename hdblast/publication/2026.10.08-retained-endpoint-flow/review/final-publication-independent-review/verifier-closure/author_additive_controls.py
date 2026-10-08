#!/usr/bin/env python3
"""Adapt controls in a new directory; manufacture no new publication claim."""
from pathlib import Path
import difflib
import hashlib

HERE = Path(__file__).absolute().parent
OLD = Path('/workspace/HDblast/hdblast/publication/2026.10.07-stored-state-comparison/test_manufactured_witnesses.py')
before = OLD.read_text(); source = before
source = source.replace('PRIOR27_CONTENT_RECEIPT.json', 'PRIOR29_CONTENT_RECEIPT.json')
source = source.replace('manufactured_29_witness_verifier', 'manufactured_31_witness_verifier')
source = source.replace("'expected_total_files_final': 29", "'expected_total_files_final': 31")
source = source.replace("seeds['23225288-record.json']", "seeds['23228395-record.json']")
source = source.replace('Explicitly fabricated 29-file verifier fixture', 'Explicitly fabricated 31-file verifier fixture')
source = source.replace("metadata['metadata']['version'] = 'FABRICATED_ONLY_29_FILE_CONTROL'", "metadata['metadata']['version'] = v.EXPECTED_VERSION\n    metadata['metadata']['additional_descriptions'][0]['description'] = '<p>Manufactured current Notes only.</p>' + v.NOTES_HISTORY_MARKER + seed_prior['metadata']['additional_descriptions'][0]['description']")
source = source.replace("parent={'id': v.PARENT},", "parent={'id': v.PARENT, 'access': {'owned_by': {'user': v.EXPECTED_OWNER}},\n                              'pids': {'doi': {'identifier': '10.5281/zenodo.' + v.PARENT}}},\n                      versions={'is_latest': True},")
source = source.replace('PASS_COMPLETE_29_FILE_EDITION_READBACK', 'PASS_COMPLETE_31_FILE_EDITION_READBACK')
source = source.replace("'file_count': 29", "'file_count': 31")
source = source.replace('PASS_OFFLINE_29_FILE_PUBLICATION_WITNESSES', 'PASS_OFFLINE_31_FILE_PUBLICATION_WITNESSES')
source = source.replace("result['fresh_complete_public_streams'] == 29", "result['fresh_complete_public_streams'] == 31")
source = source.replace('all29', 'all31').replace('prior27', 'prior29')
source = source.replace("expected_total_files_final', 29.0", "expected_total_files_final', 31.0")
source = source.replace('manufactured-29-publication-', 'manufactured-31-publication-')
source = source.replace('PASS_MANUFACTURED_OFFLINE_29_FILE_PUBLICATION_CONTROLS', 'PASS_MANUFACTURED_OFFLINE_31_FILE_PUBLICATION_CONTROLS')
needle = "        return [copy.deepcopy(inventory), copy.deepcopy(metadata), record, files, copy.deepcopy(record),\n                verified, copy.deepcopy(prior), journal, preserved, ID, inventory_sha, metadata_sha]"
addition = '''        draft = copy.deepcopy(record)
        draft['is_draft'] = True; draft['is_published'] = False
        draft['versions'] = {'is_latest_draft': True}; draft['pids'] = {}
        draft['revision_id'] = 1
        draft['links'] = {'self': PREFIX + '/draft', 'files': PREFIX + '/draft/files'}
        draft_listing = copy.deepcopy(files)
        for item in draft_listing['entries']:
            item['links']['content'] = PREFIX + '/draft/files/' + v.quote(item['key'], safe='') + '/content'
        for item in draft['files']['entries'].values():
            item['links']['content'] = PREFIX + '/draft/files/' + v.quote(item['key'], safe='') + '/content'
        draft_verified = copy.deepcopy(verified); draft_verified['published'] = False; draft_verified['etag'] = '"1"'
        draft_journal = []
        for round_index in range(2):
            draft_journal.append({'kind': 'GET', 'url': PREFIX + '/draft/files', 'status': 200,
                                  'authenticated_request': True})
            for pin in old + new:
                draft_journal.append({'kind': 'FULL_CONTENT_STREAM_VERIFIED', 'filename': pin['filename'],
                                     'observed': {key: pin[key] for key in ('bytes', 'md5', 'sha256')}})
        draft_journal.append({'kind': 'WRITE_INTENT', 'operation': 'publish', 'url': PREFIX + '/draft/actions/publish'})
        return [copy.deepcopy(inventory), copy.deepcopy(metadata), record, files, copy.deepcopy(record),
                verified, copy.deepcopy(prior), draft_journal + journal, preserved, ID, inventory_sha, metadata_sha,
                draft, draft_listing, draft_verified]'''
if source.count(needle) != 1: raise ValueError('Fixture adaptation needle differs')
source = source.replace(needle, addition)
source = source.replace("    def embedded(values, record_index=2):", "    def public_marker(values):\n        return next(i for i, e in enumerate(values[7]) if e.get('url') == PREFIX + '/files')\n    def embedded(values, record_index=2):")
source = source.replace("'terminal public GET absent', lambda x: x[7].pop(0)", "'terminal public GET absent', lambda x: x[7].pop(public_marker(x))")
source = source.replace("'terminal GET authenticated', lambda x: x[7][0].__setitem__('authenticated_request', True)", "'terminal GET authenticated', lambda x: x[7][public_marker(x)].__setitem__('authenticated_request', True)")
source = source.replace("'stream appears before latest listing', lambda x: x[7].append(copy.deepcopy(x[7][0]))", "'stream appears before latest listing', lambda x: x[7].append(copy.deepcopy(x[7][public_marker(x)]))")
source = source.replace("'one inherited stream missing', lambda x: x[7].pop(1)", "'one inherited stream missing', lambda x: x[7].pop(public_marker(x) + 1)")
for label, field, value in [('stream wrong byte count','bytes','1'), ('stream wrong MD5','md5',"'0' * 32"), ('stream wrong SHA256','sha256',"'0' * 64")]:
    source = source.replace("'"+label+"', lambda x: x[7][1]['observed'].__setitem__('"+field+"', "+value+")", "'"+label+"', lambda x: x[7][public_marker(x) + 1]['observed'].__setitem__('"+field+"', "+value+")")
source = source.replace("        name = 'journal.jsonl'; raw =", "        for role, value in zip(('new_draft_record', 'new_draft_files', 'verified_draft'), values[12:]):\n            name = role + '.json'; raw = canonical(value)\n            (root / name).write_bytes(raw); roles[role] = name\n            files[name] = {'bytes': len(raw), 'sha256': sha(raw)}\n        name = 'journal.jsonl'; raw =")
new_controls = '''    for index in (2, 4, 12):
        reject('strict owner wrong at role' + str(index), lambda x, i=index: x[i]['parent']['access']['owned_by'].__setitem__('user', '1386318'))
        reject('strict owner boolean at role' + str(index), lambda x, i=index: x[i]['parent']['access']['owned_by'].__setitem__('user', True))
        reject('strict owner omitted at role' + str(index), lambda x, i=index: x[i]['parent'].pop('access'))
        reject('strict concept DOI wrong at role' + str(index), lambda x, i=index: x[i]['parent']['pids']['doi'].__setitem__('identifier', '10.5281/zenodo.22922927'))
        reject('strict latest flag false at role' + str(index), lambda x, i=index: x[i]['versions'].__setitem__('is_latest_draft' if i == 12 else 'is_latest', False))
        reject('strict latest flag integer at role' + str(index), lambda x, i=index: x[i]['versions'].__setitem__('is_latest_draft' if i == 12 else 'is_latest', 1))
    reject('wrong planned semantic version', lambda x: x[1]['metadata'].__setitem__('version', 'FABRICATED_WRONG_VERSION'))
    reject('prior Notes exact suffix omitted', lambda x: x[1]['metadata']['additional_descriptions'][0].__setitem__('description', '<p>No historical Notes.</p>'))
    reject('prior Notes suffix replaced', lambda x: x[1]['metadata']['additional_descriptions'][0].__setitem__('description', x[1]['metadata']['additional_descriptions'][0]['description'] + ' '))
    reject('protected creator changed', lambda x: x[1]['metadata']['creators'][0]['person_or_org'].__setitem__('name', 'Other owner'))
    reject('protected inherited relation changed', lambda x: x[1]['metadata']['related_identifiers'][0].__setitem__('identifier', '10.5281/zenodo.1'))
    reject('unreviewed related URL appended', lambda x: x[1]['metadata']['related_identifiers'].append({'identifier': 'https://example.org', 'scheme': 'url', 'relation_type': {'id': 'references'}, 'resource_type': {'id': 'software'}}))
    reject('draft ID wrong', lambda x: x[12].__setitem__('id', v.PRIOR))
    reject('draft marked public', lambda x: x[12].__setitem__('is_published', True))
    reject('draft URL wrong', lambda x: x[12]['links'].__setitem__('self', PREFIX))
    reject('draft assigned DOI wrong', lambda x: x[12].__setitem__('pids', {'doi': {'identifier': '10.5281/zenodo.1'}}))
    reject('draft revision floating', lambda x: x[12].__setitem__('revision_id', 1.0))
    reject('draft revision boolean', lambda x: x[12].__setitem__('revision_id', True))
    reject('draft ETag weak', lambda x: x[14].__setitem__('etag', 'W/"1"'))
    reject('draft ETag padded', lambda x: x[14].__setitem__('etag', '"01"'))
    reject('draft ETag revision disagrees', lambda x: x[14].__setitem__('etag', '"2"'))
    reject('draft verification marked public', lambda x: x[14].__setitem__('published', True))
    reject('draft verification wrong inventory pin', lambda x: x[14].__setitem__('inventory_sha256', '0' * 64))
    reject('draft complete content row omitted', lambda x: x[14]['content'].pop())
    reject('draft inherited stream reused', lambda x: x[14]['content'][0].__setitem__('basis', 'SAME_IMMUTABLE_FILE_AND_CONTENT_VERSION_AS_PRIOR_FULL_STREAM'))
    reject('draft fresh row wrong version', lambda x: x[14]['content'][0].__setitem__('remote_version_id', 'wrong'))
    reject('draft file listing wrong content URL', lambda x: x[13]['entries'][0]['links'].__setitem__('content', PREFIX + '/files/wrong/content'))
    reject('draft file listing omitted entry', lambda x: x[13]['entries'].pop())
    reject('draft first full round absent', lambda x: x[7].__delitem__(slice(0, 32)))
    reject('draft second full round absent', lambda x: x[7].__delitem__(slice(32, 64)))
    reject('draft first inherited stream absent', lambda x: x[7].pop(1))
    reject('draft first complete stream wrong SHA', lambda x: x[7][1]['observed'].__setitem__('sha256', '0' * 64))
    reject('publish intent absent', lambda x: x[7].__delitem__(64))
    reject('publish intent repeated', lambda x: x[7].insert(64, copy.deepcopy(x[7][64])))
    reject('publish intent wrong URL', lambda x: x[7][64].__setitem__('url', PREFIX + '/actions/publish'))
    reject('mutation after last complete draft round', lambda x: x[7].insert(64, {'kind': 'WRITE_INTENT', 'operation': 'save_metadata', 'url': PREFIX + '/draft'}))
    def later_partial_draft(x):
        x[7][64:64] = [copy.deepcopy(x[7][0]), copy.deepcopy(x[7][1])]
    reject('partial draft round after last complete round', later_partial_draft)
    def missing_one_in_later_draft(x):
        x[7][64:64] = copy.deepcopy(x[7][:31])
    reject('latest draft round missing one after two older complete rounds', missing_one_in_later_draft)
    empty_listing = fixture(); empty_listing[7].insert(64, copy.deepcopy(empty_listing[7][0]))
    checked(empty_listing); positives.append('v2publish final identity-only draft listing after the full stream round')
    assigned = fixture(); assigned[12]['pids'] = {'doi': {'identifier': '10.5281/zenodo.' + ID}}
    checked(assigned); positives.append('draft permits matching assigned DOI or unassigned PID object')

'''
source = source.replace("    with tempfile.TemporaryDirectory(prefix='manufactured-31-publication-'", new_controls + "    with tempfile.TemporaryDirectory(prefix='manufactured-31-publication-'", 1)
compile(source, str(HERE/'test_manufactured_witnesses.py'), 'exec')
(HERE/'test_manufactured_witnesses.py').write_text(source)
(HERE/'CONTROLS_ADAPTATION.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True), source.splitlines(True), fromfile='previous/test_manufactured_witnesses.py', tofile='additive31/test_manufactured_witnesses.py')))
