#!/usr/bin/env python3
"""Independent offline typed/rich-text regression controls; never run the live CLI."""
from pathlib import Path
import argparse
import copy
import hashlib
import importlib.util
import json
import socket
import ssl  # Load SSLSocket before installing the socket constructor sentinel.
import sys
import unittest

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
BASE = HERE.parent
SUBJECT_SHA = 'aa7e6cfe271559731ed3b45c5339f8eb612ebe6f017ddafa5fdadad9579e084b'
OLD_SHA = '27fd844cc903a445688a177339b7af633c0a3581f71eeedf7e20bb015bee36a7'

def sha(data):
    return hashlib.sha256(data).hexdigest()

def load(name, path, pin=None):
    data = path.read_bytes()
    if pin is not None and sha(data) != pin:
        raise RuntimeError('REVIEW_SOURCE_PIN_DIFFERS: ' + str(path))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

def forbidden(*args, **kwargs):
    raise RuntimeError('REAL_NETWORK_FORBIDDEN_BY_INDEPENDENT_REVIEW')

socket.socket = forbidden
socket.create_connection = forbidden
M = load('strict_controller_subject', BASE / 'new_edition_controller.py', SUBJECT_SHA)
OLD = load('strict_controller_old_subject', BASE / 'controller-revisions' / OLD_SHA / 'new_edition_controller.py', OLD_SHA)
M.request.urlopen = forbidden
F = load('strict_fixture_basis', BASE / 'controller-review' / 'offline_tests.py')
if F.M.sha(F.SOURCE.read_bytes()) != SUBJECT_SHA:
    raise RuntimeError('FIXTURE_SUBJECT_PIN_DIFFERS')
F.M.request.urlopen = forbidden

def baseline():
    return {
        'metadata': {
            'title': 'Independent manufactured title',
            'description': '<b>Alpha</b> <b>Beta &amp; gamma.</b>',
            'resource_type': {'id': 'software'},
            'rights': [{'id': 'cc-by-4.0'}],
            'languages': [{'id': 'eng'}],
            'creators': [{'person_or_org': {'type': 'personal', 'name': 'Independent fixture'}, 'role': {'id': 'author'}}],
            'contributors': [{'person_or_org': {'type': 'personal', 'name': 'Independent contributor'}, 'role': {'id': 'other'}}],
            'additional_descriptions': [{'description': '<p>Preserved old note &amp; precision.</p>', 'type': {'id': 'notes'}, 'lang': {'id': 'eng'}}],
            'related_identifiers': [{'identifier': '10.5281/zenodo.123', 'scheme': 'doi', 'relation_type': {'id': 'references'}, 'resource_type': {'id': 'software'}}],
            'subjects': [{'id': 'test', 'subject': 'physics'}],
        },
        'custom_fields': {'integer': 1, 'real': 1.5, 'boolean': False, 'nullable': None,
                          'deep': [{'kept': {'integer': 1, 'boolean': True, 'nullable': None}}],
                          'code:programmingLanguage': [{'id': 'python'}], 'code:developmentStatus': {'id': 'active'}},
        'access': {'record': 'public', 'files': 'public', 'embargo': {'active': False, 'reason': None}},
    }

def assign(obj, path, value):
    node = obj
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = value

def delete(obj, path):
    node = obj
    for key in path[:-1]:
        node = node[key]
    del node[path[-1]]

ROWS = []
OLD_ROWS = []

def check(label, fn):
    try:
        details = fn()
        ROWS.append({'label': label, 'status': 'PASS', 'details': details})
    except Exception as exc:
        ROWS.append({'label': label, 'status': 'FAIL', 'exception': type(exc).__name__, 'message': str(exc)})

def metadata_result(module, actual, expected):
    try:
        return module.metadata_equal(actual, expected) is True
    except Exception:
        return False

def metadata_case(label, mutation=None, expect=True, expected=None, old_must_accept=False):
    expected = copy.deepcopy(expected if expected is not None else baseline())
    actual = copy.deepcopy(expected)
    if mutation:
        mutation(actual)
    def execute():
        observed = metadata_result(M, actual, expected)
        if observed is not expect:
            raise RuntimeError('new comparator unexpectedly ' + ('accepted' if observed else 'rejected'))
        old_accepts = metadata_result(OLD, actual, expected)
        if old_must_accept and not old_accepts:
            raise RuntimeError('archived vulnerability did not reproduce')
        if old_must_accept:
            OLD_ROWS.append({'label': label, 'old_accepted': True, 'new_rejected': not observed})
        return {'expected': 'ACCEPT' if expect else 'REJECT', 'old_accepted': old_accepts}
    check(label, execute)

def file_case(label, document, expect, old_must_accept=False):
    def accepted(module):
        try:
            result = module.file_map(copy.deepcopy(document))
            return isinstance(result, dict)
        except Exception:
            return False
    def execute():
        new_accepts = accepted(M)
        if new_accepts is not expect:
            raise RuntimeError('new file map accepted=' + str(new_accepts))
        old_accepts = accepted(OLD)
        if old_must_accept and not old_accepts:
            raise RuntimeError('old file-map contradiction did not reproduce')
        if old_must_accept:
            OLD_ROWS.append({'label': label, 'old_accepted': True, 'new_rejected': not new_accepts})
        return {'expected': 'ACCEPT' if expect else 'REJECT', 'old_accepted': old_accepts}
    check(label, execute)

def run_metadata_controls():
    metadata_case('Complete typed metadata identity')
    metadata_case('Entity spelling with visible inline space preserved', lambda a: assign(a, ('metadata', 'description'), '<b>Alpha</b> <b>Beta &#38; gamma.</b>'))
    metadata_case('ASCII document-edge whitespace only', lambda a: assign(a, ('metadata', 'description'), ' \t\r\n\f' + a['metadata']['description'] + '\n\t '))
    metadata_case('Server access status decoration only', lambda a: a['access'].__setitem__('status', 'open'))
    metadata_case('Metadata readback allows unrelated record wrapper fields', lambda a: a.__setitem__('id', '24000000'))
    e = baseline(); e['metadata']['description'] = '<p>First.</p>\n<p>Second.</p>'
    metadata_case('Paragraph newline retained with equivalent entities', lambda a: assign(a, ('metadata', 'description'), '<p>First&#46;</p>\n<p>Second.</p>'), expected=e)
    metadata_case('Nested note entity spelling preserved', lambda a: assign(a, ('metadata', 'additional_descriptions', 0, 'description'), '<p>Preserved old note &#38; precision.</p>'))
    e = baseline(); e['metadata']['description'] = '<b>Alpha</b>&nbsp;<b>Beta</b>'
    metadata_case('NBSP literal and character reference equivalent', lambda a: assign(a, ('metadata', 'description'), '<b>Alpha</b>\u00a0<b>Beta</b>'), expected=e)

    vocab_paths = [
        ('metadata','resource_type'), ('metadata','rights',0), ('metadata','languages',0),
        ('metadata','creators',0,'role'), ('metadata','contributors',0,'role'),
        ('metadata','related_identifiers',0,'relation_type'), ('metadata','related_identifiers',0,'resource_type'),
        ('metadata','additional_descriptions',0,'type'), ('metadata','additional_descriptions',0,'lang'),
        ('custom_fields','code:programmingLanguage',0), ('custom_fields','code:developmentStatus'),
    ]
    for path in vocab_paths:
        def decorate(a, path=path):
            node = a
            for key in path:
                node = node[key]
            node.update({'title': {'en': 'Server presentation'}, 'description': {'en': 'Vocabulary prose'}, 'icon': 'server', 'props': {'invisible': False}})
        metadata_case('Listed vocabulary decorations at ' + '/'.join(map(str,path)), decorate)
        metadata_case('Listed vocabulary still rejects unlisted key at ' + '/'.join(map(str,path)), lambda a,p=path: assign(a,p+('unknown',),'unreviewed'), expect=False)

    event_collision = json.loads(M.canonical_bytes(M.normalized(baseline()['metadata']['description'], baseline()['metadata']['description'], ('metadata','description'))))
    metadata_case('Top description string/event-list collision', lambda a: assign(a, ('metadata','description'), event_collision), expect=False, old_must_accept=True)
    note = baseline()['metadata']['additional_descriptions'][0]['description']
    note_collision = json.loads(M.canonical_bytes(M.normalized(note,note,('metadata','additional_descriptions',0,'description'))))
    metadata_case('Nested Notes string/event-list collision', lambda a: assign(a, ('metadata','additional_descriptions',0,'description'), note_collision), expect=False, old_must_accept=True)
    metadata_case('Inline visible interword space removed', lambda a: assign(a, ('metadata','description'), '<b>Alpha</b><b>Beta &amp; gamma.</b>'), expect=False, old_must_accept=True)
    metadata_case('Inline visible interword space doubled', lambda a: assign(a, ('metadata','description'), '<b>Alpha</b>  <b>Beta &amp; gamma.</b>'), expect=False, old_must_accept=True)
    e = baseline(); e['metadata']['description'] = '<p style="display:inline">Alpha</p> <p style="display:inline">Beta</p>'
    metadata_case('Inline-styled paragraph space deletion', lambda a: assign(a, ('metadata','description'), '<p style="display:inline">Alpha</p><p style="display:inline">Beta</p>'), expect=False, expected=e, old_must_accept=True)
    e = baseline(); e['metadata']['description'] = '<p>First.</p>\n<p>Second.</p>'
    metadata_case('Interior paragraph newline deletion', lambda a: assign(a, ('metadata','description'), '<p>First.</p><p>Second.</p>'), expect=False, expected=e, old_must_accept=True)
    metadata_case('Paragraph newline replacement with space', lambda a: assign(a, ('metadata','description'), '<p>First.</p> <p>Second.</p>'), expect=False, expected=e, old_must_accept=True)
    e = baseline(); e['metadata']['description'] = '\u00a0<b>Alpha</b>'
    metadata_case('Root leading NBSP deletion', lambda a: assign(a, ('metadata','description'), '<b>Alpha</b>'), expect=False, expected=e, old_must_accept=True)
    e = baseline(); e['metadata']['description'] = '<b>Alpha</b>\u00a0'
    metadata_case('Root trailing NBSP deletion', lambda a: assign(a, ('metadata','description'), '<b>Alpha</b>'), expect=False, expected=e, old_must_accept=True)
    e = baseline(); e['metadata']['description'] = '<b>Alpha</b>&nbsp;<b>Beta</b>'
    metadata_case('Interior NBSP replaced with ASCII space', lambda a: assign(a, ('metadata','description'), '<b>Alpha</b> <b>Beta</b>'), expect=False, expected=e, old_must_accept=True)

    for insertion in ('<?injected instructions?>', '<!DOCTYPE invented>', '<![CDATA[injected text]]>'):
        metadata_case('Added declaration or PI ' + insertion, lambda a,x=insertion: assign(a,('metadata','description'),x+a['metadata']['description']), expect=False, old_must_accept=True)
    for closing in ('</>', '</b junk>', '</ b>', '</b/ >'):
        e = baseline(); e['metadata']['description'] = '<b>Alpha</b>'
        metadata_case('Malformed HTML closer ' + closing, lambda a,x=closing: assign(a,('metadata','description'), '<b>Alpha'+x), expect=False, expected=e)
    metadata_case('Comment insertion remains visible to comparator', lambda a: assign(a,('metadata','description'), '<!-- injected -->'+a['metadata']['description']), expect=False)
    metadata_case('Interior paragraph text precision changes', lambda a: assign(a,('metadata','additional_descriptions',0,'description'), '<p>Preserved old note &amp; imprecision.</p>'), expect=False)
    metadata_case('Inline element attribute mutation', lambda a: assign(a,('metadata','description'), '<b class="injected">Alpha</b> <b>Beta &amp; gamma.</b>'), expect=False)

    typed_cases = [
        (('custom_fields','integer'), True), (('custom_fields','integer'),1.0), (('custom_fields','boolean'),0),
        (('custom_fields','real'),1), (('custom_fields','nullable'),'null'), (('custom_fields','deep'),{}),
        (('custom_fields','deep',0,'kept','integer'),True), (('custom_fields','deep',0,'kept','boolean'),1),
        (('access','embargo','active'),0), (('access','embargo','reason'),False),
        (('metadata','title'),None), (('metadata','title'),['Independent manufactured title']),
        (('metadata','resource_type','id'),True), (('metadata','additional_descriptions'),{}),
        (('metadata','additional_descriptions',0,'description'),None),
    ]
    for path,value in typed_cases:
        metadata_case('Exact JSON type ' + '/'.join(map(str,path)) + ' -> ' + type(value).__name__, lambda a,p=path,v=value: assign(a,p,v), expect=False)
    for path in [('metadata','additional_descriptions'),('metadata','creators'),('custom_fields','deep')]:
        metadata_case('List truncation ' + '/'.join(path), lambda a,p=path: assign(a,p,[]), expect=False)
        def append(a, path=path):
            node = a
            for key in path:
                node = node[key]
            node.append(copy.deepcopy(node[0]))
        metadata_case('List expansion ' + '/'.join(path), append, expect=False)
    for path in [('metadata','title'),('metadata','additional_descriptions',0,'description'),('custom_fields','nullable'),('custom_fields','deep',0,'kept','nullable')]:
        metadata_case('Missing nested key ' + '/'.join(map(str,path)), lambda a,p=path: delete(a,p), expect=False)
    metadata_case('Subject ID decoration outside vocabulary rejected', lambda a: assign(a,('metadata','subjects',0,'title'),'Unreviewed'), expect=False)
    metadata_case('Vocabulary ID cannot be replaced', lambda a: assign(a,('metadata','rights',0,'id'),'cc0-1.0'), expect=False)
    metadata_case('Unknown access policy key rejected', lambda a: assign(a,('access','unknown'),'private'), expect=False)
    for key in ('metadata','custom_fields','access'):
        e = baseline()
        if key in ('metadata','custom_fields'):
            e[key] = {}
        else:
            e[key] = {}
        metadata_case('Required outer key omission ' + key, lambda a,k=key:a.pop(k), expect=False, expected=e, old_must_accept=True)
    e = baseline(); e['access']['nullpolicy'] = None
    metadata_case('Required null access policy omitted', lambda a:a['access'].pop('nullpolicy'), expect=False, expected=e, old_must_accept=True)
    e = baseline(); e['access']['embargo'] = None
    metadata_case('Required null embargo omitted', lambda a:a['access'].pop('embargo'), expect=False, expected=e, old_must_accept=True)

    def structure_precedes_parser():
        expected = baseline(); actual = copy.deepcopy(expected)
        actual['metadata']['description'] = event_collision
        original = M.HtmlEvents
        class ForbiddenParser:
            def __init__(self):
                raise RuntimeError('HTML parser instantiated before shape rejection')
        M.HtmlEvents = ForbiddenParser
        try:
            try:
                M.metadata_equal(actual,expected)
            except M.Stop as exc:
                if str(exc) != 'COMPLETE_METADATA_READBACK_DIFFERS':
                    raise
                return {'parser_instantiated':False}
            raise RuntimeError('shape collision accepted')
        finally:
            M.HtmlEvents = original
    check('Recursive shape rejection occurs before rich-text parsing', structure_precedes_parser)

    class NotJsonStr(str): pass
    class NotJsonDict(dict): pass
    for value in (NotJsonStr('Independent manufactured title'), NotJsonDict({'id':'software'})):
        path = ('metadata','title') if isinstance(value,str) else ('metadata','resource_type')
        metadata_case('Builtin-only JSON type rejects ' + type(value).__name__, lambda a,p=path,v=value: assign(a,p,v), expect=False)
    e = baseline(); e['custom_fields']['nullable'] = {'string': 1}
    metadata_case('Nonstring object key rejects', lambda a:assign(a,('custom_fields','nullable'),{1:1}), expect=False, expected=e)

def run_file_controls():
    file_case('Mapping key equals nested key and filename', {'files':{'entries':{'a.txt':{'key':'a.txt','filename':'a.txt'}}}}, True)
    file_case('Mapping key equals sole nested filename', {'entries':{'a.txt':{'filename':'a.txt'}}}, True)
    file_case('List representation keys and filenames consistent', {'files':[{'key':'a.txt','filename':'a.txt'}]}, True)
    file_case('Mapping key contradicts nested key', {'files':{'entries':{'a.txt':{'key':'b.txt'}}}}, False, True)
    file_case('Mapping key contradicts nested filename', {'entries':{'a.txt':{'filename':'b.txt'}}}, False, True)
    file_case('Nested key contradicts filename under dict', {'entries':{'a.txt':{'key':'a.txt','filename':'b.txt'}}}, False, True)
    file_case('Nested key contradicts filename under list', {'files':[{'key':'a.txt','filename':'b.txt'}]}, False, True)
    file_case('Missing nested identifier under matching mapping', {'entries':{'a.txt':{}}}, False)
    file_case('Nonstring mapping key rejected', {'entries':{1:{'key':'a.txt'}}}, False, True)
    file_case('Nonobject mapping row rejected', {'entries':{'a.txt':['a.txt']}}, False)
    file_case('Nonobject list row rejected', {'entries':['a.txt']}, False)
    file_case('Duplicate filename rejected', {'entries':[{'key':'a.txt'},{'filename':'a.txt'}]}, False)

class IndependentPublishControls(F.ControllerFixtures):
    def test_description_collision_blocks_publish_before_latch(self):
        self.prepared()
        original = self.server.document
        def collision(prior=False,public=False):
            obj=original(prior=prior,public=public)
            if not prior:
                description=obj['metadata']['description']
                obj['metadata']['description']=json.loads(F.M.canonical_bytes(F.M.normalized(description,description,('metadata','description'))))
            return obj
        self.server.document=collision
        self.stop('COMPLETE_METADATA_READBACK_DIFFERS',self.controller.publish)
        self.assertFalse(self.controller.state.get('publish_post_attempted',False))
        self.assertFalse(self.posts('/actions/publish'))

    def test_empty_customfields_omission_blocks_publish_before_latch(self):
        self.prepared()
        original=self.server.document
        def missing(prior=False,public=False):
            obj=original(prior=prior,public=public)
            if not prior:obj.pop('custom_fields')
            return obj
        self.server.document=missing
        self.stop('METADATA_READBACK_REQUIRED_FIELDS_MISSING',self.controller.publish)
        self.assertFalse(self.controller.state.get('publish_post_attempted',False))
        self.assertFalse(self.posts('/actions/publish'))

    def test_mapping_contradiction_blocks_stream_and_publish_latch(self):
        self.prepared()
        count=len(self.server.streams)
        name=next(iter(self.server.files));row=self.server.files.pop(name)
        self.server.files['contradictory-mapping-name']=row
        self.stop('REMOTE_ENTRY_MAPPING_KEY_DIFFERS',self.controller.publish)
        self.assertEqual(len(self.server.streams),count)
        self.assertFalse(self.controller.state.get('publish_post_attempted',False))
        self.assertFalse(self.posts('/actions/publish'))

def run_flow_controls():
    names=[n for n in IndependentPublishControls.__dict__ if n.startswith('test_')]
    names += [
        'test_latest_only_owner_listing_complete_create_prepare_publish_flow',
        'test_raw_verified_etag_preserved_and_publish_header_is_decimal',
        'test_invalid_publish_etag_stops_before_attempt_latch_and_post',
    ]
    suite=unittest.TestSuite(IndependentPublishControls(n) for n in names)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    return {'tests_run':result.testsRun,'status':'PASS' if result.wasSuccessful() else 'FAIL',
            'failures':[str(x[0])+': '+x[1] for x in result.failures],
            'errors':[str(x[0])+': '+x[1] for x in result.errors],'skipped':len(result.skipped)}

def genuine_prior_control():
    prior_path=BASE/'read-only-observation/CURRENT_PRIOR.raw.json'
    expected_path=Path('/workspace/hdblast-research-work/continuation-next-20261007/release-planning/FINAL_METADATA_MODERN.json')
    prior_raw=prior_path.read_bytes();expected_raw=expected_path.read_bytes()
    if sha(prior_raw)!='261d81302396c1683f5c0e1185b071189cf83bbb5fde5fb4a098ca64e8fe47ca':
        raise RuntimeError('GENUINE_PRIOR_PIN_DIFFERS')
    if sha(expected_raw)!='72415e84db788b8dec0c421012c1d264fc669e6f8e5d08d93e80f0075e3485ba':
        raise RuntimeError('GENUINE_EXPECTED_METADATA_PIN_DIFFERS')
    prior=json.loads(prior_raw);expected=json.loads(expected_raw)
    if M.metadata_equal(prior,expected) is not True:
        raise RuntimeError('genuine baseline rejected')
    if prior['metadata']['description'] != expected['metadata']['description']:
        raise RuntimeError('genuine raw description differs')
    if [x['description'] for x in prior['metadata']['additional_descriptions']] != [x['description'] for x in expected['metadata']['additional_descriptions']]:
        raise RuntimeError('genuine raw Notes differs')
    return {'prior_sha256':sha(prior_raw),'expected_metadata_sha256':sha(expected_raw),
            'raw_description':'EXACT_EQUAL','raw_all_notes':'EXACT_EQUAL','metadata_equal':True,
            'network_calls':0}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    run_metadata_controls();run_file_controls()
    check('Genuine fixed prior complete metadata remains accepted',genuine_prior_control)
    flows=run_flow_controls()
    current_sha=sha((BASE/'new_edition_controller.py').read_bytes())
    if current_sha!=SUBJECT_SHA:
        raise RuntimeError('SOURCE_CHANGED_DURING_REVIEW')
    status='PASS' if all(x['status']=='PASS' for x in ROWS) and flows['status']=='PASS' else 'FAIL'
    report={'status':status,'optimized':bool(sys.flags.optimize),'controller_sha256':current_sha,
            'old_controller_sha256':OLD_SHA,'suite_sha256':sha(Path(__file__).read_bytes()),
            'targeted_controls':len(ROWS),'targeted_passed':sum(x['status']=='PASS' for x in ROWS),
            'targeted_failed':sum(x['status']=='FAIL' for x in ROWS),'controls':ROWS,
            'archived_vulnerabilities_reproduced':OLD_ROWS,'archived_vulnerability_count':len(OLD_ROWS),
            'flow_controls':flows,'real_network':'FORBIDDEN_BY_SOCKET_AND_URLOPEN_SENTINELS',
            'real_remote_mutations':0,'live_cli_invocations':0,'scientific_callbacks':0,
            'basis':'Typed manufactured metadata; 27 small old + 2 small new fake payloads; cached pinned genuine metadata only'}
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:report[k] for k in ['status','optimized','targeted_controls','targeted_passed','targeted_failed','archived_vulnerability_count','flow_controls']},indent=2))
    return 0 if status=='PASS' else 1

if __name__=='__main__':
    raise SystemExit(main())
