#!/usr/bin/env python3
"""Offline reviewer controls: accepted Transport, fabricated opener and files only."""
import argparse
import ast
import copy
import hashlib
import importlib.util
import io
import json
import socket
import ssl
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch
from urllib import error

sys.dont_write_bytecode = True
ROOT = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
BASE = ROOT / 'publication-final/empty-archive-cleanup'
SOURCE = BASE / 'delete_empty_archive_once.py'
PIN = 'a4685ff65eab0dca73ecf8066a7710359979eeffea5dffcf8e2f3c4bdb2a1d60'
OLD = BASE / 'workflow-independent-review/candidate_493c1486_snapshot.py'
OLD_PIN = '493c1486f90d6bbde3fa2a9588ef9b1a275e21334da47df1bd07b7c3ab6d80ab'
TOKEN = 'FABRICATED_REVIEW_TOKEN_ONLY'

def check(value, reason):
    if not value:
        raise AssertionError(reason)

def sha(value):
    return hashlib.sha256(value).hexdigest()

def load(path, name, pin):
    raw = path.read_bytes()
    check(sha(raw) == pin, 'reviewed source changed')
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, loader=None, origin=str(path)))
    m.__file__ = str(path)
    sys.modules[name] = m
    exec(compile(raw, str(path), 'exec'), m.__dict__)
    return m

def forbidden(*args, **kwargs):
    raise AssertionError('NETWORK_FORBIDDEN')

def rejected(fn, contains=None):
    try:
        fn()
    except Exception as exc:
        check(not isinstance(exc, AssertionError), 'fixture assertion: ' + str(exc))
        if contains:
            check(contains in str(exc), 'unexpected rejection: ' + str(exc))
        return type(exc).__name__, str(exc)
    raise AssertionError('invalid fixture accepted')

class BrokenRead:
    def read(self, size=-1):
        raise OSError('manufactured body read failure')
    def close(self):
        pass

class RawResponse:
    def __init__(self, url, status=204, body=b'', headers=None):
        self.url, self.status, self.body = url, status, body
        self.headers = headers or {}
    def geturl(self):
        return self.url
    def read(self, size=-1):
        return self.body if size < 0 else self.body[:size]
    def __enter__(self):
        return self
    def __exit__(self, *unused):
        return False

def missing(h, status=400):
    return json.dumps({'status': status, 'message': 'File with key "' + h.NAME + '" is not available.'}).encode()

class Opener:
    def __init__(self, c):
        self.c = c
        self.calls = []
        self.probe_hook = None
        self.delete_kind = 'success'
        self.delete_count = 0
    def open(self, req, timeout):
        c = self.c
        method, url = req.get_method(), req.full_url
        check(timeout == 60 and req.data is None, 'request timeout/body changed')
        check(req.get_header('Authorization') == 'Bearer ' + TOKEN, 'fabricated bearer absent')
        self.calls.append((method, url))
        if method == 'GET':
            check(url in (c.h.URL, c.m.NATIVE_URL), 'unreviewed GET destination')
            if self.probe_hook:
                response = self.probe_hook(url)
                if response is not None:
                    if isinstance(response, Exception):
                        raise response
                    return response
            raise error.HTTPError(url, 400, 'manufactured unavailable', {'Content-Type': 'application/json'}, io.BytesIO(missing(c.h)))
        check(method == 'DELETE' and url == c.m.DELETE_URL, 'unreviewed write destination')
        self.delete_count += 1
        latch = c.side / 'ONE_DELETE_ATTEMPT_LATCH.json'
        check(latch.is_file(), 'DELETE began before persistent latch')
        intent = json.loads(latch.read_text())
        check(intent['url'] == url and intent['file_id'] == c.h.FILE_ID
              and intent['version_id'] == c.h.VERSION_ID and intent['original_pending_cleared'] is False,
              'latch not bound to exact initialized target')
        if self.delete_kind == 'unknown_before':
            raise OSError('manufactured unknown before effect')
        if self.delete_kind == 'unknown_after':
            del c.files[c.h.NAME]
            raise OSError('manufactured unknown after effect')
        if self.delete_kind in ('401', '500'):
            raise error.HTTPError(url, int(self.delete_kind), 'manufactured rejected delete', {}, io.BytesIO(b'rejected'))
        if self.delete_kind == 'read_failure':
            raise error.HTTPError(url, 204, 'manufactured failed body read', {}, BrokenRead())
        if self.delete_kind != 'still_present':
            del c.files[c.h.NAME]
        return RawResponse(url)

class Controller:
    def __init__(self, m, h, module, side):
        self.m, self.h, self.module, self.side, self.token = m, h, module, side, TOKEN
        self.old = [{'filename': 'manufactured-inherited-' + str(i), 'bytes': i+1,
                     'md5': hashlib.md5(('manufactured-' + str(i)).encode()).hexdigest()} for i in range(27)]
        self.files = {p['filename']: {'key': p['filename'], 'file_id': 'manufactured-file-' + str(i),
                      'version_id': 'manufactured-version-' + str(i), 'bucket_id': h.BUCKET_ID,
                      'size': p['bytes'], 'checksum': 'md5:' + p['md5'], 'status': 'completed'}
                      for i, p in enumerate(self.old)}
        self.files[h.NAME] = {'key': h.NAME, 'file_id': h.FILE_ID, 'version_id': h.VERSION_ID,
                             'bucket_id': h.BUCKET_ID, 'status': 'pending',
                             'links': {'self': m.DELETE_URL, 'content': h.URL,
                                       'commit': h.URL.removesuffix('/content') + '/commit'}}
        self.initialized = copy.deepcopy(self.files)
        self.expected = {'metadata': {'title': 'Fabricated review only', 'description': '<p>Two words</p>'},
                         'custom_fields': {}, 'access': {'record': 'public', 'files': 'public'}}
        self.document = copy.deepcopy(self.expected)
        self.state = {'draft_id': h.DRAFT, 'published': False,
                      'continuation_binding': {'controller_sha256': h.CONTROLLER_SHA},
                      'inventory_sha256': h.INVENTORY_SHA, 'metadata_sha256': h.METADATA_SHA,
                      'create_post_attempted': True, 'publish_post_attempted': False,
                      'content_verified_files': [], 'pending': {'kind': 'put_content', 'url': h.URL,
                      'body_sha256': h.ASSET_PIN['sha256'], 'body_bytes': h.ASSET_PIN['bytes'],
                      'context': {'filename': h.NAME}, 'utc': 'MANUFACTURED'}}
        self.original_state = copy.deepcopy(self.state)
        self.guard_calls = []
        self.current_hook = None
        self.owner = {'id': h.DRAFT}
        self.transport = module.Transport(TOKEN)
        self.opener = Opener(self)
        self.transport.opener = self.opener
    def current(self):
        self.guard_calls.append('current')
        if self.current_hook:
            self.current_hook(self)
    def owner_list(self):
        self.guard_calls.append('owner')
        return self.owner
    def draft(self, draft_id):
        check(draft_id == self.h.DRAFT, 'wrong owned draft')
        self.guard_calls.append('draft')
        return self.document, {}, copy.deepcopy(self.files)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = []
    def test(name, fn):
        try:
            details = fn()
        except Exception as exc:
            rows.append({'name': name, 'status': 'FAIL', 'exception': type(exc).__name__, 'reason': str(exc)})
        else:
            rows.append({'name': name, 'status': 'PASS', 'details': details})
    with patch.object(socket, 'socket', forbidden), patch.object(socket, 'create_connection', forbidden):
        m = load(SOURCE, '_reviewer_cleanup_current', PIN)
        old = load(OLD, '_reviewer_cleanup_archived', OLD_PIN)
        h = m.load_reviewed()
        module = h.load_controller()
        def case(fn):
            with tempfile.TemporaryDirectory(prefix='offline-cleanup-review-') as tmp:
                base = Path(tmp); side = base / 'side'; side.mkdir()
                evidence = base / 'evidence'; evidence.mkdir()
                for name in ('STATE.json', 'JOURNAL.jsonl'):
                    (evidence / name).write_bytes(('MANUFACTURED ' + name + '\n').encode())
                before = {name: sha((evidence / name).read_bytes()) for name in ('STATE.json', 'JOURNAL.jsonl')}
                c = Controller(m, h, module, side)
                with patch.object(m, 'EVIDENCE', evidence):
                    result = fn(c, side, before)
                    check(c.state == c.original_state, 'logical original state/pending changed')
                    m.original_unchanged(before)
                    return result
        def probe(c, side, url=None):
            return m.unavailable(c, h, side, TOKEN, 'MANUFACTURED', m.NATIVE_URL if url is None else url)
        def native_case(response=None, expect=None):
            def run(c, side, before):
                if response is not None:
                    c.opener.probe_hook = lambda url: response(url)
                if expect is None:
                    probe(c, side)
                else:
                    rejected(lambda: probe(c, side), expect)
                check(c.opener.calls == [('GET', m.NATIVE_URL)], 'native GET repeated or changed destination')
                check(c.opener.delete_count == 0, 'GET probe mutated')
            return case(run)
        test('archived_native_route_rejected_by_real_transport_before_opener', lambda: case(lambda c,s,b:
             (rejected(lambda: old.unavailable(c,h,s,TOKEN,'ARCHIVED',old.NATIVE_URL), 'UNSAFE_QUERY'),
              check(c.opener.calls == [], 'archived rejected route reached opener'))))
        test('accepted_generic_transport_still_rejects_version_query', lambda: case(lambda c,s,b:
             (rejected(lambda: c.transport.request('GET', m.NATIVE_URL), 'UNSAFE_QUERY'),
              check(c.opener.calls == [], 'generic policy relaxed'))))
        test('fixed_native_route_exact_400_with_real_response', lambda: native_case())
        test('fixed_native_route_normal_400_matching_response', lambda: native_case(lambda u: RawResponse(u,400,missing(h))))
        test('native_redirect_302_not_followed_or_accepted', lambda: native_case(lambda u:
             error.HTTPError(u,302,'manufactured redirect',{'Location':'https://example.invalid/'},io.BytesIO()), 'MISSING_CONTENT_GET_NOT_DEFINITIVE_400'))
        test('native_changed_final_url_rejected', lambda: native_case(lambda u: RawResponse(u+'-wrong',400,missing(h)), 'NATIVE_GET_FINAL_URL_DIFFERS'))
        test('native_changed_error_url_rejected', lambda: native_case(lambda u:
             error.HTTPError(u+'-wrong',400,'manufactured',{},io.BytesIO(missing(h))), 'NATIVE_GET_ERROR_URL_DIFFERS'))
        test('native_oversized_error_body_rejected', lambda: native_case(lambda u:
             error.HTTPError(u,400,'manufactured',{},io.BytesIO(b'x'*(h.MAX_JSON+1))), 'MISSING_CONTENT_GET_NOT_DEFINITIVE_400'))
        test('native_oversized_normal_body_rejected', lambda: native_case(lambda u:
             RawResponse(u,400,b'x'*(h.MAX_JSON+1)), 'NATIVE_GET_RESPONSE_TOO_LARGE'))
        test('native_error_body_read_failure_rejected', lambda: native_case(lambda u:
             error.HTTPError(u,400,'manufactured',{},BrokenRead()), 'MISSING_CONTENT_GET_NOT_DEFINITIVE_400'))
        test('native_200_even_matching_error_body_rejected', lambda: native_case(lambda u:
             RawResponse(u,200,missing(h)), 'MISSING_CONTENT_GET_NOT_DEFINITIVE_400'))
        test('native_float_status_body_rejected', lambda: native_case(lambda u:
             error.HTTPError(u,400,'manufactured',{},io.BytesIO(missing(h,400.0))), 'MISSING_CONTENT_ERROR_DIFFERS'))
        def malformed():
            def run(c,s,b):
                c.opener.probe_hook=lambda u:error.HTTPError(u,400,'manufactured',{},io.BytesIO(b'{'))
                rejected(lambda:probe(c,s),'')
                check(c.opener.calls == [('GET',m.NATIVE_URL)],'malformed GET repeated')
            return case(run)
        test('native_malformed_json_rejected',malformed)
        test('unreviewed_probe_url_refused_before_opener', lambda: case(lambda c,s,b:
             (rejected(lambda: probe(c,s,m.NATIVE_URL+'&other=1'),'UNREVIEWED_UNAVAILABLE_PROBE_URL'),
              check(c.opener.calls == [],'unreviewed URL reached opener'))))
        test('accepted_no_redirect_handler_declines_redirect', lambda:
             check(module.NoRedirect().redirect_request(None,None,307,'manufactured',{},'https://example.invalid/') is None,
                   'accepted redirect handler follows redirect'))
        test('both_real_401_failure_chains_and_source_pins_verify_read_only', lambda:
             check(m.verify_two_failures(h)['histories_preserved'] is True, '401 history proof differs'))
        def flow(c, side, before):
            return m.delete_once(c,h,side,TOKEN,c.initialized,c.expected,{'fabricated_workflow':True},before)
        def success(c,s,b):
            result=flow(c,s,b)
            check(result['inherited_files']==27 and result['original_pending_cleared'] is False,'successful result scope')
            check(c.opener.calls==[('GET',h.URL),('GET',m.NATIVE_URL),('DELETE',m.DELETE_URL)],'one exact workflow transport sequence changed')
            check(c.guard_calls==['current','owner','draft']*3,'complete fresh guards absent')
            check((s/'ONE_DELETE_ATTEMPT_LATCH.json').is_file(),'latch missing after success')
            check(set(c.files)=={p['filename'] for p in c.old},'inherited membership changed')
        test('real_transport_one_delete_after_two_probes_and_latch_preserves_pending',lambda:case(success))
        def failed_flow(kind, code):
            def run(c,s,b):
                c.opener.delete_kind=kind
                rejected(lambda:flow(c,s,b),code)
                check(c.opener.delete_count==1 and (s/'ONE_DELETE_ATTEMPT_LATCH.json').is_file(),'failed DELETE latch/count differs')
                check(c.guard_calls==['current','owner','draft']*2,'failed DELETE advanced or retried')
                if kind.startswith('unknown'):
                    check((s/'DELETE_UNKNOWN.json').is_file(),'unknown outcome history missing')
                return {'wire_delete_calls':1,'latch_preserved':True}
            return case(run)
        test('unknown_before_effect_latched_no_retry',lambda:failed_flow('unknown_before','DELETE_OUTCOME_UNKNOWN_GET_ONLY_NO_REPEAT'))
        test('unknown_after_effect_latched_no_retry',lambda:failed_flow('unknown_after','DELETE_OUTCOME_UNKNOWN_GET_ONLY_NO_REPEAT'))
        test('delete_401_latched_no_advance',lambda:failed_flow('401','DELETE_NOT_CONFIRMED_GET_ONLY_NO_REPEAT'))
        test('delete_500_latched_no_advance',lambda:failed_flow('500','DELETE_NOT_CONFIRMED_GET_ONLY_NO_REPEAT'))
        test('delete_body_read_failure_latched_no_advance',lambda:failed_flow('read_failure','DELETE_NOT_CONFIRMED_GET_ONLY_NO_REPEAT'))
        def unknown_reconcile(c,s,b):
            c.opener.delete_kind='unknown_after'
            rejected(lambda:flow(c,s,b),'DELETE_OUTCOME_UNKNOWN_GET_ONLY_NO_REPEAT')
            before_calls=list(c.opener.calls)
            result=m.reconcile(c,h,c.initialized,c.expected)
            check(result['writes']==0 and result['repeat_delete_authorized'] is False,'reconcile scope changed')
            check(c.opener.calls==before_calls and c.opener.delete_count==1,'GET-only reconciliation repeated DELETE')
            rejected(lambda:flow(c,s,b),'EXACT_INHERITED_AND_TARGET_MEMBERSHIP_DIFFERS')
            check(c.opener.delete_count==1,'repeat DELETE after reconciliation')
        test('unknown_after_effect_get_only_reconcile_preserves_latch_and_no_repeat',lambda:case(unknown_reconcile))
        def still_present(c,s,b):
            c.opener.delete_kind='still_present'
            rejected(lambda:flow(c,s,b),'EXACT_INHERITED_AND_TARGET_MEMBERSHIP_DIFFERS')
            check(c.opener.delete_count==1 and (s/'ONE_DELETE_ATTEMPT_LATCH.json').is_file(),'unproved absence lost latch/retried')
        test('delete_204_without_absence_not_accepted',lambda:case(still_present))
        def second_guard(c,s,b):
            def change(obj):
                if len(obj.guard_calls)==4:
                    obj.files[obj.old[0]['filename']]['version_id']='manufactured-changed-version'
            c.current_hook=change
            rejected(lambda:flow(c,s,b),'INHERITED_DRAFT_IDENTITY_OR_BUCKET_DIFFERS')
            check(c.opener.delete_count==0 and not (s/'ONE_DELETE_ATTEMPT_LATCH.json').exists(),'changed inherited identity reached DELETE latch')
        test('fresh_guard_after_probes_blocks_changed_inherited_version',lambda:case(second_guard))
        def blocked_target(c,s,b,key,value):
            c.files[h.NAME][key]=value
            rejected(lambda:flow(c,s,b),'NEW_FILE_NOT_SAME_UNAVAILABLE_PENDING_ENTRY')
            check(c.opener.calls==[] and not (s/'ONE_DELETE_ATTEMPT_LATCH.json').exists(),'stored target reached probes/delete')
        test('non_null_zero_size_blocks_delete_before_probes',lambda:case(lambda c,s,b:blocked_target(c,s,b,'size',0)))
        test('non_null_checksum_blocks_delete_before_probes',lambda:case(lambda c,s,b:blocked_target(c,s,b,'checksum','')))
        def available_probe(native):
            def run(c,s,b):
                c.opener.probe_hook=lambda u:RawResponse(u,200,b'present') if (u==m.NATIVE_URL)==native else None
                rejected(lambda:flow(c,s,b),'MISSING_CONTENT_GET_NOT_DEFINITIVE_400')
                check(c.opener.delete_count==0 and not (s/'ONE_DELETE_ATTEMPT_LATCH.json').exists(),'available body reached DELETE')
            return case(run)
        test('modern_available_body_blocks_delete',lambda:available_probe(False))
        test('native_available_body_blocks_delete',lambda:available_probe(True))
        def latch_collision(c,s,b):
            raw=b'MANUFACTURED_EXISTING_LATCH\n'
            (s/'ONE_DELETE_ATTEMPT_LATCH.json').write_bytes(raw)
            rejected(lambda:flow(c,s,b))
            check(c.opener.delete_count==0 and (s/'ONE_DELETE_ATTEMPT_LATCH.json').read_bytes()==raw,'latch collision overwritten or DELETE repeated')
        test('existing_exclusive_latch_blocks_delete_and_is_not_overwritten',lambda:case(latch_collision))
        def diff_check():
            a=ast.parse(OLD.read_bytes());b=ast.parse(SOURCE.read_bytes())
            definitions=lambda tree:{n.name:ast.dump(n,include_attributes=False) for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
            ad,bd=definitions(a),definitions(b)
            check(set(bd)-set(ad)=={'native_get'} and set(ad)<=set(bd),'unexpected added/deleted definitions')
            changed=[key for key in ad if ad[key]!=bd[key]]
            check(changed==['unavailable'],'change outside native dispatch')
            top=lambda tree:[ast.dump(n,include_attributes=False) for n in tree.body if not isinstance(n,(ast.FunctionDef,ast.ClassDef,ast.Import,ast.ImportFrom))]
            check(top(a)==top(b),'source pins/entrypoint/constants changed')
            check('native_get' in bd,'fixed route absent')
            return {'added':['native_get'],'changed':['unavailable'],'remaining_definitions_unchanged':len(ad)-1}
        test('narrow_source_diff_preserves_all_other_cleanup_definitions_and_constants',diff_check)
    report={'source':str(SOURCE),'source_sha256':PIN,'archived_source_sha256':OLD_PIN,
            'test_source_sha256':sha(Path(__file__).read_bytes()),'optimized':not __debug__,
            'scope':'actual accepted aa7 Transport plus mocked opener and manufactured 27-file guards; no helper main/live network',
            'tests':rows,'passed':sum(r['status']=='PASS' for r in rows),'failed':sum(r['status']=='FAIL' for r in rows)}
    with args.output.open('x') as output:
        json.dump(report,output,indent=2,sort_keys=True);output.write('\n')
    print(json.dumps({'source_sha256':PIN,'passed':report['passed'],'failed':report['failed'],'optimized':report['optimized']}))
    return 1 if report['failed'] else 0

if __name__=='__main__':
    raise SystemExit(main())
