#!/usr/bin/env python3
"""Manufactured cleanup API controls. No helper main or live request."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import socket
import tempfile
import threading
import types
import uuid
from unittest.mock import patch


class FixtureStop(RuntimeError):
    pass


def check(ok, reason):
    if not ok:
        raise AssertionError(reason)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def deny(*args, **kwargs):
    raise AssertionError('NETWORK_FORBIDDEN_MANUFACTURED_CONTROLS')


def rejected(m, h, fn):
    try:
        fn()
    except (m.Stop, h.Stop, FixtureStop, OSError, ValueError):
        return
    raise AssertionError('invalid manufactured cleanup fixture accepted')


class Module:
    @staticmethod
    def immutable_file_identity(f):
        fields = (f.get('file_id'), f.get('version_id'))
        if not all(type(x) is str and x for x in fields):
            raise FixtureStop('invalid fixture file identity')
        return fields
    @staticmethod
    def metadata_equal(draft, expected):
        if any(draft.get(k) != value for k, value in expected.items()):
            raise FixtureStop('fixture editable metadata changed')
    @staticmethod
    def assert_pins(files, pins, exact=True, completed=True):
        if exact and set(files) != {p['filename'] for p in pins}:
            raise FixtureStop('fixture membership')
        for p in pins:
            f = files[p['filename']]
            if type(f.get('size')) is not int or f['size'] != p['bytes'] or f.get('checksum') != 'md5:'+p['md5']:
                raise FixtureStop('fixture inherited pins changed')
            if completed and f.get('status') != 'completed':
                raise FixtureStop('fixture inherited not completed')


class Response:
    def __init__(self, status, body=b'', headers=None, truncated=False, failed=False):
        self.status, self.body = status, body
        self.headers = headers or {}
        self.body_truncated, self.body_read_failed = truncated, failed
    def json(self):
        return json.loads(self.body)


class Controller:
    def __init__(self, m, h, side):
        self.m, self.h, self.side = m, h, side
        self.module = Module
        self.old = [{'filename':'manufactured-old-'+str(i), 'bytes':i+1, 'md5':hashlib.md5(str(i).encode()).hexdigest()} for i in range(27)]
        self.files = {p['filename']:{'key':p['filename'], 'file_id':str(uuid.UUID(int=i+1)), 'version_id':str(uuid.UUID(int=i+101)), 'bucket_id':h.BUCKET_ID, 'size':p['bytes'], 'checksum':'md5:'+p['md5'], 'status':'completed'} for i,p in enumerate(self.old)}
        self.files[h.NAME] = {'key':h.NAME, 'file_id':h.FILE_ID, 'version_id':h.VERSION_ID, 'bucket_id':h.BUCKET_ID, 'status':'pending', 'links':{'self':m.DELETE_URL, 'content':h.URL, 'commit':h.URL.removesuffix('/content')+'/commit'}}
        self.initialized = copy.deepcopy(self.files)
        self.expected = {'metadata':{'title':'manufactured editable metadata'}, 'custom_fields':{}, 'access':{'record':'public'}}
        self.document = copy.deepcopy(self.expected)
        self.state = {'draft_id':h.DRAFT, 'published':False, 'continuation_binding':{'controller_sha256':h.CONTROLLER_SHA}, 'inventory_sha256':h.INVENTORY_SHA, 'metadata_sha256':h.METADATA_SHA, 'create_post_attempted':True, 'publish_post_attempted':False, 'content_verified_files':[], 'pending':{'kind':'put_content','url':h.URL,'body_sha256':h.ASSET_PIN['sha256'],'body_bytes':h.ASSET_PIN['bytes'],'context':{'filename':h.NAME},'utc':'manufactured'}}
        self.owner = {'id':h.DRAFT}
        self.calls = []
        self.guard_count = 0
        self.current_hook = None
        self.request_hook = None
        self.transport = types.SimpleNamespace(request=self.request)
    def current(self):
        self.calls.append(('GET_CURRENT',))
        self.guard_count += 1
        if self.current_hook: self.current_hook(self)
    def owner_list(self):
        self.calls.append(('GET_OWNER',))
        return self.owner
    def draft(self, identifier):
        self.calls.append(('GET_DRAFT',identifier))
        check(identifier == self.h.DRAFT, 'wrong fixture draft queried')
        return self.document, {}, copy.deepcopy(self.files)
    def request(self, method, url, data=None, headers=None):
        self.calls.append((method,url))
        if self.request_hook:
            result = self.request_hook(self, method, url, data, headers)
            if result is not None: return result
        if method == 'GET':
            check(url in {self.h.URL,self.m.NATIVE_URL}, 'unexpected fixture GET')
            return Response(400,json.dumps({'status':400,'message':'File with key "'+self.h.NAME+'" is not available.'}).encode())
        check(method == 'DELETE' and url == self.m.DELETE_URL and data is None, 'unreviewed mutation')
        check((self.side/'ONE_DELETE_ATTEMPT_LATCH.json').is_file(), 'DELETE before durable latch')
        del self.files[self.h.NAME]
        return Response(204)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--primitives',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    rows=[]
    def test(name,fn):
        try: fn()
        except Exception as error: rows.append({'name':name,'status':'FAIL','exception':type(error).__name__,'reason':str(error)})
        else: rows.append({'name':name,'status':'PASS'})
    with patch.object(socket,'socket',deny),patch.object(socket,'create_connection',deny):
        m=load(args.source,'_manufactured_cleanup')
        h=load(args.primitives,'_manufactured_cleanup_primitives')
        def case(fn):
            with tempfile.TemporaryDirectory(prefix='manufactured-empty-cleanup-') as temporary:
                base=Path(temporary); side=base/'side';side.mkdir()
                evidence=base/'evidence';evidence.mkdir()
                (evidence/'STATE.json').write_bytes(b'MANUFACTURED_STATE_ONLY\n')
                (evidence/'JOURNAL.jsonl').write_bytes(b'MANUFACTURED_JOURNAL_ONLY\n')
                before={n:hashlib.sha256((evidence/n).read_bytes()).hexdigest() for n in ['STATE.json','JOURNAL.jsonl']}
                c=Controller(m,h,side)
                with patch.object(m,'EVIDENCE',evidence),patch.object(m,'REVIEWED',args.primitives):
                    return fn(c,side,before,evidence)
        def guard_case(mutator=None,present=True,accept=False):
            def run(c,side,before,evidence):
                if not present: del c.files[h.NAME]
                if mutator: mutator(c)
                call=lambda:m.guard(c,h,c.initialized,c.expected,present)
                if accept: call()
                else: rejected(m,h,call)
                check(not any(x[0]=='DELETE' for x in c.calls),'guard mutated')
            return case(run)
        test('guard_accepts_legitimate_omitted_size_checksum',lambda:guard_case(accept=True))
        test('guard_accepts_explicit_null_size_checksum',lambda:guard_case(lambda c:c.files[h.NAME].update({'size':None,'checksum':None}),accept=True))
        for key,value in [('file_id','wrong'),('version_id','wrong'),('bucket_id','wrong'),('status','completed'),('size',0),('size',1),('checksum',''),('checksum','md5:0'),('key','wrong')]:
            test('guard_reject_target_'+key+'_'+repr(value),lambda key=key,value=value:guard_case(lambda c:c.files[h.NAME].update({key:value})))
        test('guard_reject_target_self_link',lambda:guard_case(lambda c:c.files[h.NAME]['links'].update({'self':m.DELETE_URL+'?wrong=1'})))
        test('guard_reject_target_content_link',lambda:guard_case(lambda c:c.files[h.NAME]['links'].update({'content':h.URL+'?wrong=1'})))
        test('guard_reject_wrong_owner_draft',lambda:guard_case(lambda c:setattr(c,'owner',{'id':'wrong'})))
        test('guard_reject_publication_after_current',lambda:guard_case(lambda c:setattr(c,'current_hook',lambda obj:obj.state.update({'published':True}))))
        test('guard_reject_editable_metadata_change',lambda:guard_case(lambda c:c.document['metadata'].update({'title':'wrong'})))
        test('guard_reject_original_pending_change',lambda:guard_case(lambda c:c.state['pending'].update({'kind':'publish'})))
        test('guard_reject_extra_user_file',lambda:guard_case(lambda c:c.files.update({'extra-user-file':{}})))
        test('guard_reject_missing_inherited_file',lambda:guard_case(lambda c:c.files.pop(c.old[0]['filename'])))
        for key,value in [('file_id','wrong'),('version_id','wrong'),('bucket_id','wrong'),('size',0),('checksum','wrong')]:
            test('guard_reject_inherited_'+key,lambda key=key,value=value:guard_case(lambda c:c.files[c.old[0]['filename']].update({key:value})))
        test('guard_reject_changed_initialized_target',lambda:guard_case(lambda c:c.initialized[h.NAME].update({'file_id':'wrong'})))
        def collision(c):
            name=c.old[0]['filename']; c.files[name]['file_id']=h.FILE_ID;c.initialized[name]['file_id']=h.FILE_ID
        test('guard_reject_inherited_target_id_collision',lambda:guard_case(collision))
        test('guard_accepts_exact27_target_absence',lambda:guard_case(present=False,accept=True))
        def probe_case(response,url=None):
            return case(lambda c,side,before,evidence: probe_run(c,side,response,url))
        def probe_run(c,side,response,url):
            c.request_hook=lambda *unused:response
            rejected(m,h,lambda:m.unavailable(c,h,side,'DUMMY_TOKEN','PROBE',h.URL if url is None else url))
            check(not any(x[0]=='DELETE' for x in c.calls),'probe mutated')
        for label,response in [('status401',Response(401)),('present200',Response(200,b'content')),('truncated400',Response(400,b'{}',truncated=True)),('failed400',Response(400,b'',failed=True)),('wrongbody400',Response(400,b'{"status":400,"message":"wrong"}'))]:
            test('probe_reject_'+label,lambda response=response:probe_case(response))
        test('probe_reject_wrong_native_version',lambda:probe_case(Response(400,b'{}'),m.NATIVE_URL+'-wrong'))
        def flow(mutator=None,accept=False):
            def run(c,side,before,evidence):
                original=copy.deepcopy(c.state)
                if mutator:mutator(c,side,evidence)
                call=lambda:m.delete_once(c,h,side,'DUMMY_TOKEN',c.initialized,c.expected,{'manufactured_only':True,'authorization_receipt':False},before)
                if accept:
                    r=call();check(r['inherited_files']==27 and r['original_pending_cleared'] is False,'result scope')
                    check([x for x in c.calls if x[0]=='DELETE']==[('DELETE',m.DELETE_URL)],'DELETE not exact once')
                    check(c.guard_count==3,'full fresh guard stages missing')
                    check([x for x in c.calls if x[0]=='GET']==[('GET',h.URL),('GET',m.NATIVE_URL)],'missing or alternate unavailable probes')
                else:rejected(m,h,call)
                check(c.state==original,'original logical pending/state changed')
                return c,side
            return case(run)
        test('valid_one_delete_after_both_probes_and_latch_preserves_state',lambda:flow(accept=True))
        def second_probe_fail(c,side,evidence):
            c.request_hook=lambda obj,method,url,data,headers:Response(401) if url==m.NATIVE_URL else None
        test('delete_blocked_by_native_probe_failure',lambda:flow(second_probe_fail))
        def after_probe_change(c,side,evidence):
            c.current_hook=lambda obj:obj.files[h.NAME].update({'size':0}) if obj.guard_count==2 else None
        test('delete_blocked_by_final_target_change',lambda:flow(after_probe_change))
        def evidence_change(c,side,evidence):
            def hook(obj,method,url,data,headers):
                if method=='GET' and url==m.NATIVE_URL:(evidence/'STATE.json').write_bytes(b'changed manufactured state')
            c.request_hook=hook
        test('delete_blocked_by_original_evidence_change',lambda:flow(evidence_change))
        for kind in ['file','directory','symlink']:
            def existing(c,side,evidence,kind=kind):
                path=side/'ONE_DELETE_ATTEMPT_LATCH.json'
                if kind=='file':path.write_bytes(b'prior attempt')
                elif kind=='directory':path.mkdir()
                else:path.symlink_to(evidence/'STATE.json')
            test('existing_'+kind+'_latch_refuses_delete',lambda existing=existing:flow(existing))
        for status in [401,429,500,200.0]:
            def bad_delete(c,side,evidence,status=status):
                def hook(obj,method,url,data,headers):
                    if method=='DELETE':return Response(status)
                c.request_hook=hook
            test('delete_reject_status_'+repr(status)+'_keeps_latch',lambda bad_delete=bad_delete:flow(bad_delete))
        def lost_delete(c,side,evidence):
            def hook(obj,method,url,data,headers):
                if method=='DELETE':raise RuntimeError('manufactured lost response DUMMY_TOKEN')
            c.request_hook=hook
        test('lost_delete_response_no_repeat',lambda:flow(lost_delete))
        def malformed_delete(c,side,evidence):
            c.request_hook=lambda obj,method,url,data,headers:Response(204,truncated=True) if method=='DELETE' else None
        test('truncated_delete_success_is_unknown',lambda:flow(malformed_delete))
        def reconcile_case(mutator=None,accept=False):
            def run(c,side,before,evidence):
                original=copy.deepcopy(c.state)
                del c.files[h.NAME]
                if mutator:mutator(c)
                call=lambda:m.reconcile(c,h,c.initialized,c.expected)
                if accept:
                    r=call();check(r['writes']==0 and r['repeat_delete_authorized'] is False,'reconcile authorized repeat')
                else:rejected(m,h,call)
                check(c.state==original,'reconcile altered original state')
                check(not any(x[0] in {'DELETE','POST','PUT','PATCH'} for x in c.calls),'reconcile mutated')
            return case(run)
        test('GET_only_reconcile_exact_absence',lambda:reconcile_case(accept=True))
        test('reconcile_rejects_target_still_present',lambda:reconcile_case(lambda c:c.files.update({h.NAME:c.initialized[h.NAME]})))
        test('reconcile_rejects_missing_inherited',lambda:reconcile_case(lambda c:c.files.pop(c.old[0]['filename'])))
        test('reconcile_rejects_extra_user_file',lambda:reconcile_case(lambda c:c.files.update({'extra-user-file':{}})))
        test('reconcile_rejects_changed_inherited_identity',lambda:reconcile_case(lambda c:c.files[c.old[0]['filename']].update({'version_id':'wrong'})))
        test('reconcile_rejects_changed_editable_metadata',lambda:reconcile_case(lambda c:c.document['metadata'].update({'title':'wrong'})))
    result={'status':'PASS' if all(x['status']=='PASS' for x in rows) else 'FAIL','manufactured_only':True,'authorization_receipt':False,'main_invoked':False,'live_transport_calls':0,'real_state_or_prior_files_read':False,'source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),'primitives_sha256':hashlib.sha256(args.primitives.read_bytes()).hexdigest(),'optimized':not __debug__,'pass':sum(x['status']=='PASS' for x in rows),'fail':sum(x['status']=='FAIL' for x in rows),'controls':rows}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ['status','source_sha256','optimized','pass','fail']}))


if __name__=='__main__':
    main()
