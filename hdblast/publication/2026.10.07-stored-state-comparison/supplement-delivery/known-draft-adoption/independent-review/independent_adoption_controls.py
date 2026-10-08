#!/usr/bin/env python3
"""Seventeen bounded local controls; manufactured state/evidence, no main."""
import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import socket
import tempfile
from unittest.mock import patch


def check(ok,reason):
    if not ok:raise AssertionError(reason)


def encode(value):return (json.dumps(value,sort_keys=True)+'\n').encode()


def deny(*args,**kwargs):raise AssertionError('NETWORK_FORBIDDEN_LOCAL_CONTROLS')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--controller-bytes',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    s=importlib.util.spec_from_file_location('_manufactured_local_adoption',args.source)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    controller_raw=args.controller_bytes.read_bytes()
    check(hashlib.sha256(controller_raw).hexdigest()==m.CONTROLLER_SHA256,'accepted controller byte pin differs')
    rows=[]
    def test(name,fn):
        try:fn()
        except Exception as error:rows.append({'name':name,'status':'FAIL','exception':type(error).__name__,'reason':str(error)})
        else:rows.append({'name':name,'status':'PASS'})
    def refused(fn):
        try:fn()
        except (m.Refused,OSError,ValueError):return
        raise AssertionError('invalid local fixture accepted')
    def case(fn):
        with tempfile.TemporaryDirectory(prefix='manufactured-known-draft-') as temporary:
            base=Path(temporary);old_dir=base/'old';new_dir=base/'new';inputs=base/'inputs'
            for folder in [old_dir,new_dir,inputs]:folder.mkdir()
            for folder in [old_dir,new_dir]:(folder/'CONTROLLER.lock').write_bytes(b'')
            binding={'controller_sha256':m.CONTROLLER_SHA256,'prior_record_id':'23225288','parent_id':'17088132'}
            inventory=b'MANUFACTURED_INVENTORY_NOT_ACTUAL_UPLOAD_INPUT\n';metadata=b'MANUFACTURED_METADATA_NOT_ACTUAL_PUBLIC_METADATA\n'
            old={'draft_id':m.DRAFT_ID,'create_post_attempted':True,'published':False,'pending':{'kind':'put_content','marker':'MANUFACTURED_NOT_ACTUAL_AUTHORIZATION'},'continuation_binding':binding}
            new={'pending':None,'draft_id':None,'published':False,'continuation_binding':copy.deepcopy(binding),'inventory_sha256':hashlib.sha256(inventory).hexdigest(),'metadata_sha256':hashlib.sha256(metadata).hexdigest()}
            raw={role:('MANUFACTURED_NONEXECUTABLE_'+role+'\n').encode() for role in m.EXPECTED_ROLES}
            raw.update({'old_state':encode(old),'old_journal':b'{"kind":"manufactured_original"}\n','new_state':encode(new),'new_journal':b'{"kind":"manufactured_new"}\n','inventory':inventory,'metadata':metadata,'controller':controller_raw,'creation':encode({'id':int(m.DRAFT_ID),'conceptrecid':'17088132','owner':1386319,'state':'unsubmitted','submitted':False}),'cleanup_RESULT.json':encode({'original_pending_cleared':False,'repeat_delete_authorized':False,'status':'PASS_GET_ONLY_EXACT_TARGET_ABSENT_27_INHERITED_PRESERVED','writes':0})})
            spec={}
            for role,value in raw.items():
                folder=old_dir if role in {'old_state','old_journal'} else new_dir if role in {'new_state','new_journal'} else inputs
                path=folder/role;path.write_bytes(value)
                spec[role]={'path':str(path),'bytes':len(value),'sha256':hashlib.sha256(value).hexdigest()}
            held_old={role:Path(spec[role]['path']).read_bytes() for role in ['old_state','old_journal']}
            recovery=base/'recovery'
            result=fn(spec,recovery,raw,new)
            check(all(Path(spec[role]['path']).read_bytes()==value for role,value in held_old.items()),'protected original bytes changed')
            return result
    def repin(spec,role,value):
        path=Path(spec[role]['path']);path.write_bytes(value)
        spec[role].update({'bytes':len(value),'sha256':hashlib.sha256(value).hexdigest()})
    def mutation_case(role,path,value):
        def run(spec,recovery,raw,new):
            obj=json.loads(raw[role]);ref=obj
            for key in path[:-1]:ref=ref[key]
            ref[path[-1]]=value
            repin(spec,role,encode(obj))
            refused(lambda:m.adopt(spec,recovery,apply=True))
            check(not recovery.exists(),'invalid inputs reserved adoption')
        return lambda:case(run)
    with patch.object(socket,'socket',deny),patch.object(socket,'create_connection',deny):
        def read_only(spec,recovery,raw,new):
            result=m.adopt(spec,recovery)
            check(result['writes']==0 and not recovery.exists(),'read-only path wrote recovery')
            check(all(Path(spec[role]['path']).read_bytes()==value for role,value in raw.items()),'read-only input changed')
        test('default_read_only_exact19_inputs_unchanged',lambda:case(read_only))
        def valid(spec,recovery,raw,new):
            result=m.adopt(spec,recovery,apply=True)
            actual=json.loads(Path(spec['new_state']['path']).read_bytes())
            expected=copy.deepcopy(new);expected['draft_id']=m.DRAFT_ID
            check(actual==expected and result['only_state_field_changed']=='draft_id','state changed outside draft_id')
            journal=Path(spec['new_journal']['path']).read_bytes()
            check(journal.startswith(raw['new_journal']),'journal prefix rewritten')
            events=[json.loads(row) for row in journal[len(raw['new_journal']):].splitlines()]
            check([x['kind'] for x in events]==['LOCAL_KNOWN_DRAFT_ADOPTION_INTENT','LOCAL_KNOWN_DRAFT_ADOPTED'],'auditable events differ')
            check((recovery/'ONE_LOCAL_ADOPTION_LATCH.json').is_file(),'durable latch absent')
            check(result['network_writes']==0 and result['original_pending_cleared'] is False,'scope changed')
        test('valid_sole_field_journal_prefix_and_protected_old_preserved',lambda:case(valid))
        def repeated(spec,recovery,raw,new):
            m.adopt(spec,recovery,apply=True)
            held=Path(spec['new_journal']['path']).read_bytes()
            refused(lambda:m.adopt(spec,recovery,apply=True))
            check(Path(spec['new_journal']['path']).read_bytes()==held,'repeat appended journal')
        test('persistent_recovery_refuses_repeat',lambda:case(repeated))
        def role_set(spec,recovery,raw,new):
            del spec['creation'];refused(lambda:m.adopt(spec,recovery,apply=True))
        test('reject_missing_input_role',lambda:case(role_set))
        def pin_diff(spec,recovery,raw,new):
            Path(spec['cleanup_source']['path']).write_bytes(b'corrupted source')
            refused(lambda:m.adopt(spec,recovery,apply=True));check(not recovery.exists(),'pin mismatch reserved recovery')
        test('reject_source_pin_mismatch_before_intent',lambda:case(pin_diff))
        test('reject_new_state_already_has_draft',mutation_case('new_state',['draft_id'],'another-draft'))
        test('reject_nonnull_new_pending',mutation_case('new_state',['pending'],{'kind':'put_content'}))
        test('reject_published_new_state',mutation_case('new_state',['published'],True))
        test('reject_continuation_binding_mismatch',mutation_case('new_state',['continuation_binding','parent_id'],'wrong'))
        test('reject_inventory_binding_mismatch',mutation_case('new_state',['inventory_sha256'],'wrong'))
        test('reject_wrong_genuine_created_draft',mutation_case('creation',['id'],1))
        test('reject_cleanup_nonabsence_proof',mutation_case('cleanup_RESULT.json',['repeat_delete_authorized'],True))
        def journal_boundary(spec,recovery,raw,new):
            repin(spec,'new_journal',raw['new_journal'].rstrip(b'\n'));refused(lambda:m.adopt(spec,recovery,apply=True))
        test('reject_incomplete_new_journal_boundary',lambda:case(journal_boundary))
        def input_alias(spec,recovery,raw,new):
            spec['creation']=copy.deepcopy(spec['cleanup_source']);refused(lambda:m.adopt(spec,recovery,apply=True))
        test('reject_aliased_input_leaf',lambda:case(input_alias))
        def symlink(spec,recovery,raw,new):
            path=Path(spec['new_state']['path']);target=path.parent/'actual_state_copy';path.rename(target);path.symlink_to(target)
            refused(lambda:m.adopt(spec,recovery,apply=True));check(not recovery.exists(),'symlink reserved recovery')
        test('reject_symlink_state',lambda:case(symlink))
        def recovery_alias(spec,recovery,raw,new):
            refused(lambda:m.adopt(spec,Path(spec['new_state']['path']).parent/'recovery',apply=True))
        test('reject_recovery_inside_execution',lambda:case(recovery_alias))
        def rename_failure(spec,recovery,raw,new):
            with patch.object(m.os,'replace',side_effect=OSError('manufactured rename failure')):
                refused(lambda:m.adopt(spec,recovery,apply=True))
            check((recovery/'ONE_LOCAL_ADOPTION_LATCH.json').is_file(),'interrupted operation lost latch')
            check(Path(spec['new_state']['path']).read_bytes()==raw['new_state'],'failed rename changed state')
            journal=Path(spec['new_journal']['path']).read_bytes();check(journal.startswith(raw['new_journal']),'interrupted append rewrote prefix')
            refused(lambda:m.adopt(spec,recovery,apply=True))
            check(Path(spec['new_journal']['path']).read_bytes()==journal,'interrupted intent automatically repeated')
        test('atomic_state_failure_retains_latch_and_refuses_repeat',lambda:case(rename_failure))
    result={'status':'PASS' if all(x['status']=='PASS' for x in rows) else 'FAIL','source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),'controller_bytes_sha256':hashlib.sha256(controller_raw).hexdigest(),'optimized':not __debug__,'pass':sum(x['status']=='PASS' for x in rows),'fail':sum(x['status']=='FAIL' for x in rows),'manufactured_state_and_proofs':True,'authorization_receipt':False,'main_invoked':False,'live_network_calls':0,'actual_state_journal_read_or_mutated':False,'controller_imported':False,'controls':rows}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ['status','source_sha256','optimized','pass','fail']}))


if __name__=='__main__':main()
