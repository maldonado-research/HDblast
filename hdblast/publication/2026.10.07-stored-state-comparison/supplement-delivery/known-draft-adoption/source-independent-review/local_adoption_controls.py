#!/usr/bin/env python3
"""Pinned-source local adoption review; only temporary manufactured writes."""
import argparse
import ast
import copy
import hashlib
import importlib.util
import json
import socket
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

sys.dont_write_bytecode = True
BASE = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007/publication-final/known-draft-adoption')
SOURCE = BASE / 'adopt_known_draft_once.py'
PIN = 'f3cd239eba834806961a6ed246ab183ee63c149221f1fca2a9dba7a2cb5db4da'

def check(ok, reason):
    if not ok:
        raise AssertionError(reason)

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def deny(*args, **kwargs):
    raise AssertionError('NETWORK_FORBIDDEN')

def reject(fn):
    try:
        fn()
    except Exception as exc:
        check(not isinstance(exc, AssertionError), 'fixture assertion: '+str(exc))
        return
    raise AssertionError('invalid fixture accepted')

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    raw=SOURCE.read_bytes();check(digest(raw)==PIN,'source pin changed')
    m=importlib.util.module_from_spec(importlib.util.spec_from_loader('_local_adoption_review',loader=None,origin=str(SOURCE)))
    m.__file__=str(SOURCE);sys.modules[m.__name__]=m
    exec(compile(raw,str(SOURCE),'exec'),m.__dict__)
    manifest=BASE/'INPUT_PINS.json';manifest_raw=manifest.read_bytes()
    check(digest(manifest_raw)==m.PINS_SHA256,'input manifest changed')
    pinned=json.loads(manifest_raw)
    before_state=(BASE/'live-local-adoption-001/NEW_STATE_BEFORE.json').read_bytes()
    prefix=Path(pinned['new_journal']['path']).read_bytes()[:pinned['new_journal']['bytes']]
    check(digest(prefix)==pinned['new_journal']['sha256'],'initial journal prefix changed')
    originals={role:Path(pinned[role]['path']).read_bytes() for role in ('old_state','old_journal')}
    for role in originals:
        check(digest(originals[role])==pinned[role]['sha256'],'original protected evidence changed')
    rows=[]
    def test(name,fn):
        try:fn()
        except Exception as exc:rows.append({'name':name,'status':'FAIL','exception':type(exc).__name__,'reason':str(exc)})
        else:rows.append({'name':name,'status':'PASS'})
    def case(fn):
        with tempfile.TemporaryDirectory(prefix='offline-local-adoption-') as tmp:
            root=Path(tmp);old=root/'old';new=root/'new';inputs=root/'inputs'
            for directory in (old,new,inputs):directory.mkdir()
            (old/'CONTROLLER.lock').write_bytes(b'');(new/'CONTROLLER.lock').write_bytes(b'')
            spec={}
            for role,entry in pinned.items():
                directory=old if role.startswith('old_') else new if role.startswith('new_') else inputs
                path=directory/role
                value=before_state if role=='new_state' else prefix if role=='new_journal' else Path(entry['path']).read_bytes()
                if role not in ('new_state','new_journal'):
                    check(digest(value)==entry['sha256'],'pinned input changed: '+role)
                path.write_bytes(value)
                spec[role]={'path':str(path),'bytes':len(value),'sha256':digest(value)}
            protected={role:Path(spec[role]['path']).read_bytes() for role in ('old_state','old_journal')}
            recovery=root/'recovery'
            result=fn(spec,recovery)
            check(all(Path(spec[role]['path']).read_bytes()==value for role,value in protected.items()),'original evidence modified')
            return result
    with patch.object(socket,'socket',deny),patch.object(socket,'create_connection',deny):
        def readonly(spec,recovery):
            initial={role:Path(entry['path']).read_bytes() for role,entry in spec.items()}
            result=m.adopt(spec,recovery,apply=False)
            check(result['writes']==0 and not recovery.exists(),'readonly adoption wrote recovery')
            check(all(Path(spec[role]['path']).read_bytes()==value for role,value in initial.items()),'readonly changed input')
        test('local_readonly_inputs_only_with_pinned_controller_bytes_never_imported',lambda:case(readonly))
        def success(spec,recovery):
            state_path=Path(spec['new_state']['path']);journal=Path(spec['new_journal']['path'])
            initial=json.loads(state_path.read_bytes());initial_journal=journal.read_bytes()
            original_append=m.append_event
            def guarded_append(path,value):
                check((recovery/'ONE_LOCAL_ADOPTION_LATCH.json').is_file(),'journal began before latch')
                return original_append(path,value)
            with patch.object(m,'append_event',guarded_append):result=m.adopt(spec,recovery,apply=True)
            expected=copy.deepcopy(initial);expected['draft_id']='23228395'
            check(json.loads(state_path.read_bytes())==expected,'state changed outside draft_id')
            check(result['only_state_field_changed']=='draft_id' and result['network_writes']==0,'result scope')
            check(journal.read_bytes().startswith(initial_journal),'journal prefix rewritten')
            suffix=journal.read_bytes()[len(initial_journal):].splitlines()
            check(len(suffix)==2 and all(type(json.loads(line)) is dict for line in suffix),'adoption event count')
            check((recovery/'RESULT.json').is_file(),'success receipt missing')
            final_state=state_path.read_bytes();final_journal=journal.read_bytes()
            reject(lambda:m.adopt(spec,recovery,apply=True))
            check(state_path.read_bytes()==final_state and journal.read_bytes()==final_journal,'repeat changed adopted state/journal')
        test('one_local_adoption_only_draft_id_and_append_only_journal_then_repeat_blocked',lambda:case(success))
        def journal_failure(spec,recovery):
            initial=Path(spec['new_state']['path']).read_bytes();journal=Path(spec['new_journal']['path']);prefix=journal.read_bytes()
            def failed(*args,**kwargs):raise OSError('manufactured journal failure')
            with patch.object(m,'append_event',failed):reject(lambda:m.adopt(spec,recovery,apply=True))
            check((recovery/'ONE_LOCAL_ADOPTION_LATCH.json').is_file(),'journal failure lost latch')
            check(Path(spec['new_state']['path']).read_bytes()==initial and journal.read_bytes()==prefix,'journal failure changed state/prefix')
            check(not (recovery/'RESULT.json').exists(),'failure reported success')
            reject(lambda:m.adopt(spec,recovery,apply=True))
        test('journal_failure_before_state_preserves_latch_state_prefix_and_blocks_repeat',lambda:case(journal_failure))
        def replace_failure(spec,recovery):
            state=Path(spec['new_state']['path']);journal=Path(spec['new_journal']['path']);initial=state.read_bytes();prefix=journal.read_bytes()
            def failed(*args,**kwargs):raise OSError('manufactured atomic replacement failure')
            with patch.object(m.os,'replace',failed):reject(lambda:m.adopt(spec,recovery,apply=True))
            check(state.read_bytes()==initial and journal.read_bytes().startswith(prefix),'replace failure changed state/prefix')
            check((recovery/'ONE_LOCAL_ADOPTION_LATCH.json').is_file() and not (recovery/'RESULT.json').exists(),'replace failure lost latch/reported success')
            reject(lambda:m.adopt(spec,recovery,apply=True))
        test('atomic_replace_failure_keeps_prior_state_and_latch_no_automatic_repeat',lambda:case(replace_failure))
        def success_receipt_failure(spec,recovery):
            original=m.exclusive_json
            def fail_result(path,value):
                if Path(path).name=='RESULT.json':raise OSError('manufactured final receipt failure')
                return original(path,value)
            with patch.object(m,'exclusive_json',fail_result):reject(lambda:m.adopt(spec,recovery,apply=True))
            check(json.loads(Path(spec['new_state']['path']).read_bytes())['draft_id']=='23228395','adoption completion stage not exercised')
            check((recovery/'ONE_LOCAL_ADOPTION_LATCH.json').is_file() and not (recovery/'RESULT.json').exists(),'incomplete result lost latch')
            reject(lambda:m.adopt(spec,recovery,apply=True))
        test('failure_after_state_and_journal_completion_still_latched_and_not_repeated',lambda:case(success_receipt_failure))
        def static():
            tree=ast.parse(raw)
            imports=[]
            for node in ast.walk(tree):
                if isinstance(node,ast.Import):imports.extend(x.name for x in node.names)
                elif isinstance(node,ast.ImportFrom):imports.append(node.module or '')
            check(not any(x.split('.')[0] in {'socket','ssl','urllib','http','requests','subprocess','importlib'} or 'controller' in x for x in imports),'network/controller import')
            check('exec(' not in raw.decode() and 'eval(' not in raw.decode(),'dynamic controller code execution')
        test('source_has_no_network_controller_import_or_dynamic_execution',static)
    report={'source_sha256':PIN,'manifest_sha256':digest(manifest_raw),'test_source_sha256':digest(Path(__file__).read_bytes()),
            'optimized':not __debug__,'passed':sum(x['status']=='PASS' for x in rows),'failed':sum(x['status']=='FAIL' for x in rows),
            'scope':'No live adoption/main; original inputs copied into temporary fixtures; saved pre-adoption state and exact journal prefix used because root already adopted.',
            'tests':rows}
    with args.output.open('x') as out:json.dump(report,out,indent=2,sort_keys=True);out.write('\n')
    print(json.dumps({k:report[k] for k in ('source_sha256','passed','failed','optimized')}))
    return 1 if report['failed'] else 0

if __name__=='__main__':raise SystemExit(main())
