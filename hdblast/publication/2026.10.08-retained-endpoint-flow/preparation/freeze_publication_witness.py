#!/usr/bin/env python3
"""Freeze explicit terminal publication witnesses, then check them offline.

This program makes no API calls or publication writes. A separately pinned
terminal input contract and fresh historical metadata/file-list readbacks are
required. Do not point it at a running publication controller's journal.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import re
import stat
import subprocess
import sys

W3 = Path('/workspace/hdblast-research-work/continuation-trajectory-20261008')
HERE = W3/'publication-delivery-candidate'
PACKET = HERE/'hdblast/publication/2026.10.08-retained-endpoint-flow'
PROOF = PACKET/'proof'
SEEDS = PACKET/'manufactured-seeds'
VERIFIER_SHA = 'a6d7e0db5f8a68cb893c37e24fde2aa2e37d12dad96ec13ad3c9983570b9e9bf'
INVENTORY_SHA = '2c6d0a3d4c4e227c3a8aa5f6e48dd3ff8494ac4e2f1b9045c603cf35617ced3c'
METADATA_SHA = '43feec5ca6bd9ab8c5c6a2cd925f6f42ff13418bd90b334de7ab604367dc385c'
PRIOR_STREAM_SHA = '8577ff0f006a2ef5ad75d4132066c2b0fa784a569d221c47a723fd1e21081e99'
HISTORICAL_IDS = {'23228395','23225288','23114217','22347452','23111008'}

def require(ok, why):
    if not ok: raise ValueError(why)
def sha(raw): return hashlib.sha256(raw).hexdigest()
def encoded(value): return (json.dumps(value,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+'\n').encode()
def read(path):
    path = Path(os.path.abspath(path))
    require(path.is_file() and all(not p.is_symlink() for p in (path,*path.parents)) and
            stat.S_ISREG(path.stat().st_mode) and path.stat().st_nlink == 1,'REAL_SINGLE_LINK_WITNESS_REQUIRED')
    return path.read_bytes()
def pinraw(pin):
    require(type(pin) is dict and set(pin) == {'path','bytes','sha256'} and type(pin['bytes']) is int and
            0 <= pin['bytes'] <= 20*1024*1024 and re.fullmatch('[0-9a-f]{64}',pin.get('sha256','')), 'EXACT_WITNESS_SOURCE_PIN_REQUIRED')
    raw = read(pin['path']); require(len(raw)==pin['bytes'] and sha(raw)==pin['sha256'],'TERMINAL_WITNESS_SOURCE_CHANGED')
    return raw

def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--inputs',required=True); parser.add_argument('--inputs-sha256',required=True)
    args = parser.parse_args(); input_raw = read(args.inputs)
    require(re.fullmatch('[0-9a-f]{64}',args.inputs_sha256 or '') and sha(input_raw)==args.inputs_sha256,'EXTERNAL_TERMINAL_INPUT_CONTRACT_PIN_DIFFERS')
    inputs = json.loads(input_raw)
    require(inputs.get('sealed_terminal_witness_inputs') is True and inputs.get('publisher_stopped') is True and
            inputs.get('published31_complete_content_verified') is True,'TERMINAL_ROOT_CONFIRMATION_REQUIRED')
    record_id = inputs['new_record_id']; require(type(record_id) is str and re.fullmatch('[1-9][0-9]*',record_id) and record_id not in HISTORICAL_IDS,'DISTINCT_NEW_PUBLICATION_ID_REQUIRED')
    require(set(inputs['historical_readbacks'])==HISTORICAL_IDS,'EXACT_FIVE_FRESH_HISTORICAL_READBACKS_REQUIRED')
    require(sha(read(PACKET/'verify_publication.py'))==VERIFIER_SHA,'ACCEPTED_VERIFIER_SOURCE_DIFFERS')
    roles = {'inventory':'inventory/FINAL_UPLOAD_INVENTORY.json','metadata':'inventory/FINAL_METADATA_MODERN.json'}
    files = {}
    def put(name,raw):
        destination = PROOF/name
        if destination.exists(): require(read(destination)==raw,'EXISTING_PROOF_LEAF_DIFFERS:'+name)
        else: destination.parent.mkdir(parents=True,exist_ok=True); destination.write_bytes(raw)
        files[name]={'bytes':len(raw),'sha256':sha(raw)}
    for role,expected in [('inventory',INVENTORY_SHA),('metadata',METADATA_SHA)]:
        raw=read(PROOF/roles[role]);require(sha(raw)==expected,'SEALED_EDITION_CONTRACT_DIFFERS');put(roles[role],raw)
    mapping = {'new_public_record':'public/NEW_PUBLIC_RECORD.json','new_public_files':'public/NEW_PUBLIC_FILES.json',
               'latest_public_record':'public/LATEST_PUBLIC_RECORD.json','verified_public':'verification/VERIFIED_PUBLIC.json',
               'controller_journal':'controller/JOURNAL.jsonl','new_draft_record':'draft/DRAFT_RECORD.json',
               'new_draft_files':'draft/DRAFT_FILES.json','verified_draft':'verification/VERIFIED_DRAFT.json'}
    require(set(inputs['terminal_witnesses'])==set(mapping),'EXACT_TERMINAL_PUBLIC_DRAFT_JOURNAL_SOURCE_ROLES_REQUIRED')
    for role,name in mapping.items(): put(name,pinraw(inputs['terminal_witnesses'][role]));roles[role]=name
    prior_raw = read(SEEDS/'PRIOR29_CONTENT_RECEIPT.json');require(sha(prior_raw)==PRIOR_STREAM_SHA,'PRIOR29_ATTESTATION_PIN_DIFFERS')
    roles['prior_full_stream']='ancestry/PRIOR29_CONTENT_RECEIPT.json';put(roles['prior_full_stream'],prior_raw)
    groups = {}
    for rid in sorted(HISTORICAL_IDS):
        snapshots=inputs['historical_readbacks'][rid]
        require(set(snapshots)=={'record','files'},'EXACT_FRESH_HISTORICAL_METADATA_FILE_LIST_ROLES_REQUIRED')
        group = {'baseline':'preserved/'+rid+'/baseline-record.json','baseline_files':'preserved/'+rid+'/baseline-files.json',
                 'record':'preserved/'+rid+'/fresh-record.json','files':'preserved/'+rid+'/fresh-files.json'}
        put(group['baseline'],read(SEEDS/(rid+'-record.json')));put(group['baseline_files'],read(SEEDS/(rid+'-files.json')))
        put(group['record'],pinraw(snapshots['record']));put(group['files'],pinraw(snapshots['files']));groups[rid]=group
    manifest = {'schema_version':1,'files':files,'witnesses':roles,'preserved_records':groups,'historical_original_literal_Notes_guard':'FAIL_PRESERVED'}
    raw=encoded(manifest);manifest_path=PROOF/'WITNESS_MANIFEST.json';require(not manifest_path.exists(),'FRESH_WITNESS_MANIFEST_REQUIRED');manifest_path.write_bytes(raw);manifest_sha=sha(raw)
    # Captured source is authenticated before every run. All arguments are
    # reconstructed from fixed trust roots and the sealed input record ID.
    bootstrap = '''from pathlib import Path
import hashlib,sys
p=Path(sys.argv[1]);raw=p.read_bytes()
if hashlib.sha256(raw).hexdigest()!=sys.argv[2]: raise SystemExit('Accepted verifier source differs')
sys.argv=[str(p),*sys.argv[3:]]
exec(compile(raw,str(p),'exec'),{'__name__':'__main__','__file__':str(p)})
'''
    checks=PACKET/'current-verification';checks.mkdir(exist_ok=True)
    receipts={}
    for mode in ('normal','optimized'):
        output=checks/('OFFLINE_PUBLICATION_'+mode.upper()+'.json');require(not output.exists(),'FRESH_OFFLINE_CHECK_OUTPUT_REQUIRED')
        command=[sys.executable,'-I','-S','-B']+(['-O'] if mode=='optimized' else [])+['-c',bootstrap,str(PACKET/'verify_publication.py'),VERIFIER_SHA,
                 '--root',str(PROOF),'--witness-manifest',str(manifest_path),'--witness-manifest-sha256',manifest_sha,
                 '--expected-record-id',record_id,'--inventory-sha256',INVENTORY_SHA,'--metadata-sha256',METADATA_SHA,'--output',str(output)]
        run=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        (checks/(mode+'.stdout')).write_bytes(run.stdout);(checks/(mode+'.stderr')).write_bytes(run.stderr)
        require(run.returncode==0,'OFFLINE_PUBLICATION_VERIFIER_FAILURE_PRESERVED:'+mode)
        receipts[mode]=read(output)
    require(receipts['normal']==receipts['optimized'],'NORMAL_OPTIMIZED_PUBLICATION_RECEIPTS_DIFFER')
    result=json.loads(receipts['normal']); require(result.get('status')=='PASS_OFFLINE_31_FILE_PUBLICATION_WITNESSES' and
            result.get('file_count')==31 and result.get('total_bytes')==453646943,'COMPLETE_REAL_SAVED_PUBLICATION_WITNESS_REQUIRED')
    receipt={'schema_version':1,'status':'PASS_FROZEN_SAVED_PUBLICATION_WITNESS_NORMAL_OPTIMIZED_INTERNAL_PARENT_REVIEW_PENDING',
             'frozen_utc':datetime.now(timezone.utc).isoformat(),'new_record_id':record_id,
             'source_input_contract_sha256':sha(input_raw),'witness_manifest_sha256':manifest_sha,'verifier_sha256':VERIFIER_SHA,
             'inventory_sha256':INVENTORY_SHA,'metadata_sha256':METADATA_SHA,'witness_file_count':len(files),
             'normal_optimized_receipts_equal':True,'offline_receipt_sha256':sha(receipts['normal']),
             'network_calls_by_this_freezer':0,'new_attachment_streams_by_this_freezer':0,'remote_writes':0,
             'numeric_imports':0,'source_callbacks':0,'array_decodes':0,'facts_and_readme_created':False}
    (HERE/'FROZEN_PUBLICATION_WITNESS_RECEIPT.json').write_bytes(encoded(receipt))
    print(json.dumps({k:receipt[k] for k in ('status','new_record_id','witness_manifest_sha256','witness_file_count')}))

if __name__=='__main__':main()
