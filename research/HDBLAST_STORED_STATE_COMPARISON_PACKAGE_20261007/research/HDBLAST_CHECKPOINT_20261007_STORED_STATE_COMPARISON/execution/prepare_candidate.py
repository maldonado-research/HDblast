"""Create local, fabricated-only registration; never write a public GO."""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata as metadata
import json
import shutil
import sys
from registration_guard import FIXED_CONTRACT,ZERO_KEYS,runtime_identity

NAME='HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON'

def pin(path):
    raw=path.read_bytes()
    return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}

def json_write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,sort_keys=True,indent=2)+'\n')

def prepare(work,prior,candidate):
    candidate.mkdir(parents=True,exist_ok=True)
    for folder,names in {
        'decoder':['exact_binary80.py','fabricated_fixtures.py'],
        'theory':['entire_target_engine.py','ENTIRE_TARGET_AND_DISCRETE_COMPARISON.md'],
        'transport':['discrete_transport.py','DISCRETE_ERROR_TRANSPORT.md'],
        'provenance':['INPUT_SPEC.json','LAYOUT_PROVENANCE_BASIS.json','ORIGINAL_ZIP_INVENTORIES.json'],
    }.items():
        (candidate/folder).mkdir(exist_ok=True)
        for name in names:shutil.copyfile(work/folder/name,candidate/folder/name)
    (candidate/'source').mkdir(exist_ok=True);(candidate/'prior').mkdir(exist_ok=True)
    sources=['route.py','source_algebra_baseline.py','generic_operator_baseline.py']
    for name in sources:shutil.copyfile(prior/'independent'/name,candidate/'source'/name)
    shutil.copyfile(prior/'evidence/actual/independent-normal/BUDGET.json',candidate/'prior/BUDGET.json')
    for method in ('primary','independent'):
        shutil.copyfile(prior/('evidence/actual/'+method+'-normal/DATA.json'),candidate/('prior/'+method+'_DATA.json'))
    prior_names=['source/'+n for n in sources]+['prior/BUDGET.json','prior/primary_DATA.json','prior/independent_DATA.json']
    json_write(candidate/'SOURCE_PINS.json',{'schema_version':1,
        'origin_public_repository':'https://github.com/maldonado-research/HDblast',
        'origin_checkpoint':'research/HDBLAST_CHECKPOINT_20261007_BD_PREHISTORY',
        'reuse':'BYTE_IDENTICAL_THREE_INDEPENDENT_SOURCE_FILES_AND_PRIOR_SUCCESSFUL_EVIDENCE',
        'files':{n:pin(candidate/n) for n in prior_names}})
    environment={'schema_version':1,'python_version':'3.12.14','platform':'linux',
                 'packages':{'python-flint':'0.9.0','sympy':'1.14.0','mpmath':'1.3.0'},
                 'native_guard':'libseccomp.so.2 fail-closed syscall filter',
                 'byte_hash_scope':'OBSERVED_PREPARATION_RUNTIME_PROVENANCE_ONLY',
                 'observed_preparation_runtime':runtime_identity()}
    json_write(candidate/'ENVIRONMENT_PINS.json',environment)
    contract={**FIXED_CONTRACT,'new_evaluations_before_public_GO':{k:0 for k in sorted(ZERO_KEYS)},
              'publication_status':'DRAFT_LOCAL_FABRICATED_REGISTRATION_NO_PUBLIC_GO'}
    json_write(candidate/'REGISTRATION_CONTRACT.json',contract)
    (candidate/'execution').mkdir(exist_ok=True)
    for path in sorted((work/'execution').glob('*.py')):
        shutil.copyfile(path,candidate/'execution'/path.name)
    for path in sorted((work/'execution').glob('*.md')):
        shutil.copyfile(path,candidate/'execution'/path.name)
    files={p.relative_to(candidate).as_posix():pin(p) for p in sorted(candidate.rglob('*'))
           if p.is_file() and p.name!='FULL_REGISTRATION.json'}
    json_write(candidate/'FULL_REGISTRATION.json',{'schema_version':1,'files':files})
    return {'candidate':str(candidate),'local_registration_sha256':pin(candidate/'FULL_REGISTRATION.json')['sha256'],
            'files':len(files),'actual_source_callbacks':0,'retained_array_decodes':0,
            'status':'LOCAL_FABRICATED_CANDIDATE_ONLY_NO_PUBLIC_GO'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--work',required=True);p.add_argument('--prior',required=True)
    p.add_argument('--candidate',required=True);a=p.parse_args()
    print(json.dumps(prepare(Path(a.work),Path(a.prior),Path(a.candidate)),sort_keys=True))

if __name__=='__main__':
    main()
