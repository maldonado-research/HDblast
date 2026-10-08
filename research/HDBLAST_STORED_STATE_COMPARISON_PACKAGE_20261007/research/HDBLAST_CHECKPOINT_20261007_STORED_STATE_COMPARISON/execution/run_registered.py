"""Offline entry: authenticate complete new scope before numeric imports."""
import argparse
from pathlib import Path
import os
import sys

sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent))
from registration_guard import (authenticate,authenticate_local,reauthenticate,
                               real_directory,require,runtime_identity,sha)
from kernel_guard import install,require_worker_limits,check_output,OUTPUT_BYTES

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--candidate',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--registration-sha256',required=True)
    p.add_argument('--public-go')
    p.add_argument('--public-go-sha256')
    p.add_argument('--repository-root')
    p.add_argument('--fabricated-only',action='store_true')
    a=p.parse_args()
    require_worker_limits()
    if a.fabricated_only:
        require(a.public_go is None and a.public_go_sha256 is None and a.repository_root is None,
                'Fabricated entry cannot accept physical authorization or input root')
        verified=authenticate_local(a.candidate,a.registration_sha256)
    else:
        require(a.public_go is not None and a.public_go_sha256 is not None and a.repository_root is not None,
                'Actual entry requires new external PUBLIC_GO and actual repository root')
        verified=authenticate(a.candidate,a.registration_sha256,a.public_go,a.public_go_sha256)
    output=real_directory(a.output)
    require(output!=verified['root'] and verified['root'] not in output.parents,
            'Output must be outside authenticated candidate')
    require(set(x.name for x in output.iterdir()) <= {'child.log'}, 'Fresh worker output required')
    # Everything reached here is authenticated source; numerical imports occur
    # only now, and the kernel restrictions precede all construction/decoding.
    from comparison_worker import import_numeric,execute,write_json
    modules=import_numeric(verified['root'])
    from flint import ctx
    ctx.threads=1
    require(ctx.threads==1,'Exactly one Flint thread required')
    denied=install()
    science=execute(verified,output,modules,a.repository_root)
    require(ctx.threads==1 and ctx.prec==512,'Fixed arithmetic precision/threads changed')
    from readback import validate_output
    review=validate_output(output,verified['root'],fabricated=a.fabricated_only,check_receipt=False)
    write_json(output/'OUTPUT_READBACK.json',review)
    reauthenticate(verified)
    names=('SOURCE_CERTIFICATE.json','DIAGNOSTIC_CROSSCHECK.json','NODE_TARGETS.jsonl.gz',
           'DATA.json','SCIENCE_SUMMARY.json','RESOURCE.json','SOURCE_ATTEMPTS.jsonl',
           'DECODE_ATTEMPTS.jsonl','OUTPUT_READBACK.json')
    files={name:{'bytes':(output/name).stat().st_size,'sha256':sha((output/name).read_bytes())} for name in names}
    write_json(output/'ENTRY_RECEIPT.json',{
        'status':'PASS_FABRICATED_STORED_COMPARISON' if a.fabricated_only else 'PASS_AUTHENTICATED_STORED_COMPARISON',
        'fabricated_only':a.fabricated_only,'registration_sha256':a.registration_sha256,
        'public_go_sha256':a.public_go_sha256,'freeze_commit':None if a.fabricated_only else verified['receipt']['freeze_commit'],
        'uid':os.getuid(),'runtime':runtime_identity(),'denied_syscalls':denied,
        'counter':science['counter'],'files':files})
    check_output(output,OUTPUT_BYTES-65536,forbidden_roots=(verified['root'],))
    print('PASS_FABRICATED_STORED_COMPARISON' if a.fabricated_only else 'PASS_AUTHENTICATED_STORED_COMPARISON',flush=True)

if __name__=='__main__':
    main()
