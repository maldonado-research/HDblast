#!/usr/bin/env python3
"""Synthetic verifier-import guards only: no NPY/physical files are created."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def sha_bytes(b):return hashlib.sha256(b).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);args=p.parse_args()
    if args.output.exists():raise ValueError('Fresh receipt required')
    source=Path(__file__).with_name('ledger_primary.py').read_bytes()
    helper=b'raise RuntimeError("AUTHENTICATED_SYNTHETIC_HELPER_IMPORTED")\n'
    checks=[]
    with tempfile.TemporaryDirectory(prefix='ledger-bootstrap-',dir='/tmp') as temporary:
        base=Path(temporary)
        for case in ('authenticated_import','wrong_registration_sha','mutated_helper','mutated_own_source',
                     'missing_helper_pin','helper_symlink','bad_commit_format','bad_sha_format'):
            root=base/case;(root/'code').mkdir(parents=True)
            own=root/'code'/'diagnostic_primary.py';own.write_bytes(source)
            target=root/'code'/'ledger_integrity.py';target.write_bytes(helper)
            reg={'schema_version':1,'files':{'code/diagnostic_primary.py':sha_bytes(source),
                                           'code/ledger_integrity.py':sha_bytes(helper)}}
            if case=='missing_helper_pin':del reg['files']['code/ledger_integrity.py']
            registration=root/'FULL_REGISTRATION.json'
            registration.write_text(json.dumps(reg,sort_keys=True)+'\n')
            digest=sha_bytes(registration.read_bytes());commit='a'*40
            if case=='wrong_registration_sha':digest='0'*64
            elif case=='mutated_helper':target.write_bytes(helper+b'# mutation\n')
            elif case=='mutated_own_source':own.write_bytes(source+b'\n# mutation\n')
            elif case=='helper_symlink':
                outside=base/'outside-helper.py';outside.write_bytes(helper);target.unlink();target.symlink_to(outside)
            elif case=='bad_commit_format':commit='a'*39
            elif case=='bad_sha_format':digest='A'*64
            cmd=[sys.executable,str(own),'--checkpoint-root',str(root),'--registration-sha256',digest,
                 '--freeze-commit',commit,'--output-dir',str(root/'output')]
            result=subprocess.run(cmd,capture_output=True,text=True,timeout=10)
            imported='AUTHENTICATED_SYNTHETIC_HELPER_IMPORTED' in result.stderr
            if result.returncode==0 or imported != (case=='authenticated_import'):
                raise ValueError('Wrong bootstrap guard behavior '+case+' '+result.stderr)
            if (root/'output').exists():raise ValueError('Output created before verified freeze '+case)
            checks.append({'case':case,'helper_imported':imported,'exit_nonzero':True,
                           'input_arrays_opened':False,'output_created':False})
    out={'status':'PASS_SYNTHETIC_BOOTSTRAP_GUARDS','checks':checks,'checks_count':len(checks),
         'physical_files_opened':0,'source_sha256':sha_bytes(source),'verifier_sha256':sha_bytes(Path(__file__).read_bytes())}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':out['status'],'checks_count':len(checks),'output':str(args.output)}))


if __name__=='__main__':main()
