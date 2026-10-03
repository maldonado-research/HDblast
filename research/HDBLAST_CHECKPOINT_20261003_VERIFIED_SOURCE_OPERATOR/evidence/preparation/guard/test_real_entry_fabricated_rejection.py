"""Explicit fake readback must fail the actual child entry before core import.

This test never supplies genuine readback metadata. It uses the guard fixture
whose flags say fabricated and never calls a numerical or registered source.
It deliberately preserves the supervisor's expected negative execution log.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--guard-fixture',type=Path,required=True)
    parser.add_argument('--python',required=True)
    parser.add_argument('--output-directory',type=Path,required=True)
    parser.add_argument('--optimized',action='store_true')
    args=parser.parse_args()
    fixture=args.guard_fixture.absolute();cp=fixture/'baseline/checkpoint'
    rp=fixture/'baseline/FABRICATED_READBACK.json';out=args.output_directory.absolute()
    receipt=json.loads(rp.read_bytes())
    registration=json.loads((cp/'FULL_REGISTRATION.json').read_bytes())
    if receipt.get('fabricated_fixture') is not True or registration.get('fabricated_fixture') is not True:
        raise ValueError('This control refuses genuine authorization metadata')
    cmd=[args.python,'-B',str(cp/'execution/execute_bounded.py'),
         '--root',str(cp),'--receipt',str(rp),'--receipt-sha256',sha(rp),
         '--registration-sha256',sha(cp/'FULL_REGISTRATION.json'),
         '--python',args.python,'--route','primary','--output-dir',str(out)]
    if args.optimized:cmd.append('--optimized')
    completed=subprocess.run(cmd,text=True,capture_output=True,check=False)
    execution=json.loads((out/'EXECUTION.json').read_bytes())
    log=(out/'child.log').read_text()
    if completed.returncode!=1 or execution['status']!='EXECUTION_FAILED':
        raise RuntimeError('Fabricated metadata unexpectedly authorized actual child entry')
    if 'Fabricated readback cannot authorize a real source' not in log:
        raise RuntimeError('Actual entry failed for a different reason')
    if (out/'OUTPUT.json').exists() or (out/'SOURCE_ATTEMPTS.jsonl').exists():
        raise RuntimeError('Actual entry reached output/source authorization')
    result={'status':'PASS_EXPECTED_FABRICATED_RECEIPT_REAL_ENTRY_REJECTION',
            'fabricated_control_only':True,'optimization':1 if args.optimized else 0,
            'actual_source_callbacks':0,'numerical_method_imports':0,
            'source_attempt_journal_absent_before_authorization':True,
            'checkpoint_registration_sha256':sha(cp/'FULL_REGISTRATION.json'),
            'execution_receipt_sha256':sha(out/'EXECUTION.json'),
            'log_sha256':sha(out/'child.log'),
            'test_source_sha256':sha(Path(__file__)),
            'limitation':'Fake metadata guard control; no genuine public freeze or source evaluation.'}
    with (out/'NEGATIVE_CONTROL_REVIEW.json').open('x') as handle:
        json.dump(result,handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
