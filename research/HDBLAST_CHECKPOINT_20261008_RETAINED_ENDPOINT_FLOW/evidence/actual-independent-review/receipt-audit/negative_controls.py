"""Source-free adversarial controls for the independent receipt auditor.

Loads only the externally SHA-pinned review script, never production modules.
Each mutation affects copied small exported JSON/journal artifacts. The node
stream and retained archives are never opened by this control suite.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import argparse
import json
import sys
import types

HERE=Path(__file__).resolve().parent
AUDIT=HERE/'audit_actual_receipts.py'
AUDIT_SHA='5ecad7cef60c49e73b1e67e29c20a92a22cb1cb1e6434e850681b76a556be060'
BASE=HERE.parents[1]
ACTUAL=BASE/'root/actual-normal005-001'
MODE='optimized' if sys.flags.optimize else 'normal'
OUT=HERE/('negative-controls-final-'+MODE)

def require(ok,message):
    if not ok:raise ValueError(message)

def pin(raw):return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}

def main():
    global ACTUAL,OUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--actual',type=Path,default=ACTUAL)
    parser.add_argument('--core',type=Path)
    parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args();ACTUAL=args.actual;OUT=args.output
    raw=AUDIT.read_bytes();require(pin(raw)['sha256']==AUDIT_SHA,'review auditor SHA mismatch')
    module=types.ModuleType('independent_receipt_auditor');module.__file__=str(AUDIT)
    exec(compile(raw,str(AUDIT),'exec'),module.__dict__)
    if args.core:module.CORE=args.core
    OUT.mkdir(exist_ok=False)
    source_names=['SOURCE_CERTIFICATE.json','TARGET_BUDGETS.json','SOURCE_ATTEMPTS.jsonl']
    auth_names=['INPUT_AUTHENTICATION.json','DECODE_ATTEMPTS.jsonl']
    original={n:(ACTUAL/n).read_bytes() for n in source_names+auth_names+['ENTRY_RECEIPT.json']}
    require(module.source_accounting(ACTUAL)['budget_vectors_verified']==4,'positive source baseline')
    require(module.authentication_accounting(ACTUAL)['planned_real_slots_derived_from_shapes']==688152,'positive auth baseline')
    module.verify_entry_bindings(json.loads(original['ENTRY_RECEIPT.json']))
    cases=[]
    def run(name,kind,mutate,expected_message):
        dest=OUT/name;dest.mkdir()
        names=source_names if kind=='source' else auth_names if kind=='auth' else ['ENTRY_RECEIPT.json']
        for filename in names:(dest/filename).write_bytes(original[filename])
        mutate(dest)
        error=None
        try:
            if kind=='source':module.source_accounting(dest)
            elif kind=='auth':module.authentication_accounting(dest)
            else:module.verify_entry_bindings(json.loads((dest/'ENTRY_RECEIPT.json').read_bytes()))
        except ValueError as exc:error=str(exc)
        require(error is not None and expected_message in error,'mutation not rejected for expected reason: '+name+' '+str(error))
        cases.append({'name':name,'status':'PASS_EXPECTED_REJECTION','error':error,
                      'artifacts':{n:pin((dest/n).read_bytes()) for n in names}})
    def change_json(dest,name,mutate):
        obj=json.loads((dest/name).read_bytes());mutate(obj)
        (dest/name).write_text(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n')
    run('source_rounding_omitted','source',
        lambda d:change_json(d,'SOURCE_CERTIFICATE.json',lambda x:x['rows'][0].update(source_error_rounding='0/1')),
        'complete source error rounded upward exactly once')
    def double_rounding(x):
        r=x['rows'][0];rounding=F(r['source_error_rounding'])
        require(rounding>0,'nonzero actual source rounding required')
        r['source_error_rounding']=module.canonical(2*rounding)
        r['source_error']=module.canonical(F(r['source_error'])+rounding)
    run('source_rounding_doubled','source',lambda d:change_json(d,'SOURCE_CERTIFICATE.json',double_rounding),
        'complete source error rounded upward exactly once')
    run('kernel_tail_erased','source',
        lambda d:change_json(d,'TARGET_BUDGETS.json',lambda x:x['positive_B']['-4'].update(kernel_tail_W='0')),
        'recomputed complete target budget')
    def missing_source_attempt(d):
        p=d/'SOURCE_ATTEMPTS.jsonl';lines=p.read_bytes().splitlines(keepends=True);del lines[63];p.write_bytes(b''.join(lines))
    run('source_journal_attempt_63_missing','source',missing_source_attempt,'128 exact ordered source attempts')
    run('authenticated_member_count_lowered','auth',
        lambda d:change_json(d,'INPUT_AUTHENTICATION.json',lambda x:x['capsules'][0].update(authenticated_members=18)),
        'input authentication canonical metadata bindings')
    run('inputspec_digest_mismatch','auth',
        lambda d:change_json(d,'INPUT_AUTHENTICATION.json',lambda x:x['capsules'][0].update(input_spec_sha256='0'*64)),
        'input authentication canonical metadata bindings')
    run('entry_go_digest_mismatch','entry',
        lambda d:change_json(d,'ENTRY_RECEIPT.json',lambda x:x.update(go_sha256='0'*64)),
        'entry external controls')
    run('entry_registration_digest_mismatch','entry',
        lambda d:change_json(d,'ENTRY_RECEIPT.json',lambda x:x.update(registration_sha256='0'*64)),
        'entry external controls')
    receipt={'status':'PASS_EIGHT_SOURCE_FREE_NEGATIVE_CONTROLS','optimized':bool(sys.flags.optimize),
             'auditor':pin(raw),'control_script':pin(Path(__file__).read_bytes()),'positive_baselines':3,
             'negative_cases':cases,'original_small_artifact_pins':{n:pin(raw) for n,raw in original.items()},
             'retained_archive_opens':0,'physical_source_reconstructions':0,'production_modules_imported':0,
             'node_stream_opens':0,'network_requests':0}
    target=OUT/'RESULT.json';target.write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'optimized':receipt['optimized'],'receipt':str(target),
                      'sha256':pin(target.read_bytes())['sha256']},sort_keys=True))

if __name__=='__main__':main()
