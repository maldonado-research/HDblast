"""Exercise actual staged stdlib guard on explicitly fabricated readback metadata.

No remote requests occur, no numerical method module is imported, and no
registered source callback or scientific array decoder is invoked. Fake public
receipt fields exist solely to test the offline custodian guard's contract.
"""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys

BASE = Path('/workspace/hdblast-research-work/verified-integration-20261003')
ROOT = Path(__file__).resolve().parent

def sha(data):
    return hashlib.sha256(data).hexdigest()

def dump(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')

def module(path):
    spec = importlib.util.spec_from_file_location('fabricated_guard_only', path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

def source_map():
    independent = next((BASE/'independent-design'/name for name in
        ('independent_route.py', 'route.py')
        if (BASE/'independent-design'/name).is_file()), None)
    if independent is None:
        raise ValueError('Actual staged independent route is not yet present')
    return {
        'execution/registration_guard.py':BASE/'execution-design/registration_guard.py',
        'execution/run_registered.py':BASE/'execution-design/run_registered.py',
        'execution/execute_bounded.py':BASE/'execution-design/execute_bounded.py',
        'execution/register_and_verify.py':BASE/'execution-design/register_and_verify.py',
        'execution/review_completed.py':BASE/'execution-design/review_completed.py',
        'execution/replay_checkpoint.py':BASE/'execution-design/replay_checkpoint.py',
        'primary/verified_moments.py':BASE/'primary-design/verified_moments.py',
        'primary/route.py':BASE/'primary-design/primary_route.py',
        'independent/generic_operator.py':BASE/'independent-design/generic_operator.py',
        'independent/source_models.py':BASE/'independent-design/source_models.py',
        'independent/route.py':independent,
        'independent/fabricated_complete_route.py':BASE/'independent-design/fabricated_complete_route.py',
        'protocol/operator_probe_contract.py':BASE/'protocol/operator_probe_contract.py',
    }

def seal(checkpoint, receipt_path):
    """Construct fake but internally bound metadata, never an execution GO."""
    contract=json.loads((checkpoint/'REGISTRATION_CONTRACT.json').read_text())
    files={}
    for path in sorted(checkpoint.rglob('*')):
        if path.is_file() and path.name!='FULL_REGISTRATION.json':
            data=path.read_bytes()
            files[path.relative_to(checkpoint).as_posix()]={'sha256':sha(data),'bytes':len(data)}
    reg={'schema_version':1,'scientific_outcome_before_freeze':'UNCOMPUTED_SOURCE_OPERATOR',
         'new_target_evaluations_before_freeze':{'primary':0,'independent':0},
         'files':files,'frozen_configuration':contract,
         'fabricated_fixture':True,'not_a_real_remote_freeze':True}
    dump(checkpoint/'FULL_REGISTRATION.json',reg)
    return write_receipt(checkpoint,receipt_path)

def write_receipt(checkpoint, receipt_path):
    regraw=(checkpoint/'FULL_REGISTRATION.json').read_bytes()
    reg=json.loads(regraw)
    remote={}
    for name in sorted(set(reg['files'])|{'FULL_REGISTRATION.json'}):
        path=checkpoint/name
        if not path.is_file():
            raise ValueError('Receipt build requires all local copied fixture bytes')
        raw=path.read_bytes()
        remote[name]={'sha256':sha(raw),'bytes':len(raw),
          'independent_public_blob_download_sha256_verified':True,
          'git_blob':hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()}
    receipt={'status':'PASS_REMOTE_REGISTERED_SOURCE_GO',
        'public_repository':'https://github.com/maldonado-research/HDblast',
        'input_scope':'NO_RETAINED_ARRAYS_SOURCE_OPERATOR_ONLY',
        'all_registered_public_blobs_independently_downloaded':True,
        'new_target_evaluations_before_public_verification':{'primary':0,'independent':0},
        'freeze_commit':'f'*40,'registration_sha256':sha(regraw),
        'remote_registered_files':remote,'remote_files_verified':len(remote),
        'fabricated_fixture':True,'not_a_real_remote_freeze':True,
        'actual_public_blob_downloads':0,'actual_source_callbacks':0,
        'warning':'OFFLINE FABRICATED CUSTODIAN GUARD FIXTURE. NEVER USE TO AUTHORIZE A REAL SOURCE.'}
    dump(receipt_path,receipt)
    return sha(receipt_path.read_bytes()),sha(regraw)

def build_fixture(parent):
    checkpoint=parent/'checkpoint';checkpoint.mkdir(parents=True)
    original_pins={}
    for relative,origin in source_map().items():
        raw=origin.read_bytes();original_pins[str(origin)]={'sha256':sha(raw),'bytes':len(raw)}
        target=checkpoint/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
    contract=json.loads((BASE/'protocol/EXPERIMENT.json').read_text())
    contract.update(freeze_authorization='NOT_AUTHORIZED_BEFORE_REMOTE_BYTE_VERIFICATION',include_Lg=True,
        fabricated_fixture=True,not_a_real_remote_freeze=True)
    contract['resources'].update(each_route_wall_seconds=120,each_route_peak_rss_kib=131072,
        status='FABRICATED_GUARD_FIXTURE_LIMITS_NOT_PRODUCTION_EVIDENCE')
    dump(checkpoint/'REGISTRATION_CONTRACT.json',contract)
    (checkpoint/'PROTOCOL.md').write_bytes((BASE/'protocol/PROTOCOL.md').read_bytes())
    reviewed={p.relative_to(checkpoint).as_posix():sha(p.read_bytes())
              for p in checkpoint.rglob('*') if p.is_file()}
    review=checkpoint/'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json'
    review.parent.mkdir(parents=True)
    dump(review,{'status':'GO_FOR_PROSPECTIVE_PUBLIC_FREEZE','reviewed_source_sha256':reviewed,
        'fabricated_fixture':True,'not_a_real_remote_freeze':True,
        'warning':'Dummy GO metadata exercises the offline guard; not an actual review authorization.'})
    receipt=parent/'FABRICATED_READBACK.json';pins=seal(checkpoint,receipt)
    return checkpoint,receipt,pins,original_pins

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-directory',type=Path,required=True)
    args=parser.parse_args();out=args.output_directory.resolve()
    if out.exists():raise ValueError('Fresh external fixture directory required')
    out.mkdir(parents=True)
    guard=module(BASE/'execution-design/registration_guard.py')
    baseline,receipt,pins,original_pins=build_fixture(out/'baseline')
    checks=[];mutations=[]
    verified=guard.authenticate(baseline,receipt,*pins)
    checks.append('actual staged guard accepts internally bound fabricated baseline')
    def rejected(label,change):
        case=out/('mutation_'+str(len(mutations)+1));shutil.copytree(out/'baseline',case)
        cp=case/'checkpoint';rp=case/'FABRICATED_READBACK.json';rs,gs=pins
        result=change(cp,rp,rs,gs)
        if result is not None:rs,gs=result
        try:guard.authenticate(cp,rp,rs,gs)
        except (ValueError,KeyError,TypeError,FileNotFoundError) as err:
            mutations.append({'name':label,'rejected':True,'diagnostic':str(err)})
        else:raise RuntimeError('Guard mutation survived: '+label)
    def receipt_edit(rp,edit):
        r=json.loads(rp.read_text());edit(r);dump(rp,r);return sha(rp.read_bytes())
    def reg_edit(cp,rp,edit):
        p=cp/'FULL_REGISTRATION.json';r=json.loads(p.read_text());edit(r);dump(p,r)
        return write_receipt(cp,rp)
    rejected('tampered core source',lambda cp,rp,rs,gs:(cp/'primary/verified_moments.py').write_text((cp/'primary/verified_moments.py').read_text()+'\n# mutation\n') and None)
    rejected('missing required member',lambda cp,rp,rs,gs:(cp/'independent/source_models.py').unlink())
    rejected('receipt boolean chronology',lambda cp,rp,rs,gs:(receipt_edit(rp,lambda r:r.update(new_target_evaluations_before_public_verification={'primary':False,'independent':0})),gs))
    rejected('receipt string chronology',lambda cp,rp,rs,gs:(receipt_edit(rp,lambda r:r.update(new_target_evaluations_before_public_verification={'primary':'0','independent':0})),gs))
    rejected('registration boolean chronology',lambda cp,rp,rs,gs:reg_edit(cp,rp,lambda r:r.update(new_target_evaluations_before_freeze={'primary':0,'independent':False})))
    rejected('registration renamed chronology',lambda cp,rp,rs,gs:reg_edit(cp,rp,lambda r:r.__setitem__('new_target_evaluations_before_freeze',{'producer':0,'independent':0})))
    rejected('incorrect receipt hash',lambda cp,rp,rs,gs:('0'*64,gs))
    rejected('incorrect registration hash',lambda cp,rp,rs,gs:(rs,'0'*64))
    rejected('remote inventory omission',lambda cp,rp,rs,gs:(receipt_edit(rp,lambda r:r['remote_registered_files'].pop('primary/route.py')),gs))
    rejected('remote Git blob mismatch',lambda cp,rp,rs,gs:(receipt_edit(rp,lambda r:r['remote_registered_files']['primary/route.py'].update(git_blob='0'*40)),gs))
    rejected('unregistered Python source',lambda cp,rp,rs,gs:(cp/'primary/unregistered.py').write_text('raise RuntimeError("MUST NEVER IMPORT")\n') and None)
    rejected('injected cached bytecode',lambda cp,rp,rs,gs:(cp/'primary/injected.pyc').write_bytes(b'not executed') and None)
    rejected('injected native module',lambda cp,rp,rs,gs:(cp/'primary/injected.so').write_bytes(b'not executed') and None)
    rejected('checkpoint symlink',lambda cp,rp,rs,gs:(cp/'primary/unregistered-link').symlink_to(cp/'primary/route.py'))
    def empty_review(cp,rp,rs,gs):
        path=cp/'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json';r=json.loads(path.read_text())
        r['reviewed_source_sha256']={};dump(path,r);return seal(cp,rp)
    rejected('empty reviewed-source coverage',empty_review)
    def changed_review(cp,rp,rs,gs):
        path=cp/'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json';r=json.loads(path.read_text())
        r['reviewed_source_sha256']['primary/route.py']='0'*64;dump(path,r);return seal(cp,rp)
    rejected('reviewed core hash mismatch',changed_review)
    for raw,label in [(' {"same":1,"same":2}', 'duplicate JSON key'),
                      ('{"value":NaN}', 'NaN literal'),('{"value":Infinity}', 'Infinity literal')]:
        try:guard.load(raw)
        except ValueError:mutations.append({'name':label,'rejected':True})
        else:raise RuntimeError(label+' survived')
    journal=out/'SIMULATED_METADATA_ATTEMPTS.jsonl'
    auth=guard.SourceAuthorization(verified,'primary',journal)
    for source in guard.SOURCES:
        for j in range(64):
            auth('real_source_taylor',{'source':source,'center':str(Fraction(-9,2)+Fraction(2*j+1,128))})
    if auth.complete()['source_bundle_constructions']!=128:
        raise RuntimeError('exact fabricated event membership failed')
    journal_raw=journal.read_bytes()
    rows=[guard.load(line) for line in journal_raw.splitlines()]
    if len(rows)!=128 or [r['index'] for r in rows]!=list(range(1,129)):
        raise RuntimeError('Persisted authorization metadata chronology failed')
    if any(r['route']!='primary' or set(r)!= {'index','route','source','center'} for r in rows):
        raise RuntimeError('Persisted metadata includes unexpected fields')
    checks.append('all 128 simulated metadata attempts persisted before any possible caller; no numerical callback')
    checks.append('all 128 simulated authorization metadata events accepted; no source callback invoked')
    for label,event,detail in [
        ('repeated event','real_source_taylor',{'source':'positive_B','center':str(Fraction(-9,2)+Fraction(1,128))}),
        ('out-of-domain event','real_source_taylor',{'source':'positive_B','center':'0'}),
        ('unknown source event','real_source_taylor',{'source':'other','center':str(Fraction(-9,2)+Fraction(1,128))}),
        ('unknown event','source_decode',{'source':'positive_B','center':'0'}),
    ]:
        try:auth(event,detail)
        except ValueError:mutations.append({'name':label,'rejected':True})
        else:raise RuntimeError(label+' survived')
    if journal.read_bytes()!=journal_raw:
        raise RuntimeError('Rejected source metadata altered persisted chronology')
    checks.append('rejected metadata does not mutate the persisted attempt journal')
    try:guard.SourceAuthorization(verified,'primary',journal)
    except FileExistsError:mutations.append({'name':'occupied authorization journal','rejected':True})
    else:raise RuntimeError('Occupied authorization journal was reused')
    broken_journal=out/'FABRICATED_MISSING_DIRECTORY'/'ATTEMPTS.jsonl'
    try:guard.SourceAuthorization(verified,'primary',broken_journal)
    except FileNotFoundError:mutations.append({'name':'unwritable authorization journal','rejected':True})
    else:raise RuntimeError('Missing journal directory survived')
    failed_append=out/'SIMULATED_FAILED_APPEND.jsonl'
    blocked_auth=guard.SourceAuthorization(verified,'primary',failed_append)
    failed_append.unlink();failed_append.mkdir()
    try:blocked_auth('real_source_taylor',{'source':'positive_B',
                'center':str(Fraction(-9,2)+Fraction(1,128))})
    except IsADirectoryError:mutations.append({'name':'attempt journal append failure','rejected':True})
    else:raise RuntimeError('Source authorization continued after journal append failure')
    if blocked_auth.events or blocked_auth.seen:
        raise RuntimeError('Failed journal append recorded a completed authorization')
    checks.append('journal append failure blocks authorization before any possible numerical callback')
    try:guard.SourceAuthorization(verified,'primary').complete()
    except ValueError:mutations.append({'name':'missing authorization events','rejected':True})
    else:raise RuntimeError('Incomplete callback universe survived')
    for name,pin in original_pins.items():
        data=Path(name).read_bytes()
        if {'sha256':sha(data),'bytes':len(data)}!=pin:
            raise RuntimeError('Staged author bytes changed during guard review; rerun fresh')
    checks.append('all copied author source hashes unchanged during review')
    if any(name in sys.modules for name in ('flint','numpy','route','verified_moments','source_models','generic_operator')):
        raise RuntimeError('A numerical or actual route module was imported')
    result={'status':'PASS_FABRICATED_STDLIB_GUARD_REVIEW_ONLY',
        'scope':'Not a genuine public freeze or scientific execution GO',
        'actual_source_callbacks':0,'archive_arrays_decoded':0,'remote_requests':0,
        'numerical_method_imports':0,'checks':checks,'mutation_controls':mutations,
        'copied_source_pins':original_pins,'guard_sha256':original_pins[str(BASE/'execution-design/registration_guard.py')]['sha256'],
        'test_source_sha256':sha(Path(__file__).read_bytes())}
    dump(out/'REVIEW_RECEIPT.json',result)
    print(json.dumps({'status':result['status'],'checks':len(checks),'mutations':len(mutations),'guard_sha256':result['guard_sha256']}))

if __name__=='__main__':
    main()
