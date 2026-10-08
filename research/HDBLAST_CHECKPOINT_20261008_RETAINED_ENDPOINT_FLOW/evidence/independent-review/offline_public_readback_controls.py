"""Offline public-readback protocol controls; every API reply is fabricated.

This creates explicitly labelled MOCK receipts, never network requests or real
GO authorization. No candidate code is imported or evaluated.
"""
from pathlib import Path
import base64
import contextlib
import hashlib
import io
import json
import shutil
import sys

BASE=Path(__file__).resolve().parent
HELPER=BASE.parent/'root/public_endpoint_readback_go.py'
def digest(raw):return hashlib.sha256(raw).hexdigest()
def jbytes(value):return (json.dumps(value,sort_keys=True)+'\n').encode()
def blob_id(raw):return hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
def require(ok,msg):
    if not ok:raise RuntimeError(msg)

def run(label):
    destination=BASE/label;destination.mkdir();source=HELPER.read_bytes()
    (destination/'reviewed_helper.py.snapshot').write_bytes(source)
    plans=('positive','registration_race','review_race','private_repository','anonymous_mismatch',
           'wrong_remote_path','wrong_blob_bytes','wrong_blob_size','support_pin_mismatch',
           'extra_support_pin','duplicate_support_name','bool_zero_counter','positive_counter',
           'missing_mandatory_source','empty_source','bool_source_size','unlisted_core_file',
           'symlink_support','hardlink_support','unsafe_support_path','output_inside_core',
           'unknown_empty_directory','hidden_execute_only_package','unreadable_registered_directory')
    results=[]
    for label_case in plans:
        root=destination/label_case;root.mkdir();shutil.copytree(BASE/'frozen-source-004',root/'core')
        registration=json.loads((BASE/'REVIEWED_REGISTRATION_004.json').read_bytes())
        # The fabricated local tree must have exactly the parent directories of files.
        for directory in sorted((root/'core').rglob('*'),key=lambda p:len(p.parts),reverse=True):
            if directory.is_dir() and not any(directory.iterdir()):directory.rmdir()
        inaccessible_directory=None
        (root/'evidence').mkdir();support=b'FABRICATED SUPPORTING REVIEW\n'
        spath=root/'evidence/support.json';spath.write_bytes(support)
        readme=b'OFFLINE FABRICATED README\n';(root/'README.md').write_bytes(readme)
        support_names=['evidence/support.json'];support_pins={'evidence/support.json':{'bytes':len(support),'sha256':digest(support)}}
        if label_case=='missing_mandatory_source':
            registration['files'].pop('engine/endpoint_engine.py');(root/'core/engine/endpoint_engine.py').unlink()
        if label_case=='empty_source':
            registration['files']={};shutil.rmtree(root/'core');(root/'core').mkdir()
        if label_case=='bool_source_size':registration['files']['PROTOCOL.md']['bytes']=True
        if label_case=='unlisted_core_file':(root/'core/unlisted.py').write_text('raise RuntimeError("MUST_NOT_IMPORT")\n')
        if label_case=='unknown_empty_directory':(root/'core/hidden').mkdir()
        if label_case=='hidden_execute_only_package':
            inaccessible_directory=root/'core/source/later_source';inaccessible_directory.mkdir()
            (inaccessible_directory/'__init__.py').write_text('raise RuntimeError("MUST_NOT_IMPORT")\n')
            inaccessible_directory.chmod(0o111)
        if label_case=='unreadable_registered_directory':
            inaccessible_directory=root/'core/source';inaccessible_directory.chmod(0o111)
        if label_case=='support_pin_mismatch':support_pins['evidence/support.json']['sha256']='0'*64
        if label_case=='extra_support_pin':support_pins['evidence/unlisted.json']={'bytes':0,'sha256':digest(b'')}
        if label_case=='duplicate_support_name':support_names.append('evidence/support.json')
        if label_case in ('symlink_support','hardlink_support'):
            outside=destination/(label_case+'-outside.txt');outside.write_bytes(support);spath.unlink()
            if label_case=='symlink_support':spath.symlink_to(outside)
            else:spath.hardlink_to(outside)
        if label_case=='unsafe_support_path':
            support_names=['../escape'];support_pins={'../escape':{'bytes':0,'sha256':digest(b'')}}
        regraw=jbytes(registration);(root/'FULL_REGISTRATION.json').write_bytes(regraw)
        review={'status':'GO_FOR_PROSPECTIVE_PUBLIC_FREEZE','registration_sha256':digest(regraw),
                'registered_file_pins':registration['files'],'independent_review_pass':True,
                'physical_execution_authorized_by_this_receipt':False,
                'new_target_evaluations_before_freeze':{key:0 for key in ('physical_source_callbacks',
                     'retained_array_decodes','stored_endpoint_comparisons','observational_likelihood_evaluations')},
                'supporting_evidence_files':support_names,'supporting_evidence_pins':support_pins}
        if label_case=='bool_zero_counter':review['new_target_evaluations_before_freeze']['physical_source_callbacks']=False
        if label_case=='positive_counter':review['new_target_evaluations_before_freeze']['physical_source_callbacks']=1
        reviewraw=jbytes(review);(root/'evidence/FINAL_PRE_FREEZE_REVIEW.json').write_bytes(reviewraw)
        namespace={'__name__':'offline_public_readback_control'}
        exec(compile(source,str(HELPER),'exec'),namespace)
        calls=[];blobs={};prefix=namespace['PREFIX']
        def fake_get(path,anonymous=False):
            calls.append({'path':path,'anonymous':anonymous})
            if path=='':
                if label_case=='registration_race':
                    changed=dict(registration);changed['UNREVIEWED_AFTER_CAPTURE']=True
                    (root/'FULL_REGISTRATION.json').write_bytes(jbytes(changed))
                if label_case=='review_race':
                    changed=dict(review);changed['UNREVIEWED_AFTER_CAPTURE']=True
                    (root/'evidence/FINAL_PRE_FREEZE_REVIEW.json').write_bytes(jbytes(changed))
                return {'full_name':'maldonado-research/HDblast','private':label_case=='private_repository'}
            if path.startswith('contents/'):
                relative=path.split('?ref=',1)[0].removeprefix('contents/'+prefix+'/')
                body=(root/relative).read_bytes();sha=blob_id(body);blobs[sha]=body
                reply={'type':'file','path':prefix+'/'+relative,'sha':sha,'size':len(body),
                       'encoding':'base64','content':base64.b64encode(body).decode()}
                if label_case=='anonymous_mismatch' and anonymous:reply['content']=base64.b64encode(b'WRONG README').decode()
                if label_case=='wrong_remote_path' and not anonymous:reply['path']='wrong/path'
                return reply
            if path.startswith('git/blobs/'):
                sha=path.removeprefix('git/blobs/');body=blobs[sha]
                if label_case=='wrong_blob_bytes':body=bytes([body[0]^1])+body[1:] if body else b'x'
                return {'sha':sha,'size':len(body)+(1 if label_case=='wrong_blob_size' else 0),
                        'encoding':'base64','content':base64.b64encode(body).decode()}
            raise RuntimeError('unexpected simulated API path:'+path)
        namespace['get']=fake_get
        output=(root/'core' if label_case=='output_inside_core' else root)/'MOCK_GO_NOT_AUTHORIZATION.json'
        sys.argv=['offline-fixture','--checkpoint',str(root),'--commit','a'*40,
                  '--registration-sha256',digest(regraw),'--review-sha256',digest(reviewraw),'--output',str(output)]
        log=io.StringIO();error=None
        try:
            with contextlib.redirect_stdout(log):namespace['main']()
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        if inaccessible_directory is not None:inaccessible_directory.chmod(0o755)
        (root/'CONTROL.log').write_text(log.getvalue()+(error or '')+'\n')
        if label_case=='positive':
            require(error is None and output.is_file(),'positive manufactured byte-readback failed:'+str(error))
            mock=json.loads(output.read_bytes())
            require(mock['registration_sha256']==digest(regraw) and mock['prospective_review_sha256']==digest(reviewraw),
                    'external pins not bound')
            require(mock['remote_files_verified']==len(registration['files'])+4,'complete source/support coverage differs')
            require(mock['remote_file_receipts']['FULL_REGISTRATION.json']['sha256']==digest(regraw) and
                    mock['remote_file_receipts']['evidence/FINAL_PRE_FREEZE_REVIEW.json']['sha256']==digest(reviewraw),
                    'downloaded control pins differ')
            passed=True
        else:passed=error is not None and not output.exists()
        results.append({'case':label_case,'passed':passed,'rejected':error is not None,'error':error,
                        'simulated_api_calls':len(calls),'mock_receipt_written':output.exists()})
    receipt={'status':'PASS_OFFLINE_PUBLIC_READBACK_CONTROLS' if all(x['passed'] for x in results) else 'FAIL_OFFLINE_PUBLIC_READBACK_CONTROLS',
             'helper_sha256':digest(source),'checks':len(results),'cases':results,
             'network_requests':0,'real_public_GO_issued':False,'retained_decodes':0,'physical_source_calls':0}
    with (destination/'RESULT.json').open('x') as out:json.dump(receipt,out,sort_keys=True,indent=2);out.write('\n')
    print(json.dumps(receipt,sort_keys=True))
    require(all(x['passed'] for x in results),'a public-readback protocol control failed')
if __name__=='__main__':run(sys.argv[1])
