"""Offline API simulation only: demonstrate mismatched externally pinned readback.

Every 'public' byte is self-written fixture data returned by fake_get. This
does not issue a network request or a real public/physical GO.
"""
from pathlib import Path
import base64
import hashlib
import json
import shutil
import sys

BASE=Path(__file__).resolve().parent
HELPER=BASE.parent/'root/public_endpoint_readback_go.py'
raw=HELPER.read_bytes();(BASE/'initial_public_readback_helper.py.snapshot').write_bytes(raw)
namespace={'__name__':'offline_public_readback_under_review'}
exec(compile(raw,str(HELPER),'exec'),namespace)
def digest(raw):return hashlib.sha256(raw).hexdigest()
def blob_id(raw):return hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
def jbytes(value):return (json.dumps(value,sort_keys=True)+'\n').encode()
root=BASE/'fabricated-public-readback-pin-race';root.mkdir()
shutil.copytree(BASE/'frozen-source-003',root/'core')
regraw=(BASE/'REVIEWED_REGISTRATION_003.json').read_bytes();registration=json.loads(regraw)
(root/'FULL_REGISTRATION.json').write_bytes(regraw);(root/'README.md').write_bytes(b'OFFLINE FABRICATED README\n')
(root/'evidence').mkdir()
review={'status':'GO_FOR_PROSPECTIVE_PUBLIC_FREEZE','registration_sha256':digest(regraw),
        'registered_file_pins':registration['files'],'independent_review_pass':True,
        'physical_execution_authorized_by_this_receipt':False,
        'new_target_evaluations_before_freeze':{key:0 for key in ('physical_source_callbacks',
             'retained_array_decodes','stored_endpoint_comparisons','observational_likelihood_evaluations')},
        'supporting_evidence_files':[],'supporting_evidence_pins':{}}
reviewraw=jbytes(review);(root/'evidence/FINAL_PRE_FREEZE_REVIEW.json').write_bytes(reviewraw)
calls=[];blobs={}
def fake_get(path,anonymous=False):
    calls.append(path)
    if path=='':
        replacement=dict(registration);replacement['unreviewed_addition_after_initial_pin_check']=True
        (root/'FULL_REGISTRATION.json').write_bytes(jbytes(replacement))
        return {'full_name':'maldonado-research/HDblast','private':False}
    if path.startswith('contents/'):
        relative=path.split('?ref=',1)[0].removeprefix('contents/'+namespace['PREFIX']+'/')
        body=(root/relative).read_bytes();sha=blob_id(body);blobs[sha]=body
        return {'type':'file','path':namespace['PREFIX']+'/'+relative,'sha':sha,'size':len(body),
                'encoding':'base64','content':base64.b64encode(body).decode()}
    if path.startswith('git/blobs/'):
        sha=path.removeprefix('git/blobs/');body=blobs[sha]
        return {'sha':sha,'size':len(body),'encoding':'base64','content':base64.b64encode(body).decode()}
    raise RuntimeError('unexpected fake API path:'+path)
namespace['get']=fake_get
output=root/'MOCK_GO_NOT_AUTHORIZATION.json'
sys.argv=['offline-fixture','--checkpoint',str(root),'--commit','a'*40,
          '--registration-sha256',digest(regraw),'--review-sha256',digest(reviewraw),'--output',str(output)]
result={'status':'INITIAL_PUBLIC_READBACK_PIN_RACE_CONTROL','helper_sha256':digest(raw),
        'network_requests':0,'retained_decodes':0,'physical_source_calls':0}
try:namespace['main']()
except Exception as error:result.update(rejected=True,error_type=type(error).__name__,error=str(error))
else:
    generated=json.loads(output.read_bytes())
    result.update(rejected=False,claimed_registration_sha256=generated['registration_sha256'],
                  downloaded_registration_sha256=generated['remote_file_receipts']['FULL_REGISTRATION.json']['sha256'])
result['simulated_api_calls']=len(calls)
with (BASE/'INITIAL_PUBLIC_READBACK_PIN_RACE.json').open('x') as out:json.dump(result,out,sort_keys=True,indent=2);out.write('\n')
print(json.dumps(result,sort_keys=True))
