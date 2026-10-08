"""Meaningful fabricated-only registration and payload controls."""
from pathlib import Path
from fractions import Fraction
import argparse
import gzip
import hashlib
import json
import shutil
import sys
import tempfile
from registration_guard import authenticate_local,authenticate,require,load,sha,SourceAuthorization,DecodeAuthorization

def expect_failure(call,message):
    try:call()
    except (ValueError,OSError) as error:
        return {'control':message,'rejected':True,'reason':str(error)}
    raise ValueError('Control unexpectedly passed: '+message)

def write_json(path,obj):
    path.write_text(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n')

def run(candidate,regsha,evidence,output=None):
    evidence.mkdir(parents=True,exist_ok=False)
    controls=[];authenticate_local(candidate,regsha)
    controls.append(expect_failure(lambda:authenticate_local(candidate,'0'*64),'external registration pin tamper'))
    for label,mutation in (
        ('registered source bytes',lambda r:(r/'source/route.py').write_bytes((r/'source/route.py').read_bytes()+b'\n# mutation\n')),
        ('unregistered helper',lambda r:(r/'execution/rogue.py').write_text('raise RuntimeError("unregistered")\n')),
        ('registered payload symlink',lambda r:((r/'prior/independent_DATA.json').unlink(),
                                              (r/'prior/independent_DATA.json').symlink_to(candidate/'prior/independent_DATA.json'))),
    ):
        clone=evidence/label.replace(' ','_');shutil.copytree(candidate,clone)
        mutation(clone)
        controls.append(expect_failure(lambda r=clone:authenticate_local(r,regsha),label))
    missing=evidence/'missing-PUBLIC_GO.json'
    controls.append(expect_failure(lambda:authenticate(candidate,regsha,missing,'0'*64),'missing external PUBLIC_GO'))
    wrong=evidence/'wrong-scope-PUBLIC_GO.json'
    write_json(wrong,{'status':'PASS_REMOTE_REGISTERED_BD_PREHISTORY_GO',
                      'input_scope':'NO_RETAINED_ARRAYS_BD_PREHISTORY_ONLY'})
    controls.append(expect_failure(lambda:authenticate(candidate,regsha,wrong,sha(wrong.read_bytes())),
                                   'old nine-probe source contract cannot authorize selected decodes'))
    verified=authenticate_local(candidate,regsha)
    controls.append(expect_failure(lambda:SourceAuthorization(verified,evidence/'blocked-source.jsonl',[Fraction(i) for i in range(11)]),
                                   'fabricated guard cannot authorize a physical callback'))
    controls.append(expect_failure(lambda:DecodeAuthorization(verified,evidence/'blocked-decode.jsonl',{}),
                                   'fabricated guard cannot authorize retained arrays'))
    if output:
        from readback import validate_output
        validate_output(output,candidate,True)
        clone=evidence/'counter-tamper';shutil.copytree(output,clone)
        data=load((clone/'DATA.json').read_bytes());data['counter']['target_node_evaluations']-=1
        for name in ('DATA.json','SCIENCE_SUMMARY.json'):write_json(clone/name,data)
        controls.append(expect_failure(lambda:validate_output(clone,candidate,True),'exported all-node counter weakened'))
        clone=evidence/'component-tamper';shutil.copytree(output,clone)
        # Preserve a complete exported stream but alter one target endpoint.
        # Scientific summaries/entry pins are held unchanged; standalone
        # recomputation must reject the altered numerical payload.
        with gzip.open(clone/'NODE_TARGETS.jsonl.gz','rb') as z:lines=z.readlines()
        row=load(lines[0]);row['target_dyadic96']['U']['real'][0]-=1<<96
        lines[0]=(json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode()
        with (clone/'NODE_TARGETS.jsonl.gz').open('wb') as raw:
            with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as z:z.writelines(lines)
        controls.append(expect_failure(lambda:validate_output(clone,candidate,True),'target component payload tamper'))
        clone=evidence/'raw-state-tamper';shutil.copytree(output,clone)
        with gzip.open(clone/'NODE_TARGETS.jsonl.gz','rb') as z:lines=z.readlines()
        row=load(lines[0]);row['stored_raw_u'][0]='0/1'
        lines[0]=(json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode()
        with (clone/'NODE_TARGETS.jsonl.gz').open('wb') as raw:
            with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as z:z.writelines(lines)
        controls.append(expect_failure(lambda:validate_output(clone,candidate,True),'raw state payload tamper'))
    result={'status':'PASS_FABRICATED_REGISTRATION_PAYLOAD_CONTROLS','controls':controls,'count':len(controls),
            'physical_source_callbacks':0,'retained_arrays_opened':0,'public_GO_created':False}
    write_json(evidence/'CONTROLS.json',result)
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--candidate',required=True);p.add_argument('--registration-sha256',required=True)
    p.add_argument('--evidence',required=True);p.add_argument('--output');a=p.parse_args()
    print(json.dumps(run(Path(a.candidate),a.registration_sha256,Path(a.evidence),Path(a.output) if a.output else None),sort_keys=True))

if __name__=='__main__':main()
