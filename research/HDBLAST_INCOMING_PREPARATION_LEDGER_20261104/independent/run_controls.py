#!/usr/bin/env python3
"""One canonical normal/optimized pair and eight real rejection children."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
PRODUCER_SHA='675fc7af7c9decd96154b617bd52b85dc799c2c223671b169c1942c4db2d256c'
PIN_SHA='d7ecbf4f0b6721b1ab7e283030cc6cb8d2b05e5f5c6a45f2ec8043baeae6124d'

class Failure(Exception):
    pass

def need(ok,witness):
    if not ok:
        raise Failure(witness)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def regular(p):
    return p.is_file() and not p.is_symlink() and all(not q.is_symlink() for q in p.absolute().parents)

def canonical(v):
    return (json.dumps(v,sort_keys=True,indent=2)+'\n').encode()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['repo-root','source-candidate','theory-candidate','output']:
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    paths=[getattr(args,key) for key in ['repo_root','source_candidate','theory_candidate','output']]
    need(all(p.is_absolute() for p in paths),'ABSOLUTE_ARGUMENT_PATHS')
    need(all(not q.is_symlink() for p in paths for q in [p,*p.parents]),'CONTROL_PATH_SYMLINK')
    repo,source,theory=map(lambda p:p.resolve(strict=True),paths[:3])
    out=paths[3].resolve()
    need(not out.exists(),'CONTROL_FRESH_OUTPUT')
    need(out.parent.is_dir(),'CONTROL_OUTPUT_PARENT')
    need(all(base != out and base not in out.parents for base in [repo,HERE]),'CONTROL_OUTPUT_OUTSIDE_SOURCE')
    need(regular(HERE/'recompute.py') and sha(HERE/'recompute.py') == PRODUCER_SHA,'CONTROL_PRODUCER_BYTES')
    need(regular(HERE/'SOURCE_PINS.json') and sha(HERE/'SOURCE_PINS.json') == PIN_SHA,'CONTROL_PIN_BYTES')
    pins=json.loads((HERE/'SOURCE_PINS.json').read_bytes())['files']
    artifacts=out.with_suffix('.artifacts')
    need(not artifacts.exists(),'CONTROL_FRESH_ARTIFACTS')
    artifacts.mkdir()
    positives=[]
    negatives=[]
    cases={'source_pin':'SOURCE_PIN_MISMATCH','candidate_source_budget':'CANDIDATE_SOURCE_LEDGER_MISMATCH','malformed_fraction':'NUMERIC_PARSE_ERROR','occupied_output':'FRESH_OUTPUT_FILE'}
    def child(root,source_path,theory_path,target,optimized):
        flags=['-I','-S','-B',*(['-O'] if optimized else [])]
        command=[sys.executable,*flags,str(HERE/'recompute.py'),'--repo-root',str(root),'--source-candidate',str(source_path),'--theory-candidate',str(theory_path),'--output',str(target)]
        completed=subprocess.run(command,capture_output=True,text=True,timeout=60)
        need(completed.stderr == '' and len(completed.stdout.splitlines()) == 1,'EXACT_CHILD_OUTPUT_PROTOCOL')
        document=json.loads(completed.stdout)
        need(type(document) is dict,'CHILD_OBJECT_WITNESS')
        return completed,document,flags
    for optimized in [False,True]:
        target=artifacts/('OPTIMIZED.json' if optimized else 'NORMAL.json')
        completed,document,flags=child(repo,source,theory,target,optimized)
        need(completed.returncode == 0 and document == {'status':'PASS','witness':'INDEPENDENT_INCOMING_RECONSTRUCTION','finite_checks':45,'output_sha256':sha(target)},'CANONICAL_CHILD_WITNESS')
        report=json.loads(target.read_bytes())
        need(report['status'] == 'PASS_INDEPENDENT_INCOMING_PREPARATION_RECONSTRUCTION','CANONICAL_REPORT_STATUS')
        positives.append({'optimized':optimized,'flags':flags,'returncode':completed.returncode,'stdout':document,'stderr_empty':True,'report_sha256':sha(target)})
        for case,witness in cases.items():
            with tempfile.TemporaryDirectory(prefix='hdblast-independent-control-',dir=out.parent) as temporary:
                tmp=Path(temporary)
                fixture=tmp/'repo'
                fixture.mkdir()
                for name in pins:
                    destination=fixture/name
                    destination.parent.mkdir(parents=True,exist_ok=True)
                    shutil.copyfile(repo/name,destination)
                sf,tf=tmp/'source-candidate.json',tmp/'theory-candidate.json'
                shutil.copyfile(source,sf)
                shutil.copyfile(theory,tf)
                output=tmp/'result.json'
                previous=None
                if case == 'source_pin':
                    p=fixture/next(iter(pins))
                    p.write_bytes(p.read_bytes()+b'\n')
                elif case == 'candidate_source_budget':
                    changed=json.loads(sf.read_bytes())
                    row=changed['budgets']['wedge_zero_retained_cell_caps']
                    row['energy_source_ledger_upper']=str(F(row['energy_source_ledger_upper'])+1)
                    sf.write_bytes(canonical(changed))
                elif case == 'malformed_fraction':
                    changed=json.loads(tf.read_bytes())
                    changed['forced_energy']['bulk_square_factor']='1/0'
                    tf.write_bytes(canonical(changed))
                elif case == 'occupied_output':
                    previous=b'preexisting output bytes must survive\n'
                    output.write_bytes(previous)
                completed,document,flags=child(fixture,sf,tf,output,optimized)
                need(completed.returncode == 2 and document == {'status':'FAIL','witness':witness},'EXACT_REJECTION_WITNESS')
                preserved=(not output.exists()) if previous is None else output.read_bytes() == previous
                need(preserved,'REJECTED_ARTIFACT_CONTRACT')
                negatives.append({'case':case,'optimized':optimized,'flags':flags,'returncode':completed.returncode,'stdout':document,'stderr_empty':True,'output_absent':previous is None and not output.exists(),'preexisting_bytes_preserved':previous is not None and preserved})
    need(positives[0]['report_sha256'] == positives[1]['report_sha256'],'NORMAL_OPTIMIZED_REPORT_IDENTITY')
    need(sha(HERE/'recompute.py') == PRODUCER_SHA and sha(HERE/'SOURCE_PINS.json') == PIN_SHA,'CONTROL_LOCAL_BYTES_UNCHANGED')
    report={'schema':'HDBLAST_INCOMING_INDEPENDENT_CONTROLS_V1','status':'PASS_STRICT_INDEPENDENT_INCOMING_CONTROLS','canonical_positive_children':len(positives),'actual_negative_children':len(negatives),'distinct_rejection_cases':len(cases),'canonical_children':positives,'negative_children':negatives,'canonical_pair_sha256':positives[0]['report_sha256'],'new_PDE_solves':0,'new_field_background_bump_or_source_queries':0,'scope':'One canonical pair and eight actual isolated rejection children; controls authenticate and compare finite exact reconstructions, not continuum or past-existence premises.'}
    with out.open('xb') as stream:
        stream.write(canonical(report))
    print(json.dumps({'status':'PASS','witness':'STRICT_INDEPENDENT_INCOMING_CONTROLS','canonical_positive_children':len(positives),'actual_negative_children':len(negatives),'output_sha256':sha(out)},sort_keys=True))

if __name__ == '__main__':
    try:
        main()
    except Failure as exc:
        print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True))
        sys.exit(2)
    except ArithmeticError:
        print(json.dumps({'status':'FAIL','witness':'NUMERIC_PARSE_ERROR'},sort_keys=True))
        sys.exit(2)
    except (ValueError,OSError,KeyError,TypeError,subprocess.SubprocessError):
        print(json.dumps({'status':'FAIL','witness':'RUNTIME_INPUT_ERROR'},sort_keys=True))
        sys.exit(2)
