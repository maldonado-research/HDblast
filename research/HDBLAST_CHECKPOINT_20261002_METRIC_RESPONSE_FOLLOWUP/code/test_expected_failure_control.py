#!/usr/bin/env python3
"""Saved-evidence classifier guards; no producer, profile or quadrature calls."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('original_failure_classifier',HERE/'reproduce_original_expected_failure.py')
M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--saved-failure',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise RuntimeError('Refusing overwrite')
    args.output.mkdir(parents=True);checks=[]
    positive=M.qualify_saved_failure(args.saved_failure,1)
    if not (positive['control_status']=='EXPECTED_FAILURE_CONFIRMED' and positive['original_scientific_status']=='ORIGINAL_EXPERIMENT_FAILED'):raise RuntimeError('Saved failure classification mismatch')
    checks.append('preserved original saved failure qualifies without producer execution')
    def reject(name,mutation,exit_code=1,timed_out=False):
        q=args.output/name;shutil.copytree(args.saved_failure,q);mutation(q)
        try:M.qualify_saved_failure(q,exit_code,timed_out)
        except RuntimeError:checks.append(name+' rejected')
        else:raise RuntimeError('Classifier mutation survived: '+name)
    def change(q,file,key,value):
        d=json.loads((q/file).read_text());d[key]=value;(q/file).write_text(json.dumps(d,indent=2)+'\n')
    reject('unexpected_success_exit',lambda q:None,0)
    reject('timeout',lambda q:None,1,True)
    reject('wrong_exception',lambda q:change(q,'failure.json','exception',"RuntimeError('different failure')"))
    reject('wrong_route',lambda q:change(q,'failure.json','traceback','finite_history failure'))
    reject('wrong_source_time',lambda q:change(q,'partial_results.json','active',{'source':'positive_B','eta':-4}))
    reject('different_completed_rows',lambda q:change(q,'partial_results.json','rows',[]))
    reject('wrong_registration',lambda q:change(q,'started.json','registration_sha256','0'*64))
    reject('wrong_producer',lambda q:change(q,'started.json','producer_sha256','0'*64))
    reject('completed_results_exist',lambda q:(q/'results.json').write_text('{}'))
    result={'status':'PASS_SAVED_EVIDENCE_CLASSIFIER_ONLY','check_count':len(checks),'checks':checks,'physical_evaluations':0,'child_producers_executed':0,'fresh_original_reproduction_performed':False,'python_optimization':sys.flags.optimize,'helper_sha256':hashlib.sha256((HERE/'reproduce_original_expected_failure.py').read_bytes()).hexdigest(),'test_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'saved_input_sha256':{name:hashlib.sha256((args.saved_failure/name).read_bytes()).hexdigest() for name in ('failure.json','partial_results.json','started.json')}}
    (args.output/'CHECKS.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({k:result[k] for k in ('status','check_count','physical_evaluations','child_producers_executed')}))

if __name__=='__main__':main()
