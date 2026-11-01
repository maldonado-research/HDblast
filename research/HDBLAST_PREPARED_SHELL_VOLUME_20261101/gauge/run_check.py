"""Authenticated CLI orchestration; scientific modules load only after checks."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE))
from authenticate import authenticate, output_path, failure_payload


def main(kind):
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo-root',required=True,type=Path)
    parser.add_argument('--output',required=True)
    parser.add_argument('--mutation',default=None)
    args=parser.parse_args()
    try:
        dest=output_path(args.output,args.repo_root)
        before,pins=authenticate(args.repo_root)
        path=BASE/'checks'/('_'+kind+'_science.py')
        spec=importlib.util.spec_from_file_location('checked_'+kind,path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        if kind=='preparation':result=module.run(args.repo_root,pins,args.mutation)
        else:result=module.run(args.mutation)
        after,_=authenticate(args.repo_root)
        if before!=after:raise ValueError('authentication changed during science')
        result['authentication']={'before':before,'after':after,'equal':True}
        result['component']='prepared_shell_volume/'+kind
        payload=json.dumps(result,indent=2)+'\n'
        if dest is None:sys.stdout.write(payload)
        else:
            with dest.open('x') as stream:stream.write(payload)
    except Exception as exc:
        sys.stdout.write(failure_payload(exc,args.mutation))
        return 2
    return 0
