#!/usr/bin/env python3
"""Validate registered numerical output without importing/evaluating sources.

Authentication of the frozen checkpoint, proof and runtime precedes this CLI
in the outer registered driver. Scientific negative results are serialized and
return zero; malformed execution artifacts return one.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

from operator_probe_contract import load_json_no_duplicates, validate_pair, validate_payload


def read_output(path):
    raw=path.read_bytes()
    return load_json_no_duplicates(raw),hashlib.sha256(raw).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path)
    parser.add_argument('--primary-output',type=Path)
    parser.add_argument('--independent-output',type=Path)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    if args.report.exists() or args.report.is_symlink():
        raise RuntimeError('Choose a fresh validator report path')
    code=0
    try:
        if args.input is not None:
            if args.primary_output is not None or args.independent_output is not None:
                raise ValueError('Choose single-output or both-route mode')
            value,pin=read_output(args.input)
            report=validate_payload(value)
            report['output_sha256']=pin
        else:
            if args.primary_output is None or args.independent_output is None:
                raise ValueError('Both route outputs are required')
            p,pin_p=read_output(args.primary_output)
            i,pin_i=read_output(args.independent_output)
            report=validate_pair(p,i)
            report['output_sha256']={'primary':pin_p,'independent':pin_i}
    except (ValueError,KeyError,TypeError,OSError) as error:
        report={'status':'FAIL_EXECUTION','reason':str(error),
                'scientific_classification':'UNCOMPUTED',
                'outer_bootstrap_proof_and_resource_checks_required':True}
        code=1
    args.report.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return code


if __name__=='__main__':
    raise SystemExit(main())
