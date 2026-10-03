#!/usr/bin/env python3
"""Symbolic coefficient bridge to the independently prepared tail source jets."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import traceback
import sympy as S
from source_jet import POLYNOMIALS


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tail-module",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise RuntimeError("Refusing to overwrite source-jet symbolic evidence")
    provenance={"verifier_sha256":sha(__file__),"primary_source_jet_sha256":sha(Path(__file__).with_name("source_jet.py")),
                "tail_source_jet_sha256":sha(args.tail_module),"python_optimization":sys.flags.optimize,"sympy":S.__version__}
    try:
        spec=importlib.util.spec_from_file_location("independent_tail_polynomial_constants",args.tail_module)
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        u=S.symbols("u"); d=1-u*u
        poly=lambda c: sum(S.Integer(x)*u**(len(c)-i-1) for i,x in enumerate(c))
        checks=[]
        for source in ("positive_B","signed_uB"):
            for n in range(6):
                primary=poly(POLYNOMIALS[n])
                if source=="signed_uB":
                    primary=u*primary+(n*d*d*poly(POLYNOMIALS[n-1]) if n else 0)
                tail=poly(module.DERIVATIVE_POLYNOMIALS[source][n])
                if S.expand(primary-tail)!=0:
                    raise RuntimeError("Independent source polynomial mismatch: "+source+" derivative "+str(n))
                checks.append({"source":source,"derivative_order":n,"passed":True})
        report={"passed":True,"check_count":len(checks),"checks":checks,"provenance":provenance,
                "scope":"Pure symbolic coefficients only; no physical source, derivative norm, root, quadrature, mode, or response evaluated."}
    except BaseException as exc:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps({"passed":False,"provenance":provenance,"exception":repr(exc),"traceback":traceback.format_exc()},indent=2)+"\n")
        raise
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("Pure symbolic independent source-jet bridge passed:",len(checks),"identities")


if __name__=="__main__":
    main()
