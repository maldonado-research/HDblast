#!/usr/bin/env python3
"""Real failed-candidate fixture: continue q only through degree two."""
import argparse,hashlib,importlib.util,json,sys
from fractions import Fraction as Q
from pathlib import Path
sys.dont_write_bytecode=True
def need(ok,why):
    if not ok:raise RuntimeError(why)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--helper',type=Path,required=True);ap.add_argument('--helper-sha256',required=True);ap.add_argument('--source-pins',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    need(not a.output.exists(),'FRESH_OUTPUT_REQUIRED');need(hashlib.sha256(a.helper.read_bytes()).hexdigest()==a.helper_sha256,'HELPER_HASH')
    spec=importlib.util.spec_from_file_location('independent_raw_fixture',a.helper);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    ev=m.Evaluator(a.source,a.source_pins);r=m.I(Q(7,4));rr=m.J([r,1/m.ELL]+[0]*4);den=(rr-1)*(2-rr);q=m.expjetminus(1/den-4);e=m.expjetminus(m.ELL*rr)
    raw_quadratic=m.J(q.c[:3]+[0]*3)
    # Shared jets0..2 cancel exactly, while raw exterior jets3..5 vanish.
    dq=m.J([0,0,0]+q.c[3:]);qp=m.P(q,raw_quadratic,[dq,m.C(0),m.C(0)])
    bg={name:m.P(z.b) for name,z in ev.background(r).items()}
    uu,_=m.registered(bg,qp,m.P(e));jump=uu[12].d[0].c[0]
    need(jump.minimum()>0,'CONTROL_MUST_EXHIBIT_NONZERO_VALUE_MISMATCH')
    need(jump.lo<=0<=jump.hi,'REQUIRED_RAW_QUADRATIC_SCALAR_VALUE_MATCH_FAILED:label=w:0:radius=7/4:left_minus_right='+str(jump.row()))
    a.output.write_text(json.dumps({'status':'UNEXPECTED_SUCCESS'})+'\n')
if __name__=='__main__':main()
