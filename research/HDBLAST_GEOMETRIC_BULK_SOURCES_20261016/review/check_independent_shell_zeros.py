#!/usr/bin/env python3
"""Independent endpoint-rational sign proof for the four shell mean zeros."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]
sys.dont_write_bytecode=True;sys.path.insert(0,str(ROOT/'numeric'))
import actual_first_jets as actual

def pair(x):return (x.lo,x.hi)
def point(x):return (F(x),F(x))
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[1],-a[0])
def mul(a,b):
    v=[x*y for x,y in product(a,b)];return (min(v),max(v))
def square(a):return (0 if a[0]<=0<=a[1] else min(x*x for x in a),max(x*x for x in a))
def inverse(a):
    if a[0]<=0<=a[1]:raise ValueError('ZERO_IN_DENOMINATOR')
    return (1/a[1],1/a[0])
def divide(a,b):return mul(a,inverse(b))
def row(a):return list(map(str,a))

def bracket(left,right):
    _,s,b=actual.actual_jet(actual.I(left,right),actual.I(0))
    v,vt,vtt=map(pair,(s['s'],s['st'],s['stt']))
    if v[0]<=0<=v[1]:raise ValueError('V_ZERO_NOT_EXCLUDED')
    N=add(add(mul(v,vtt),neg(square(vt))),mul(pair(b['theta']),mul(v,vt)))
    if N[1]>=0:raise ValueError('H_STRICT_MONOTONICITY_NOT_ESTABLISHED')
    if b['c'].lo<=0 or b['rho'].lo<=0:raise ValueError('POSITIVE_GEOMETRIC_SCALES_NOT_ESTABLISHED')
    ends=[]
    for j,x in enumerate((left,right)):
        _,es,eb=actual.actual_jet(actual.I(x),actual.I(0))
        ev,et=map(pair,(es['s'],es['st']))
        H=divide(mul(pair(eb['c']),et),ev)
        if j==0 and not(H[0]>0 and H[0]**2>8):raise ValueError('LEFT_SQRT8_THRESHOLD')
        if j==1 and not(H[1]<0 and H[1]**2>8):raise ValueError('RIGHT_MINUS_SQRT8_THRESHOLD')
        mean=divide(add(neg(square(et)),divide(mul(point(8),square(ev)),square(pair(eb['c'])))),mul(point(3),square(pair(eb['rho']))))
        if mean[1]>=0:raise ValueError('ENDPOINT_MEAN_NOT_NEGATIVE')
        ends.append({'x':str(x),'H':row(H),'mean_second_derivative':row(mean)})
    return {'x_bracket':list(map(str,(left,right))),'v':row(v),'vtt':row(vtt),'N':row(N),'endpoints':ends,
        'theorem':'v nonzero; H_t=c*N/v^2<0 throughout bracket; H runs above sqrt(8) to below -sqrt(8); mean=v^2(8-H^2)/(3*rho^2*c^2) has exactly two simple zeros, with positive mean between and negative mean outside within the bracket.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    if args.output.exists():raise ValueError('OUTPUT_ALREADY_EXISTS')
    brackets=[bracket(F(5,4),F(11,8)),bracket(F(13,8),F(7,4))]
    paths=[Path(__file__),ROOT/'numeric/actual_first_jets.py',ROOT/'numeric/interval_core.py',ROOT/'numeric/exterior_error_caps.py',ROOT/'shell-invariant/intrinsic_scalar.py',ROOT/'shell-invariant/SHELL_INTRINSIC_SCALAR.md',ROOT/'shell-invariant/check_shell_invariant.py']
    out={'status':'PASS_INDEPENDENT_FOUR_SIMPLE_SHELL_MEAN_ZEROS','brackets':brackets,
        'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'arithmetic':'Exact Fraction endpoints; independent interval products, powers, reciprocals and sums, no producer interval arithmetic after input jets.',
        'scope':'Conditional inherited actual first PDE/jet premise; intrinsic quadratic coefficient mean on fixed background foliation; not a finite-amplitude force or nonlinear claim.'}
    args.output.write_text(json.dumps(out,indent=2)+'\n');print(out['status'])

if __name__=='__main__':main()
