#!/usr/bin/env python3
"""Exact one-sided polynomial derivative audit at the formerly nonsmooth point."""
from fractions import Fraction as F
from math import comb
import argparse, contextlib, io, json, runpy, sys
from pathlib import Path

def require(x,m):
    if not x: raise RuntimeError(m)
p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True);args=p.parse_args()
sys.argv=[str(args.candidate),'--repo-root',str(args.repo_root)]
with contextlib.redirect_stdout(io.StringIO()): d=runpy.run_path(str(args.candidate))
ridge,obj,free=d['ridge'],d['obj'],d['free']; M=d['p']['matrices']['M'];delta=d['D']
x,r,h=F(1,1000),F(1499,1000),F(1,65536);n=11
require(n*h<delta,'stencil must stay inside one smoothed knot window')
require(all(abs((r-(x+j*h))-F(383,256))>delta and abs((r-(x+j*h))-F(3,2))>delta for j in range(-n,n+1)),'other endpoint crosses smoothing window')
def f(x):
    return (free(x,r)['dx'][0]+obj.feedback(x,r,M)['dx'][0])
f0=f(x)
weights=[-sum(F(1,k) for k in range(1,n+1))]+[F((-1)**(k+1)*comb(n,k),k) for k in range(1,n+1)]
right=sum(weights[k]*f(x+k*h) for k in range(n+1))/h
left=-sum(weights[k]*f(x-k*h) for k in range(n+1))/h
require(right==left,'old knot has a residual second-derivative jump')
print(json.dumps({'status':'PASS','former_kink':['1/1000','1499/1000'],'method':'exact 12-point one-sided polynomial derivative stencils','right_minus_left_smoothed_candidate_xx':str(right-left),'degree_at_most':11,'scope':'one local consistency check on nominal candidate; not an actual-field accuracy or PDE residual certificate'},sort_keys=True))
