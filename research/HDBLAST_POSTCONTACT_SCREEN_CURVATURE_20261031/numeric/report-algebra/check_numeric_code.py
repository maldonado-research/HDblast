#!/usr/bin/env python3
"""Focused arithmetic and finite cover audit of the independent evaluator."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
if args.output.exists() or args.output.is_symlink():raise ValueError('FRESH_OUTPUT')
TARGET=HERE.parent/'independent_numeric.py'
spec=importlib.util.spec_from_file_location('reviewed_numeric_code',TARGET)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
G=1<<112

def bounds(lo,hi):return F((lo*G).__floor__(),G),F((hi*G).__ceil__(),G)
def need(ok,name):
    if not ok:raise ValueError(name)

checks=0;negative_inverses=0
for n in range(-20,21):
 for k in range(n,21):
  x=m.Box(F(n,3),F(k,7)) if F(n,3)<=F(k,7) else m.Box(F(k,7),F(n,3))
  a,b=x.lo,x.hi
  need((x*x).row()==list(map(str,bounds(min(a*a,a*b,b*b),max(a*a,a*b,b*b)))),'MULTIPLICATION');checks+=1
  for p in (2,4,6):
   lo=F(0) if a<=0<=b else min(a**p,b**p);hi=max(a**p,b**p)
   y=x**p;need(y.lo<=lo<=hi<=y.hi,'EVEN_POWER');checks+=1
  if a>0 or b<0:
   z=x.inv();need(z.row()==list(map(str,bounds(1/b,1/a))),'INVERSE');checks+=1
   negative_inverses+=b<0
T,R=F(5,4),F(1,8);U,L=T+R,T-R
count=0;omitted=0
for i in range(8):
 a,b=T*i/8,T*(i+1)/8
 for j in range(8):
  c,d=U*j/8,U*(j+1)/8
  xmax=min(b,U-c,L+d);rlo=max(c,F(0),a-L);rhi=min(d,U-a)
  # At s=a, the nested physical radial slices attain the projection extremes.
  exists=(a<=b and max(c,F(0),a-L)<=min(d,U-a))
  algorithm_keeps=not(xmax<a or rhi<rlo)
  need(exists==algorithm_keeps,'CELL_OMISSION_EXACT');checks+=1
  if not algorithm_keeps:omitted+=1;continue
  count+=1
  # Endpoints of all relevant affine pieces suffice for these inequalities.
  for s in sorted({a,xmax,max(a,min(xmax,L)),max(a,min(xmax,U-c)),max(a,min(xmax,L+d))}):
   lower=max(c,F(0),s-L);upper=min(d,U-s)
   if lower<=upper:
    need(a<=s<=xmax and rlo<=lower<=upper<=rhi,'PROJECTED_RECTANGLE_CONTAINS_SLICE');checks+=1
need(count==43 and omitted==21,'FIXED_COVER_COUNT')
report={'status':'PASS_FOCUSED_NUMERIC_CODE_REVIEW_CHECKS','checks':checks,'negative_inverse_cases':negative_inverses,'clipped_cover_cells':count,'omitted_empty_cells':omitted,'reviewed_source_sha256':hashlib.sha256(TARGET.read_bytes()).hexdigest(),'scope':'Imports only the reviewed module definitions; does not call scientific model input functions or rerun numerical source evaluations.'}
args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
print(json.dumps(report))
