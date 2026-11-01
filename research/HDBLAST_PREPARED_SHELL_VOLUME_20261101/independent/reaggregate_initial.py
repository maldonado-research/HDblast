#!/usr/bin/env python3
"""Independent quadrature and primitive reaggregation from retained raw cells."""
from fractions import Fraction as F
from math import factorial
import json,sys
from pathlib import Path
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent))
from portable_io import arguments,authenticate,finish
from reconstruct_response import Box,point,readbox,load,need
ELL=F(1,10**6)
def rounded(b):
    g=1<<160
    return Box(F((b.lo*g).__floor__(),g),F((b.hi*g).__ceil__(),g))
def exponential_minus(x):
    need(x>=0,'EXP_ARGUMENT')
    if x==0:return point(1)
    y=x;n=0
    while y>1:y/=2;n+=1
    upper=sum(((-y)**j/F(factorial(j)) for j in range(41)),F(0))
    lower=upper-y**41/F(factorial(41))
    need(lower>0,'EXP_POSITIVITY')
    v=rounded(Box(lower,upper))
    for _ in range(n):v=rounded(v*v)
    return v
def bump(lo,hi,tail):
    if tail:
        edge=F(1,64);s=1/(edge*(1-edge))
        need(s>=F(4096,63),'TAIL_DOMAIN')
        return Box(F(0),exponential_minus(s-4).hi)
    offset=Box(lo-F(3,2),hi-F(3,2))
    distance=point(F(1,4))-offset**2
    s=1/distance
    return Box(exponential_minus(s.hi-4).lo,exponential_minus(s.lo-4).hi)
def contains(outer,inner,label):
    need(outer.lo<=inner.lo<=inner.hi<=outer.hi,label)
def hull_integral(b):return Box(min(F(0),b.lo),max(F(0),b.hi))
def run(a):
    raw=load(a.initial_report);rows=raw['profile_cells']
    need(len(rows)==994,'CELL_COUNT')
    I=point(0);J=point(0);reportedI=point(0);reportedJ=point(0)
    rebuilt=[];cursor=F(1);tails=0
    for row in rows:
        lo,hi=map(F,row['r']);need(lo==cursor and hi>lo,'PARTITION');cursor=hi
        tail=row['tail'];tails+=int(tail)
        if tail:need((lo,hi) in [(F(1),F(65,64)),(F(127,64),F(2))],'TAIL_PLACEMENT')
        p={k:readbox(v) for k,v in row['profiles'].items()}
        b=bump(lo,hi,tail);width=ELL*(hi-lo)
        di=width*p['rho']**5*p['QH'];dj=width*p['rho']**3*b
        contains(readbox(row['integral_I']),di,'CELL_I_ENCLOSURE')
        contains(readbox(row['integral_J']),dj,'CELL_J_ENCLOSURE')
        I+=di;J+=dj;reportedI+=readbox(row['integral_I']);reportedJ+=readbox(row['integral_J'])
        rebuilt.append((row,p,b,di,dj))
    need(cursor==2 and tails==2,'SUPPORT_AND_TAILS')
    need(reportedI==readbox(raw['integrals']['I']),'RAW_I_REAGGREGATION')
    need(reportedJ==readbox(raw['integrals']['J']),'RAW_J_REAGGREGATION')
    contains(readbox(raw['integrals']['I']),I,'GLOBAL_I_ENCLOSURE')
    contains(readbox(raw['integrals']['J']),J,'GLOBAL_J_ENCLOSURE')
    need(I.hi<0 and J.lo>0,'GLOBAL_INTEGRAL_SIGNS')
    ratio=I/J;c=F(2,9)*ratio
    moment=F(2,3)*I
    if a.mutant=='wrong_moment':moment=I/3
    contains(readbox(raw['integrals']['coefficient_2I_over_9J']),c,'F_COEFFICIENT')
    contains(readbox(raw['integrals']['metric_moment_2I_over_3']),moment,'METRIC_MOMENT')
    KH=point(0);Kb=point(0);primitives=0
    def xi_value(p,b,kh,kb):
        f=c*b
        factor=F(-3,2) if a.mutant=='wrong_primitive_sign' else F(3,2)
        return factor*(p['fy']-p['a']*f)/p['p']**2-(kh-ratio*kb)/(p['rho']**4*p['p']**2)
    for row,p,b,di,dj in reversed(rebuilt):
        kh=KH+hull_integral(di);kb=Kb+hull_integral(dj)
        contains(p['KH'],kh,'CELL_KH_PREFIX')
        contains(p['Kb'],kb,'CELL_KB_PREFIX')
        contains(p['f'],c*b,'CELL_F')
        contains(p['xi'],xi_value(p,b,kh,kb),'CELL_XI')
        primitives+=1;KH+=di;Kb+=dj
    points=0
    for target in raw['targets']:
        r=F(target['r']);p={k:readbox(v) for k,v in target['profiles'].items()}
        kh=sum((di for row,_,_,di,_ in rebuilt if F(row['r'][0])>=r),point(0))
        kb=sum((dj for row,_,_,_,dj in rebuilt if F(row['r'][0])>=r),point(0))
        contains(p['KH'],kh,'TARGET_KH');contains(p['Kb'],kb,'TARGET_KB')
        contains(p['xi'],xi_value(p,bump(r,r,False),kh,kb),'TARGET_XI');points+=1
    return {'status':'PASS_INDEPENDENT_INITIAL_QUADRATURE_AND_PRIMITIVE','counts':{'cells':len(rows),'endpoint_tails':tails,'cell_integral_inclusions':2*len(rows),'cell_primitive_inclusions':primitives,'point_primitive_inclusions':points,'raw_global_integral_equalities':2},
        'independent_integrals':{k:[str(v.lo),str(v.hi)] for k,v in [('I',I),('J',J),('moment_2I_over_3',moment)]},
        'premises':['Retained raw background, QH and f_y boxes are shared validated input enclosures; they are not independently re-oracled.',
        'Bump values and flat endpoint-tail caps are re-enclosed independently with degree-40 alternating rational exponentials on a finer 2^160 grid.',
        'All 994 quadrature contributions, global moments, suffix primitives and both point primitives are independently reconstructed; exact unrounded intervals must lie inside retained outward bounds.']}
if __name__=='__main__':
    a=arguments(mutants=['wrong_moment','wrong_primitive_sign']);authenticate(a)
    try:finish(a,run(a))
    except ValueError as error:
        print(json.dumps({'status':'REJECT','stage':'SCIENCE','reason':str(error)},sort_keys=True))
        raise SystemExit(2)
