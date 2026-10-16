#!/usr/bin/env python3
"""Actual common-gauge one-copy physical shell affine Hessian values.

Uses the authenticated on-shell first physical rows, only at r=0.  These
rows eliminate normal derivatives in the affine source; they are never
differentiated normally.  The fixed tension has sigma'=-2p,
sigma''=4(1+eta_background), sigma'''=4.
"""
import sys
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import actual_first_jets as j

SCALAR=(F(1,6),F(1,2),F(1))
VECTOR=(F(1,4),F(1,2))
YYK=(F(1,6),F(1,2),F(1),0,0,0,0)
YHESS=(F(-4,9),F(-4,3),F(-8,3),F(-1,4),F(1,6),F(2,3),F(2,3))

def raw_fields(H,D,E,L,M,s,a,p,sigma2,mutation=None):
 sigma1=-2*p
 A=a*H**2/4-sigma1*H*s/6-sigma2*s**2/6
 B=a*H+sigma1*s/3
 P=p*H**2/4+sigma2*H*s/2+2*s**2
 if mutation=='omit_normal_normalization':
  # Difference between the true unit-normal affine Hessian and one that
  # replaces N^-1 by its background value before variation.
  A-=F(3,4)*a*H**2;P-=F(3,4)*p*H**2
 return {'00':-A-B*D,'0a':-B*E,'ab_trace':A-B*L,'ab_hessian':-B*M,'scalar':P}

def projected(raw):
 return {'00':tuple(v*raw['00'] for v in SCALAR),
         '0a':tuple(v*raw['0a'] for v in VECTOR),
         'ab':tuple(v*raw['ab_trace']+w*raw['ab_hessian'] for v,w in zip(YYK,YHESS)),
         'scalar':tuple(v*raw['scalar'] for v in SCALAR)}

def actual_shell_sources(x,mutation=None):
 m,s,b=j.actual_jet(x,j.I(0));data,_,_=j.inputs()
 eta0=j.I(*data['shell']['boxes']['eta_b']['exact']);sigma2=4*(1+eta0)
 raw=raw_fields(m['A'],m['D'],m['E'],m['L'],m['M'],s['s'],b['a'],b['p'],sigma2,mutation)
 Q=projected(raw)
 return {'one_copy_outward_hessian':Q,
         'two_copy_outward_hessian':{key:tuple(2*v for v in row) for key,row in Q.items()},
         'two_copy_inward_hessian':{key:tuple(-2*v for v in row) for key,row in Q.items()},
         'coordinate_physical_forcing':{key:tuple((-1 if key=='scalar' else -2)*v for v in row) for key,row in Q.items()}}

def algebra_check(mutation=None):
 H,D,E,L,M,s,a,p,sigma2=map(F,(2,3,5,7,11,13,17,19,23));sigma1=-2*p
 raw=raw_fields(H,D,E,L,M,s,a,p,sigma2,mutation)
 # Independent pre-elimination affine unit-normal formulas.  m_y below
 # differentiates COORDINATE tangential components in the parallel frame
 # and therefore includes the explicit 2 a m warp-factor term.
 expected={}
 for key,m,metric in (('00',D,-1),('0a',E,0),('ab_trace',L,1),('ab_hessian',M,0)):
  my=2*a*m+(a*H+sigma1*s/3)*metric
  expected[key]=F(3,4)*a*H**2*metric-H*my/2-(sigma2*s*s*metric+2*sigma1*s*m)/6
 sy=p*H/2-sigma2*s/2
 expected['scalar']=F(3,4)*p*H**2-H*sy+2*s*s
 j.need(raw==expected,'physical shell unit-normal affine Hessian and first-row elimination')
 return {'status':'PASS_EXACT_PHYSICAL_SHELL_SOURCE_ALGEBRA','independent_component_families':5,
         'normalization':'one-copy outward affine Hessian; two-copy outward=2Q; coordinate physical metric forcing=-2QK, scalar forcing=-Qp'}
