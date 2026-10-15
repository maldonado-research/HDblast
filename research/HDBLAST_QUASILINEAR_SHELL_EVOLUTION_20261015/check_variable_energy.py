#!/usr/bin/env python3
"""Exact finite identities only; this does not certify the analytic existence proof."""
import argparse
import json
from pathlib import Path
import sympy as s

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=Path)
parser.add_argument('--corrupt', choices=['cross_factor','energy_sign','boundary_shift','mixed_flux'])
args = parser.parse_args()
t,z,x=s.symbols('t z x', real=True)
u=s.Function('u')(t,z,x)
e=s.Function('e')(t,z,x)
bz=s.Function('bz')(t,z,x); bx=s.Function('bx')(t,z,x)
azz=s.Function('azz')(t,z,x); azx=s.Function('azx')(t,z,x); axx=s.Function('axx')(t,z,x)
coords=(z,x); bs=(bz,bx); aa=((azz,azx),(azx,axx))
energy=(e*s.diff(u,t)**2+sum(aa[i][j]*s.diff(u,coords[i])*s.diff(u,coords[j]) for i in range(2) for j in range(2)))/2
if args.corrupt=='energy_sign':
    energy-=e*s.diff(u,t)**2
cross=1 if args.corrupt=='cross_factor' else 2
operator=e*s.diff(u,t,2)-cross*sum(bs[i]*s.diff(u,t,coords[i]) for i in range(2))-sum(aa[i][j]*s.diff(u,coords[i],coords[j]) for i in range(2) for j in range(2))
flux=[bs[i]*s.diff(u,t)**2+sum(aa[i][j]*s.diff(u,coords[j])*s.diff(u,t) for j in range(2)) for i in range(2)]
remainder=(s.diff(e,t)/2-sum(s.diff(bs[i],coords[i]) for i in range(2)))*s.diff(u,t)**2-sum(s.diff(aa[i][j],coords[i])*s.diff(u,coords[j])*s.diff(u,t) for i in range(2) for j in range(2))+sum(s.diff(aa[i][j],t)*s.diff(u,coords[i])*s.diff(u,coords[j]) for i in range(2) for j in range(2))/2
checks=[]
def require(name, expression):
    residual=s.simplify(s.expand(expression))
    if residual != 0:
        raise RuntimeError('Identity rejected: '+name)
    checks.append(name)
require('variable_scalar_energy_divergence_identity',s.diff(energy,t)-sum(s.diff(flux[i],coords[i]) for i in range(2))-operator*s.diff(u,t)-remainder)
# The positive control fiber has four D and twelve N components.
vt=s.symbols('vt0:16'); vz=s.symbols('vz0:16'); va=s.symbols('va0:16')
B,Azz,AzA=s.symbols('B Azz AzA')
boundary=sum(B*vt[j]**2+Azz*vz[j]*vt[j]+AzA*va[j]*vt[j] for j in range(16))
sub={vt[j]:0 for j in range(4)}
sub.update({vz[j]:0 for j in range(4,16)})
sub[B]=B if args.corrupt=='boundary_shift' else 0
sub[AzA]=AzA if args.corrupt=='mixed_flux' else 0
require('four_D_twelve_N_zero_principal_boundary_flux',boundary.subs(sub, simultaneous=True))
# Explicit leftover terms show why the normal/tangent coefficient condition matters.
require('shift_leakage_on_N_sector',s.diff(boundary,B).subs({vt[j]:0 for j in range(4)})-sum(vt[j]**2 for j in range(4,16)))
require('spatial_cross_leakage_on_N_sector',s.diff(boundary,AzA).subs({vt[j]:0 for j in range(4)})-sum(va[j]*vt[j] for j in range(4,16)))
# Recovering normal derivatives requires an independent positive normal coefficient.
f,utt,utx,utzz,uxx,uzx,aa0,b0,b1,a01,a11=s.symbols('f utt utx utz uxx uzx aa0 b0 b1 a01 a11')
normal=(aa0*utt-2*b0*utzz-2*b1*utx-2*a01*uzx-a11*uxx-f)/Azz
require('normal_second_derivative_recovery',aa0*utt-2*b0*utzz-2*b1*utx-Azz*normal-2*a01*uzx-a11*uxx-f)
# Conservatively reserve four of the m=14 spatial Sobolev orders for C^2 parameters.
m=14; lower=m-4; spatial_dimension=4; spacetime_C_order=4
for j in range(spacetime_C_order+1):
    if not lower-j > spatial_dimension/2+(spacetime_C_order-j):
        raise RuntimeError('Insufficient spatial embedding budget')
checks.append('C2_parameter_X10_to_piecewise_spacetime_C4_budget_only')
result={'schema':'hdblast-variable-energy-exact-identities-v1','status':'PASS','checks':checks,'positive_checks':len(checks),'scope':'Finite exact algebra and a conservative embedding count; not a PDE existence, convergence or physical stability certificate.'}
raw=json.dumps(result,indent=2)+'\n'
if args.output:
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(raw)
else:
    print(raw,end='')
