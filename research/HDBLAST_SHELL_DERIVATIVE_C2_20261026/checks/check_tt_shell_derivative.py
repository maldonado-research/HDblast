#!/usr/bin/env python3
"""Bounds registered TT affine shell forcing and its x derivative.
Uses only an authenticated inherited first-field report. It does not solve an evolution.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def jet_add(a,b): return a[0]+b[0],a[1]+b[1]
def jet_scale(a,k): return a[0]*k,a[1]*k
def jet_mul(a,b): return a[0]*b[0],a[1]*b[0]+a[0]*b[1]

# Test the actual physical-to-scaled derivative and both product terms,
# independent of the bound computation and with all factors nonzero.
ell,rho,a,sigma=F(2,5),F(7,3),F(11,4),F(-13,5)
J,Jt,v,vt,d,dt=F(3,2),F(17,6),F(-4,7),F(19,8),F(23,9),F(-29,10)
# Coordinate jets contain x derivatives ell times physical-time derivatives.
j=J,ell*Jt; z=v,ell*vt; f=d,ell*dt
actual=jet_scale(jet_mul(jet_add(jet_scale(j,a),jet_scale(z,sigma/3)),f),-F(8,3)*ell*rho)
expected=-F(8,3)*ell**2*rho*((a*Jt+sigma*vt/3)*d+(a*J+sigma*v/3)*dt)
need(actual[1]==expected,'product rule and dx conversion')
wrong_drop_dt=-F(8,3)*ell**2*rho*((a*Jt+sigma*vt/3)*d)
wrong_scaling=-F(8,3)*ell*rho*((a*Jt+sigma*vt/3)*d+(a*J+sigma*v/3)*dt)
wrong_thirds=-F(8,3)*ell**2*rho*((a*Jt+sigma*vt)*d+(a*J+sigma*v)*dt)
need(all(actual[1]!=v for v in (wrong_drop_dt,wrong_scaling,wrong_thirds)), 'wrong-formula controls')

parser=argparse.ArgumentParser()
parser.add_argument('--repo-root',type=Path,required=True)
args=parser.parse_args()
path=args.repo_root/'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/shell-producer/COMMON_HARMONIC_SHELL_NORMAL.json'
raw=path.read_bytes()
pin='559f7587feca262cc68990f2ee79a3aeea22df6955e5aafb698f3326fdacd94a'
need(hashlib.sha256(raw).hexdigest()==pin,'inherited report hash')
report=json.loads(raw)
need(report['status']=='PASS_CONTINUOUS_FULL_COMMON_FIRST_HARMONIC_PHYSICAL_SHELL_FORCING','report status')
need(report['domain']['scaled_shell']=='0<=x<=3,r=0; entire unit S3','report domain')
need(report['scope']['actual_registered_first_pair'] is True and report['scope']['continuous_values_and_first_tangential_covariant_derivatives'] is True,'scope')
B=report['bounds']; A=B['actual_field_amplitude_caps']; P=B['parameters']
get=lambda k:F(A[k]['exact_upper'])
par=lambda k:F(P[k])
ell=par('ell'); rho_lo,rho_hi=map(F,P['rho_b_interval'])
a,sigma=par('a_abs'),par('sigma_prime_abs')
need(ell==F(1,10**6) and 0<rho_lo<=rho_hi and a>=0 and sigma>=0,'parameters')
J,Jt,v,vt,d,dt=map(get,['J','J_t','v','v_t','d','d_t'])
need(all(x>=0 for x in (J,Jt,v,vt,d,dt)),'amplitude caps')
beta=a*J+sigma*v/3
betat=a*Jt+sigma*vt/3
value=F(8,3)*ell*rho_hi*beta*d
slope=F(8,3)*ell**2*rho_hi*(betat*d+beta*dt)
value_cap=F(5,10**21)
slope_cap=F(301,10**21)
need(value<value_cap,'value claimed cap')
need(slope<slope_cap,'x derivative claimed cap')
print(json.dumps({'status':'PASS','inherited_report_sha256':pin,'scaled_shell_domain':'0<=x<=3,r=0','same_bound_for':'Tplus and Tminus','exact_positive_bound_G':str(value),'exact_positive_bound_G_x':str(slope),'strict_readable_bound_G':'5/10^21','strict_readable_bound_G_x':'301/10^21','formula':'G_T=-(8/3)*ell*rho*(a*J+sigma_prime*v/3)*d; d/dx=ell*d/dt','controls':['omitted d_t term','omitted physical/scaled ell factor','omitted scalar factor 1/3'],'scope':'conditional inherited first-field caps; affine shell source only, no all-row residual or nonlinear solution'},sort_keys=True))
