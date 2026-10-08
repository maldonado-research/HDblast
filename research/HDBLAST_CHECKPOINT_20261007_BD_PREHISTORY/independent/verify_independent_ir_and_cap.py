#!/usr/bin/env python3
"""Exact independent algebra and fabricated controls; no physical callbacks/data."""
from __future__ import annotations
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
import sympy as S


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checks, controls = [], []

    def eq(name, value, target=0):
        residual = S.simplify(S.expand(value-target))
        if residual != 0:
            raise RuntimeError(name+': '+str(residual))
        checks.append(name)

    def reject(name, value):
        residual = S.simplify(S.expand(value))
        if residual == 0:
            raise RuntimeError('Missed negative control: '+name)
        controls.append({'name': name, 'nonzero_witness': str(residual)})

    k, tau, L, K, Pi = S.symbols('k tau L K Pi', positive=True)
    x, y, p, q, C, T, f, d = S.symbols('x y p q C T f d', real=True)
    I = S.I
    E = S.cos(2*k*tau)+I*S.sin(2*k*tau)
    Phi = (E-1)/(2*I*k)
    eq('phase derivative in lag', S.diff(E,tau),2*I*k*E)
    eq('drift derivative in lag', S.diff(Phi,tau),E)
    eq('drift continuous at zero momentum', S.limit(Phi,k,0),tau)
    eq('real drift equals imaginary phase over2k',S.re(Phi),S.im(E)/(2*k))
    eq('imaginary drift equals1−realphase over2k',S.im(Phi),(1-S.re(E))/(2*k))

    # W=-∫E r and U=-∫Phi r: each real quadrature atom has c=0.
    Uatom, Watom = -f*Phi, -f*E
    eq('real source atom canonical invariant',S.re(Uatom)-S.im(Watom)/(2*k))
    eq('real source atom phase real part',S.re(Watom/(2*I*k)),-f*S.sin(2*k*tau)/(2*k))
    eq('real source atom phase imaginary part',S.im(Watom/(2*I*k)),f*S.cos(2*k*tau)/(2*k))
    eq('source atom p finite limit',S.limit(S.re(Watom/(2*I*k)),k,0),-f*tau)
    eq('source atom imaginary residue',S.limit(k*S.im(Watom/(2*I*k)),k,0),f/2)
    reject('generic source atom is not bounded A at k0',S.limit(k*Watom/(2*I*k),k,0))

    # No unit-modulus premise in this scalar normalization invariant.
    arbitrary_phase=C+I*T
    arbitrary_drift=(arbitrary_phase-1)/(2*I*k)
    oldU, oldW=x+I*y,p+I*q
    newU=oldU+arbitrary_drift*oldW
    newW=arbitrary_phase*oldW
    eq('canonical invariant for arbitrary paired propagator',S.re(newU)-S.im(newW)/(2*k),x-q/(2*k))
    reject('normalization does not impose phase unitarity',S.expand(arbitrary_phase*S.conjugate(arbitrary_phase)-1))
    eq('wrong zero phase still preserves canonical invariant',(S.re(newU)-S.im(newW)/(2*k)-(x-q/(2*k))).subs({C:0,T:0}))
    reject('wrong zero phase fails exact propagation',arbitrary_phase.subs({C:0,T:0})-E)

    # Entire measure-weighted state stress kernels for a real source atom.
    pa,qa=S.re(Watom/(2*I*k)),S.im(Watom/(2*I*k))
    R=3*L**2*pa/(2*k)+L*qa
    P=-(2*k/S.Integer(3)+L**2/(2*k))*pa+L*qa
    RK=k**2*R/f
    PK=k**2*P/f
    eq('weighted density kernel',RK,-3*L**2*S.sin(2*k*tau)/4+L*k*S.cos(2*k*tau)/2)
    eq('weighted pressure kernel',PK,(k**2/S.Integer(3)+L**2/4)*S.sin(2*k*tau)+L*k*S.cos(2*k*tau)/2)
    eq('weighted density vanishes at k0',S.limit(RK,k,0))
    eq('weighted pressure vanishes at k0',S.limit(PK,k,0))
    eq('weighted density linear coefficient',S.limit(RK/k,k,0),(L-3*L**2*tau)/2)
    eq('weighted pressure linear coefficient',S.limit(PK/k,k,0),(L+L**2*tau)/2)
    D=lambda expression:S.diff(expression,tau)+L**2*S.diff(expression,L)
    eq('direct homogeneous Ward source-atom identity',D(R),L*(R-3*P))
    eq('direct pressure primitive source-atom identity',D(-R/(3*L)),P)
    reject('omitting real prehistory gives nonzero direct density',RK)
    reject('omitting real prehistory gives nonzero direct pressure',PK)

    M0,M1=S.symbols('M0 M1',nonnegative=True)
    majorR=(3*L**2*M1+L*M0)/(2*k)
    majorP=2*k*M1/3+(L**2*M1+L*M0)/(2*k)
    eq('integrated correlated density majorant',S.integrate(k**2*majorR/(2*Pi**2),(k,0,K)),K**2*(3*L**2*M1+L*M0)/(8*Pi**2))
    eq('integrated correlated pressure majorant',S.integrate(k**2*majorP/(2*Pi**2),(k,0,K)),M1*K**4/(12*Pi**2)+K**2*(L**2*M1+L*M0)/(8*Pi**2))

    # Twice-integrated endpoint formula: no evaluation of either physical source.
    s=S.symbols('s',real=True,nonzero=True)
    h=S.Function('h')(s)
    lag=S.symbols('a',real=True)-s
    Es=S.exp(2*I*k*lag)
    Ls=-1/s
    g=4*Ls**2*h-2*Ls*S.diff(h,s)-S.diff(h,s,2)
    boundary=Es*(-S.diff(h,s)-2*Ls*h-2*I*k*h)
    interior=Es*(4*k**2-4*I*k*Ls+6*Ls**2)*h
    eq('twice integrated source identity',S.diff(boundary,s)+interior,Es*g)
    Phis=(Es-1)/(2*I*k)
    boundary_U=Phis*(-S.diff(h,s)-2*Ls*h)-Es*h
    interior_U=(6*Ls**2*Phis-(2*I*k+2*Ls)*Es)*h
    eq('twice integrated nested drift source identity',S.diff(boundary_U,s)+interior_U,Phis*g)
    reject('wrong imaginary interior sign',S.diff(boundary,s)+Es*(4*k**2+4*I*k*Ls+6*Ls**2)*h-Es*g)
    reject('omitted boundary derivative',interior-Es*g)

    z=S.symbols('z',real=True)
    B=S.exp(1-1/(1-z**2))
    Dz=1-z**2
    eq('bump first derivative algebra',S.diff(B,z)/B,-2*z/Dz**2)
    eq('bump second derivative numerator',S.diff(B,z,2)/B,(6*z**4-2)/Dz**4)
    delta=F(1,128); Ddelta=delta*(2-delta); lc=1/(5-delta)
    exponent=F(1)-1/Ddelta
    if not (exponent < -63 and F(65,24)>F(8,3)):
        raise RuntimeError('Rational exponential cap argument failed')
    checks.append('rational e lower bound and endpoint exponent')
    H=F(3,8)**63
    pd=2*(1-delta)/Ddelta**2
    cap_rows=[]
    for cutoff in (64,128,256):
        Ccap=H*(pd+1+2*lc+2*cutoff+delta*(4*cutoff**2+4*cutoff*lc+6*lc**2))
        rb=Ccap*(F(cutoff**2,252)+F(cutoff,294))
        pb=Ccap*(F(cutoff**3,162)+F(5*cutoff,1764))
        cap_rows.append({'K':cutoff,'common_signed_source_Ccap':str(Ccap),'density_cap_upper':str(rb),'pressure_cap_upper':str(pb),'pressure_cap_less_than_2e_minus18':pb<F(2,10**18)})
    if not all(row['pressure_cap_less_than_2e_minus18'] for row in cap_rows):
        raise RuntimeError('Claimed cap threshold false')
    checks.append('three exact rational cap pressure thresholds')
    cap=S.symbols('cap',positive=True)
    density_triangle=cap*(L/(2*k)+3*L**2/(4*k**2))
    pressure_triangle=cap*(S.Rational(1,3)+5*L**2/(8*k**2))
    eq('candidate cap density coefficient',S.integrate(k**2*density_triangle/(2*Pi**2),(k,0,K)).subs({L:S.Rational(2,7),Pi:3})/cap,K**2/252+K/294)
    eq('candidate cap pressure coefficient',S.integrate(k**2*pressure_triangle/(2*Pi**2),(k,0,K)).subs({L:S.Rational(2,7),Pi:3})/cap,K**3/162+5*K/1764)

    # Nonlinear normalized family confirms normalization leaves arbitrary tangent beta.
    eps,beta=S.symbols('eps beta',real=True)
    alpha=S.sqrt(1+eps**2*beta**2)
    eq('exact canonical Bogoliubov identity',alpha**2-eps**2*beta**2,1)
    eq('arbitrary beta is a normalized tangent',S.diff(alpha,eps).subs(eps,0))
    reject('normalization does not set beta to zero',S.diff(eps*beta,eps).subs(eps,0))

    result={'status':'PASS','checks':checks,'check_count':len(checks),'negative_controls':controls,'negative_control_count':len(controls),'cap_rows':cap_rows,'scope':'Exact symbolic identities and rational analytic cap bounds only. No physical callback, no saved quantum array decoding, no runtime physical mode or source evaluation. Mathematical proof accompanies these checks.','actual_incoming_state':'NOT_ENCLOSED','full_twelve_case_certificate':'UNRESOLVED','metric_calibration':'FAIL','higher_dimensional_Big_Bang_origin':'NOT_ESTABLISHED'}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),'negative_controls':len(controls),'output':str(args.output)}))


if __name__=='__main__':
    main()
