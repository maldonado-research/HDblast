#!/usr/bin/env python3
"""Exact symbolic identities only; no sources, stored modes or responses loaded."""
from __future__ import annotations
import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import sys
import sympy as s


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def audit():
    checks, mutations = [], []
    def zero(name, expression):
        require(s.simplify(s.expand(expression)) == 0, name)
        checks.append(name)
    def reject(name, expression):
        require(s.simplify(s.expand(expression)) != 0, 'Mutation survived: '+name)
        mutations.append(name)
    L,k,e = s.symbols('L k epsilon', positive=True)
    X,Y,Z,T = s.symbols('X Y Z T', real=True)
    h = s.symbols('h0:5', real=True)
    g = 4*L**2*h[0]-2*L*h[1]-h[2]
    def D(expression):
        return s.expand(L**2*s.diff(expression,L)+Z*s.diff(expression,X)+T*s.diff(expression,Y)
            +(-2*k*T-e*g)*s.diff(expression,Z)+2*k*Z*s.diff(expression,T)
            +sum(h[j+1]*s.diff(expression,h[j]) for j in range(4)))
    C = 2*L**2
    for n in range(7):
        zero('correct C derivative '+str(n), C-2*factorial(n+1)*L**(n+2))
        C = L**2*s.diff(C,L)
    reject('wrong reader factorial at n1', factorial(3)*L**3-2*factorial(2)*L**3)
    R0 = (2*k**2+3*L**2)/(4*k)
    P0 = (2*k**2/3-L**2)/(4*k)
    Rmode = ((2*k**2+3*L**2)*X/k-T-L*Z/k)/(2*e)
    Pmode = ((2*k**2/3-L**2)*X/k-T-L*Z/k)/(2*e)
    R = Rmode+(L*h[1]+2*L**2*h[0])/(2*k)-4*h[0]*R0
    P = Pmode+(L*h[1]-2*L**2*h[0])/(2*k)-4*h[0]*P0
    zero('bare background Ward',D(R0)-L*R0+3*L*P0)
    zero('full bare metric Ward without Wronskian projection',D(R)-L*R+3*L*P+3*h[1]*(R0+P0))
    zero('unprojected Wronskian is constant',D(2*X-T/k))
    reject('drop fixed K baseline term',3*h[1]*(R0+P0))
    reject('wrong force sign',2*g*s.diff(R,Z)*e)
    dt,theta = s.symbols('Delta theta', positive=True)
    f = s.symbols('F0:5')
    rr = s.symbols('R0 R2 R4')
    cell0 = rr[1]-rr[0]-dt*(f[0]+4*f[1]+f[2])/3
    cell1 = rr[2]-rr[1]-dt*(f[2]+4*f[3]+f[4])/3
    total = rr[2]-rr[0]-dt*(f[0]+4*f[1]+2*f[2]+4*f[3]+f[4])/3
    zero('signed discrete Ward telescopes exactly',cell0+cell1-total)
    harmonic = 2*s.I*(s.sin(theta)-theta*(2+s.cos(theta))/3)
    explicit = (s.exp(s.I*theta)-s.exp(-s.I*theta)
                -s.I*theta*(s.exp(-s.I*theta)+4+s.exp(s.I*theta))/3)
    zero('exact Simpson harmonic defect',s.expand_complex(explicit)-harmonic)
    zero('Simpson leading harmonic error',s.limit(harmonic/theta**5,theta,0)+s.I/90)
    transfer = theta*(2+s.cos(theta))/(3*s.sin(theta))
    zero('Simpson transfer leading fourth order',s.limit((transfer-1)/theta**4,theta,0)-s.Rational(1,180))
    zero('Nyquist Simpson nonzero alias',harmonic.subs(theta,s.pi)+2*s.I*s.pi/3)
    reject('assume Nyquist pair integrates derivative exactly',harmonic.subs(theta,s.pi))
    # Arbitrary Re(c), including any retained normalization/Wronskian defect.
    c,dr,di,co,si = s.symbols('c_R d_R d_I cos_phase sin_phase', real=True)
    x = c+dr*co-di*si
    z = -2*k*(di*co+dr*si)
    t = 2*k*(dr*co-di*si)
    phase_real,phase_imag = dr*co-di*si,di*co+dr*si
    Rfree = ((2*k**2+3*L**2)*c+3*L**2*phase_real+2*k*L*phase_imag)/(2*k*e)
    Pfree = ((2*k**2/3-L**2)*c-(4*k**2/3+L**2)*phase_real+2*k*L*phase_imag)/(2*k*e)
    zero('free density from unrestricted complex amplitudes',Rmode.subs({X:x,Z:z,T:t})-Rfree)
    zero('free pressure from unrestricted complex amplitudes',Pmode.subs({X:x,Z:z,T:t})-Pfree)
    Ffree = 3*L**3*c/(k*e)+((2*k**2*L+3*L**3)*phase_real-2*k*L**2*phase_imag)/(k*e)
    zero('free ledger integrand independently expanded',L*(Rfree-3*Pfree)-Ffree)
    Dfree = lambda q:s.expand(L**2*s.diff(q,L)-2*k*si*s.diff(q,co)+2*k*co*s.diff(q,si))
    zero('free analytic ledger differentiates to direct integrand',Dfree(Rfree)-Ffree)
    J1,B1,B2,d = s.symbols('J1 boundary_LE boundary_L2E d')
    J2 = B1-2*s.I*k*J1
    J3 = (B2-2*s.I*k*J2)/2
    exact = d*(2*k**2*J1+3*J3)+(2*s.I*k*d)*J2
    closed = d*(3*B2-2*s.I*k*B1)/2
    zero('oscillatory moment recurrence cancels master integral',exact-closed)
    reject('drop retained real constant amplitude',3*L**3*c/(k*e))
    dX,dZ,dT = s.symbols('delta_X delta_Z delta_T', real=True)
    flow_projection = ((2*k**2+3*L**2)*dX/k-dT-L*dZ/k)/(2*e)
    zero('exact linear flow density projection',Rmode.subs({X:X+dX,Z:Z+dZ,T:T+dT})-Rmode-flow_projection)
    endpoint,simpson,exact_ledger = s.symbols('endpoint Simpson exact_ledger')
    zero('exact ledger attribution identity',(endpoint-simpson)-((endpoint-exact_ledger)+(exact_ledger-simpson)))
    reject('reverse quadrature attribution sign',2*(exact_ledger-simpson))
    return {'status':'PASS_PURE_MATHEMATICS','identities':checks,'identity_count':len(checks),
            'mutations':mutations,'mutation_count':len(mutations),'physical_evaluations':0,
            'saved_scientific_values_loaded':0,'scope':'Generic action/ODE and discrete quadrature identities only.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(not args.output.exists(),'Fresh pure proof receipt required')
    result=audit()
    result.update(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  sympy_version=s.__version__,python_optimization=sys.flags.optimize)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','identity_count','mutation_count','physical_evaluations')}))


if __name__=='__main__':
    main()
