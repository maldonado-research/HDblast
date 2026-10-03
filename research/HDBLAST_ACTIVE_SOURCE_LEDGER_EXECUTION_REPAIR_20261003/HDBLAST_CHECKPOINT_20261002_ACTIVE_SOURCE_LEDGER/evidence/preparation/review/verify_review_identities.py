#!/usr/bin/env python3
"""Independent exact-symbol review; never imports sources or array readers."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import sympy as s


def audit():
    checks, rejected = [], []
    def zero(name, expression):
        if s.factor(s.expand(expression)) != 0:
            raise RuntimeError(name)
        checks.append(name)
    def nonzero(name, expression):
        if s.factor(s.expand(expression)) == 0:
            raise RuntimeError('Undetected mutation: '+name)
        rejected.append(name)

    L,k,e = s.symbols('L k epsilon', positive=True)
    X,Y,Z,T,f = s.symbols('X Y Z T real_force', real=True)
    h = s.symbols('h0:5', real=True)
    g = 4*L**2*h[0]-2*L*h[1]-h[2]
    def D(expression):
        return s.expand(L**2*s.diff(expression,L)+Z*s.diff(expression,X)
                        +T*s.diff(expression,Y)+(-2*k*T-e*f)*s.diff(expression,Z)
                        +2*k*Z*s.diff(expression,T)
                        +sum(h[j+1]*s.diff(expression,h[j]) for j in range(4)))

    R0 = (2*k**2+3*L**2)/(4*k)
    P0 = (2*k**2/3-L**2)/(4*k)
    Rm = ((2*k**2+3*L**2)*X/k-T-L*Z/k)/(2*e)
    Pm = ((2*k**2/3-L**2)*X/k-T-L*Z/k)/(2*e)
    CR = (L*h[1]+2*L**2*h[0])/(2*k)-4*h[0]*R0
    CP = (L*h[1]-2*L**2*h[0])/(2*k)-4*h[0]*P0
    R,P = Rm+CR,Pm+CP
    Fm = L*(Rm-3*Pm)
    FC = L*(CR-3*CP)-3*h[1]*(R0+P0)
    F = Fm+FC
    cR = X-T/(2*k)
    G = (3*L**2*X-L*Z)/(2*k*e)

    zero('real forcing preserves unrestricted c_R',D(cR))
    zero('direct mode density decomposes into constant plus G',Rm-k*cR/e-G)
    zero('mode primitive derivative includes real forcing work',D(G)-Fm-L*f/(2*k))
    zero('direct metric contact derivative cancels canonical g work',D(CR)-FC+L*g/(2*k))
    zero('source mismatch Ward identity',D(R)-F-L*(f-g)/(2*k))
    zero('consistent source obeys full bare Ward', (D(R)-F).subs(f,g))
    zero('bare finite-band background obeys Ward',D(R0)-L*R0+3*L*P0)
    zero('primitive derivative equals separately expanded full integrand', (D(G+CR)-F).subs(f,g))
    nonzero('baseline contact omission',3*h[1]*(R0+P0))
    nonzero('wrong source sign', (D(R)-F).subs(f,-g))
    nonzero('independent source jets change the target',L*(f-g)/(2*k))
    nonzero('Wronskian projection changes unrestricted density',k*cR/e)

    dX,dZ,dT = s.symbols('deltaX deltaZ deltaT', real=True)
    projection = ((2*k**2+3*L**2)*dX/k-dT-L*dZ/k)/(2*e)
    zero('density flow projection uses measured unrestricted defects',
         Rm.subs({X:X+dX,Z:Z+dZ,T:T+dT})-Rm-projection)
    zero('equal contact target cancels from endpoint density defect',
         R.subs({X:X+dX,Z:Z+dZ,T:T+dT})-R-projection)

    # Duhamel map, with arbitrary initial complex modes and a real M0. The
    # exp-moment is arbitrary; no physical source/phase quadrature occurs.
    co,si,mr,mi,m0 = s.symbols('cos sin Mexp_real Mexp_imag M0', real=True)
    # W_b = E W_a-epsilon Mexp; U_b=U_a+(E-1)W_a/(i2k)
    #       -epsilon(Mexp-M0)/(i2k).
    Ub_re = X+(si*Z+(co-1)*T)/(2*k)-e*mi/(2*k)
    Wb_im = si*Z+co*T-e*mi
    zero('Duhamel endpoint retains c_R without projection',Ub_re-Wb_im/(2*k)-cR)

    deltaR,S,I,Eflow = s.symbols('DeltaR Simpson Integral Eflow', real=True)
    zero('signed decomposition is algebraic not evidence',
         deltaR-S-(deltaR-I)-(I-S))
    zero('signed operator-flow closure is a separate reported diagnostic',
         (deltaR-I)-Eflow-(deltaR-I-Eflow))

    return {'status':'PASS_EXACT_SYMBOLIC_REVIEW','checks':checks,
            'check_count':len(checks),'detected_mutations':rejected,
            'mutation_count':len(rejected),'physical_evaluations':0,
            'physical_arrays_loaded':0,
            'scope':'Unrestricted real-forced modes, direct bare operators and arbitrary symbolic source jets only; renormalized contact inventory and quadrature still require separate review.'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():
        raise RuntimeError('Fresh review receipt required')
    result=audit()
    result.update(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  sympy_version=s.__version__,python_optimization=sys.flags.optimize)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('status','check_count','mutation_count','physical_evaluations')}))


if __name__=='__main__':
    main()
