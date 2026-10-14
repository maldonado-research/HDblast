#!/usr/bin/env python3
"""Direct 5D Ricci linearization for the actual l=2 S3 Weyl tensor harmonic."""
import argparse
import json
from pathlib import Path
import sympy as s


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    checks={}

    def clean(v):
        v=s.factor(s.trigsimp(s.simplify(s.expand(v))))
        if v != 0:
            v=s.simplify(s.trigsimp(s.expand_trig(v).rewrite(s.exp)))
        return v

    def zero(name,v):
        v=clean(v)
        if v != 0:
            raise RuntimeError(f'{name}: {v}')
        checks[name]='passed'

    y,t,th,ph,ps=coords=s.symbols('y t th ph ps',real=True)
    rho=s.Function('rho')(y)
    H=s.Function('H')(y,t)
    a=s.diff(rho,y)/rho
    c=s.cosh(t)
    sn,cs=s.sin(th),s.cos(th)
    # n=(sin th cos ph,sin th sin ph,cos th cos ps,cos th sin ps).
    # T_ab=C_ikjl n^k n^l partial_a n^i partial_b n^j, same C as the source note.
    TT=[s.Rational(1,6),sn**2*(s.Rational(1,6)-sn**2/2),
        cs**2*(s.Rational(1,6)-cs**2/2)]
    g=[s.Integer(1),-rho**2,rho**2*c**2,rho**2*c**2*sn**2,rho**2*c**2*cs**2]
    h=[s.Integer(0),s.Integer(0)]+[rho**2*c**2*H*v for v in TT]
    G,dG={},{}
    for A in range(5):
        for B in range(5):
            for C in range(5):
                dg=((s.diff(g[A],coords[B]) if A==C else 0)
                    +(s.diff(g[A],coords[C]) if A==B else 0)
                    -(s.diff(g[B],coords[A]) if B==C else 0))
                dh=((s.diff(h[A],coords[B]) if A==C else 0)
                    +(s.diff(h[A],coords[C]) if A==B else 0)
                    -(s.diff(h[B],coords[A]) if B==C else 0))
                G[A,B,C]=s.simplify(dg/(2*g[A]))
                dG[A,B,C]=s.simplify(dh/(2*g[A])-h[A]*G[A,B,C]/g[A])

    def dRic(A,B):
        v=sum(s.diff(dG[C,A,B],coords[C])-s.diff(dG[C,A,C],coords[B]) for C in range(5))
        for C in range(5):
            for D in range(5):
                v+=(dG[C,C,D]*G[D,A,B]+G[C,C,D]*dG[D,A,B]
                    -dG[C,B,D]*G[D,A,C]-G[C,B,D]*dG[D,A,C])
        return v

    Ufactor=3/rho**2-s.diff(rho,y,2)/rho-3*a**2
    Lop=s.diff(H,y,2)+4*a*s.diff(H,y)-(s.diff(H,t,2)+3*s.tanh(t)*s.diff(H,t)+8*H/c**2)/rho**2
    zero('S3_trace',TT[0]+TT[1]/sn**2+TT[2]/cs**2)
    zero('S3_divergence_theta',s.diff(TT[0],th)+s.cot(th)*(TT[0]-TT[1]/sn**2)-s.tan(th)*(TT[0]-TT[2]/cs**2))
    for A in range(5):
        for B in range(A,5):
            target=-rho**2*c**2*TT[A-2]*Lop/2 if A==B and A>=2 else 0
            dE=dRic(A,B)-(Ufactor*h[A] if A==B else 0)
            zero(f'linear_Einstein_TT_{A}{B}',dE-target)
    result={'status':'passed','checks':checks,'tensor_harmonic':'Weyl(B wedge B), B=diag(1,-1,0,0)',
            'operator':'H_yy+4 a H_y-rho^-2(H_tt+3 tanh(t) H_t+8 sech(t)^2 H)',
            'scope':'Exact direct Christoffel variation, all 15 independent residual components. Only the tangential background Einstein equation is substituted.'}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
