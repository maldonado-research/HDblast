#!/usr/bin/env python3
"""Independent rational affine-jet boundary audit; no scientific imports."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path


class Failure(Exception):
    pass


class D:
    """A rational scalar and its four first intrinsic coordinate derivatives."""
    def __init__(self, v=0, g=None):
        self.v = F(v)
        self.g = tuple(F(x) for x in (g or (0,)*4))
    def __add__(self, b):
        b = dual(b)
        return D(self.v+b.v, [x+y for x,y in zip(self.g,b.g)])
    __radd__ = __add__
    def __neg__(self):
        return D(-self.v, [-x for x in self.g])
    def __sub__(self,b):
        return self+-dual(b)
    def __rsub__(self,b):
        return dual(b)+-self
    def __mul__(self,b):
        b=dual(b)
        return D(self.v*b.v,[x*b.v+self.v*y for x,y in zip(self.g,b.g)])
    __rmul__=__mul__
    def __truediv__(self,b):
        b=dual(b)
        return self*D(1/b.v,[-x/b.v**2 for x in b.g])
    def __eq__(self,b):
        b=dual(b)
        return self.v==b.v and self.g==b.g


def dual(x):
    return x if isinstance(x,D) else D(x)


class P:
    """Taylor COEFFICIENTS (not derivatives), truncated after epsilon squared."""
    def __init__(self,a=0,b=0,c=0):
        self.c=[dual(a),dual(b),dual(c)]
    def __add__(self,b):
        b=poly(b)
        return P(*(x+y for x,y in zip(self.c,b.c)))
    __radd__=__add__
    def __neg__(self):
        return P(*(-x for x in self.c))
    def __sub__(self,b):
        return self+-poly(b)
    def __rsub__(self,b):
        return poly(b)+-self
    def __mul__(self,b):
        b=poly(b)
        return P(*(sum((self.c[j]*b.c[k-j] for j in range(k+1)),D()) for k in range(3)))
    __rmul__=__mul__
    def inv(self):
        a,b,c=self.c
        return P(D(1)/a,-b/(a*a),b*b/(a*a*a)-c/(a*a))
    def __truediv__(self,b):
        return self*poly(b).inv()
    def sqrt1(self):
        a,b,c=self.c
        if a != 1:
            raise Failure('sqrt base is not one')
        return P(1,b/2,c/2-b*b/8)
    def deriv(self,a):
        return P(*(D(x.g[a]) for x in self.c))


def poly(x):
    return x if isinstance(x,P) else P(x)


def matmul(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(4)),P()) for j in range(4)] for i in range(4)]


def trace(a):
    return sum((a[i][i] for i in range(4)),P())


def invseries(g, diag):
    """Independent recurrence for inverse-matrix Taylor coefficients."""
    out=[[[D(int(i==j))/diag[i] for j in range(4)] for i in range(4)]]
    for n in (1,2):
        out.append([[-sum((g[i][k].c[j]*out[n-j][k][l] for j in range(1,n+1) for k in range(4)),D())/diag[i]
                     for l in range(4)] for i in range(4)])
    return [[P(*(out[n][i][j] for n in range(3))) for j in range(4)] for i in range(4)]


def check(injection=None):
    count=0
    def equal(actual, expected, label):
        nonlocal count
        if actual != expected:
            raise Failure(label)
        count+=1
    # Independently construct rational tangent frames by ambient reflection.
    for xyz in ((F(1),F(0),F(0)),(F(1,2),F(1,3),F(1,5)),
                (F(2,3),F(1,7),F(-2,5)),(F(1),F(1),F(1))):
        norm=sum(x*x for x in xyz)
        n=[2*x/(1+norm) for x in xyz]+[(1-norm)/(1+norm)]
        rv=[F(int(i==3))-n[i] for i in range(4)]
        rvnorm=sum(x*x for x in rv)
        e=[[F(int(i==j))-2*rv[i]*rv[j]/rvnorm for i in range(4)] for j in range(3)]
        Y=n[0]*n[0]-n[1]*n[1]; S=n[0]*n[0]+n[1]*n[1]
        grad=[2*(v[0]*n[0]-v[1]*n[1]) for v in e]
        hess=[[2*(e[i][0]*e[j][0]-e[i][1]*e[j][1])-2*Y*int(i==j)
               for j in range(3)] for i in range(3)]
        gg=sum(x*x for x in grad); hh=sum(x*x for row in hess for x in row)
        equal(gg,4*(S-Y*Y),'actual harmonic gradient contraction')
        equal(hh,8-8*S+24*Y*Y,'actual harmonic Hessian contraction')
        U0=F(7,13); L=F(-5,17); J=F(11,19); dd=F(3,23)
        mf=[[F(0) for j in range(4)] for i in range(4)]
        mf[0][0]=-2*U0*Y
        for i in range(3):
            mf[0][i+1]=mf[i+1][0]=J*grad[i]
            for j in range(3):
                mf[i+1][j+1]=2*(L*Y*int(i==j)+dd*hess[i][j])
        lorentz=(-1,1,1,1)
        square=sum(lorentz[i]*lorentz[j]*mf[i][j]*mf[j][i]
                   for i in range(4) for j in range(4))
        equal(square,(4*U0*U0+12*L*L-64*L*dd)*Y*Y+4*dd*dd*hh-2*J*J*gg,
              'actual Lorentzian full tangential square')
    for seed in (1,3,7):
        a=F(2+seed,11); p=F(-3-seed,13)
        sp=-2*p; spp=F(4*(7+seed),9); sppp=F(4)
        if injection=='scalar-sign':
            sp=2*p
        diag=[F(-4),F(9),F(16),F(25)]
        base=[[diag[i] if i==j else F(0) for j in range(4)] for i in range(4)]
        def item(i,j,offset):
            lo,hi=sorted((i,j))
            return D(F((lo+1)*(hi+2)+seed+offset,17+offset),
                     [F((lo+2)*(hi+3)+(q+1)*seed+offset,31+q+offset) for q in range(4)])
        m=[[item(i,j,0) for j in range(4)] for i in range(4)]
        k=[[item(i,j,5) for j in range(4)] for i in range(4)]
        H=D(F(4+seed,19),[F(seed+j+1,37+j) for j in range(4)])
        eta=D(F(5+seed,23),[F(2*seed+j,41+j) for j in range(4)])
        kk=D(F(7+seed,29)); w=D(F(9+seed,31))
        tau=sum((m[i][i]/diag[i] for i in range(4)),D())
        aa2=sum((m[i][j]*m[j][i]/(diag[i]*diag[j]) for i in range(4) for j in range(4)),D())
        my=[[2*a*m[i][j]+(a*H+sp*eta/3)*base[i][j] for j in range(4)] for i in range(4)]
        etay=p*H/2-spp*eta/2
        Hy=-4*a*H+2*a*tau+4*sp*eta/3
        qr=[[(a*H*H/4-sp*H*eta/6-spp*eta*eta/6)*base[i][j]
             -(a*H+sp*eta/3)*m[i][j] for j in range(4)] for i in range(4)]
        qp=p*H*H/4+spp*H*eta/2+sppp*eta*eta/2
        qc=4*a*H*H-3*a*H*tau+2*a*aa2-4*sp*H*eta/3+sp*eta*tau/3
        fm=[[-2*qr[i][j] for j in range(4)] for i in range(4)]
        if injection=='factorial':
            fm=[[-qr[i][j] for j in range(4)] for i in range(4)]
        if injection=='isotropic-source':
            fm=[[fm[i][j] if i==j else D() for j in range(4)] for i in range(4)]
        fp=-qp
        if injection=='scalar-third':
            fp=fp+sppp*eta*eta/4
        gn=4*spp*eta*eta/3+4*sp*H*eta-10*a*H*H+8*a*H*tau-4*a*aa2
        if injection=='omit-matrix-square':
            gn=gn+4*a*aa2
        ky=[[2*a*k[i][j]+(a*kk+sp*w/3)*base[i][j]+fm[i][j] for j in range(4)] for i in range(4)]
        wy=-spp*w/2-sp*kk/4+fp
        kky=-4*a*kk+2*a*sum((k[i][i]/diag[i] for i in range(4)),D())+4*sp*w/3+gn
        # Direct exact laws, first on affine paths, then on completed second paths.
        for affine in (True,False):
            g=[[P(base[i][j],m[i][j],0 if affine else k[i][j]/2) for j in range(4)] for i in range(4)]
            gy=[[P(2*a*base[i][j],my[i][j],0 if affine else ky[i][j]/2) for j in range(4)] for i in range(4)]
            gi=invseries(g,diag)
            lapse2=P(1,H,0 if affine else kk/2)
            lapse=lapse2.sqrt1()
            ly=P(0,Hy,0 if affine else kky/2)
            phiy=P(p,etay,0 if affine else wy/2)
            sig=P(6*a,sp*eta,(spp*eta*eta+(0 if affine else sp*w))/2)
            sigp=P(sp,spp*eta,(sppp*eta*eta+(0 if affine else spp*w))/2)
            extr=[[gy[i][j]/(2*lapse) for j in range(4)] for i in range(4)]
            residual=[[extr[i][j]-sig*g[i][j]/6 for j in range(4)] for i in range(4)]
            scalar=phiy/lapse+sigp/2
            normal_gauge=ly/(2*lapse2*lapse2)-trace(matmul(gi,gy))/(2*lapse2)
            normal_gauge=normal_gauge+a*sum((gi[i][i]*diag[i] for i in range(4)),P())
            for i in range(4):
                for j in range(i,4):
                    equal(residual[i][j].c[1],D(),'first physical metric')
                    equal(2*residual[i][j].c[2],qr[i][j] if affine else D(),'physical metric Hessian/completion')
            equal(scalar.c[0],D(),'scalar convention')
            equal(scalar.c[1],D(),'first physical scalar')
            equal(2*scalar.c[2],qp if affine else D(),'physical scalar Hessian/completion')
            equal(normal_gauge.c[1],D(),'first normal gauge')
            equal(2*normal_gauge.c[2],qc if affine else D(),'normal gauge Hessian/completion')
            if not affine:
                n2=kk/2-H*H/4
                if injection=='lapse-conversion':
                    n2=kk/2
                equal(2*lapse.c[2],n2,'metric versus lapse coefficient')
            if affine:
                # Direct Codazzi contraction on an anisotropic rational boundary jet.
                mixed=matmul(gi,extr)
                trK=trace(mixed)
                conn=[[[sum((gi[d][e]*(g[c][e].deriv(b)+g[b][e].deriv(c)-g[b][c].deriv(e))/2
                              for e in range(4)),P()) for c in range(4)] for b in range(4)] for d in range(4)]
                for aa in range(4):
                    ena=sum((mixed[b][aa].deriv(b) for b in range(4)),P())-trK.deriv(aa)
                    ena+=sum((conn[b][b][d]*mixed[d][aa]-conn[d][b][aa]*mixed[b][d]
                              for b in range(4) for d in range(4)),P())
                    ena-=phiy/lapse*P(0,D(eta.g[aa]))
                    expected=sum((qr[b][aa].g[b]/diag[b] for b in range(4)),F())
                    expected-=sum((qr[b][b].g[aa]/diag[b] for b in range(4)),F())
                    equal(ena.c[1].v,F(0),'first Codazzi matter cancellation')
                    equal(2*ena.c[2].v,expected,'quadratic tangential source compatibility')
        equal(sum((fm[i][i]/diag[i] for i in range(4)),D())-2*qc,gn,'gauge/physical trace join')
    return count


CONTROLS={
    'scalar-sign':'scalar convention',
    'factorial':'physical metric Hessian/completion',
    'isotropic-source':'physical metric Hessian/completion',
    'scalar-third':'physical scalar Hessian/completion',
    'omit-matrix-square':'normal gauge Hessian/completion',
    'lapse-conversion':'metric versus lapse coefficient',
}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',required=True)
    ap.add_argument('--run-controls',action='store_true')
    ap.add_argument('--inject-error',choices=CONTROLS)
    args=ap.parse_args()
    output=Path(args.output)
    if output.exists():
        raise SystemExit('output already exists')
    try:
        count=check(args.inject_error)
        controls={}
        if args.run_controls:
            for name,reason in CONTROLS.items():
                try:
                    check(name)
                except Failure as exc:
                    if str(exc)!=reason:
                        raise Failure('wrong control reason: '+name+': '+str(exc))
                    controls[name]=reason
                else:
                    raise Failure('control unexpectedly passed: '+name)
        report={'status':'PASS','exact_comparisons':count,'rational_anisotropic_boundary_jets':3,
                'controls_rejected_for_named_reason':controls,
                'checks':['affine metric/scalar Hessians','all ten physical metric completion components',
                          'lapse versus yy-metric conversion','direct relative-normal-gauge Hessian',
                          'normal-gauge/physical-trace join','Codazzi quadratic source compatibility',
                          'actual harmonic gradient/Hessian and Lorentzian tensor contractions'],
                'scope':'Exact finite witnesses supplement displayed all-point identities; no actual-field numeric caps, U2 evolution or IBVP theorem.'}
        output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    except Failure as exc:
        raise SystemExit(str(exc))


if __name__=='__main__':
    main()
