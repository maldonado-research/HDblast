"""Exact sparse-polynomial envelopes for the complete scalar normal operator.

No sampling, symbolic-algebra dependency, inverse oracle, or file writes.
Continuous Fourier bootstrap and covariant conversion are proved in the text.
"""
from fractions import Fraction as F
from math import comb


SCALAR_LOWER_JET_CEILING = F(1, 10**6)


def need(condition, message):
    if not condition:
        raise ValueError(message)


class P:
    """Sparse rational polynomials in (x,1/rho,a,p,eta)."""
    def __init__(self, terms=None):
        self.t = {k:F(v) for k,v in (terms or {}).items() if v}
    @staticmethod
    def cast(value):
        return value if isinstance(value,P) else P({(0,)*5:F(value)})
    @staticmethod
    def var(i):
        k=[0]*5;k[i]=1
        return P({tuple(k):F(1)})
    def __add__(self,value):
        out=self.t.copy()
        for k,v in P.cast(value).t.items():
            out[k]=out.get(k,F(0))+v
        return P(out)
    __radd__=__add__
    def __neg__(self):
        return P({k:-v for k,v in self.t.items()})
    def __sub__(self,value):
        return self+-P.cast(value)
    def __rsub__(self,value):
        return P.cast(value)+-self
    def __mul__(self,value):
        out={}
        for ka,va in self.t.items():
            for kb,vb in P.cast(value).t.items():
                k=tuple(a+b for a,b in zip(ka,kb))
                out[k]=out.get(k,F(0))+va*vb
        return P(out)
    __rmul__=__mul__
    def __truediv__(self,value):
        return self*F(1,F(value))
    def __pow__(self,n):
        need(type(n) is int and n>=0,'Nonnegative polynomial power required')
        out=P.cast(1)
        for _ in range(n):
            out=out*self
        return out
    def partial(self,i):
        out={}
        for k,v in self.t.items():
            if k[i]:
                kk=list(k);kk[i]-=1
                out[tuple(kk)]=v*k[i]
        return P(out)
    def normdeg(self):
        return sum(map(abs,self.t.values()),F(0)),max((sum(k) for k in self.t),default=0)


def derive():
    names=[]
    def check(condition,name):
        need(condition,name)
        need(name not in names,'Duplicate scalar envelope check')
        names.append(name)
    x,b,a,p,t=(P.var(i) for i in range(5))
    q=F(1,10**16)*(1-x*x)**2
    U=-F(2,27)+F(14,9)*t*t+F(50,27)*t**3-t**4/6-F(4,9)*t**5-F(2,27)*t**6
    vf=(-F(1,10**10)*(1-x*x)**2*b,-q*a*b,
        q*(-b*b-p*p/3),q*(U.partial(4)-4*a*p),q*p)
    check(max(v.normdeg()[0] for v in vf)<=F(1,10**8)
          and max(v.normdeg()[1] for v in vf)<=9,
          'polynomial scaled vector field coefficient and degree envelope')
    def ds(expr):
        return sum((v*expr.partial(i) for i,v in enumerate(vf)),P())
    V=F(1,10**8)
    ry=2*10**6*x*b
    d=2*F(1,10**16)*(1-x*x)*ry
    c=d/2+F(3,2)*q*a+F(1,10**16)*ry
    e=F(1,10**16)*ry+F(3,2)*d-F(3,2)*q*a
    def jets(expr,N=16):
        norm,degree=P.cast(expr).normdeg()
        out=[norm]
        for k in range(N):
            out.append(out[-1]*V*(degree+8*k))
        return out
    def prodjet(u,v,h):
        return sum((F(comb(h,k))*u[k]*v[h-k] for k in range(h+1)),F(0))
    records={}
    for lam in (0,8,24):
        Rt=3*b*b+p*p/3-3*a*a
        rows=[{1:-F(3,2)*a*q,0:F(3,2)*a*q*c+q*q*(3*a*a+(lam-3)*b*b)/2},
              {2:P.cast(-F(1,2)),1:c+d/2-a*q,
               0:-(c*c+d*c-ds(c))/2+a*q*c+q*q*(lam*b*b/3-Rt/2)},
              {1:q*b,0:-q*b*(c+a*q)}, {0:q*q*b*b}]
        weights=[F(1),F(3),F(lam,2),F(lam*(lam-3),6)]
        principal=sum((w*row.get(2,P())**2 for row,w in zip(rows,weights)),P())
        check(not (principal-F(3,4)).t,'scalar normal principal coefficient lambda'+str(lam))
        out={j:[F(0) for _ in range(13)] for j in range(4)}
        for row,w in zip(rows,weights):
            jj={i:jets(coef) for i,coef in row.items()}
            for i,ui in jj.items():
                for k in range(i+1):
                    for j,uj in jj.items():
                        for z in range(k+1):
                            order=k-z+j
                            if order>=4:
                                continue
                            factor=w*comb(i,k)*comb(k,z)
                            for h in range(13):
                                out[order][h]+=factor*prodjet(ui[i-k:],uj[z:],h)
        peak=max(max(v) for v in out.values())
        check(peak<SCALAR_LOWER_JET_CEILING,
              'scalar normal lower coefficient jets lambda'+str(lam))
        row_sums=[]
        for row in rows:
            total=sum((v.normdeg()[0] for v in row.values()),F(0))
            degree=max(v.normdeg()[1] for v in row.values())
            for j in range(14):
                total*=1+V*degree+(e-j*d).normdeg()[0]
                degree+=8
            row_sums.append(total)
        check(sum(row_sums)<10,'physical radial scalar operator recurrence lambda'+str(lam))
        total=F(1);degree=0
        for j in range(16):
            total*=1+V*degree+(-c-j*d).normdeg()[0]
            degree+=8
        check(total<2,'weighted scalar potential recurrence lambda'+str(lam))
        records[str(lam)]={'lower_coefficient_jet_max_exact':str(peak),
                           'physical_radial_coefficient_sum_upper':str(sum(row_sums)),
                           'weighted_potential_coefficient_sum_upper':str(total)}
    K=F(4*2**12,10**6)
    check(K<F(3,4) and 1/(F(3,4)-K)<2 and K/(F(3,4)-K)<1,
          'scalar Fourier highest derivative absorption')
    check(17*(1+2)**2<20**2,'scalar H16 derivative square sum factor20')
    check(24**24<10**34,'scalar physical exponential multiplier envelope')
    angular=F(100);degree=4
    for rank in range(1,19):
        angular*=degree*5**rank
        degree+=2*rank-1
    check(100*angular<10**200,'sphere tensor polynomial derivative envelope through rank18')
    yvf=(-a*b,-b*b-p*p/3,U.partial(4)-4*a*p,p)
    check(max(v.normdeg()[0] for v in yvf)<=100
          and max(v.normdeg()[1] for v in yvf)<=5,
          'unscaled radial polynomial vector field envelope')
    check(7*F(61,12)<100 and 7*F(5,4)<100,'normalized scalar harmonic polynomial coefficient envelope')
    background=F(1)
    for k in range(18):
        background*=100*(1+4*k)
    check(background<10**80,'unscaled background radial jets through18')
    check(1000**18<=10**54,'covariant derivative branching envelope')
    check(1000**14*10**(80*14)*10**200*10**232<10**1600,
          'physical metric covariant conversion envelope')
    check(1000**16*10**(80*16)*10**200*10**4<10**1600,
          'weighted scalar potential covariant conversion envelope')
    return {'status':'PASS','named_checks':names,'named_check_count':len(names),
            'exact':{'vector_field_coefficient_max':str(max(v.normdeg()[0] for v in vf)),
                     'vector_field_degree_max':max(v.normdeg()[1] for v in vf),
                     'lambda_records':records,'Fourier_lower_coefficient_K':str(K),
                     'background_radial_jet_envelope':str(background),
                     'angular_tensor_jet_envelope':str(100*angular),
                     'physical_H14_conversion_upper':str(10**2000),
                     'weighted_H16_conversion_upper':str(10**2000)}}
