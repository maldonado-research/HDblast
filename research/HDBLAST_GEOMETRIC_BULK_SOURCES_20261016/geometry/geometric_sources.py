"""Exact registered-frame geometric bulk source adapter (stdlib only).

Indices: 0=y, 1=timelike frame vector, 2..5=ambient tangent slots on S^3.
The redundant sphere slots are contracted with delta; all are tangent.
Default angular backend is sibling projection/angular_generic.py. All angular
operations route through the explicit backend argument. No sampled evaluation
or field equation enters the tensor construction.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import sys

from indexed_geometry import indexed_geometry, natural_sources

MULTI = ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))


class Jet:
    """Ordinary two-variable derivatives, truncated to total degree <= n."""
    __slots__ = ('v', 'n')
    def __init__(self, values, n=2):
        self.n=n; self.v=tuple(values[:(1,3,6)[n]])
        if len(self.v)!=(1,3,6)[n]: raise ValueError('JET_LENGTH')
    def coerce(self, x):
        return x if isinstance(x,Jet) else Jet((x,)+(0,)*5,self.n)
    def truncate(self,n):return Jet(self.v,min(n,self.n))
    def __bool__(self):return any(bool(x) for x in self.v)
    def __add__(self,x):
        x=self.coerce(x); n=min(self.n,x.n)
        return Jet(tuple(a+b for a,b in zip(self.v,x.v)),n)
    __radd__=__add__
    def __neg__(self):return Jet(tuple(-x for x in self.v),self.n)
    def __sub__(self,x):return self+-self.coerce(x)
    def __rsub__(self,x):return self.coerce(x)+-self
    def __mul__(self,x):
        x=self.coerce(x); n=min(self.n,x.n); out=[]
        for k in MULTI[:(1,3,6)[n]]:
            out.append(sum((comb(k[0],i[0])*comb(k[1],i[1])*self.v[j]*x.v[MULTI.index((k[0]-i[0],k[1]-i[1]))]
                       for j,i in enumerate(MULTI[:(1,3,6)[n]]) if i[0]<=k[0] and i[1]<=k[1]),0))
        return Jet(out,n)
    __rmul__=__mul__
    def inverse(self):
        v=self.v; inv=F(1)/v[0]
        vals=[inv]
        if self.n>=1:vals.extend((-v[1]*inv*inv,-v[2]*inv*inv))
        if self.n>=2:vals.extend((2*v[1]*v[1]*inv**3-v[3]*inv*inv,
                                 2*v[1]*v[2]*inv**3-v[4]*inv*inv,
                                 2*v[2]*v[2]*inv**3-v[5]*inv*inv))
        return Jet(vals,self.n)
    def __truediv__(self,x):return self*self.coerce(x).inverse()
    def __rtruediv__(self,x):return self.coerce(x)*self.inverse()
    def derivative(self,d):
        if not self.n:raise ValueError('JET_DERIVATIVE_ORDER')
        return Jet(tuple(self.v[MULTI.index((i[0]+(d==0),i[1]+(d==1)))]
                         for i in MULTI[:(1,3,6)[self.n-1]]),self.n-1)


class Angular:
    """Small arithmetic wrapper around the selected angular polynomial ring."""
    __slots__=('p','ring')
    def __init__(self,p,ring):self.p=p;self.ring=ring
    def coerce(self,x):
        if isinstance(x,Angular):return x
        return Angular({(0,0,0,0):x} if bool(x) else {},self.ring)
    def __bool__(self):return bool(self.p)
    def __add__(self,x):
        x=self.coerce(x);return Angular(self.ring.add(self.p,x.p),self.ring)
    __radd__=__add__
    def __neg__(self):return Angular(self.ring.scale(-1,self.p),self.ring)
    def __sub__(self,x):return self+-self.coerce(x)
    def __rsub__(self,x):return self.coerce(x)+-self
    def __mul__(self,x):
        x=self.coerce(x);return Angular(self.ring.mul(self.p,x.p),self.ring)
    __rmul__=__mul__
    def derivative(self,d):
        return Angular(self.ring.clean({e:v.derivative(d) for e,v in self.p.items()}),self.ring)
    def angular_diff(self,d):return Angular(self.ring.diff(self.p,d),self.ring)
    def truncate(self,n):
        return Angular(self.ring.clean({e:v.truncate(n) for e,v in self.p.items()}),self.ring)
    def value(self):
        return Angular(self.ring.clean({e:v.v[0] for e,v in self.p.items()}),self.ring)


def _backend():
    path=str(Path(__file__).resolve().parent.parent/'projection')
    if path not in sys.path:sys.path.insert(0,path)
    import angular_generic
    return angular_generic


def frame_two_jet(m,s,b,angular_backend=None,mutation=None):
    """Build q,dq,ddq and eta derivatives, including every connection slot."""
    for mapping in (m,s,b):
        if any(isinstance(value,float) for value in mapping.values()):
            raise ValueError('FLOAT_INPUT_REQUIRES_EXACT_COEFFICIENT_RING')
    R=angular_backend or _backend(); I=range(6)
    wrap=lambda p:Angular(p,R)
    const=lambda c:wrap({(0,0,0,0):c} if bool(c) else {})
    coeff=lambda key:Jet(tuple(m[key+suffix] for suffix in ('','y','t','yy','yt','tt')))
    rho,a,ay,c,th=(b[k] for k in ('rho','a','ay','c','theta'))
    # Only first derivatives of connection coefficients enter nabla²q. Values
    # assigned to unused second derivatives below are deliberately zero.
    rj=Jet((rho,a*rho,0,(ay+a*a)*rho,0,0))
    cj=Jet((c,0,c*th,0,0,c))
    aj=Jet((a,ay,0,0,0,0))
    tj=Jet((th,0,1-th*th,0,0,0))
    ir=const(rj.inverse()); j=const((rj*cj).inverse())
    wa=const(aj); wh=const(tj/rj)
    X=[wrap(p) for p in R.X]; K={(i,k):wrap(v) for (i,k),v in R.P.items()}
    Y=wrap(R.Y); dY={k:wrap(v) for k,v in R.grad(R.Y).items()}
    hY={k:wrap(v) for k,v in R.hess(R.Y).items()}
    A,B,C,D,E,L,M=(const(coeff(k)) for k in ('A','B','C','D','E','L','M'))
    zero=const(Jet((0,)*6)); q={(i,k):zero for i,k in product(I,repeat=2)}
    q[0,0]=A*Y;q[0,1]=q[1,0]=B*Y;q[1,1]=D*Y
    for i in range(4):
        q[0,i+2]=q[i+2,0]=C*dY[i,]
        q[1,i+2]=q[i+2,1]=E*dY[i,]
        for k in range(4):q[i+2,k+2]=L*Y*K[i,k]+M*hY[i,k]

    def cov_component(t,idx,d,order):
        value=t[idx]
        if d==0:return value.derivative(0).truncate(order)
        if d==1:out=ir*value.derivative(1)
        else:
            i=d-2
            intrinsic=sum((K[i,k]*value.angular_diff(k) for k in range(4)),zero)
            # Intrinsic sphere connection acts on every angular tensor slot,
            # including a prior derivative slot in t=dq.
            for slot,v in enumerate(idx):
                if v>=2:
                    repl=list(idx);repl[slot]=d
                    intrinsic+=X[v-2]*t[tuple(repl)]
            out=j*intrinsic
        for slot,v in enumerate(idx):
            if mutation=='omit_derivative_slot' and len(idx)==3 and slot==0:continue
            repl=list(idx)
            if d==1:
                if v in (0,1):
                    repl[slot]=1-v;out-=wa*t[tuple(repl)]
            elif v==0:
                repl[slot]=d;out-=wa*t[tuple(repl)]
            elif v==1:
                repl[slot]=d;out-=wh*t[tuple(repl)]
            else:
                repl[slot]=0;out+=wa*K[d-2,v-2]*t[tuple(repl)]
                repl[slot]=1;out-=wh*K[d-2,v-2]*t[tuple(repl)]
        return out.truncate(order)

    dqj={(d,a,b0):cov_component(q,(a,b0),d,1) for d,a,b0 in product(I,repeat=3)}
    ddq={(e,d,a,b0):cov_component(dqj,(d,a,b0),e,0).value()
         for e,d,a,b0 in product(I,repeat=4)}
    sj=const(Jet(tuple(s[k] for k in ('s','sy','st','syy','syt','stt'))))
    eta=sj*Y; etaj={():eta}
    detaj={(d,):cov_component(etaj,(),d,1) for d in I}
    ddeta={(e,d):cov_component(detaj,(d,),e,0).value() for e,d in product(I,repeat=2)}
    g0={(a,b0):const(0) for a,b0 in product(I,repeat=2)}
    g0[0,0]=const(1);g0[1,1]=const(-1)
    for i,k in product(range(4),repeat=2):g0[i+2,k+2]=K[i,k]
    bggrad={(a,):const(b['p'] if a==0 else 0) for a in I}
    bghess={(a,b0):g0[a,b0]*(b['a']*b['p']) for a,b0 in product(I,repeat=2)}
    bghess[0,0]=const(b['py'])
    return dict(q={k:v.value() for k,v in q.items()},dq={k:v.value() for k,v in dqj.items()},
                ddq=ddq,eta=eta.value(),deta={k:v.value() for k,v in detaj.items()},
                ddeta=ddeta,g0=g0,background_gradient=bggrad,background_hessian=bghess,
                constant=const)


def evaluate_sources(m,s,b,angular_backend=None,harmonic=True,mutation=None,include_jets=False):
    """Return full exact angular source tensors for supplied ordinary y,t jets.

    harmonic=True selects registered actual-first-field formulas; it does NOT
    certify the input satisfies C1=0. C1 is always returned for diagnosis.
    general QEred and Qphi are also returned to support off-shell audits.
    """
    data=frame_two_jet(m,s,b,angular_backend,mutation)
    sig=(1,-1,1,1,1,1)
    geo=indexed_geometry(data['q'],data['dq'],data['ddq'],sig,mutation)
    c=data['constant']
    kw=dict(q=data['q'],eta=data['eta'],deta=data['deta'],ddeta=data['ddeta'],
            background_hessian=data['background_hessian'],background_gradient=data['background_gradient'],
            g0=data['g0'],signature=sig,geometry=geo,Vp=c(b['Vp']),Vpp=c(b['Vpp']),Vppp=c(b['Vppp']))
    selected=natural_sources(**kw,harmonic=harmonic,mutation=mutation)
    general=natural_sources(**kw,harmonic=False,mutation=mutation)
    sources={key:geo[key] for key in ('QRic','QC','QC_actual','QC_background','C1','dQC')}
    sources['QCup']={(i,):sig[i]*geo['QC_background'][i,] for i in range(6)}
    sources.update(selected)
    sources['QEred_general']=general['QEred'];sources['Qphi_general']=general['Qphi']
    sources['gauge_gradient']=geo['harmonic_reduction']
    sources['PC1']={k:geo['harmonic_reduction'][k]-geo['general_reduction'][k] for k in geo['QRic']}
    factor=c(b['ell']*b['ell']*b['rho']*b['rho'])
    sources['metric_forcing']={k:-2*factor*v for k,v in sources['QEred'].items()}
    sources['scalar_forcing']=factor*sources['Qphi']
    if include_jets:sources.update({k:data[k] for k in ('q','dq','ddq','deta','ddeta')})
    return {key:({k:v.p for k,v in val.items()} if isinstance(val,dict) else val.p)
            for key,val in sources.items()}
