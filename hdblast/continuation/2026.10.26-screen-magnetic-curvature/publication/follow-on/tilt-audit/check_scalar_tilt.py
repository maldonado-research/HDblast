#!/usr/bin/env python3
"""Independent exact pointwise tilt algebra and formal gauge substitutions.

No frozen imports, actual field values, interval evaluation, or sign claim.
"""
from fractions import Fraction as F
from itertools import product,permutations
import json

I=range(3)
def eps(a,b,c):
    p=(a,b,c)
    return 0 if len(set(p))!=3 else (-1)**sum(p[i]>p[j] for i in I for j in range(i+1,3))
def sym(t):return {(a,b):(t[a,b]+t[b,a])/2 for a,b in product(I,repeat=2)}
def axial(v,S):return sym({(a,b):sum((F(eps(a,c,d))*v[c]*S[b,d] for c,d in product(I,repeat=2)),F(0)) for a,b in product(I,repeat=2)})
def wedge(S,b,e,c,d):return (int(b==c)*S[e,d]-int(b==d)*S[e,c]-int(e==c)*S[b,d]+int(e==d)*S[b,c])
def basis_sym(i,j):return {(a,b):F((a==i and b==j) or (a==j and b==i)) for a,b in product(I,repeat=2)}

# Tiny formal Laurent ring; variables are independent and denominators rho,c.
NAMES=('rho','c','theta','b','bt','by','ft','dt','dtt','dy','dyt')
class P:
    def __init__(self,d=None):self.d={k:F(v) for k,v in (d or {}).items() if v}
    @staticmethod
    def cv(x):return x if isinstance(x,P) else P({(0,)*len(NAMES):F(x)})
    def __add__(self,x):
        z=self.d.copy()
        for k,v in P.cv(x).d.items():z[k]=z.get(k,F(0))+v
        return P(z)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.d.items()})
    def __sub__(self,x):return self+-P.cv(x)
    def __rsub__(self,x):return P.cv(x)+-self
    def __mul__(self,x):
        z={}
        for k,v in self.d.items():
            for l,w in P.cv(x).d.items():
                h=tuple(a+b for a,b in zip(k,l));z[h]=z.get(h,F(0))+v*w
        return P(z)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=P.cv(x)
        if len(x.d)!=1:raise ValueError('MONOMIAL_DIVISOR_REQUIRED')
        k,v=next(iter(x.d.items()))
        return self*P({tuple(-i for i in k):1/v})
    def __pow__(self,n):
        if n<0:return P.cv(1)/(self**(-n))
        out=P.cv(1)
        for _ in range(n):out=out*self
        return out

def variable(name):
    e=[0]*len(NAMES);e[NAMES.index(name)]=1
    return P({tuple(e):1})

def main():
    angular_count=0
    for k in I:
        v={a:F(a==k) for a in I}
        for i in I:
            for j in range(i,3):
                S=basis_sym(i,j)
                direct=sym({(a,b):sum((F(eps(a,c,d))*v[e]*wedge(S,b,e,c,d) for e,c,d in product(I,repeat=3)),F(0)) for a,b in product(I,repeat=2)})
                A=axial(v,S)
                if any(direct[key]!=-2*A[key] for key in A):raise ValueError('ANGULAR_RIEMANN_TILT_SIGN')
                angular_count+=1
        # Scalar radial-boost block C_bncd=delta_bc v_d-delta_bd v_c.
        radial=sym({(a,b):sum((F(eps(a,c,d))*(int(b==c)*v[d]-int(b==d)*v[c]) for c,d in product(I,repeat=2)),F(0)) for a,b in product(I,repeat=2)})
        if any(radial.values()):raise ValueError('RADIAL_BOOST_NOT_ZERO')
        trace=axial(v,{(a,b):F(a==b) for a,b in product(I,repeat=2)})
        if any(trace.values()):raise ValueError('TRACE_TILT_NOT_ZERO')
    rho,c,theta,b,bt,by,ft,dt,dtt,dy,dyt=map(variable,NAMES)
    R=rho*c
    # Independent differentiation of E=c*d_t-b/c and C=f/R+R*d_y.
    Et=c*dtt+c*theta*dt-bt/c+theta*b/c
    Ey=c*dyt-by/c
    Ct_plus_thetaC=ft/R+2*theta*R*dy+R*dyt
    B=ft/rho-rho*by
    D_plus_L=-2*bt+2*theta*b
    Pactual=(2*dtt+2*theta*dt)/(2*rho**2)-Et/(rho*R)+D_plus_L/(2*R**2)
    Qactual=(2*dyt+2*theta*dy)/(2*rho)-Ey/(2*R)-Ct_plus_thetaC/(2*rho*R)+B/(2*R**2)
    if Pactual.d:raise ValueError('REGISTERED_GAUGE_P_NOT_ZERO')
    if Qactual.d:raise ValueError('REGISTERED_GAUGE_Q_NOT_ZERO')
    print(json.dumps({'status':'PASS_EXACT_SCALAR_TILT_ALGEBRA','all_symmetric_tensor_vector_basis_pairs':angular_count,'angular_block_factor':'-2','radial_scalar_boost_block':'zero','trace_blocks':'zero','compact_correction':'Delta H2=2(V P-W Q) A=4(V P-W Q)(Tplus-Tminus)','registered_amplitude_P':'0','registered_amplitude_Q':'0','formal_substitution':'registered metric = Newtonian scalar plus arbitrary f,b,d gauge amplitudes','field_equations_used':False,'actual_values_evaluated':False,'general_observer_extension_certified':False,'scope':'pointwise algebra and registered scalar-tilt coefficient only'},indent=2))
if __name__=='__main__':main()
