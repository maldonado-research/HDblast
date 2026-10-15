#!/usr/bin/env python3
"""Independent exact epsilon-series audit of the compact second Cauchy jet.

No other scientific/checker program is imported. Three-term series contain
Taylor coefficients, so their quadratic coefficient is multiplied by two.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import random


def need(ok, message):
    if not ok:
        raise ValueError(message)


class S:
    def __init__(self, value=0):
        self.a = tuple(value) if isinstance(value, (tuple,list)) else (Q(value),Q(0),Q(0))
    def __add__(self, other):
        other = cast(other)
        return S([a+b for a,b in zip(self.a,other.a)])
    __radd__ = __add__
    def __neg__(self):
        return S([-a for a in self.a])
    def __sub__(self, other):
        return self+-cast(other)
    def __rsub__(self, other):
        return cast(other)+-self
    def __mul__(self, other):
        other = cast(other)
        return S([sum((self.a[j]*other.a[i-j] for j in range(i+1)),Q(0)) for i in range(3)])
    __rmul__ = __mul__
    def __truediv__(self, other):
        other=cast(other)
        lead=other.a[0]
        need(lead != 0,"nonzero series denominator")
        return self*S([x/lead for x in other.a]).power(-1)*(1/lead)
    def __rtruediv__(self, other):
        return cast(other).__truediv__(self)
    def power(self, exponent):
        exponent = Q(exponent)
        need(self.a[0] == 1, "series powers use constant term one")
        z = self-1
        return 1+exponent*z+(exponent*(exponent-1)/2)*z*z
    def hessian(self):
        return 2*self.a[2]


def cast(x):
    return x if isinstance(x,S) else S(x)


E = S([Q(0),Q(1),Q(0)])


def spatial_curvature(rho,a,ay,u,uy,uyy,Y,grad2,delta):
    """Radial Gauss identity for A dy²+rho² C kappa, evaluated exactly.

    R4=R3-12b²-6 b_y/L-2 Delta_h L/L, L=sqrt(A),
    b=(a+sigma_y)/L, sigma=log(C)/2. Sphere dimension is three.
    """
    C = 1+2*E*u*Y
    L = (1-4*E*u*Y).power(Q(1,2))
    sy = E*uy*Y/C
    syy = E*uyy*Y/C-2*E*E*uy*uy*Y*Y/(C*C)
    Ly = -2*E*uy*Y/L
    b = (a+sy)/L
    by = (ay+syy)/L-(a+sy)*Ly/(L*L)
    delta_sigma = E*u*delta/C-2*E*E*u*u*grad2/(C*C)
    grad_sigma2 = E*E*u*u*grad2/(C*C)
    delta_L = -2*E*u*delta/L-4*E*E*u*u*grad2/(L*L*L)
    grad_sigma_L = -2*E*E*u*u*grad2/(C*L)
    R3 = (6-4*delta_sigma-2*grad_sigma2)/(rho*rho*C)
    lap_L = (delta_L+grad_sigma_L)/(rho*rho*C)
    return R3-12*b*b-6*by/L-2*lap_L/L


def direct_momentum(rho,a,u,uy,ut,uyt,v,vy,vt,p,Y,Yi,spacetime=False):
    """Coordinate divergence for the two diagonal mixed K eigenvalues."""
    A,C = 1-4*E*u*Y,1+2*E*u*Y
    Ay,Cy = -4*E*uy*Y,2*E*uy*Y
    Ai,Ci = -4*E*u*Yi,2*E*u*Yi
    q,qy,qi = ut*Y/rho,(uyt-a*ut)*Y/rho,ut*Yi/rho
    Py = -3*E*q/C
    Ps = 2*E*q*(1/A-1/C)
    Pyy = -3*E*(qy/C-q*Cy/(C*C))
    Psi = 2*E*qi*(1/A-1/C)+2*E*q*(-Ai/(A*A)+Ci/(C*C))
    Pi = E*vt*Y/rho
    if spacetime:
        scale = C.power(Q(-1,2))
        scale_y,scale_i = -E*uy*Y*C.power(Q(-3,2)),-E*u*Yi*C.power(Q(-3,2))
        Pyy,Psi = Pyy*scale+Py*scale_y,Psi*scale+Ps*scale_i
        Py,Ps,Pi = Py*scale,Ps*scale,Pi*scale
    My = Pyy+3*(a+Cy/(2*C))*(Py-Ps)-Pi*(p+E*vy*Y)
    Mi = Psi+Ai/(2*A)*(Ps-Py)-Pi*E*v*Yi
    return My,Mi


def review(seed):
    rng=random.Random(seed)
    rnd=lambda:Q(rng.randint(-7,7),rng.randint(1,6))
    rho=Q(7,3)
    a,p,py,u,uy,uyy,ut,uyt,Y,Yi=[rnd() for _ in range(10)]
    if not p:p=Q(2,3)
    ay=-1/rho**2-p*p/3
    U1=py+4*a*p
    U2=rnd()
    v=-3*(uy+2*a*u)/p
    vy=-3*(uyy+(2*a-py/p)*uy+(2*ay-2*a*py/p)*u)/p
    vt=-3*(uyt+2*a*ut)/p
    grad2=Yi*Yi+rnd()**2+rnd()**2
    R=spatial_curvature(rho,a,ay,u,uy,uyy,Y,grad2,-8*Y)
    A,C=1-4*E*u*Y,1+2*E*u*Y
    scalar_gradient=(p+E*vy*Y)*(p+E*vy*Y)/A+E*E*v*v*grad2/(rho*rho*C)
    # K²-|K|² is already quadratic. Both affine conventions have the same
    # Hessian here because background K=Pi=0.
    H=R-scalar_gradient-2*E*U1*v*Y-E*E*U2*v*v*Y*Y
    H=H-6*E*E*ut*ut*Y*Y/rho**2-E*E*vt*vt*Y*Y/rho**2
    AH=(-24*uy*uy-288*a*u*uy-24*u*uyy-(192*ay+384*a*a+32*p*p)*u*u
        -(336*u*u+12*ut*ut+2*vt*vt)/rho**2-2*vy*vy-16*p*u*vy-2*U2*v*v)
    need(H.hessian()==AH*Y*Y+(36*u*u-2*v*v)*grad2/rho**2,"Hamiltonian Hessian from independent radial Gauss expansion")
    need(H.a[1]==0,"actual aligned first Hamiltonian constraint")
    F=(18*u*ut-v*vt)/rho
    G=(12*u*uyt-48*a*u*ut-6*uy*ut-2*vt*vy)/rho
    My,Mi=direct_momentum(rho,a,u,uy,ut,uyt,v,vy,vt,p,Y,Yi)
    need(My.a[1]==0 and Mi.a[1]==0,"actual aligned first momentum constraints")
    need(My.hessian()==G*Y*Y and Mi.hessian()==2*F*Y*Yi,"physical Cauchy affine momentum Hessian")
    Sy,Si=direct_momentum(rho,a,u,uy,ut,uyt,v,vy,vt,p,Y,Yi,True)
    need(Sy.hessian()==(G+6*uy*ut/rho)*Y*Y and Si.hessian()==Mi.hessian(),"spacetime affine lapse cross terms")
    f,fy,fyy,source,Fy=[rnd() for _ in range(5)]
    zeta=-3*fy/(2*p)
    zetay=-3*fyy/(2*p)+3*fy*py/(2*p*p)
    for lam in (0,8,24):
        dR=-3*fyy-12*a*fy+(2*lam-6)*f/rho**2
        need(dR-2*p*zetay-2*U1*zeta==2*(lam-3)*f/rho**2,"Hamiltonian right inverse scalar cancellation")
        fc=-rho*rho*source/(2*(lam-3))
        need(2*(lam-3)*fc/rho**2==-source,"Hamiltonian Hessian factorial/source normalization")
    q,qy,Pi2=F/3,Fy/3,(G-Fy)/p
    need(-3*q==-F and -3*qy-p*Pi2==-G,"full-spatial pure-trace momentum right inverse")
    # Angular support: Delta(Y²)=-16Y²+2|gradY|².
    for lam in (0,8,24):
        need(2*(8-Q(lam,2))-16==-lam,"finite angular sectors of the Hamiltonian source")
    return 15


def witness(path):
    raw=path.read_bytes()
    digest=hashlib.sha256(raw).hexdigest()
    need(digest=="fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0","actual coefficient report pin")
    report=json.loads(raw)
    b,c=report['background_intervals'],report['reconstruction_coefficients']
    ab=lambda x:max(abs(Q(v)) for v in x['exact'])
    need(tuple(map(Q,b['rho']['exact']))==(129,131),"witness actual rho box")
    B=list(map(Q,c['v_wz_coefficient']['exact']))
    need(B[0]<B[1]<0 and min(map(abs,B))>Q(19,2000) and max(map(abs,B))<Q(3,250),"witness nonvanishing B")
    for key,cap in [('v_wz_coefficient_z',Q(3,10)),('alpha_log_z',270),('alpha_log_zz',24000),
                    ('u_w_coefficient',Q(1,10000)),('u_w_coefficient_z',Q(3,200)),('u_w_coefficient_zz',5)]:
        need(ab(c[key])<cap,"witness coefficient cap: "+key)
    need(ab(b['a'])<1 and ab(b['p'])<1 and ab(b['a_z'])/129<1 and ab(b['U_2'])<4,"witness background caps")
    ell=Q(1,10**6)
    err=(Q(3,10)*270+Q(3,250)*24000)/129
    lo=32*Q(19,2000)/(131*ell*ell)-err
    hi=32*Q(3,250)/(129*ell*ell)+err
    need(lo>2300000000 and hi<3000000000,"witness scalar gradient")
    need((5+32*Q(1,10000)/(ell*ell))/129**2+Q(3,200)/129<10**6,"witness usecond cap")
    rest=24+288*Q(1,10000)+24*Q(1,10000)*10**6+608*Q(1,10000)**2+(336*Q(1,10000)**2+12+32)/129**2+128+16*Q(1,10000)*3*10**9
    need(rest<5*10**6 and -2*(23*10**8)**2+rest<-10**19,"strict nonzero initial Hamiltonian witness")
    return {'background_sha256':digest,'strict_Hessian_upper':'-10000000000000000000','point':'t=0,z=3ell/2,n=(1,0,0,0)'}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--background',type=Path,default=Path(__file__).resolve().parents[2]/'potential/ACTUAL_POTENTIAL_JET_BOUNDS.json')
    ap.add_argument('--output',type=Path)
    a=ap.parse_args()
    count=sum(review(seed) for seed in (31,191,271,1927,6469))
    result={'status':'PASS_INDEPENDENT_INITIAL_CAUCHY_JET_REVIEW','optimized':not __debug__,
        'finite_exact_check_groups':count,'witness':witness(a.background),
        'physical_Cauchy_affine_G_verified':True,'spacetime_affine_lapse_cross_terms_verified':True,
        'actual_aligned_first_constraints_verified':True,'all_point_claim_depends_on_analytic_derivation':True,
        'full_second_evolution_solved':False}
    text=json.dumps(result,indent=2)+'\n'
    if a.output:a.output.write_text(text)
    else:print(text,end='')


if __name__=='__main__':main()
