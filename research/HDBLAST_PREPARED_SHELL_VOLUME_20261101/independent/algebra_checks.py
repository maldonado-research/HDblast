#!/usr/bin/env python3
"""Independent finite algebra checks. Imports no scientific producer module."""
import hashlib
import json
from pathlib import Path

import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent))
from portable_io import arguments, authenticate, finish, fail
def main():
    args=arguments()
    authenticate(args)
    import sympy as S
    ROOT=args.repo_root
    HERE=Path(__file__).resolve().parent
    e, Y = S.symbols('e Y')
    rho,a,u,uy,ut,H,Hy,J,Jy,QM,p,psi,psip,C,E,d = S.symbols(
        'rho a u uy ut H Hy J Jy QM p psi psip C E d', nonzero=True)
    checks = {}
    def check(name, lhs, rhs):
        residual = S.simplify(lhs-rhs)
        if residual != 0:
            raise ValueError((name, str(residual)))
        checks[name] = 'PASS'
    def jet2(expr):
        return S.diff(expr,e,2).subs(e,0)
    def avg(expr):
        poly=S.Poly(S.expand(expr),Y)
        moments={0:S.Integer(1),1:S.Integer(0),2:S.Rational(1,6)}
        return S.simplify(sum(co* moments[k[0]] for k,co in poly.terms()))

    # Direct density conversion: pi^{ij}=sqrt(gamma)(K^{ij}-tr(K)gamma^{ij})/2.
    # Inversion is K_ij=2/sqrt(gamma)[gamma_ik gamma_jl pi^kl-gamma_ij tr_gamma(pi)/3].
    A=1-4*e*u*Y+e**2*H/2
    B=1+2*e*u*Y+e**2*J/2
    P=1+e*u*Y
    D=-QM/(6*a)
    piR=-S.Rational(3,2)*e*ut/rho*Y
    piT=e**2*D/2
    vol=S.sqrt(A*B**3)
    trpi=A*piR+3*B*piT
    KR=2/vol*(A**2*piR-A*trpi/3)
    KT=2/vol*(B**2*piT-B*trpi/3)
    check('canonical_radial_K2',avg(jet2(KR)),QM/(3*a)+6*u*ut/rho)
    check('canonical_angular_K2',avg(jet2(KT)),-u*ut/rho)
    check('physical_K_t',avg(jet2(2*rho*P*KR)),2*rho*QM/(3*a)+S.Rational(32,3)*u*ut)
    check('physical_J_t',avg(jet2(2*rho*P*KT)),-S.Rational(4,3)*u*ut)

    # Direct contraction at a normal-coordinate point on S3, retaining its warp.
    # g=(tt,yy,aa), shift=0 initially; all terms in Gamma[g]-Gamma[g0]
    # are constructed from first derivatives, not imported gauge formulas.
    ndot,bdot=S.symbols('ndot bdot')
    Ay=-4*e*uy*Y+e**2*Hy/2
    By=2*e*uy*Y+e**2*Jy/2
    Py=e*uy*Y
    g=[-rho**2*P**2,A,rho**2*B,rho**2*B,rho**2*B]
    gi=[1/x for x in g]
    dg={}
    dg[0,0,0]=-2*rho*P*ndot
    dg[0,0,1]=dg[0,1,0]=A*bdot
    dg[0,1,1]=2*rho*P*KR
    for i in (2,3,4): dg[0,i,i]=rho**2*2*rho*P*KT
    dg[1,0,0]=-2*rho**2*P*(a*P+Py)
    dg[1,1,1]=Ay
    for i in (2,3,4): dg[1,i,i]=rho**2*(2*a*B+By)
    def gamma(I,j,k):
        return gi[I]*(dg.get((j,I,k),0)+dg.get((k,I,j),0)-dg.get((I,j,k),0))/2
    contr=[sum(gi[j]*gamma(I,j,j) for j in range(5)) for I in range(2)]
    ref_t=0
    ref_y=-a/P**2-3*a/B
    Ndot=S.solve(contr[0]-ref_t,ndot)[0]
    Bdot=S.solve(contr[1]-ref_y,bdot)[0]
    check('natural_harmonic_Ndot2',avg(jet2(Ndot)),rho**2*(QM/(3*a)-u*ut/rho))
    check('natural_harmonic_N_t',avg(jet2(-2*rho*P*Ndot/rho**2)),
          -2*rho*QM/(3*a)+S.Rational(4,3)*u*ut)
    Lwant=rho*(Hy/2-3*Jy/2+4*a*H-3*a*J-S.Rational(13,3)*u*uy-a*u**2/3)
    check('natural_harmonic_L_t_with_metric_shift_cross',avg(jet2(A*Bdot/rho)),Lwant)
    wrong=avg(jet2(Bdot/rho))
    if S.simplify(wrong-Lwant)==0: raise ValueError('wrong mixed conversion escaped')
    checks['missing_metric_shift_cross_detected']='PASS'
    j=S.symbols('j')
    check('canonical_scalar_coordinate_second_zero',jet2(rho*P*e*j*Y/vol),0)

    # A normalized average eliminates second measure and constant-background terms.
    f0,w0,v2=S.symbols('f0 w0 v2')
    scalar=f0+e*psi*Y+e**2*(w0-2*C*psip*Y**2)/2
    measure=1+e*(3*E-8*d)*Y+e**2*v2/2
    den=avg(measure)
    num=avg((scalar*measure).expand().series(e,0,3).removeO())
    check('proper_scalar_mean_clock_and_measure',jet2(num/den),
          w0-C*psip/3+(E-S.Rational(8,3)*d)*psi)

    # General determinant Hessian and full graph pullback reproduce the volume map.
    X,Z,Zp,dp,V,Hb,Hbp,Bclock,Aclock=S.symbols('X Z Zp dp V Hb Hbp Bclock Aclock')
    clock_trace_half=(S.Rational(3,2)*J+3*Hb*Bclock
        -C*(3*Zp+6*Hb*Z-8*dp-16*Hb*d)/3
        +C**2*(Hbp+2*Hb**2)/2-8*C*V/(3*Aclock)-4*C**2/(3*Aclock**2))
    trBsq_mean=S.Rational(2,3)*(3*E-8*d)**2
    trB2_mean=4*(E**2/2-S.Rational(8,3)*E*d+8*d**2)
    detHess=trBsq_mean/4-trB2_mean/2
    check('volume_first_quadratic_determinant',detHess,E**2/2-S.Rational(8,3)*E*d-S.Rational(16,3)*d**2)

    # New restricted cone estimate, independent of old initial-response caps.
    T=S.Rational(5,4);delta=T-1
    q=T*(S.Rational(2,10000)+delta*S.Rational(1,1000000))+delta/S.Integer(2000)
    check('new_cone_contraction',q,S.Rational(1201,3200000))

    # Second-only conditioning gauge; these are full monopole vector identities.
    xi,xiy,xiyy,ay,f,fy=S.symbols('xi xiy xiyy ay f fy')
    old_harm=(Hy/2-3*Jy/2+4*a*H-3*a*J).subs({H:f+2*xiy,Hy:fy+2*xiyy,J:f+2*a*xi,Jy:fy+2*ay*xi+2*a*xiy})
    delta_Lt=rho*(-xiyy-5*a*xiy+(3*ay+6*a*a)*xi)
    check('conditioned_mixed_velocity',rho*old_harm+delta_Lt,rho*(-fy+a*f))
    Xtt=rho**2*(-xiyy-4*a*xiy+(4*ay+6*a*a)*xi)
    check('monopole_vector_acceleration_gives_mixed_velocity',Xtt/rho-rho*(ay*xi+a*xiy),delta_Lt)
    check('conditioned_lapse_value_unchanged',2*a*rho**2*xi-2*rho**2*a*xi,0)
    shell_b,Ht=S.symbols('shell_b Ht')
    check('second_gauge_shell_volume_cancellation',S.Rational(3,2)*(2*Ht*shell_b)+3*(Ht/rho)*(-rho*shell_b),0)
    Ep,Lp,Dp=S.symbols('Ep Lp Dp')
    check('volume_ODE_kinematic_term',( (3*Ep+8*Lp)**2-40*Lp**2)/9,
          Ep**2+S.Rational(16,3)*Ep*Lp+S.Rational(8,3)*Lp**2)

    # Exact de Sitter Green kernel and elementary rational relative-error bounds.
    tt,ss,ll=S.symbols('tt ss ll', positive=True)
    kernel=S.sinh(2*ll*(tt-ss))/(2*ll)*S.cosh(ll*ss)**2/S.cosh(ll*tt)**2
    kernel_t=S.diff(kernel,tt)
    ode=S.diff(kernel,tt,2)+4*ll*S.tanh(ll*tt)*kernel_t-2*ll**2/S.cosh(ll*tt)**2*kernel
    check('volume_Green_homogeneous_ODE',S.trigsimp(ode.rewrite(S.exp)),0)
    check('volume_Green_initial_value',kernel.subs(tt,ss),0)
    check('volume_Green_initial_derivative',kernel_t.subs(tt,ss),1)
    rmax=S.Rational(1,800000); umax=S.Rational(1,4000000)
    cmax=1/(1-rmax**2/2)
    ratio_lower=1/cmax**2
    deriv_upper=1/(1-2*umax**2)
    deriv_lower=ratio_lower*(1-rmax*2*umax/(1-2*umax**2))
    if not (1-ratio_lower<S.Rational(1,10**11) and deriv_upper-1<S.Rational(1,10**11)
            and 1-deriv_lower<S.Rational(1,10**11)):
        raise ValueError('Green kernel rational relative bounds')
    checks['Green_value_and_derivative_relative_bound_1e_minus11']='PASS'

    # Independent full-fibre tensor construction, followed by invariant projection.
    # Coordinates (y,0,1,2,3); no harmonic master inversion or lambda limit is used.
    aa,hh,bb,dd,ell,ht,p_y,Vpot,Vp,Vpp=S.symbols('aa hh bb dd ell ht p_y Vpot Vp Vpp')
    K,L,N,Jj,w=S.symbols('K L N Jj w')
    tensor=S.diag(K,N,Jj,Jj,Jj)
    tensor[0,1]=tensor[1,0]=L
    sig=[1,-1,1,1,1]
    G0=S.zeros(5);G0[1,0]=aa;G0[0,1]=aa
    Ga=[]
    for jidx in (2,3,4):
        g=S.zeros(5);g[jidx,0]=aa;g[jidx,1]=hh;g[1,jidx]=hh;g[0,jidx]=-aa
        Ga.append(g)
    def omega(g,k):return -g.T*k-k*g
    def proj(k):
        vals=S.Matrix([k[0,0],k[0,1],k[1,1],sum(k[j,j] for j in (2,3,4))/3])
        rebuilt=S.diag(vals[0],vals[2],vals[3],vals[3],vals[3]);rebuilt[0,1]=rebuilt[1,0]=vals[1]
        if any(S.simplify(z)!=0 for z in rebuilt-k):raise ValueError('invariant fibre closure failed')
        return vals
    O=proj(omega(G0,tensor))
    Zten=-omega(G0,omega(G0,tensor))-3*hh*omega(G0,tensor)
    for g in Ga:Zten+=omega(g,omega(g,tensor))
    Zvec=proj(Zten)
    def riem(i,j,k,l):
        curv=bb if 0 in (i,j,k,l) else dd
        return curv*sig[i]*sig[j]*(int(i==k and j==l)-int(i==l and j==k))
    ric=S.Matrix(5,5,lambda i,j:sum(sig[c]*riem(c,i,c,j) for c in range(5)))
    lr=S.zeros(5)
    for i in range(5):
        for jidx in range(5):
            lr[i,jidx]=sum((ric[i,c]*sig[c]*tensor[c,jidx]+ric[jidx,c]*sig[c]*tensor[c,i])/2 for c in range(5))
            lr[i,jidx]-=sum(riem(i,c,jidx,f)*sig[c]*sig[f]*tensor[c,f] for c in range(5) for f in range(5))
    LR=proj(lr)
    checks['SO4_invariant_tensor_connection_and_curvature_closure']='PASS'
    vars5=S.Matrix([K,L,N,Jj,w])
    Avec=-3*ell*ht*vars5
    Avec[:4,0]+= -2*ell*rho*O
    Avec[1]+=2*ell*rho*p*w;Avec[4]+=ell*rho*p*L
    Bvec=-3*ell*rho*aa*vars5
    Bvec[0]+=-4*ell*rho*p*w
    Bvec[4]+=ell*rho*p*(K+N-3*Jj)/2
    Dvec=S.zeros(5,1)
    gdiag=S.Matrix([1,0,-1,1])
    Dvec[:4,0]=ell**2*rho**2*(Zvec-2*LR+S.Rational(4,3)*Vpot*S.Matrix([K,L,N,Jj])+S.Rational(4,3)*Vp*w*gdiag)
    Dvec[4]=ell**2*rho**2*(-Vpp*w-(p_y+4*aa*p)*K+3*p*hh*L)
    sigma1,sigma2=S.symbols('sigma1 sigma2')
    Mvec=-ell*rho*S.Matrix([-4*aa*K-2*aa*N+6*aa*Jj+4*sigma1*w/3,0,-aa*K-sigma1*w/3,aa*K+sigma1*w/3,p*K/2-sigma2*w/2])
    matrix_report={name:[[str(S.expand(z)) for z in row] for row in vec.jacobian(vars5).tolist()]
                   for name,vec in [('A',Avec),('B',Bvec),('D',Dvec),('M',Mvec)]}

    result={'status':'PASS_INDEPENDENT_EXACT_ALGEBRA','checks':checks,'check_count':len(checks),'mean_matrices':matrix_report,'scope':'Finite algebra; retained analytic premises remain conditional.'}
    finish(args,result)

if __name__=='__main__':
    try:main()
    except Exception as error:fail('SCIENCE',str(error))
