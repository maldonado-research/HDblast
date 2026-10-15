#!/usr/bin/env python3
"""Separate follow-on: exact initial-constraint obstruction and local inverses.

This file does not modify or belong to the frozen PR73 capsule.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import random
import subprocess
import sys

def require(ok,label):
    if not ok:raise ValueError(label)

def eq(a,b,label):require(a==b,label+': '+str(a-b))

def arr(shape):return [arr(shape[1:]) for _ in range(shape[0])] if shape else F(0)

def momentum_check(seed,injection):
    """Direct coordinate expansion of M=div(K-trK gamma)-Pi dphi."""
    n=4;I=range(n);S=range(1,n);rng=random.Random(seed)
    rnd=lambda:F(rng.randint(-7,7),rng.randint(1,6))
    rho=F(7,3);a,u,uy,ut,uyt,v,vy,vt,p,Y=[rnd() for _ in range(10)]
    if not p:p=F(2,3)
    Yi=[F(0)]+[rnd() for _ in S]
    g=arr([n,n]);g[0][0]=F(1)
    for i in S:g[i][i]=rho*rho
    inv=[[1/g[i][i] if i==j else F(0) for j in I] for i in I]
    dg,h,dh=arr([n,n,n]),arr([n,n]),arr([n,n,n])
    K,dK=arr([n,n]),arr([n,n,n])
    h[0][0]=-4*u*Y;K[0][0]=-2*ut*Y/rho
    dh[0][0][0]=-4*uy*Y;dK[0][0][0]=-2*(uyt-a*ut)*Y/rho
    for i in S:
        dg[0][i][i]=2*a*g[i][i]
        h[i][i]=2*u*Y*g[i][i]
        dh[0][i][i]=2*(uy+2*a*u)*Y*g[i][i]
        K[i][i]=ut*Y*g[i][i]/rho
        dK[0][i][i]=(uyt+a*ut)*Y*g[i][i]/rho
    for c in S:
        dh[c][0][0]=-4*u*Yi[c];dK[c][0][0]=-2*ut*Yi[c]/rho
        for i in S:
            dh[c][i][i]=2*u*Yi[c]*g[i][i]
            dK[c][i][i]=ut*Yi[c]*g[i][i]/rho
    i1=[[-sum(inv[i][k]*h[k][l]*inv[l][j] for k in I for l in I) for j in I] for i in I]
    di0,di1=arr([n,n,n]),arr([n,n,n])
    for c,i,j in itertools.product(I,I,I):
        di0[c][i][j]=-sum(inv[i][k]*dg[c][k][l]*inv[l][j] for k in I for l in I)
        di1[c][i][j]=-sum(di0[c][i][k]*h[k][l]*inv[l][j]+inv[i][k]*dh[c][k][l]*inv[l][j]
                         +inv[i][k]*h[k][l]*di0[c][l][j] for k in I for l in I)
    gam0,gam1=arr([n,n,n]),arr([n,n,n])
    for i,j,k in itertools.product(I,I,I):
        gam0[i][j][k]=sum(inv[i][l]*(dg[j][k][l]+dg[k][j][l]-dg[l][j][k]) for l in I)/2
        gam1[i][j][k]=sum(i1[i][l]*(dg[j][k][l]+dg[k][j][l]-dg[l][j][k])
                            +inv[i][l]*(dh[j][k][l]+dh[k][j][l]-dh[l][j][k]) for l in I)/2
    P1,P2,dP2=arr([n,n]),arr([n,n]),arr([n,n,n])
    T1=sum(inv[i][j]*K[j][i] for i in I for j in I)
    T2=sum(i1[i][j]*K[j][i] for i in I for j in I)
    for i,j in itertools.product(I,I):
        P1[i][j]=sum(inv[i][l]*K[l][j] for l in I)-(T1 if i==j else 0)
        P2[i][j]=sum(i1[i][l]*K[l][j] for l in I)-(T2 if i==j else 0)
        for c in I:
            dt2=sum(di1[c][k][l]*K[l][k]+i1[k][l]*dK[c][l][k] for k in I for l in I)
            dP2[c][i][j]=sum(di1[c][i][l]*K[l][j]+i1[i][l]*dK[c][l][j] for l in I)-(dt2 if i==j else 0)
    gradchi=[vy*Y]+[v*Yi[i] for i in S]
    direct=[]
    for i in I:
        direct.append(2*(sum(dP2[j][j][i] for j in I)
            +sum(gam0[j][j][l]*P2[l][i]-gam0[l][j][i]*P2[j][l]
                 +gam1[j][j][l]*P1[l][i]-gam1[l][j][i]*P1[j][l] for j in I for l in I)
            -vt*Y*gradchi[i]/rho))
    G=(12*u*uyt-48*a*u*ut-6*uy*ut-2*vt*vy)/rho
    V=(18*u*ut-v*vt)/rho
    if injection=='drop_radial_momentum_cross':G+=6*uy*ut/rho
    if injection=='half_angular_momentum':V/=2
    eq(direct[0],G*Y*Y,'actual affine-data radial momentum Hessian')
    for i in S:eq(direct[i],2*V*Y*Yi[i],'actual affine-data tangential momentum Hessian')
    # Direct linearized pure-trace momentum completion, source GY²+d(VY²).
    V_y=rnd();q=V/3;q_y=V_y/3;pi2=(G-V_y)/p if p else F(0)
    if p:
        eq(-3*q,-V,'pure-trace angular momentum inverse')
        eq(-3*q_y-p*pi2,-G,'scalar-velocity radial momentum inverse')
    # The actual first momentum identity relates vt to uyt,ut.
    actual_vt=-3*(uyt+2*a*ut)/p
    induced=(6*uy*ut+6*u*uyt+12*a*u*ut+2*u*actual_vt*p)/rho
    eq(induced,6*uy*ut/rho,'affine spacetime versus physical Cauchy source conversion')
    return n+3

def hamiltonian_inverse_check(seed,injection):
    rng=random.Random(seed);rnd=lambda:F(rng.randint(-7,7),rng.randint(1,6))
    rho=F(7,3);a,p,py,f,fy,fyy=[rnd() for _ in range(6)]
    if not p:p=F(2,3)
    U1=py+4*a*p
    zeta=-3*fy/(2*p);zetay=-3*fyy/(2*p)+3*fy*py/(2*p*p)
    if injection=='flip_scalar_compensator':zeta=-zeta;zetay=-zetay
    for lam in (0,8,24):
        dR=-3*fyy-12*a*fy+(2*lam-6)*f/rho**2
        eq(dR-2*p*zetay-2*U1*zeta,2*(lam-3)*f/rho**2,'compact Hamiltonian inverse')
        src=rnd();fc=-rho**2*src/(2*(lam-3))
        if injection=='wrong_sector_denominator':fc=-rho**2*src/(2*(lam+3))
        eq(2*(lam-3)*fc/rho**2,-src,'sector Hamiltonian source cancellation')
    return 6

def scalar_curvature_hessian(g,h,dg,dh,ddg,ddh):
    """Direct inverse/Christoffel/Ricci epsilon expansion in four dimensions."""
    n=4;I=range(n);mm=lambda a,b:[[sum(a[i][k]*b[k][j] for k in I) for j in I] for i in I]
    inv0=[[1/g[i][i] if i==j else F(0) for j in I] for i in I]
    ih=mm(mm(inv0,h),inv0)
    inv=[inv0,[[-v for v in row] for row in ih],mm(mm(ih,h),inv0)]
    gd=[dg,dh,arr([n,n,n])];gdd=[ddg,ddh,arr([n,n,n,n])]
    dinv=arr([3,n,n,n]);gam=arr([3,n,n,n]);dgam=arr([3,n,n,n,n])
    for k,l,a,b in itertools.product(range(3),I,I,I):
        dinv[k][l][a][b]=-sum(inv[r][a][c]*gd[s][l][c][d]*inv[k-r-s][d][b]
            for r in range(k+1) for s in range(k-r+1) for c in I for d in I)
    for k,a,b,c in itertools.product(range(3),I,I,I):
        gam[k][a][b][c]=sum(inv[r][a][d]*(gd[k-r][b][c][d]+gd[k-r][c][b][d]-gd[k-r][d][b][c])
            for r in range(k+1) for d in I)/2
        for l in I:
            dgam[k][l][a][b][c]=sum(
                dinv[r][l][a][d]*(gd[k-r][b][c][d]+gd[k-r][c][b][d]-gd[k-r][d][b][c])
                +inv[r][a][d]*(gdd[k-r][l][b][c][d]+gdd[k-r][l][c][b][d]-gdd[k-r][l][d][b][c])
                for r in range(k+1) for d in I)/2
    Ric=arr([3,n,n])
    for k,a,b in itertools.product(range(3),I,I):
        Ric[k][a][b]=sum(dgam[k][c][c][a][b]-dgam[k][b][c][a][c] for c in I)
        Ric[k][a][b]+=sum(gam[r][c][c][d]*gam[k-r][d][a][b]-gam[r][c][b][d]*gam[k-r][d][a][c]
                         for r in range(k+1) for c in I for d in I)
    R0=sum(inv[0][a][b]*Ric[0][a][b] for a in I for b in I)
    R2=2*sum(inv[r][a][b]*Ric[2-r][a][b] for r in range(3) for a in I for b in I)
    return R0,R2,inv

def hamiltonian_source_check(seed,injection):
    rng=random.Random(seed);rnd=lambda:F(rng.randint(-7,7),rng.randint(1,6))
    n=4;I=range(n);S=range(1,n);rho=F(7,3)
    a,ay,u,uy,uyy,ut,v,vy,vt,p,U2,Y=[rnd() for _ in range(12)]
    Ys=[F(0)]+[rnd() for _ in S];Yss=arr([n,n])
    for i in S:
        for j in range(i,n):Yss[i][j]=Yss[j][i]=rnd()
    Yss[3][3]=-8*Y-Yss[1][1]-Yss[2][2]
    g=arr([n,n]);g[0][0]=F(1)
    for i in S:g[i][i]=rho*rho
    dg,ddg=arr([n,n,n]),arr([n,n,n,n])
    for i in S:
        dg[0][i][i]=2*a*g[i][i];ddg[0][0][i][i]=2*(ay+2*a*a)*g[i][i]
    delta=lambda i,j:F(i==j)
    for k,l,i,j in itertools.product(S,S,S,S):
        ddg[k][l][i][j]=-rho*rho*(2*delta(i,j)*delta(k,l)-delta(i,l)*delta(k,j)-delta(i,k)*delta(l,j))/3
    psi=u*Y;ps1=[uy*Y]+[u*Ys[i] for i in S];ps2=arr([n,n])
    ps2[0][0]=uyy*Y
    for i in S:ps2[0][i]=ps2[i][0]=uy*Ys[i]
    for i,j in itertools.product(S,S):ps2[i][j]=u*Yss[i][j]
    h,dh,ddh=arr([n,n]),arr([n,n,n]),arr([n,n,n,n])
    h[0][0]=-4*psi
    for k in I:
        dh[k][0][0]=-4*ps1[k]
        for l in I:ddh[k][l][0][0]=-4*ps2[k][l]
    for i,j in itertools.product(S,S):
        h[i][j]=2*psi*g[i][j]
        for k in I:
            dh[k][i][j]=2*(ps1[k]*g[i][j]+psi*dg[k][i][j])
            for l in I:ddh[k][l][i][j]=2*(ps2[k][l]*g[i][j]+ps1[k]*dg[l][i][j]+ps1[l]*dg[k][i][j]+psi*ddg[k][l][i][j])
    R0,R2,inv=scalar_curvature_hessian(g,h,dg,dh,ddg,ddh)
    eq(R0,6/rho**2-6*ay-12*a*a,'normal-coordinate sphere background curvature')
    G=sum(Ys[i]**2 for i in S)
    A=(-24*uy**2-288*a*u*uy-24*u*uyy-(192*ay+384*a*a+32*p*p)*u*u
       -(336*u*u+12*ut*ut+2*vt*vt)/rho**2-2*vy*vy-16*p*u*vy-2*U2*v*v)
    D=(36*u*u-2*v*v)/rho**2
    ph0=[p,F(0),F(0),F(0)];ph1=[vy*Y]+[v*Ys[i] for i in S]
    gradient2=2*sum(inv[2][i][j]*ph0[i]*ph0[j]+inv[1][i][j]*(ph0[i]*ph1[j]+ph1[i]*ph0[j])
                          +inv[0][i][j]*ph1[i]*ph1[j] for i in I for j in I)
    direct=R2-12*ut*ut*Y*Y/rho**2-2*vt*vt*Y*Y/rho**2-gradient2-2*U2*v*v*Y*Y
    candidate=A*Y*Y+D*G
    if injection=='half_Hamiltonian_source':candidate/=2
    if injection=='drop_angular_gradient_source':candidate-=D*G
    eq(candidate,direct,'full actual-state Hamiltonian source formula')
    # Scalar harmonic projections tested as independent polynomial identities.
    # N=n1²+n2²; G=4(N-Y²), Y²=H4+Z2/2+1/6.
    N=rnd();Z=N-F(1,2);H4=Y*Y-N/2+F(1,12)
    S0=A/6+4*D/3;S2=A/2+2*D;S4=A-4*D
    eq(S0+S2*Z+S4*H4,A*Y*Y+4*D*(N-Y*Y),'all three Hamiltonian sector projections')
    return 3

def witness(background,injection):
    b=background['background_intervals'];c=background['reconstruction_coefficients']
    ab=lambda q:max(abs(F(x)) for x in q['exact'])
    rr=[F(x) for x in b['rho']['exact']]
    require(rr==[129,131],'witness radius certificate mismatch')
    Braw=[F(x) for x in c['v_wz_coefficient']['exact']]
    require(max(Braw)<0,'actual B interval must stay negative')
    B=[abs(x) for x in Braw]
    require(min(B)>F(19,2000) and max(B)<F(3,250),'B witness box')
    for key,cap in [('v_wz_coefficient_z',F(3,10)),('alpha_log_z',270),('alpha_log_zz',24000),
                    ('u_w_coefficient',F(1,10000)),('u_w_coefficient_z',F(3,200)),('u_w_coefficient_zz',5)]:
        require(ab(c[key])<cap,'witness coefficient cap: '+key)
    for key in ('a','p'):require(ab(b[key])<1,'witness state cap: '+key)
    require(ab(b['a_z'])/129<1 and ab(b['U_2'])<4,'witness background derivative caps')
    ell=F(1,1000000)
    h=F(1,4);hp=F(0);hpp=F(-2)
    exponent=4-1/h;exponentp=hp/h**2;exponentpp=hpp/h**2-2*hp**2/h**3
    eq(exponent,F(0),'packet maximum exponent')
    eq(exponentp,F(0),'packet maximum first derivative')
    eq(exponentpp+exponentp**2,F(-32),'packet maximum scaled second derivative')
    bump_second=F(32)
    if injection=='omit_bump_second_jet':bump_second=F(0)
    err=(F(3,10)*270+F(3,250)*24000)/129
    vylo=bump_second*F(19,2000)/(131*ell**2)-err
    vyhi=bump_second*F(3,250)/(129*ell**2)+err
    require(vylo>2300000000 and vyhi<3000000000,'nonzero initial scalar-gradient witness')
    require(F(3,200)/129<1,'packet-maximum u_y coarse cap')
    require((5+32*F(1,10000)/ell**2)/129**2+F(3,200)/129<1000000,
            'packet-maximum u_yy coarse cap')
    require(F(1,10000)<1 and F(3,250)*270<4 and F(3,250)*271<4,
            'packet-maximum velocity/scalar coarse caps')
    # Conservative derivative bounds at q maximum. Exponential e^-z is <=1.
    u=F(1,10000);uy=F(1);uyy=F(1000000);ut=F(1);v=F(4);vt=F(4)
    remainder=(24*uy**2+288*u*uy+24*u*uyy+(192+384+32)*u**2
               +(336*u**2+12*ut**2+2*vt**2)/129**2+16*u*3000000000+2*4*v**2)
    require(remainder<5000000,'Hamiltonian witness remainder bound')
    upper=-2*F(2300000000)**2+5000000
    require(upper<-10**19,'Hamiltonian witness must exclude zero second data')
    return {'point':{'t':'0','z':'3/2000000','n':['1','0','0','0']},
            'exact_initial_master_at_point':{'w':'1','w_z':'0','w_zz':str(-32/ell**2)},
            'actual_abs_v_y_strict_lower':'2300000000','actual_abs_v_y_strict_upper':'3000000000',
            'other_Hessian_terms_abs_upper':str(remainder),'Hamiltonian_Hessian_strict_upper':str(upper)}

def main():
    ap=argparse.ArgumentParser();root=Path(__file__).resolve().parent
    ap.add_argument('--background',type=Path,default=root/'inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json')
    ap.add_argument('--output',type=Path)
    mutations={'drop_radial_momentum_cross':'actual affine-data radial momentum Hessian',
               'half_angular_momentum':'actual affine-data tangential momentum Hessian',
               'flip_scalar_compensator':'compact Hamiltonian inverse',
               'wrong_sector_denominator':'sector Hamiltonian source cancellation',
               'half_Hamiltonian_source':'full actual-state Hamiltonian source formula',
               'drop_angular_gradient_source':'full actual-state Hamiltonian source formula',
               'omit_bump_second_jet':'nonzero initial scalar-gradient witness'}
    ap.add_argument('--inject-error',choices=['none']+list(mutations),default='none')
    ap.add_argument('--run-controls',action='store_true')
    a=ap.parse_args();raw=a.background.read_bytes()
    require(hashlib.sha256(raw).hexdigest()=='fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0','retained actual background pin mismatch')
    checks=sum(momentum_check(seed,a.inject_error)+hamiltonian_inverse_check(seed,a.inject_error)
               +hamiltonian_source_check(seed,a.inject_error) for seed in (271,1927,6469))
    report={'status':'PASS_SEPARATE_EXPLICIT_INITIAL_CONSTRAINT_COMPLETION',
            'finite_exact_checks':checks,'background_sha256':hashlib.sha256(raw).hexdigest(),
            'physical_Cauchy_affine_variables':['spatial_metric gamma','covariant K=+Lie_n(gamma)/2','scalar phi','unit-normal scalar velocity Pi'],
            'angular_sectors':{'degrees':[0,2,4],'eigenvalues':[0,8,24],'Hamiltonian_inverse_denominators':[-6,10,42]},
            'witness':witness(json.loads(raw),a.inject_error),
            'scope':{'separate_follow_on':True,'frozen_capsule_modified':False,
                     'zero_second_initial_data_compatible':False,'explicit_second_constraint_jet':True,
                     'chosen_scalar_second_data_zero':False,'chosen_scalar_velocity_second_data_zero':False,
                     'second_initial_norms_numerically_enclosed':False,'gauge_completion_computed':False,
                     'U2_evolution_computed':False}}
    if a.run_controls:
        import tempfile
        controls=[]
        with tempfile.TemporaryDirectory(prefix='initial-jet-controls-') as td:
            for mutation,reason in mutations.items():
                output=Path(td)/(mutation+'.json');cmd=[sys.executable]
                if sys.flags.optimize:cmd.append('-O')
                cmd += [str(Path(__file__).resolve()),'--background',str(a.background),
                        '--inject-error',mutation,'--output',str(output)]
                result=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
                require(result.returncode!=0 and ('ValueError: '+reason) in result.stderr and not output.exists(),
                        'control did not reject for the intended reason: '+mutation)
                controls.append({'mutation':mutation,'rejection_prefix':reason,'returncode':result.returncode,
                                 'no_output_report':True})
        report['controls']={'status':'PASS','optimized':bool(sys.flags.optimize),'cases':controls}
    text=json.dumps(report,indent=2)+'\n'
    if a.output:a.output.write_text(text)
    print(text,end='')

if __name__=='__main__':main()
