#!/usr/bin/env python3
"""All-point algebra independent of numerical curvature fixtures."""
import hashlib
import argparse
import itertools
import json
from pathlib import Path
import sympy as s


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def zero(expr,msg):
    need(s.factor(expr)==0,msg)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',required=True,type=Path,help='accepted for uniform CLI; no production sources are read')
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(),'FRESH_OUTPUT')
    need(Path(__file__).resolve().parent not in args.output.resolve().parents,'OUTPUT_OUTSIDE_CODE')
    A,B,C,D,E,L,M,My,Ly,Mt,Lt,Ct,Ey,rho,R,theta,a=s.symbols('A B C D E L M My Ly Mt Lt Ct Ey rho R theta a')
    H=theta/rho
    mr=My/2-C/R;mt=Mt/(2*rho)-E/R
    lr=Ly/2-a*A/2+H*B;betar=mr+a*M;betat=mt+H*M
    eta=R*(Ct-theta*C)/(2*rho)-R*Ey/2-B/2
    eta_r=R*Ey/2-R*(Ct-theta*C)/(2*rho)-B/2
    radial_four=(-A*betar/2+B*betat-C*(L-2*M)/R
                 -(L/2-M)*mr+M*(-lr+mr)+mt*eta_r)
    Fy=((-A+4*M-L)*My-2*M*Ly)/4+C*(A-L)/(2*R)+(E/R-Mt/(2*rho))*eta
    zero(radial_four-Fy,'RADIAL_FOUR_CONTRIBUTION_IDENTITY')
    Mtt,Mty,Myy,Et,Cy=s.symbols('Mtt Mty Myy Et Cy')
    mt_t=Mtt/(2*rho)-Et/R+theta*E/R
    mt_y=Mty/(2*rho)-a*Mt/(2*rho)-Ey/R+a*E/R
    mr_y=Myy/2-Cy/R+a*C/R
    g=-L/(2*R**2)+H*mt-a*mr
    et=-mt_t/rho-2*H*mt-D/(2*R**2)+a*mr
    ey=-mr_y-2*a*mr-A/(2*R**2)
    k=eta/R**2-mt_y-a*mt-H*mr
    At=-(Mtt+theta*Mt)/(2*rho**2)+Et/(rho*R)-(D+L)/(2*R**2)
    Ay=-(Myy+a*My)/2+Cy/R-theta*Mt/(2*rho**2)+theta*E/(rho*R)+(L-A)/(2*R**2)
    K=-(Mty+theta*My)/(2*rho)+(Ct+theta*C)/(2*rho*R)+Ey/(2*R)-B/(2*R**2)
    for z,name in ((g+et-At,'AT'),(ey-g-Ay,'AY'),(k-K,'K')):zero(z,'FIRST_WEYL_INVARIANT_'+name)
    # Arbitrary pointwise scalar curvature coefficients and Hessian/gradient.
    gg,ee,kk,p,q=s.symbols('g e k p q')
    v=s.symbols('v0:3');h=s.MutableDenseMatrix(3,3,[0]*9)
    for i in range(3):
        for j in range(i,3):h[i,j]=h[j,i]=s.Symbol('h'+str(i)+str(j))
    def ep(i,j,k):return s.LeviCivita(i,j,k)
    aa=s.zeros(3);direct=s.zeros(3)
    for i,j in itertools.product(range(3),repeat=2):
        for c,d in itertools.product(range(3),repeat=2):
            aa[i,j]+=(ep(i,c,d)*v[c]*h[j,d]+ep(j,c,d)*v[c]*h[i,d])/2
            for b in range(3):
                # Observer angular variation contracts the screen Riemann.
                cur=gg*(int(j==c)*h[b,d]+int(b==d)*h[j,c]-int(j==d)*h[b,c]-int(b==c)*h[j,d])
                direct[i,j]+=ep(i,c,d)*p*v[b]*cur/2
            direct[i,j]+=ep(i,c,d)*(-p*ee*v[c]*h[j,d]+p*ee*v[d]*h[j,c]-q*kk*v[c]*h[j,d]+q*kk*v[d]*h[j,c])/2
    direct=(direct+direct.T)/2
    for i,j in itertools.product(range(3),repeat=2):zero(direct[i,j]+(p*(gg+ee)+q*kk)*aa[i,j],'LORENTZ_ALL_POINT')
    # Hyperbolic reduction with c=(z+z^-1)/2,s=(z-z^-1)/2.
    z,alpha,beta=s.symbols('z alpha beta',nonzero=True);co=(z+1/z)/2;si=(z-1/z)/2
    pp=(rho*co*alpha-si*beta)/R;qq=(-rho*si*alpha+co*beta)/R
    Echi=co**2*At+si**2*Ay+2*co*si*K
    Kchi=co*si*(At+Ay)+(co**2+si**2)*K
    zero(pp*Echi+qq*Kchi-(rho*alpha*(co*At+si*K)+beta*(si*Ay+co*K))/R,'GAUSSIAN_TILT_REDUCTION')
    # Coordinate connection variation independently yields exact source slots.
    Dt,Dy,Atm,Aym,Bt,By=s.symbols('Dt Dy Atmetric Aymetric Bt By')
    dGttt=-Dt/2+a*rho*B;dGtty=-Dy/2
    dGtyy=-By/rho-a*B/rho+Atm/(2*rho**2)
    dGytt=rho*Bt-rho**2*(a*(A+D)+Dy/2)
    dGyty=Atm/2-a*rho*B;dGyyy=Aym/2
    SP=a*B-(co**2*Dt-si**2*Atm)/(2*rho)-co*si*Dy-si**2*By
    SQ=co**2*(Bt/rho-a*(A+D)-Dy/2)+co*si*(Atm/rho-2*a*B)+si**2*Aym/2
    zero(rho*(dGttt*(co/rho)**2+2*dGtty*co/rho*si+dGtyy*si**2)-SP,'GEODESIC_SP')
    zero(dGytt*(co/rho)**2+2*dGyty*co/rho*si+dGyyy*si**2-SQ,'GEODESIC_SQ')
    report={'status':'PASS_ALL_POINT_SYMBOLIC_IDENTITIES','checks':['radial_four_contributions','first_Weyl_At_Ay_K','full_Lorentz_tilt','Gaussian_alpha_beta_reduction','geodesic_SP_SQ'], 'assert_statements_used':False,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    report['arithmetic']='Exact Sympy rational algebra; no production scientific modules or source values read'
    with args.output.open('x')as handle:
        json.dump(report,handle,sort_keys=True,indent=2);handle.write('\n')
    print(json.dumps({'status':report['status'],'checks':len(report['checks'])},sort_keys=True))


if __name__=='__main__':main()
