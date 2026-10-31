#!/usr/bin/env python3
"""Independent local-coordinate Weyl two-jet and Lorentz-tilt audit.

No production scientific modules are imported. Angular derivatives are exact
rational stereographic derivatives; parameter differentiation uses degree-two
polynomial convolution, not finite differences or symbolic formula fitting.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parent


def need(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def angular_jets(point):
    w = s.symbols('w0:3')
    den = 1 + sum(z*z for z in w)
    lam = 2/den
    ns = [(1-sum(z*z for z in w))/den] + [2*z/den for z in w]
    Y = ns[0]**2-ns[1]**2
    kap = s.eye(3)*lam**2
    vv = [s.diff(Y,z) for z in w]
    ll = [s.diff(lam,z)/lam for z in w]
    hh = s.Matrix(3,3,lambda a,b:s.diff(Y,w[a],w[b])-
                  ll[a]*vv[b]-ll[b]*vv[a]+int(a==b)*sum(ll[c]*vv[c] for c in range(3)))
    at = dict(zip(w,point))
    def jet(expr):
        out = np.zeros((4,4))
        out[0,0] = float(expr.subs(at))
        for a in range(3):
            out[a+1,0] = out[0,a+1] = float(s.diff(expr,w[a]).subs(at))
            for b in range(3):
                out[a+1,b+1] = float(s.diff(expr,w[a],w[b]).subs(at))
        return out
    return {'Y':jet(Y), 'v':[jet(v) for v in vv],
            'k':[[jet(kap[a,b]) for b in range(3)] for a in range(3)],
            'h':[[jet(hh[a,b]) for b in range(3)] for a in range(3)],
            'lam':float(lam.subs(at))}


def outer_jet(f, ang):
    """f=(f,ft,fy,ftt,fty,fyy), times angular function two-jet."""
    g = f[0]*ang[0,0]
    d = np.zeros(5); dd = np.zeros((5,5))
    d[:2] = np.array(f[1:3])*ang[0,0]
    dd[0,0] = f[3]*ang[0,0]; dd[0,1] = dd[1,0] = f[4]*ang[0,0]
    dd[1,1] = f[5]*ang[0,0]
    for a in range(3):
        d[a+2] = f[0]*ang[a+1,0]
        dd[0,a+2] = dd[a+2,0] = f[1]*ang[a+1,0]
        dd[1,a+2] = dd[a+2,1] = f[2]*ang[a+1,0]
        for b in range(3):
            dd[a+2,b+2] = f[0]*ang[a+1,b+1]
    return g,d,dd


def factor_jet(rhopow,cpow,m,b):
    rho,c,ct,a = (b[k] for k in ('rho','c','ct','a'))
    v,vt,vy = m[:3]
    vtt,vty,vyy=m[3:] if len(m)==6 else (0,0,0)
    p = rho**rhopow*c**cpow
    pt = cpow*ct/c
    py = rhopow*a
    ptt = cpow+cpow*(cpow-1)*(ct/c)**2
    return (p*v,p*(vt+pt*v),p*(vy+py*v),
            p*(vtt+2*pt*vt+ptt*v),p*(vty+pt*vy+py*vt+pt*py*v),
            p*(vyy+2*py*vy+py*py*v))


def metric_jets(m,b,ang):
    g = np.zeros((3,5,5)); d = np.zeros((3,5,5,5)); dd = np.zeros((3,5,5,5,5))
    one = np.zeros((4,4));one[0,0]=1
    def add(deg,i,j,rp,cp,val,aa):
        z,x,y=outer_jet(factor_jet(rp,cp,val,b),aa)
        g[deg,i,j]+=z;d[deg,:,i,j]+=x;dd[deg,:,:,i,j]+=y
        if i!=j:
            g[deg,j,i]+=z;d[deg,:,j,i]+=x;dd[deg,:,:,j,i]+=y
    add(0,0,0,2,0,(-1,0,0),one)
    add(0,1,1,0,0,(1,0,0),one)
    add(1,0,0,2,0,m['D'],ang['Y'])
    add(1,0,1,1,0,m['B'],ang['Y'])
    add(1,1,1,0,0,m['A'],ang['Y'])
    for i in range(3):
        add(1,0,i+2,2,1,m['E'],ang['v'][i])
        add(1,1,i+2,1,1,m['C'],ang['v'][i])
        for j in range(i,3):
            add(0,i+2,j+2,2,2,(1,0,0),ang['k'][i][j])
            # Product Y*k is diagonal, with derivatives computed directly
            # by the angular jet product rule below.
            yk=angular_product(ang['Y'],ang['k'][i][j])
            add(1,i+2,j+2,2,2,m['L'],yk)
            add(1,i+2,j+2,2,2,m['M'],ang['h'][i][j])
            if 'N' in m:
                add(1,i+2,j+2,2,2,m['N'],angular_product(ang['Y'],ang['h'][i][j]))
    return g,d,dd


def angular_product(a,b):
    z=np.zeros((4,4));z[0,0]=a[0,0]*b[0,0]
    for q in range(1,4):
        z[q,0]=z[0,q]=a[q,0]*b[0,0]+a[0,0]*b[q,0]
        for r in range(1,4):
            z[q,r]=a[q,r]*b[0,0]+a[q,0]*b[r,0]+a[r,0]*b[q,0]+a[0,0]*b[q,r]
    return z


def conv_einsum(subs,*args):
    out=None
    for degs in itertools.product(range(3),repeat=len(args)):
        deg=sum(degs)
        if deg>2:continue
        term=np.einsum(subs,*(a[k] for a,k in zip(args,degs)),optimize=True)
        if out is None:out=np.zeros((3,)+np.shape(term))
        out[deg]+=term
    return out


def weyl(g,d,dd):
    gi=np.zeros_like(g);gi[0]=np.linalg.inv(g[0]);gi[1]=-gi[0]@g[1]@gi[0]
    gi[2]=gi[0]@g[1]@gi[0]@g[1]@gi[0]
    gam=np.zeros((3,5,5,5)) # first-kind gamma_a,bc
    for a,b,c in itertools.product(range(5),repeat=3):
        gam[:,a,b,c]=(d[:,b,a,c]+d[:,c,a,b]-d[:,a,b,c])/2
    rr=np.zeros((3,5,5,5,5))
    for a,b,c,e in itertools.product(range(5),repeat=4):
        rr[:,a,b,c,e]=(dd[:,b,c,a,e]+dd[:,a,e,b,c]-dd[:,b,e,a,c]-dd[:,a,c,b,e])/2
    rr+=conv_einsum('ef,ead,fbc->abcd',gi,gam,gam)
    rr-=conv_einsum('ef,eac,fbd->abcd',gi,gam,gam)
    ric=conv_einsum('cd,cadb->ab',gi,rr)
    scalar=conv_einsum('ab,ab->',gi,ric)
    sch=(ric-conv_einsum(',ab->ab',scalar,g)/8)/3
    cc=rr-conv_einsum('ac,bd->abcd',g,sch)+conv_einsum('ad,bc->abcd',g,sch)
    cc+=conv_einsum('bc,ad->abcd',g,sch)-conv_einsum('bd,ac->abcd',g,sch)
    need(np.max(np.abs(cc[0]))<2e-12,'BACKGROUND_WEYL_ZERO')
    return cc


def epsilon():
    ee=np.zeros((3,3,3))
    for p in itertools.permutations(range(3)):
        ee[p]=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
    return ee


EPS=epsilon()


def H(cc,U,screen):
    x=conv_einsum('ai,j,bk,cl,ijkl->abc',screen,U,screen,screen,cc)
    h=np.einsum('acd,kbcd->kab',EPS,x)/2
    return (h+h.swapaxes(1,2))/2


def frames(m,b,ang,ch,p=0.,q=0.):
    rho,c=(b[k] for k in ('rho','c'));R=rho*c;lam=ang['lam']
    eu=np.eye(5)[0]/rho;en=np.eye(5)[1];ea=np.eye(5)[2:]/(R*lam)
    vv=np.array([z[0,0] for z in ang['v']])/lam
    hh=np.array([[z[0,0] for z in row] for row in ang['h']])/lam**2
    Y=ang['Y'][0,0]
    A,B,C,D,E,L,M=(m[k][0] for k in ('A','B','C','D','E','L','M'))
    Ur=np.zeros((3,5));nr=np.zeros_like(Ur);screen=np.zeros((3,3,5))
    Ur[0]=eu;nr[0]=en;screen[0]=ea
    Ur[1]=D*Y*eu/2-B*Y*en-E*np.einsum('a,ai->i',vv,ea)
    nr[1]=-A*Y*en/2-C*np.einsum('a,ai->i',vv,ea)
    screen[1]=-np.einsum('ab,bi->ai',L*Y*np.eye(3)+M*hh,ea)/2
    U=np.cosh(ch)*Ur+np.sinh(ch)*nr
    n=np.sinh(ch)*Ur+np.cosh(ch)*nr
    # Relative first-order Lorentz transformation, both screen tilts retained.
    U[1]+=p*np.einsum('a,ai->i',vv,ea)
    screen[1]+=np.einsum('a,i->ai',vv,p*U[0]+q*n[0])
    aa=(np.einsum('acd,c,bd->ab',EPS,vv,hh)+np.einsum('bcd,c,ad->ab',EPS,vv,hh))/2
    return U,screen,aa


def oldF(m,b):
    rho,c,ct=(b[k] for k in ('rho','c','ct'));R=rho*c
    B,C,D,E,L,M=(m[k][0] for k in ('B','C','D','E','L','M'))
    Mt,Lt,Ct,Ey,My=m['M'][1],m['L'][1],m['C'][1],m['E'][2],m['M'][2]
    return ((D+4*M-L)*Mt-2*M*Lt)/(4*rho)-E*(D+L)/(2*R)+(C/R-My/2)*(R*(Ct-ct/c*C)/(2*rho)-R*Ey/2+B/2)


def radialF(m,b):
    rho,c,ct=(b[k] for k in ('rho','c','ct'));R=rho*c
    A,B,C,E,L,M=(m[k][0] for k in ('A','B','C','E','L','M'))
    Mt,My,Ly,Ct,Ey=m['M'][1],m['M'][2],m['L'][2],m['C'][1],m['E'][2]
    eta=R*(Ct-ct/c*C)/(2*rho)-R*Ey/2-B/2
    return ((-A+4*M-L)*My-2*M*Ly)/4+C*(A-L)/(2*R)+(E/R-Mt/(2*rho))*eta


def expected_invariants(m,b):
    rho,c,ct,a=(b[k] for k in ('rho','c','ct','a'));R=rho*c;th=ct/c
    A,B,C,D,E,L=(m[k][0] for k in ('A','B','C','D','E','L'))
    Mt,My,Mtt,Mty,Myy=m['M'][1:]
    Ct,Cy,Et,Ey=m['C'][1],m['C'][2],m['E'][1],m['E'][2]
    At=-(Mtt+th*Mt)/(2*rho*rho)+Et/(rho*R)-(D+L)/(2*R*R)
    Ay=-(Myy+a*My)/2+Cy/R-th*Mt/(2*rho*rho)+th*E/(rho*R)+(L-A)/(2*R*R)
    K=-(Mty+th*My)/(2*rho)+(Ct+th*C)/(2*rho*R)+Ey/(2*R)-B/(2*R*R)
    return np.array([At,Ay,K])


def curvature_invariants(cc,b,ang):
    rho,c=(b[k] for k in ('rho','c'));R=rho*c;lam=ang['lam']
    eu=np.eye(5)[0]/rho;en=np.eye(5)[1];ea=np.eye(5)[2:]/(R*lam)
    hh=np.array([[z[0,0] for z in row] for row in ang['h']])/lam**2
    ht=hh-np.trace(hh)*np.eye(3)/3
    d=np.eye(3)
    hk=(np.einsum('ac,bd->abcd',d,ht)+np.einsum('bd,ac->abcd',d,ht)
        -np.einsum('ad,bc->abcd',d,ht)-np.einsum('bc,ad->abcd',d,ht))
    angular=np.einsum('ai,bj,ck,dl,ijkl->abcd',ea,ea,ea,ea,cc[1],optimize=True)
    g=np.sum(angular*hk)/np.sum(hk*hk)
    e0=np.einsum('ai,j,bk,l,ijkl->ab',ea,eu,ea,eu,cc[1],optimize=True)
    ey=np.einsum('ai,j,bk,l,ijkl->ab',ea,en,ea,en,cc[1],optimize=True)
    em=np.einsum('ai,j,bk,l,ijkl->ab',ea,eu,ea,en,cc[1],optimize=True)
    return np.array([g+np.sum(e0*ht)/np.sum(ht*ht),np.sum(ey*ht)/np.sum(ht*ht)-g,np.sum(em*ht)/np.sum(ht*ht)])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',required=True,type=Path,help='accepted for uniform CLI; no production sources are read')
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--raw-diagnostic-output',type=Path)
    args=parser.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(),'FRESH_OUTPUT')
    need(ROOT not in args.output.resolve().parents,'OUTPUT_OUTSIDE_CODE')
    points=[(s.Rational(1,3),s.Rational(1,4),s.Rational(1,5)),
            (s.Rational(-1,5),s.Rational(2,7),s.Rational(1,4))]
    report=[];linear=[]
    for j,point in enumerate(points):
        ang=angular_jets(point)
        for i in range(1,5):
            m={k:tuple((((i+3)*(n+2)*(v+3)*7)%37-18)/(v+n+3) for v in range(6)) for n,k in enumerate('ABC DELM'.replace(' ',''))}
            b={'rho':1.5,'c':5/3,'ct':4/3,'a':(i-2)/7}
            cc=weyl(*metric_jets(m,b,ang))
            U,S,A=frames(m,b,ang,0)
            hh=H(cc,U,S)
            need(np.max(np.abs(hh[1]))<2e-11,'FIRST_TIME_TENSOR_ZERO')
            exact=4*oldF(m,b)/(b['rho']*b['c'])
            got=4*np.sum(hh[2]*A)/np.sum(A*A)
            need(abs(got-exact)<2e-10*(1+abs(exact)),'OLD_COORDINATE_CONVENTION')
            Ur,Sr,A=frames(m,b,ang,0)
            # Pure radial frame through first order: obtain from fixed boosts.
            Up,Sp,_=frames(m,b,ang,.7)
            hr=(H(cc,Up,Sp)-np.cosh(.7)*hh)/np.sinh(.7)
            need(np.max(np.abs(hr[1]))<2e-11,'FIRST_RADIAL_TENSOR_ZERO')
            radial=4*np.sum(hr[2]*A)/np.sum(A*A)
            expect_radial=4*radialF(m,b)/(b['rho']*b['c'])
            need(abs(radial-expect_radial)<2e-10*(1+abs(radial)),'RADIAL_COORDINATE_QUADRATIC')
            # Check all-point span at the fixture independently of projection.
            need(np.max(np.abs(hr[2]-radial*A/4))<2e-10,'RADIAL_PURE_A')
            p,q=.23,-.41
            Ut,St,_=frames(m,b,ang,.7,p,q)
            tilt=4*np.sum((H(cc,Ut,St)[2]-H(cc,Up,Sp)[2])*A)/np.sum(A*A)
            inv=curvature_invariants(cc,b,ang)
            expected_inv=expected_invariants(m,b)
            need(np.max(np.abs(inv-expected_inv))<2e-10,'WEYL_HESSIAN_INVARIANTS')
            ch=.7;co=np.cosh(ch);si=np.sinh(ch)
            At,Ay,K=inv
            pred_tilt=-4*(p*(co*co*At+si*si*Ay+2*co*si*K)+q*(co*si*(At+Ay)+(co*co+si*si)*K))
            need(abs(tilt-pred_tilt)<2e-10*(1+abs(tilt)),'LORENTZ_TILT_MAP')
            report.append({'point':j,'case':i,'time_quadratic':got,'radial_quadratic':radial,'tilt_second_scalar':tilt,'invariants_At_Ay_K':list(inv),'m':m,'b':b})
        m={k:(0,0,0,0,0,0) for k in 'ABCDELM'}
        m['N']=(2,3,-5,7,-11,13)
        b={'rho':1.5,'c':5/3,'ct':4/3,'a':3/7}
        cc=weyl(*metric_jets(m,b,ang));R=b['rho']*b['c']
        for ch in (0.,.7,-.9):
            U,S,A=frames(m,b,ang,ch)
            h1=H(cc,U,S)[1]
            expected=(np.cosh(ch)*m['N'][1]/b['rho']+np.sinh(ch)*m['N'][2])/(2*R)*A
            need(np.max(np.abs(h1-expected))<2e-11,'LINEAR_NORMAL_DERIVATIVE_SIGN_AND_NO_AU')
            linear.append({'point':j,'chi':ch,'max_tensor_residual':float(np.max(np.abs(h1-expected)))})
    rejected={k:0. for k in ('drop_radial_quadratic','radial_lapse_wrong_sign','radial_eta_B_wrong_sign','drop_time_epsilon_volume','drop_mixed_curvature_in_tilt','reverse_radial_screen_tilt','half_second_derivative_tilt')}
    for case in report:
        m,b=case['m'],case['b'];rho,c,ct,a=(b[k] for k in ('rho','c','ct','a'));R=rho*c
        A,B,C,D,E,L,M=(m[k][0] for k in 'ABCDELM')
        Mt,My,Lt,Ly=m['M'][1],m['M'][2],m['L'][1],m['L'][2]
        expected=case['radial_quadratic']
        tests={'drop_radial_quadratic':0.,
               'radial_lapse_wrong_sign':4*(radialF(m,b)+A*My/2)/R,
               'radial_eta_B_wrong_sign':4*(radialF(m,b)+(E/R-Mt/(2*rho))*B)/R}
        ell=Lt/(2*rho)+ct/c/rho*D/2-a*B;mm=Mt/(2*rho)-E/R
        tests['drop_time_epsilon_volume']=case['time_quadratic']-4*M*(-ell+mm)/R
        co,si=np.cosh(.7),np.sinh(.7);p,q=.23,-.41
        At,Ay,K=case['invariants_At_Ay_K']
        ep=co*co*At+si*si*Ay+2*co*si*K
        kp=co*si*(At+Ay)+(co*co+si*si)*K
        tests['drop_mixed_curvature_in_tilt']=-4*(p*(co*co*At+si*si*Ay)+q*co*si*(At+Ay))
        tests['reverse_radial_screen_tilt']=-4*(p*ep-q*kp)
        tests['half_second_derivative_tilt']=case['tilt_second_scalar']/2
        for name,value in tests.items():
            target=case['time_quadratic'] if name=='drop_time_epsilon_volume' else case['tilt_second_scalar'] if name.endswith('tilt') else expected
            rejected[name]=max(rejected[name],abs(value-target))
    for name,gap in rejected.items():need(gap>1e-3,'MUTANT_NOT_REJECTED:'+name)
    out={'status':'PASS_SOURCE_FREE_COORDINATE_WEYL_DIAGNOSTIC',
         'arithmetic':'NumPy floating arrays with Sympy-derived rational stereographic angular inputs; diagnostic checks only',
         'scope':'Independent finite-coordinate diagnostic; no exact identity or validated enclosure follows from floating residual thresholds',
         'production_source_reads':0,'parameter_polynomial_degree':2,
         'scalar_fixtures':[{'point':v['point'],'case':v['case'],'passed':True}for v in report],
         'linear_normal_derivative_fixture_count':len(linear),
         'floating_tolerances':{'first_tensor_zero':'2e-11','linear_TT_tensor':'2e-11','curvature_absolute_or_relative':'2e-10','wrong_formula_gap_minimum':'1e-3'},
         'wrong_formula_controls':{name:{'rejected':True,'minimum_floating_gap_exceeded':'1e-3'}for name in rejected},
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.raw_diagnostic_output:
        need(not args.raw_diagnostic_output.exists() and ROOT not in args.raw_diagnostic_output.resolve().parents,'FRESH_EXTERNAL_RAW_DIAGNOSTIC_OUTPUT')
        with args.raw_diagnostic_output.open('x')as handle:
            json.dump({'fixtures':report,'linear_normal_derivative_checks':linear,'wrong_formula_max_residuals':rejected},handle,sort_keys=True,indent=2);handle.write('\n')
    with args.output.open('x')as handle:
        json.dump(out,handle,sort_keys=True,indent=2);handle.write('\n')
    print(json.dumps({'status':out['status'],'scalar_fixtures':len(report),'wrong_formula_controls':len(rejected)},sort_keys=True))


if __name__=='__main__':
    main()
