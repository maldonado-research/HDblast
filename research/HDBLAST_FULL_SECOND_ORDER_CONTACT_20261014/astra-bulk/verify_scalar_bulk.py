#!/usr/bin/env python3
"""Independent exact local-jet checks of the 2+3 scalar bulk equations.

No scientific imports from another capsule. The Ricci implementation below
uses the defining metric/Christoffel formula, not warped-product reductions.
All base jets are algebraically independent symbols. Angular dependence is
checked separately as an exact rational identity on stereographic S^3.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import sys
import sympy as s

CHECKS=[]

def zero(name, expression):
    value=s.factor(s.cancel(expression))
    if value != 0:
        raise RuntimeError(f"{name}: nonzero residual {value}")
    CHECKS.append(name)

def equal(name, actual, expected):
    zero(name, actual-expected)

def jet(name):
    return s.symbols(name+' '+name+'_y '+name+'_t '+name+'_yy '+name+'_yt '+name+'_tt')

def basejet(j):
    v,vy,vt,vyy,vyt,vtt=j
    return v,[vy,vt],[[vyy,vyt],[vyt,vtt]]

def multiply(a,b):
    av,ad,add=a; bv,bd,bdd=b
    dim=len(ad)
    return av*bv,[ad[i]*bv+av*bd[i] for i in range(dim)],[[add[i][j]*bv+ad[i]*bd[j]+ad[j]*bd[i]+av*bdd[i][j] for j in range(dim)] for i in range(dim)]

def scale(a,k):
    v,d,dd=a
    return k*v,[k*x for x in d],[[k*x for x in row] for row in dd]

def lift_base(a):
    v,d,dd=a
    return v,d+[s.S.Zero]*3,[[dd[i][j] if i<2 and j<2 else s.S.Zero for j in range(5)] for i in range(5)]

def lift_angular(v,d,dd):
    return v,[s.S.Zero]*2+d,[[dd[i-2][j-2] if i>=2 and j>=2 else s.S.Zero for j in range(5)] for i in range(5)]

def const(v,dim=5):
    return s.sympify(v),[s.S.Zero]*dim,[[s.S.Zero]*dim for _ in range(dim)]

def metric_from_entries(entries):
    g=s.zeros(5); dg=[s.zeros(5) for _ in range(5)]; ddg=[[s.zeros(5) for _ in range(5)] for _ in range(5)]
    for (i,j),(v,d,dd) in entries.items():
        g[i,j]=g[j,i]=v
        for a in range(5):
            dg[a][i,j]=dg[a][j,i]=d[a]
            for b in range(5):
                ddg[a][b][i,j]=ddg[a][b][j,i]=dd[a][b]
    return g,dg,ddg

def direct_geometry(background,perturbation,order):
    """Coefficient arrays of Ricci using g+epsilon*h, without ansatz formulas."""
    g,dg,ddg=background; h,dh,ddh=perturbation
    inv=[g.inv(),-g.inv()*h*g.inv()]
    if order==2: inv.append(g.inv()*h*g.inv()*h*g.inv())
    gd=[dg,dh]+([[s.zeros(5) for _ in range(5)]] if order==2 else [])
    gdd=[ddg,ddh]+([[[s.zeros(5) for _ in range(5)] for _ in range(5)]] if order==2 else [])
    dinv=[]
    for k in range(order+1):
        dinv.append([])
        for a in range(5):
            total=s.zeros(5)
            for i in range(k+1):
                for j in range(k-i+1):
                    ell=k-i-j
                    if ell<=1: total-=inv[i]*gd[ell][a]*inv[j]
            dinv[k].append(total)
    G=[[[[s.S.Zero for b in range(5)] for a in range(5)] for c in range(5)] for k in range(order+1)]
    dG=[[[[[s.S.Zero for b in range(5)] for a in range(5)] for c in range(5)] for d in range(5)] for k in range(order+1)]
    for k in range(order+1):
        for c in range(5):
            for a in range(5):
                for b in range(a,5):
                    val=0
                    for i in range(k+1):
                        j=k-i
                        if j>1:continue
                        for d in range(5):
                            val+=inv[i][c,d]*(gd[j][a][d,b]+gd[j][b][d,a]-gd[j][d][a,b])/2
                    G[k][c][a][b]=G[k][c][b][a]=s.expand(val)
                    for e in range(5):
                        val=0
                        for i in range(k+1):
                            j=k-i
                            if j>1:continue
                            for d in range(5):
                                val+=(dinv[i][e][c,d]*(gd[j][a][d,b]+gd[j][b][d,a]-gd[j][d][a,b])+inv[i][c,d]*(gdd[j][e][a][d,b]+gdd[j][e][b][d,a]-gdd[j][e][d][a,b]))/2
                        dG[k][e][c][a][b]=dG[k][e][c][b][a]=s.expand(val)
    ric=[s.zeros(5) for _ in range(order+1)]
    for k in range(order+1):
        for a in range(5):
            for b in range(a,5):
                val=sum(dG[k][c][c][a][b]-dG[k][b][c][a][c] for c in range(5))
                for i in range(k+1):
                    j=k-i
                    val+=sum(G[i][c][c][d]*G[j][d][a][b]-G[i][c][b][d]*G[j][d][a][c] for c in range(5) for d in range(5))
                ric[k][a,b]=ric[k][b,a]=s.factor(s.expand(val))
    return ric,inv,G


def sphere_geometry(point):
    x=s.symbols('x0:3'); den=1+sum(z*z for z in x)
    kappa=4/den**2
    ambient=[2*z/den for z in x]+[(1-sum(z*z for z in x))/den]
    Y=ambient[0]**2-ambient[1]**2
    F2=ambient[0]**2+ambient[1]**2-s.Rational(1,2)
    Q=Y**2-(ambient[0]**2+ambient[1]**2)/2+s.Rational(1,12)
    Gamma=[[[s.factor(((s.diff(kappa,x[j]) if i==k else 0)+(s.diff(kappa,x[k]) if i==j else 0)-(s.diff(kappa,x[i]) if j==k else 0))/(2*kappa)) for k in range(3)] for j in range(3)] for i in range(3)]
    def tensor(f):
        df=[s.diff(f,z) for z in x]
        hess=[[s.diff(df[i],x[j])-sum(Gamma[k][i][j]*df[k] for k in range(3)) for j in range(3)] for i in range(3)]
        return df,hess
    dy,hy=tensor(Y); dq,hq=tensor(Q); df2,hf2=tensor(F2)
    zero('sphere Y eigenvalue',sum(hy[i][i] for i in range(3))/kappa+8*Y)
    zero('sphere Q eigenvalue',sum(hq[i][i] for i in range(3))/kappa+24*Q)
    zero('sphere F2 eigenvalue',sum(hf2[i][i] for i in range(3))/kappa+8*F2)
    zero('sphere scalar product decomposition',Y**2-Q-F2/2-s.Rational(1,6))
    zero('sphere gradient product decomposition',sum(z*z for z in dy)/kappa+4*Q-2*F2-s.Rational(4,3))
    # Ambient Weyl tensor, restricted directly to the tangent bundle.
    sectionals={(0,1):-s.Rational(1,3),(2,3):-s.Rational(1,3),
                (0,2):s.Rational(1,6),(0,3):s.Rational(1,6),
                (1,2):s.Rational(1,6),(1,3):s.Rational(1,6)}
    TT=s.zeros(3)
    for (a,b),weight in sectionals.items():
        oneform=[s.diff(ambient[a],z)*ambient[b]-s.diff(ambient[b],z)*ambient[a] for z in x]
        for i in range(3):
            for j in range(3):TT[i,j]+=weight*oneform[i]*oneform[j]
    normDY=sum(z*z for z in dy)/kappa
    zero('sphere ambient Weyl trace',sum(TT[i,i] for i in range(3))/kappa)
    for i in range(3):
        for j in range(i,3):
            trace=kappa if i==j else 0
            actual=dy[i]*dy[j]-trace*normDY/3
            expected=(hf2[i][j]+trace*8*F2/3)/2+(hq[i][j]+trace*8*Q)/3-s.Rational(4,3)*TT[i,j]
            zero(f'sphere full scalar-TT product decomposition {i}{j}',actual-expected)
    sub=dict(zip(x,point))
    def ev(e):return s.factor(e.subs(sub))
    def aj(f):return lift_angular(ev(f),[ev(s.diff(f,z)) for z in x],[[ev(s.diff(f,z,w)) for w in x] for z in x])
    return {'k':ev(kappa),'kj':aj(kappa),'Y':ev(Y),'Q':ev(Q),'Yj':aj(Y),'Qj':aj(Q),'dY':[ev(z) for z in dy],'hY':[[ev(z) for z in row] for row in hy],'dQ':[ev(z) for z in dq],'hQ':[[ev(z) for z in row] for row in hq]}


def setup(point):
    sphere=sphere_geometry(point)
    r,ry,ryy=s.symbols('rho rho_y rho_yy',nonzero=True)
    c=s.Rational(5,4); sh=s.Rational(3,4)
    rho2=(r*r,[2*r*ry,0],[[2*(ry*ry+r*ryy),0],[0,0]])
    cosh2=(c*c,[0,2*c*sh],[[0,0],[0,2*(c*c+sh*sh)]])
    R2=multiply(rho2,cosh2)
    bg={(0,0):const(1),(1,1):scale(lift_base(rho2),-1)}
    for i in range(2,5):bg[i,i]=multiply(lift_base(R2),sphere['kj'])
    return sphere,(r,ry,ryy,c,sh),rho2,R2,metric_from_entries(bg)


def quadratic_check(data, corrupt):
    sp,(r,ry,ryy,c,sh),rho2,R2,bg=data
    u,uy,ut,uyy,uyt,utt=ju=jet('u')
    v,vy,vt,vyy,vyt,vtt=jv=jet('v')
    p,py,U0,U1,U2,U3=s.symbols('p p_y U U1 U2 U3')
    psi=multiply(lift_base(basejet(ju)),sp['Yj'])
    hp={(0,0):scale(psi,-4),(1,1):scale(multiply(lift_base(rho2),psi),-2)}
    for i in range(2,5):hp[i,i]=scale(multiply(multiply(lift_base(R2),sp['kj']),psi),2)
    pert=metric_from_entries(hp)
    ric,inv,G=direct_geometry(bg,pert,2)
    y=sp['Y'];dy=sp['dY'];hy=sp['hY'];k=sp['k'];a=ry/r;ap=ryy/r-a*a;hh=sh/c
    boxu=-utt-3*hh*ut-8*u/c**2
    boxv=-vtt-3*hh*vt-8*v/c**2
    grad2=-ut**2*y**2+u**2*sum(z*z for z in dy)/(c*c*k)
    trace2=2*u*y*boxu*y+2*grad2-r*r*(4*uy**2+4*u*uyy+52*a*u*uy+24*(ap+4*a*a)*u*u)*y*y
    expected=s.zeros(5)
    expected[0,0]=(-4*uy**2+8*u*uyy-16*a*u*uy)*y*y+(-4*u*boxu*y*y+8*grad2)/r**2
    expected[0,1]=expected[1,0]=(6*u*uyt-24*a*u*ut)*y*y
    for i in range(3):expected[0,2+i]=expected[2+i,0]=(6*u*uy-24*a*u*u)*y*dy[i]
    expected[1,1]=(12*u*utt+6*ut*ut)*y*y-trace2
    for i in range(3):expected[1,2+i]=expected[2+i,1]=(18*u*ut-12*hh*u*u)*y*dy[i]
    for i in range(3):
        for j in range(3):
            expected[2+i,2+j]=12*u*u*y*hy[i][j]+6*u*u*dy[i]*dy[j]
            if i==j:expected[2+i,2+j]+=(-12*c*c*hh*u*ut*y*y+c*c*trace2)*k
    if corrupt:expected[0,1]+=u*uyt*y*y
    for i in range(5):
        for j in range(i,5):equal(f'direct quadratic Ricci {i}{j}',ric[2][i,j],expected[i,j])
    # Scalar operator from g^{-1}(partial partial phi - Gamma partial phi).
    phi0d=[p,0,0,0,0];phi0dd=s.zeros(5);phi0dd[0,0]=py
    _,chid,chidd=multiply(lift_base(basejet(jv)),sp['Yj'])
    box2=0
    for i in range(5):
        for j in range(5):
            box2+=inv[2][i,j]*(phi0dd[i,j]-sum(G[0][k0][i][j]*phi0d[k0] for k0 in range(5)))
            box2+=inv[1][i,j]*(chidd[i][j]-sum(G[0][k0][i][j]*chid[k0]+G[1][k0][i][j]*phi0d[k0] for k0 in range(5)))
            box2-=inv[0][i,j]*sum(G[1][k0][i][j]*chid[k0]+G[2][k0][i][j]*phi0d[k0] for k0 in range(5))
    expectbox=(6*uy*vy+24*p*u*uy+4*u*(vyy+4*a*vy)+16*u*u*(py+4*a*p)-2*u*boxv/r**2)*y*y
    equal('direct quadratic scalar dAlembertian',box2,expectbox)
    # Coefficient algebra for complete l=4 matter terms and Taylor normalization.
    q=s.symbols('Q'); dQ=s.symbols('dQ'); hQ=s.symbols('hQ');kQ=s.symbols('kQ')
    # Independent scalar product rules: [Y²]_4=Q, [|DY|²]_4=-4Q,
    # [Y DY]_4=dQ/2, [(DY DY)^TF]_4=hQ/3.
    L=2*u*boxu+2*(-ut**2-4*u*u/c**2)-r*r*(4*uy**2+4*u*uyy+52*a*u*uy+24*(ap+4*a*a)*u*u+U2*v*v/3+4*U1*u*v/3)
    Syy=-8*uy**2+16*u*uyy-32*a*u*uy+(8*u*utt+24*hh*u*ut-16*ut**2)/r**2-2*vy**2-2*U2*v*v/3+16*U1*u*v/3
    equal('l4 yy complete source',Syy,2*((-4*uy**2+8*u*uyy-16*a*u*uy)+(-4*u*boxu+8*(-ut**2-4*u*u/c**2))/r**2-vy**2-U2*v*v/3+8*U1*u*v/3))
    equal('l4 yt complete source',12*u*uyt-48*a*u*ut-2*vy*vt,2*(6*u*uyt-24*a*u*ut-vy*vt))
    equal('l4 yi complete source',6*u*uy-24*a*u*u-v*vy,2*(6*u*uy-24*a*u*u-v*vy)/2)
    equal('l4 ti complete source',18*u*ut-12*hh*u*u-v*vt,2*(18*u*ut-12*hh*u*u-v*vt)/2)
    equal('l4 angular trace complete source',c*c*(2*L-24*hh*u*ut)-80*u*u+s.Rational(8,3)*v*v,2*c*c*L-24*c*c*hh*u*ut+(24*u*u*(-8)+(12*u*u-2*v*v)*(-4))/3)
    equal('l4 angular TF complete source',8*u*u-s.Rational(2,3)*v*v,12*u*u+(-12*u*u-2*v*v)/3)
    equal('l4 scalar complete source',12*uy*vy+48*p*u*uy+8*u*(vyy+4*a*vy)+32*u*u*(py+4*a*p)-4*u*boxv/r**2-U3*v*v,2*(expectbox/y**2-U3*v*v/2))
    return {'source_jet_variables':list(map(str,ju+jv+(p,py,U0,U1,U2,U3,r,ry,ryy)))}


def linear_check(data):
    sp,(r,ry,ryy,c,sh),rho2,R2,bg=data
    jf=[jet('N'),jet('B'),jet('C')];jk=jet('K');jw=jet('w')
    hp={(0,0):multiply(lift_base(basejet(jf[0])),sp['Qj']),(0,1):multiply(lift_base(basejet(jf[1])),sp['Qj']),(1,1):multiply(lift_base(basejet(jf[2])),sp['Qj'])}
    for i in range(2,5):hp[i,i]=multiply(multiply(multiply(lift_base(R2),sp['kj']),lift_base(basejet(jk))),sp['Qj'])
    pert=metric_from_entries(hp)
    ric,inv,G=direct_geometry(bg,pert,1)
    # Build the proposed covariant 2D formula independently by differentiating
    # its printed base expressions before substituting arbitrary base jets.
    y,t=s.symbols('y t');coords=[y,t]
    rr=s.Function('r')(y);cc=s.cosh(t)*c+s.sinh(t)*sh
    qq=s.diag(1,-rr**2);qi=qq.inv();rad=rr*cc
    fc=s.Matrix([[s.Function('N')(y,t),s.Function('B')(y,t)],[s.Function('B')(y,t),s.Function('C')(y,t)]])
    KK=s.Function('K')(y,t);ww=s.Function('w')(y,t);phi=s.Function('phi')(y)
    gam=[[[sum(qi[a,d]*(s.diff(qq[d,b],coords[c0])+s.diff(qq[d,c0],coords[b])-s.diff(qq[b,c0],coords[d]))/2 for d in range(2)) for c0 in range(2)] for b in range(2)] for a in range(2)]
    rab=[s.diff(rad,z)/rad for z in coords];rup=list(qi*s.Matrix(rab));fup=qi*fc;ftrace=s.trace(fup)
    df=[[[s.diff(fc[a,b],coords[d])-sum(gam[e][d][a]*fc[e,b]+gam[e][d][b]*fc[a,e] for e in range(2)) for b in range(2)] for a in range(2)] for d in range(2)]
    AA=[[[sum(qi[a,d]*(df[b][d][c0]+df[c0][d][b]-df[d][b][c0])/2 for d in range(2)) for c0 in range(2)] for b in range(2)] for a in range(2)]
    def hess(z,a,b):return s.diff(z,coords[a],coords[b])-sum(gam[d][a][b]*s.diff(z,coords[d]) for d in range(2))
    def divA(a,b):
        return sum(s.diff(AA[d][a][b],coords[d])+sum(gam[d][d][e]*AA[e][a][b]-gam[e][d][a]*AA[d][e][b]-gam[e][d][b]*AA[d][a][e] for e in range(2)) for d in range(2))
    ricab=s.Matrix(2,2,lambda a,b:divA(a,b)-hess(ftrace,a,b)/2+3*sum(rab[d]*AA[d][a][b] for d in range(2))-s.Rational(3,2)*(hess(KK,a,b)+rab[a]*s.diff(KK,coords[b])+rab[b]*s.diff(KK,coords[a]))+12*fc[a,b]/rad**2)
    divf=[]
    for a in range(2):
        divf.append(sum(s.diff(fup[b,a],coords[b])+sum(gam[b][b][d]*fup[d,a]-gam[d][b][a]*fup[b,d] for d in range(2)) for b in range(2)))
    ricai=[(divf[a]-s.diff(ftrace,coords[a])-2*s.diff(KK,coords[a])+sum(rab[b]*fup[b,a] for b in range(2))+rab[a]*ftrace)/2 for a in range(2)]
    BB=[sum(rup[b]*fup[a,b] for b in range(2))-rup[a]*KK-sum(qi[a,b]*s.diff(KK,coords[b])/2 for b in range(2)) for a in range(2)]
    divB=sum(s.diff(BB[a],coords[a])+sum(gam[a][a][b]*BB[b] for b in range(2)) for a in range(2))
    J=ftrace+KK
    iso=rad**2*(divB+3*sum(rab[a]*BB[a] for a in range(2))-sum(rup[a]*s.diff(J,coords[a])/2 for a in range(2)))+12*KK+4*J
    sub={rr:r,s.diff(rr,y):ry,s.diff(rr,y,2):ryy}
    for fn,j in zip([fc[0,0],fc[0,1],fc[1,1],KK,ww],jf+[jk,jw]):
        sub.update({fn:j[0],s.diff(fn,y):j[1],s.diff(fn,t):j[2],s.diff(fn,y,2):j[3],s.diff(fn,y,t):j[4],s.diff(fn,t,2):j[5]})
    def ev(z):return s.factor(s.expand(z.doit().subs(sub).subs({y:0,t:0})))
    for a in range(2):
        for b in range(a,2):equal(f'direct linear Ricci base {a}{b}',ric[1][a,b],ev(ricab[a,b])*sp['Q'])
    for a in range(2):
        for i in range(3):equal(f'direct linear Ricci mixed {a}{i}',ric[1][a,2+i],ev(ricai[a])*sp['dQ'][i])
    for i in range(3):
        for j in range(i,3):
            htf=sp['hQ'][i][j]+(8*sp['k']*sp['Q'] if i==j else 0)
            expected=-ev(J)*htf/2+(ev(iso)*sp['k']*sp['Q'] if i==j else 0)
            equal(f'direct linear Ricci angular {i}{j}',ric[1][2+i,2+j],expected)
    # Direct scalar linearization versus printed covariant formula.
    p,py=s.symbols('p p_y');phi0d=[p,0,0,0,0];phi0dd=s.zeros(5);phi0dd[0,0]=py
    _,wd,wdd=multiply(lift_base(basejet(jw)),sp['Qj'])
    box1=sum(inv[1][a,b]*(phi0dd[a,b]-sum(G[0][d][a][b]*phi0d[d] for d in range(5)))+inv[0][a,b]*(wdd[a][b]-sum(G[0][d][a][b]*wd[d]+G[1][d][a][b]*phi0d[d] for d in range(5))) for a in range(5) for b in range(5))
    fuu=qi*fc*qi
    divfuu=[sum(s.diff(fuu[a,b],coords[a])+sum(gam[a][a][d]*fuu[d,b]+gam[b][a][d]*fuu[a,d] for d in range(2)) for a in range(2)) for b in range(2)]
    scal=sum(qi[a,b]*hess(ww,a,b) for a in range(2) for b in range(2))+3*sum(rup[a]*s.diff(ww,coords[a]) for a in range(2))-24*ww/rad**2
    scal-=sum(fuu[a,b]*hess(phi,a,b) for a in range(2) for b in range(2))
    scal-=sum((divfuu[b]+3*sum(rab[a]*fuu[a,b] for a in range(2))-sum(qi[b,a]*s.diff(ftrace+3*KK,coords[a])/2 for a in range(2)))*s.diff(phi,coords[b]) for b in range(2))
    sub.update({s.diff(phi,y):p,s.diff(phi,y,2):py})
    equal('direct scalar linear dAlembertian',box1,ev(scal)*sp['Q'])
    # General connection decomposition gives gauge fixing and TF relation.
    N,B,C=jf[0][0],jf[1][0],jf[2][0];K=jk[0]
    equal('RW scalar trace combination',ev(J),N-C/r**2+K)
    return {'linear_jet_variables':list(map(str,sum(jf,())+jk+jw+(p,py,r,ry,ryy)))}


def run(args):
    point=(s.Rational(1,3),s.Rational(1,2),s.Rational(2,5))
    data=setup(point)
    info=quadratic_check(data,args.inject_corruption)
    info.update(linear_check(data))
    if args.model_root:
        model_root=Path(args.model_root)
        input_names=[
            'HDBLAST_REGULAR_CONE_EVOLUTION_20261013/manuscript/REGULAR_CONE_EVOLUTION.md',
            'HDBLAST_MATERIAL_JUNCTIONS_20261014/MATERIAL_SECOND_ORDER_JUNCTIONS.md',
            'HDBLAST_QUADRATIC_TENSOR_CONTACT_20261014/bulk/FULL_QUADRATIC_TT_BULK_SOURCE.md',
            'HDBLAST_QUADRATIC_TENSOR_CONTACT_20261014/boundary/FORCED_TENSOR_CONTACT_IBVP.md']
        info['model_source_sha256']={name:hashlib.sha256((model_root/name).read_bytes()).hexdigest() for name in input_names}
    info.update({'status':'passed','check_count':len(CHECKS),'checks':CHECKS,'python':platform.python_version(),'sympy':s.__version__,'optimized':not __debug__,'angular_point':list(map(str,point)),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Exact arbitrary base-jet identities at one nonspecial rational sphere point, plus globally exact stereographic scalar projection identities; analytic connection derivation supplies the all-point operator proof.'})
    if args.output:
        target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps({k:info[k] for k in ['status','check_count','python','sympy','optimized','script_sha256']},sort_keys=True))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    parser.add_argument('--model-root',help='Research root used for recording the four controlling model/source SHA-256 hashes.')
    parser.add_argument('--inject-corruption',action='store_true',help='Change the actual proposed quadratic yt formula; the exact geometry check must reject it.')
    try:run(parser.parse_args())
    except Exception as exc:
        print(f'VERIFICATION FAILED: {exc}',file=sys.stderr)
        sys.exit(1)
