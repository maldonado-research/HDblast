"""Exact angular adapters for the registered common-gauge first pair.

All inputs may be Fraction or symbolic numbers. Sources are second parameter
variations; no eps**2 coefficient or forcing-sign conversion is implicit.
Qphi/QC formulas simplify using the actual first harmonic condition DC[q]=0.
"""
from fractions import Fraction as Q
Y2=(Q(1,6),Q(1,2),Q(1))
GRAD2=(Q(4,3),Q(2),Q(-4))
HESS2=(Q(8),Q(4),Q(24))
# Tensor order: kappa, Z*kappa, F*kappa, HessTF(Z), HessTF(F), T+, T-.
DYDY=(Q(4,9),Q(2,3),Q(-4,3),Q(1,2),Q(1,3),Q(-2,3),Q(-2,3))
YHESS=(Q(-4,9),Q(-4,3),Q(-8,3),Q(-1,4),Q(1,6),Q(2,3),Q(2,3))
HESS_SQUARE=(Q(8,3),Q(4,3),Q(8),Q(5,2),Q(-1),Q(-2),Q(-2))

def add(*vectors):return tuple(sum(x) for x in zip(*vectors))
def scale(a,v):return tuple(a*x for x in v)
def scalar_project(y2,g2=0,h2=0):return add(scale(y2,Y2),scale(g2,GRAD2),scale(h2,HESS2))
def vector_project(ydy,hdy=0):return (ydy/4+hdy,ydy/2-2*hdy)

def connection_contractions(m,b):
    """q^{AB} P_{iAB}; lower i=y,0,a. y/0 output in (Y²,|DY|²,|HessY|²), a in (YDY,HessY·DY).
    m has metric amplitudes A,B,C,D,E,L,M and their ordinary y,t derivatives.
    b has rho, c=cosh(t), theta=tanh(t), a=rho_y/rho.
    """
    A,B,C,D,E,L,M=(m[k] for k in ('A','B','C','D','E','L','M'))
    Ay,At,By,Bt,Cy,Ct,Dy,Dt,Ey,Et,Ly,Lt,My,Mt=(m[k] for k in ('Ay','At','By','Bt','Cy','Ct','Dy','Dt','Ey','Et','Ly','Lt','My','Mt'))
    rho,c,a,h=b['rho'],b['c'],b['a'],b['theta']/b['rho']
    j=1/(rho*c)
    Ry=a*(A-L)-h*B-Ly/2
    Ty=j*C-a*M-My/2
    vy=(A*Ay/2-B*(At/rho-2*a*B)+D*(Bt/rho-a*(A+D)-Dy/2)+(3*L-8*M)*Ry-8*L*Ty,
        C*(j*A-2*a*C)-E*(Ct/rho+j*B-h*C-2*a*E-Ey),M*Ty)
    R0=a*B-h*(D+L)-Lt/(2*rho)
    T0=j*E-h*M-Mt/(2*rho)
    v0=(A*(By-At/(2*rho)+a*B)-B*Dy+D*(Dt/(2*rho)-a*B)+(3*L-8*M)*R0-8*L*T0,
        C*(Ey+j*B-h*C-Ct/rho)-E*(j*D-2*h*E),M*T0)
    va=(A*(Cy-j*A/2+a*C)-B*(Ey+Ct/rho-j*B+h*C)+D*(Et/rho-a*C-j*D/2+h*E)
        +C*Ly-E*Lt/rho+(a*C-h*E-j*L/2)*(3*L-8*M)+j*L*(L-2*M),
        C*My-E*Mt/rho+j*M*(L-2*M))
    return {'y':vy,'0':v0,'a':va}

def gauge_source(m,b):
    """Lower QC components in scalar (1,Z,F), vector (DZ,DF) bases."""
    p=connection_contractions(m,b)
    return {'y':scale(-2,scalar_project(*p['y'])),
            '0':scale(-2,scalar_project(*p['0'])),
            'a':scale(-2,vector_project(*p['a']))}

def scalar_bulk_source(m,s,b):
    """Qphi projections (1,Z,F), with actual DC[q]=0.
    s gives s,sy,st,syy,syt,stt for eta=sY; b additionally gives p,py,Vppp.
    """
    A,B,C,D,E,L,M=(m[k] for k in ('A','B','C','D','E','L','M'))
    rho,c,a,theta=(b[k] for k in ('rho','c','a','theta'))
    p,py,Vppp=(b[k] for k in ('p','py','Vppp'))
    z,zy,zt,zyy,zyt,ztt=(s[k] for k in ('s','sy','st','syy','syt','stt'))
    qsq=(py*(A*A-B*B)+a*p*(D*D-B*B+3*L*L-16*L*M),
         py*C*C+a*p*(C*C-2*E*E),a*p*M*M)
    radial_angular=a*zy-theta*zt/(rho*rho)
    qhess=(A*zyy-2*B*(zyt-a*zt)/rho+D*(ztt/(rho*rho)-a*zy)
           +(3*L-8*M)*radial_angular-8*L*z/(rho*rho*c*c),
           2*C*(zy-a*z)/(rho*c)-2*E*(zt-theta*z)/(rho*rho*c),
           M*z/(rho*rho*c*c))
    qp=connection_contractions(m,b)['y']
    primitive=tuple(2*(x-y+p*zv) for x,y,zv in zip(qsq,qhess,qp))
    primitive=(primitive[0]-Vppp*z*z,primitive[1],primitive[2])
    return scalar_project(*primitive)

def shell_sources(*,a,p,sigmap,sigmapp,sigmappp,h,X,V,L,d,v):
    """One-copy outward physical Hessians and normal gauge source.
    H=hY, eta=vY, m00=-2XY,m0a=V DY,mab=2(LY kappa+d HessY).
    """
    a0=a*h*h/4-sigmap*h*v/6-sigmapp*v*v/6
    l0=a*h+sigmap*v/3
    p0=p*h*h/4+sigmapp*h*v/2+sigmappp*v*v/2
    tau=2*X+6*L-16*d
    asquare=(4*X*X+12*L*L-64*L*d,-2*V*V,4*d*d)
    qc=(4*a*h*h-3*a*h*tau-4*sigmap*h*v/3+sigmap*v*tau/3,0,0)
    qc=add(qc,scale(2*a,asquare))
    gn=(4*sigmapp*v*v/3+4*sigmap*h*v-10*a*h*h+8*a*h*tau,0,0)
    gn=add(gn,scale(-4*a,asquare))
    tensor=add(scale(a0-2*l0*L,(*Y2,0,0,0,0)),scale(-2*l0*d,YHESS))
    return {'QK00':scalar_project(-a0+2*l0*X),'QK0a':vector_project(-l0*V),
            'QKab':tensor,'Qp':scalar_project(p0),'QC_y':scalar_project(*qc),
            'GN':scalar_project(*gn)}

def initial_frame_data(*,rho,a,p,u,uy,ut,v,vt,Fmom,Gmom,Fmom_y,f,fy):
    """Frame metric/scalar initial values and ordinary t velocities at t=0.
    f,fy are coefficient triples in (1,Z,F). Fmom,Gmom are source profiles
    from EXPLICIT_COMPACT_SECOND_INITIAL_JET (its F,G, not angular F).
    A matching PDE in x=t/ell must multiply every velocity by ell.
    """
    zero=(0,0,0)
    zeta=scale(-3/(2*p),fy)
    q=Fmom/3
    values={'yy':zero,'y0':zero,'ya':(0,0),'00':scale(-2*u*u,Y2),
            '0a':(0,0),'ab_trace_amplitude':tuple(f),'ab_HessTF':(0,0),
            'ab_TT':(0,0),'scalar':zeta}
    mixed=add(scale(rho*(-26*u*uy-2*a*u*u),Y2),scale(-3*rho/2,fy),scale(-3*rho*a,f))
    angular=vector_project(22*u*u)
    angular=(angular[0]-f[1]/2,angular[1]-f[2]/2)
    velocities={'yy':scale(2*rho*q-8*u*ut,Y2),'y0':mixed,'ya':(0,0),
                '00':scale(-8*rho*q+44*u*ut,Y2),'0a':angular,
                'ab_trace_amplitude':scale(2*rho*q+4*u*ut,Y2),
                'ab_HessTF':(0,0),'ab_TT':(0,0),
                'scalar':scale(rho*(Gmom-Fmom_y)/p+2*u*vt,Y2)}
    return {'values':values,'t_velocities':velocities}

def einstein_matter_source(m,s,b):
    """Non-Ricci terms in QE, in scalar/vector/full angular tensor bases.
    b additionally supplies Vp,Vpp. All first scalar derivative products and
    the potential product are included, with the second-variation factor 2.
    """
    A,B,C,D,E,L,M=(m[k] for k in ('A','B','C','D','E','L','M'))
    z,zy,zt=(s[k] for k in ('s','sy','st'))
    rho,c,Vp,Vpp=(b[k] for k in ('rho','c','Vp','Vpp'))
    pot=Q(2,3)*Vpp*z*z
    cross=Q(4,3)*Vp*z
    angular=add(scale(-2*z*z/(rho*rho*c*c),DYDY),
                scale(-pot-cross*L,(*Y2,0,0,0,0)),scale(-cross*M,YHESS))
    return {'yy':scalar_project(-2*zy*zy-pot-cross*A),
            'y0':scalar_project(-2*zy*zt/rho-cross*B),
            'ya':vector_project(-2*zy*z/(rho*c)-cross*C),
            '00':scalar_project(-2*zt*zt/(rho*rho)+pot-cross*D),
            '0a':vector_project(-2*zt*z/(rho*rho*c)-cross*E),'ab':angular}
