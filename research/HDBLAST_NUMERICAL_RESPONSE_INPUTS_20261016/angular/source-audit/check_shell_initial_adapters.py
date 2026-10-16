"""Symbolic coefficient identities for shell, Einstein matter and initial data."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from fractions import Fraction as Q
import angular_exact as e
import source_adapters as a
import source_laurent as ring
from check_scalar_gauge_adapters import rank
POINTS=((Q(1),Q(0),Q(0),Q(0)),(Q(0),Q(0),Q(1),Q(0)),(Q(3,5),Q(0),Q(4,5),Q(0)),(Q(1,3),Q(2,3),Q(2,3),Q(0)))
def peval(poly,point):return sum((v*product(point[i]**power for i,power in enumerate(exp)) for exp,v in poly.items()),Q(0))
def product(values):
    value=Q(1)
    for item in values:value*=item
    return value
def teval(t,point):return {key:peval(poly,point) for key,poly in t.items()}
def combine(coeff,basis,key):return sum((x*y[key] for x,y in zip(coeff,basis)),Q(0))
def require_equal(x,y,label):
    if x-y!=0:raise ValueError(label)

def validate():
    sym=dict(zip(ring.NAMES,ring.variables()))
    rho,aa,p,sp,spp,sppp=(sym[k] for k in ('rho','a','p','sigma1','sigma2','sigma3'))
    hs,X,V,L,d,v=(sym[k] for k in ('A','B','C','D','E','v'))
    shell=a.shell_sources(a=aa,p=p,sigmap=sp,sigmapp=spp,sigmappp=sppp,h=hs,X=X,V=V,L=L,d=d,v=v)
    m={k:sym[k] for k in ('A','B','C','D','E','L','M','Ay','At','By','Bt','Cy','Ct','Dy','Dt','Ey','Et','Ly','Lt','My','Mt')}
    b={k:sym[k] for k in ('rho','c','a','theta','p','py','Vppp','Vp','Vpp')}
    sc={k:sym[k] for k in ('s','sy','st','syy','syt','stt')}
    matter=a.einstein_matter_source(m,sc,b)
    ikeys=('rho','a','p','u','uy','ut','v','vt','Fmom','Gmom','Fmom_y')
    ini={k:sym[k] for k in ikeys};ini['f']=tuple(sym[k] for k in ('f0','f2','f4'));ini['fy']=tuple(sym[k] for k in ('fy0','fy2','fy4'))
    initial=a.initial_frame_data(**ini)
    tensor_basis=[e.P,e.tmul(e.Z,e.P),e.tmul(e.F,e.P),*e.HESS,*e.make_tt()]
    tensor_rows=[]
    for point in POINTS:
        Y=peval(e.Y,point);P=teval(e.P,point);DY=teval(e.grad(e.Y),point);HY=teval(e.hess(e.Y),point)
        SB=[teval(t,point) for t in e.SCALARS];VB=[teval(t,point) for t in e.GRADS];TB=[teval(t,point) for t in tensor_basis]
        tensor_rows.extend([t[i,j] for t in TB] for i in range(4) for j in range(4))
        H=hs*Y;eta=v*Y
        # Direct affine physical row from N=sqrt(1+eps H), before eliminating m_y.
        gamma={(0,0):Q(-1)};mt={(0,0):-2*X*Y}
        for i in range(4):
            gamma[0,i+1]=gamma[i+1,0]=Q(0);mt[0,i+1]=mt[i+1,0]=V*DY[i,]
            for j in range(4):gamma[i+1,j+1]=P[i,j];mt[i+1,j+1]=2*(L*Y*P[i,j]+d*HY[i,j])
        my={key:2*aa*mt[key]+(aa*H+sp*eta/3)*g for key,g in gamma.items()}
        QR={key:(Q(3,4)*aa*H*H-spp*eta*eta/6)*g-H*my[key]/2-sp*eta*mt[key]/3 for key,g in gamma.items()}
        require_equal(QR[0,0],combine(shell['QK00'],SB,()),'SOURCE_SHELL_00')
        for i in range(4):
            require_equal(QR[0,i+1],combine(shell['QK0a'],VB,(i,)),'SOURCE_SHELL_0A')
            for j in range(4):require_equal(QR[i+1,j+1],combine(shell['QKab'],TB,(i,j)),'SOURCE_SHELL_TENSOR')
        eta_y=p*H/2-spp*eta/2
        QP=Q(3,4)*p*H*H-H*eta_y+sppp*eta*eta/2
        require_equal(QP,combine(shell['Qp'],SB,()),'SOURCE_SHELL_SCALAR')
        # Join to full qP_y formula, using radially parallel shell jet relations.
        tau0=2*X+6*L-16*d;mu=aa*hs+sp*v/3
        shelljet={k:Q(0) for k in m}
        shelljet.update(A=hs,D=-2*X,E=V,L=2*L,M=2*d,Ay=-4*aa*hs+2*aa*tau0+4*sp*v/3,Dy=-mu,Ly=mu)
        qcy=a.gauge_source(shelljet,b)['y']
        require_equal(combine(qcy,SB,()),combine(shell['QC_y'],SB,()),'SOURCE_SHELL_GAUGE_JOIN')
        trace=-QR[0,0]+sum(QR[i+1,i+1] for i in range(4))
        require_equal(-2*trace-2*combine(qcy,SB,()),combine(shell['GN'],SB,()),'SOURCE_NORMAL_GAUGE_JOIN')
        # Direct Einstein stress/potential Hessian in all frame slots.
        qq={(0,0):m['A']*Y,(0,1):m['B']*Y,(1,0):m['B']*Y,(1,1):m['D']*Y}
        gg={(0,0):Q(1),(0,1):Q(0),(1,0):Q(0),(1,1):Q(-1)}
        grad={0:sc['sy']*Y,1:sc['st']*Y/rho}
        for i in range(4):
            grad[i+2]=sc['s']*DY[i,]/(rho*b['c'])
            qq[0,i+2]=qq[i+2,0]=m['C']*DY[i,];qq[1,i+2]=qq[i+2,1]=m['E']*DY[i,]
            gg[0,i+2]=gg[i+2,0]=gg[1,i+2]=gg[i+2,1]=Q(0)
            for j in range(4):qq[i+2,j+2]=m['L']*Y*P[i,j]+m['M']*HY[i,j];gg[i+2,j+2]=P[i,j]
        mat={key:-2*grad[key[0]]*grad[key[1]]-Q(2,3)*b['Vpp']*sc['s']**2*Y**2*g-Q(4,3)*b['Vp']*sc['s']*Y*qq[key] for key,g in gg.items()}
        for name,key in [('yy',(0,0)),('y0',(0,1)),('00',(1,1))]:require_equal(mat[key],combine(matter[name],SB,()),'SOURCE_MATTER_'+name.upper())
        for i in range(4):
            require_equal(mat[0,i+2],combine(matter['ya'],VB,(i,)),'SOURCE_MATTER_YA')
            require_equal(mat[1,i+2],combine(matter['0a'],VB,(i,)),'SOURCE_MATTER_0A')
            for j in range(4):require_equal(mat[i+2,j+2],combine(matter['ab'],TB,(i,j)),'SOURCE_MATTER_AB')
        # Independent ADM parameter expansion and frame conversion at t=0.
        U=sym['u']*Y;Uy=sym['uy']*Y;T=sym['ut']*Y;Vt=sym['vt']*Y
        f=combine(ini['f'],SB,());fy=combine(ini['fy'],SB,());q=sym['Fmom']*Y**2/3
        N0=rho;N1=rho*U;Ndot1=rho*T;Ndot2=4*rho*rho*q-24*rho*U*T
        betadot1=-6*rho*rho*(Uy+4*aa*U)
        betadot2=rho*rho*(-74*U*Uy-Q(3,2)*fy-194*aa*U*U-3*aa*f)
        val=initial['values'];vel=initial['t_velocities']
        require_equal(-2*N1*N1/(rho*rho),combine(val['00'],SB,()),'SOURCE_INITIAL_LAPSE_VALUE')
        require_equal(-3*fy/(2*p),combine(val['scalar'],SB,()),'SOURCE_INITIAL_SCALAR_VALUE')
        require_equal((-2*N0*Ndot2-4*N1*Ndot1)/(rho*rho),combine(vel['00'],SB,()),'SOURCE_INITIAL_LAPSE_VELOCITY')
        require_equal((betadot2+2*(-4*U)*betadot1)/rho,combine(vel['y0'],SB,()),'SOURCE_INITIAL_MIXED_VELOCITY')
        require_equal(2*N0*q+4*N1*(-2*T/rho),combine(vel['yy'],SB,()),'SOURCE_INITIAL_YY_VELOCITY')
        require_equal(2*N0*q+4*N1*T/rho,combine(vel['ab_trace_amplitude'],SB,()),'SOURCE_INITIAL_ANGULAR_VELOCITY')
        Pi1=Vt/rho;Pi2=(sym['Gmom']-sym['Fmom_y'])*Y*Y/p
        require_equal(N0*Pi2+2*N1*Pi1,combine(vel['scalar'],SB,()),'SOURCE_INITIAL_SCALAR_VELOCITY')
        for i in range(4):require_equal(22*U*sym['u']*DY[i,]-combine(ini['f'][1:],VB,(i,))/2,combine(vel['0a'],VB,(i,)),'SOURCE_INITIAL_GRADIENT_VELOCITY')
    if rank(tensor_rows)!=7:raise ValueError('SOURCE_TENSOR_WITNESSES_NOT_COMPLETE')
    print('PASS: universal shell, matter and initial symbolic identities on rank-complete angular witnesses (tensor rank 7).')
    return {'angular_tensor_witness_rank':7,'symbolic_background_and_jets':len(ring.NAMES)}
if __name__=='__main__':validate()
