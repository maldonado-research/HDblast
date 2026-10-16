#!/usr/bin/env python3
"""23-amplitude exact rational angular matrices for the full linear system.

Background parameters multiply fixed rational matrices.  No symbolic or numeric
angular quadrature is used.  Actual source amplitudes are supplied separately.
"""
from fractions import Fraction as Q
from itertools import product

N=23
LABELS=[f'{field}:{l}' for field in ('K','L','N','J','w') for l in (0,2,4)] + [f'{field}:{l}' for field in ('v','u') for l in (2,4)] + ['H:2','H:4','T:+','T:-']
IDX={name:i for i,name in enumerate(LABELS)}
def zeros(n=N,m=None):return [[Q(0) for _ in range(n if m is None else m)] for _ in range(n)]
def ident(n=N):return [[Q(i==j) for j in range(n)] for i in range(n)]
def add(*ms):return [[sum((m[i][j] for m in ms),Q(0)) for j in range(len(ms[0][0]))] for i in range(len(ms[0]))]
def scale(c,m):return [[c*v for v in row] for row in m]
def mm(a,b):return [[sum((a[i][k]*b[k][j] for k in range(len(b))),Q(0)) for j in range(len(b[0]))] for i in range(len(a))]
def put(m,out,inc,c,l):
    if f'{out}:{l}' in IDX and f'{inc}:{l}' in IDX:m[IDX[f'{out}:{l}']][IDX[f'{inc}:{l}']]+=Q(c)

FIBRES=[(i,j) for i in range(5) for j in range(i,5)]
SIG=[1,-1,1,1,1]
def omega(which,a,h):
    """Direct covariant connection on the 15 symmetric fibre components."""
    gamma=[[Q(0) for _ in range(5)] for _ in range(5)]
    if which==1:gamma[0][1]=a;gamma[1][0]=a
    else:
        gamma[0][which]=a;gamma[1][which]=h
        gamma[which][0]=-a;gamma[which][1]=h
    out=zeros(15)
    for col,(c,d) in enumerate(FIBRES):
        k=[[Q((i,j)==(c,d) or (j,i)==(c,d)) for j in range(5)] for i in range(5)]
        for row,(i,j) in enumerate(FIBRES):
            out[row][col]=-sum((gamma[i][b]*k[b][j]+gamma[j][b]*k[i][b] for b in range(5)),Q(0))
    return out
def curvature(b,d):
    """L_R from the five-dimensional Riemann/Ricci contractions."""
    def R(i,j,k,l):
        coeff=b if 0 in (i,j,k,l) else d
        return coeff*SIG[i]*SIG[j]*Q((i==k and j==l)-(i==l and j==k))
    ric=[[sum((SIG[c]*R(c,i,c,j) for c in range(5)),Q(0)) for j in range(5)] for i in range(5)]
    out=zeros(15)
    for col,(c,d0) in enumerate(FIBRES):
        k=[[Q((i,j)==(c,d0) or (j,i)==(c,d0)) for j in range(5)] for i in range(5)]
        for row,(i,j) in enumerate(FIBRES):
            v=sum((ric[i][e]*SIG[e]*k[e][j]+ric[j][e]*SIG[e]*k[e][i] for e in range(5)),Q(0))/2
            v-=sum((R(i,e,j,f)*SIG[e]*SIG[f]*k[e][f] for e,f in product(range(5),repeat=2)),Q(0))
            out[row][col]=v
    return out

CHANNELS=['K','L','N','J','v','u','H']
def fibre_embedding(channel):
    vals={k:Q(0) for k in FIBRES}
    if channel=='K':vals[0,0]=1
    if channel=='L':vals[0,1]=1
    if channel=='N':vals[1,1]=1
    if channel=='J':
        for i in (2,3,4):vals[i,i]=1
    if channel=='v':vals[0,2]=1
    if channel=='u':vals[1,2]=1
    if channel=='H':vals[2,2]=1;vals[3,3]=-1
    return [[vals[k]] for k in FIBRES]
def lift_fibre(m):
    result=zeros()
    for inc in CHANNELS:
        ans=mm(m,fibre_embedding(inc));v={k:ans[i][0] for i,k in enumerate(FIBRES)}
        vals={'K':v[0,0],'L':v[0,1],'N':v[1,1],'J':sum(v[i,i] for i in (2,3,4))/3,'v':v[0,2],'u':v[1,2],'H':(v[2,2]-v[3,3])/2}
        recon=add(*(scale(vals[ch],fibre_embedding(ch)) for ch in CHANNELS))
        if recon!=ans:raise ValueError('ISOTROPIC_FIBRE_RECONSTRUCTION')
        for out,c in vals.items():
            for l in (0,2,4):put(result,out,inc,c,l)
            if inc==out=='H':
                for name in ('T:+','T:-'):result[IDX[name]][IDX[name]]=c
    return result

def matrices():
    out={'I':ident()}
    o=omega(1,Q(1),Q(0));out['O']=lift_fibre(o)
    aa=[omega(i,Q(1),Q(0)) for i in (2,3,4)]
    hh=[omega(i,Q(0),Q(1)) for i in (2,3,4)]
    out['Sa']=lift_fibre(add(*(mm(a,a) for a in aa)))
    out['Sh']=lift_fibre(add(*(mm(h,h) for h in hh)))
    out['Sah']=lift_fibre(add(*(add(mm(a,h),mm(h,a)) for a,h in zip(aa,hh))))
    out['Rb']=lift_fibre(curvature(Q(1),Q(0)))
    out['Rd']=lift_fibre(curvature(Q(0),Q(1)))
    out['O2']=mm(out['O'],out['O'])
    for key in ('Delta','Ca','Ch','Aphi','Bphi','Ephi','Pmet','Pg','Pww','W_K','W_L','By_a','By_p','By_sigma1','By_sigma2'):out[key]=zeros()
    for name,i in IDX.items():
        if name.startswith('T:'):lam=6
        else:
            ch,l=name.split(':');l=int(l);lam=l*(l+2)
            if ch in ('v','u'):lam-=2
            if ch=='H':lam-=6
        out['Delta'][i][i]=-Q(lam)
        if not name.startswith('w:'):out['Pmet'][i][i]=1
    for l in (0,2,4):
        lam=l*(l+2);beta=Q(2)-Q(2,3)*lam
        for target,src,c in [('K','v',2*lam),('L','u',lam),('v','K',1),('v','J',-1),('v','H',-beta),('u','L',1),('J','v',-Q(2,3)*lam),('H','v',2)]:put(out['Ca'],target,src,c,l)
        for target,src,c in [('L','v',lam),('N','u',2*lam),('v','L',-1),('u','N',-1),('u','J',-1),('u','H',-beta),('J','u',Q(2,3)*lam),('H','u',-2)]:put(out['Ch'],target,src,c,l)
        for target,src,c in [('L','w',2),('w','L',1)]:put(out['Aphi'],target,src,c,l)
        for target,src,c in [('K','w',-4),('w','K',Q(1,2)),('w','N',Q(1,2)),('w','J',-Q(3,2))]:put(out['Bphi'],target,src,c,l)
        for target,src,c in [('v','w',2),('w','v',lam)]:put(out['Ephi'],target,src,c,l)
        for target,c in [('K',1),('N',-1),('J',1)]:put(out['Pg'],target,'w',c,l)
        put(out['Pww'],'w','w',1,l);put(out['W_K'],'w','K',1,l);put(out['W_L'],'w','L',1,l)
        for target,src,c in [('N','K',-1),('J','K',1),('K','K',-4),('K','N',-2),('K','J',6)]:put(out['By_a'],target,src,c,l)
        put(out['By_p'],'w','K',Q(1,2),l)
        for target,c in [('N',-Q(1,3)),('J',Q(1,3)),('K',Q(4,3))]:put(out['By_sigma1'],target,'w',c,l)
        put(out['By_sigma2'],'w','w',-Q(1,2),l)
    out['PD']=zeros();out['PN']=ident()
    for name,i in IDX.items():
        if name.startswith(('L:','v:')):out['PD'][i][i]=1;out['PN'][i][i]=0
    return out

def assemble(ell,rho,a,h,c,Ht,p,py,b,d,V,V1,V2,sigma1,sigma2):
    """Accept Fraction for exact arithmetic or ordinary floats for solver use.

    t=ell*x, c=cosh(t), Ht=tanh(t), h=Ht/rho,
    a=rho_y/rho, b=-(a_y+a*a), d=rho**-2-a*a.
    U_xx-U_rr=A U_x+B U_r+D U+F; PN U_r=M PN U+G.
    """
    m=matrices()
    A=add(scale(-3*ell*Ht,m['I']),scale(-2*ell*rho*a,m['O']),scale(ell*rho*p,m['Aphi']))
    B=add(scale(-3*ell*rho*a,m['I']),scale(ell*rho*p,m['Bphi']))
    Z=add(scale(a*a,add(m['Sa'],scale(-1,m['O2']))),scale(a*h,add(m['Sah'],scale(-3,m['O']))),scale(h*h,m['Sh']))
    C=add(Z,scale(-2*b,m['Rb']),scale(-2*d,m['Rd']),scale(Q(4,3)*V,m['Pmet']),scale(Q(4,3)*V1,m['Pg']),scale(-V2,m['Pww']),scale(-(py+4*a*p),m['W_K']),scale(3*p*h,m['W_L']))
    D=scale(ell*ell,add(scale(1/(c*c),m['Delta']),scale(2*rho/c,add(scale(a,m['Ca']),scale(h,m['Ch']))),scale(rho*rho,C),scale(rho*p/c,m['Ephi'])))
    M=scale(-ell*rho,add(scale(a,m['By_a']),scale(p,m['By_p']),scale(sigma1,m['By_sigma1']),scale(sigma2,m['By_sigma2'])))
    return {'A':A,'B':B,'D':D,'M':M,'PD':m['PD'],'PN':m['PN']}

def gram():
    norms={0:Q(1),2:Q(1,12),4:Q(19,720)}
    out=[]
    for name in LABELS:
        if name.startswith('T:'):g=Q(1,6)
        else:
            ch,l=name.split(':');l=int(l);g=norms[l]
            if ch=='L':g*=2
            if ch=='J':g*=3
            if ch in ('v','u'):g*=2*l*(l+2)
            if ch=='H':g*=Q(2,3)*l*(l+2)*(l*(l+2)-3)
        out.append(g)
    return out

def pack_fields(fields):
    """Pack separated field dictionaries without discarding either TT slot.

    Scalar slots yy,y0,00,scalar are triples; ya,0a are pairs.
    ab is the seven-component angular tensor basis used in angular_exact.
    Missing fields are zero. Supplied shapes are checked explicitly.
    """
    out=[Q(0)]*23
    mapping={'yy':('K',(0,2,4)),'y0':('L',(0,2,4)),
             '00':('N',(0,2,4)),'scalar':('w',(0,2,4)),
             'ya':('v',(2,4)),'0a':('u',(2,4))}
    for key,(ch,degrees) in mapping.items():
        if key in fields:
            if len(fields[key])!=len(degrees):raise ValueError('FIELD_SHAPE_'+key)
            for degree,value in zip(degrees,fields[key]):out[IDX[f'{ch}:{degree}']]=value
    if 'ab' in fields:
        if len(fields['ab'])!=7:raise ValueError('FIELD_SHAPE_ab')
        for name,value in zip(['J:0','J:2','J:4','H:2','H:4','T:+','T:-'],fields['ab']):out[IDX[name]]=value
    if set(fields)-set(mapping)-{'ab'}:raise ValueError('UNKNOWN_FIELD_KEY')
    return out

def bulk_forcing(ell,rho,QEred,Qphi):
    """Natural affine second variations -> equation (3) forcing signs.

    The complete QEred must be supplied by the caller; the current checkpoint
    supplies its matter part only and does not synthesize missing Ricci terms.
    """
    if set(QEred)!={'yy','y0','00','ya','0a','ab'}:raise ValueError('COMPLETE_QERED_REQUIRED')
    out=pack_fields(QEred)
    for i in range(23):out[i]*=-2*ell*ell*rho*rho
    scalar=pack_fields({'scalar':Qphi})
    return [x+ell*ell*rho*rho*y for x,y in zip(out,scalar)]

def boundary_forcing(ell,rho,sources):
    """All physical and normal-gauge sources -> the 18 Robin amplitudes."""
    fields={'yy':tuple(-ell*rho*x for x in sources['GN']),
            '00':tuple(2*ell*rho*x for x in sources['QK00']),
            '0a':tuple(2*ell*rho*x for x in sources['QK0a']),
            'ab':tuple(2*ell*rho*x for x in sources['QKab']),
            'scalar':tuple(ell*rho*x for x in sources['Qp'])}
    return pack_fields(fields)

def pack_initial(data,ell):
    """Return U0,V0 from source_adapters.initial_frame_data's t velocities."""
    def one(fields):
        combined=dict(fields)
        combined['ab']=tuple(combined.pop('ab_trace_amplitude'))+tuple(combined.pop('ab_HessTF'))+tuple(combined.pop('ab_TT'))
        return pack_fields(combined)
    return one(data['values']),[ell*x for x in one(data['t_velocities'])]
