#!/usr/bin/env python3
"""Independent exact frame-index audit of the direct gauge-constraint formulas."""
from fractions import Fraction as Q
import json


def need(ok,message):
    if not ok:raise RuntimeError(message)

SIG=(1,-1,1,1,1)

def zero():return [[Q(0)]*5 for _ in range(5)]

def sphere_jet(l):
    # At the north pole of X0=sqrt(1-|z|^2), Xi=zi.  These are actual
    # restrictions of harmonic ambient polynomials:
    # l2: X0^2-X1^2+X0X1;
    # l4: Re(X0+iX1)^4+Im(X0+iX1)^4/4.
    if l==0:return Q(1),[Q(0)]*3,[[Q(0)]*3 for _ in range(3)],[[[Q(0)]*3 for _ in range(3)] for _ in range(3)]
    Y=Q(1);grad=[Q(1),Q(0),Q(0)]
    hdiag=(-4,-2,-2) if l==2 else (-16,-4,-4)
    hess=[[Q(hdiag[i] if i==j else 0) for j in range(3)] for i in range(3)]
    third=[[[Q(0)]*3 for _ in range(3)] for _ in range(3)]
    c111,c122=(-3,-1) if l==2 else (-15,-3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ix=sorted((i,j,k))
                partial=Q(c111 if ix==[0,0,0] else c122 if ix in ([0,1,1],[0,2,2]) else 0)
                # d_i Hess_jk=Y_ijk-(d_i Gamma^c_jk)Y_c,
                # with d_i Gamma^c_jk=delta_ic delta_jk at north pole.
                third[i][j][k]=partial-(grad[i] if j==k else 0)
    return Y,grad,hess,third


def assemble_metric(l,values,dy,dt,rho,invc):
    lam=l*(l+2);Y,g,h,T=sphere_jet(l)
    k=zero();partial=[zero() for _ in range(5)]
    def put(a,b,value,y,t,angular):
        k[a][b]=k[b][a]=value
        partial[0][a][b]=partial[0][b][a]=y
        partial[1][a][b]=partial[1][b][a]=t/rho
        for i in range(3):partial[i+2][a][b]=partial[i+2][b][a]=angular[i]*invc/rho
    for a,b,ch in ((0,0,'K'),(0,1,'L'),(1,1,'N')):
        put(a,b,values[ch]*Y,dy[ch]*Y,dt[ch]*Y,[values[ch]*z for z in g])
    for i in range(3):
        for a,ch in ((0,'v'),(1,'u')):
            put(a,i+2,values[ch]*g[i],dy[ch]*g[i],dt[ch]*g[i],[values[ch]*h[j][i] for j in range(3)])
    for i in range(3):
        for j in range(i,3):
            tensor=h[i][j]+(Q(lam,3)*Y if i==j else 0)
            value=(values['J']*Y if i==j else 0)+values['H']*tensor
            yd=(dy['J']*Y if i==j else 0)+dy['H']*tensor
            td=(dt['J']*Y if i==j else 0)+dt['H']*tensor
            angular=[(values['J']*g[a] if i==j else 0)+values['H']*(T[a][i][j]+(Q(lam,3)*g[a] if i==j else 0)) for a in range(3)]
            put(i+2,j+2,value,yd,td,angular)
    return k,partial


def index_DC(k,partial,a,h):
    # Gamma[A][B][C] is coefficient C of nabla_(E_A) E_B.
    gamma=[zero() for _ in range(5)]
    gamma[1][0][1]=gamma[1][1][0]=a
    for b in range(2,5):
        gamma[b][0][b]=a;gamma[b][1][b]=h
        gamma[b][b][0]=-a;gamma[b][b][1]=h
    cov=[zero() for _ in range(5)]
    for A in range(5):
        for B in range(5):
            for C in range(5):
                cov[A][B][C]=partial[A][B][C]-sum(gamma[A][B][d]*k[d][C]+gamma[A][C][d]*k[B][d] for d in range(5))
    return [sum(SIG[A]*cov[A][A][B] for A in range(5))-sum(SIG[A]*partial[B][A][A] for A in range(5))/2 for B in range(5)]


def modal_DC(l,value,dy,dt,rho,invc,a,h,drop_PD_normal=False):
    K,N,J,L,v,u,H=[value[ch] for ch in ('K','N','J','L','v','u','H')]
    lam=l*(l+2);beta=Q(2)-Q(2,3)*lam
    R=(dy['K']+dy['N']-3*dy['J'])/2-dt['L']/rho-lam*invc*v/rho+4*a*K+a*N-3*a*J-3*h*L
    T=(0 if drop_PD_normal else dy['L'])-(dt['K']+dt['N']+3*dt['J'])/(2*rho)-lam*invc*u/rho+5*a*L-3*h*(N+J)
    A=(0 if drop_PD_normal else dy['v'])-dt['u']/rho+5*a*v-4*h*u+invc/rho*((-K+N-J)/2+beta*H)
    Y,g,_,_=sphere_jet(l)
    return [R*Y,T*Y]+[A*z for z in g]


channels=('K','L','N','J','v','u','H')
rho=Q(7,3);invc=Q(8,9);a=Q(5,11);h=Q(13,17)
groups=[];cases=0
for l in (0,2,4):
    for family in range(3):
        for channel in channels:
            if not l and channel in ('v','u','H'):continue
            data=[{ch:Q(0) for ch in channels} for _ in range(3)]
            data[family][channel]=Q(5,7)
            k,dp=assemble_metric(l,*data,rho,invc)
            need(index_DC(k,dp,a,h)==modal_DC(l,*data,rho,invc,a,h),'full five-index contraction agrees:'+str((l,family,channel)))
            cases+=1
    groups.append('five_index_contraction_degree_'+str(l))

# Actual spherical Hessian divergence constants at both nonconstant degrees.
for l in (2,4):
    _,grad,_,third=sphere_jet(l);lam=l*(l+2)
    need([sum(third[i][i][b] for i in range(3)) for b in range(3)]==[(2-lam)*z for z in grad],'Ricci-correct Hessian divergence')
    beta=Q(2)-Q(2,3)*lam
    need([(2-lam)*z+Q(lam,3)*z for z in grad]==[beta*z for z in grad],'tracefree Hessian divergence')
groups.append('actual_sphere_harmonic_jet_Ricci_commutation')

# PD value/tangential zero does not force PD normal derivative zero.
v={ch:Q(0) for ch in channels};y=v.copy();t=v.copy();y['L']=Q(3);y['v']=Q(5)
proper=modal_DC(2,v,y,t,rho,invc,a,h)
wrong=modal_DC(2,v,y,t,rho,invc,a,h,drop_PD_normal=True)
need(proper!=wrong and proper[1]!=0 and proper[2]!=0,'dropped Dirichlet-normal corruption rejected')
groups.append('Dirichlet_normal_derivatives_retained')

# Generic homogeneous shell normal identity, including the scalar tension row.
K,N,J,w=map(Q,(2,3,5,7));ell=Q(1,1000);sigma1=Q(11);alpha=ell*rho*a
Kr=alpha*(4*K+2*N-6*J)-Q(4,3)*ell*rho*sigma1*w
Nr=alpha*K+ell*rho*sigma1*w/3
Jr=-alpha*K-ell*rho*sigma1*w/3
need(-(Kr+Nr-3*Jr)/2+alpha*(4*K+N-3*J)==0,'homogeneous normal gauge boundary identity')
need(-(Kr+Nr-3*Jr)/2-alpha*(4*K+N-3*J)!=0,'reversed radial orientation rejected')
groups.append('normal_boundary_row_cancellation_with_scalar_coupling')

# Exact second variation of the full metric-lowered covector.  The Christoffel
# numerator S is linear in the covariant metric first jet; only one inverse
# metric remains after lowering C^A by the FULL metric.
q=[[Q((i+2)*(j+2),13) for j in range(5)] for i in range(5)]
dq=[zero() for _ in range(5)]
dq[0][1][1]=dq[0][2][2]=Q(1) # nontrivial first jet with DC[q]=0
S=lambda B,C,D:(dq[C][B][D]+dq[D][B][C]-dq[B][C][D])/2
first=[sum(SIG[C]*S(B,C,C) for C in range(5)) for B in range(5)]
need(first==[0]*5,'manufactured first field is harmonic')
QC=[-2*sum(SIG[C]*SIG[D]*q[C][D]*S(B,C,D) for C in range(5) for D in range(5)) for B in range(5)]
need(any(QC),'nontrivial second-parameter quadratic source')
QCvec=[]
for A in range(5):
    cross=-2*sum(q[A][D]*SIG[A]*SIG[D]*first[D] for D in range(5))
    nonlinear=-2*sum(SIG[C]*SIG[D]*q[C][D]*SIG[A]*S(A,C,D) for C in range(5) for D in range(5))
    QCvec.append(cross+nonlinear)
need([SIG[A]*QCvec[A] for A in range(5)]==QC,'background lowering equals full lowering when first constraint vanishes')
need([z/2 for z in QC]!=QC,'eps-squared instead of second-variation corruption rejected')
k1=[Q(3),Q(5),Q(7),Q(11),Q(13)];k2=[Q(17),Q(19),Q(23),Q(29),Q(31)]
need([(a+q0)-(b+q0) for a,b,q0 in zip(k1,k2,QC)]==[a-b for a,b in zip(k1,k2)],'same fixed QC cancellation')
groups.append('full_lowered_second_variation_and_common_source_cancellation')
print(json.dumps({'status':'PASS','named_groups':groups,'group_count':len(groups),'exact_field_jet_cases':cases,
 'wrong_formula_controls':['dropped_Dirichlet_normal','reversed_radial_orientation','eps_squared_instead_of_second_variation'],
 'scope':'independent exact five-frame-index, actual sphere polynomial jet, common quadratic source and homogeneous shell-normal algebra; not a source point evaluation or full Einstein constraint certificate'},sort_keys=True))
