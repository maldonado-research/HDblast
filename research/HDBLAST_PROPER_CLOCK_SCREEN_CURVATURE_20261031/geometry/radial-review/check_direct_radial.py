#!/usr/bin/env python3
"""Direct radial, boosted and Lorentz-tilt curvature oracle; no project imports.

Derived from the retained independent coordinate-curvature witness, adapted
for the radial normal, full radial boost, and all four Lorentz-generator slots.

One rational, curved, nonnormal point with all seven scalar metric amplitudes
and all their independent ordinary two-jets populated. Metric and curvature
are built directly in stereographic coordinates. This is a finite witness,
not the universal angular/radial proof in INDEPENDENT_SCREEN_AUDIT.md.
"""
from fractions import Fraction as F
from itertools import product, permutations
from pathlib import Path
import hashlib
import json
import sys
sys.dont_write_bytecode = True

def need(ok, label):
    if not ok:
        raise ValueError(label)

DIM = 5
ZERO = (0,) * DIM

class J:
    """Ordinary coordinate Taylor coefficients, total degree at most four."""
    def __init__(self, p):
        self.p = {m:F(v) for m,v in p.items() if v}
    @staticmethod
    def const(v): return J({ZERO:v})
    def coerce(self, v): return v if isinstance(v,J) else J.const(v)
    def __add__(self, v):
        v=self.coerce(v); p=dict(self.p)
        for m,c in v.p.items(): p[m]=p.get(m,0)+c
        return J(p)
    __radd__=__add__
    def __neg__(self): return J({m:-c for m,c in self.p.items()})
    def __sub__(self,v): return self+-self.coerce(v)
    def __rsub__(self,v): return self.coerce(v)+-self
    def __mul__(self,v):
        v=self.coerce(v); p={}
        for a,c in self.p.items():
            for b,d in v.p.items():
                m=tuple(x+y for x,y in zip(a,b))
                if sum(m)<=4: p[m]=p.get(m,0)+c*d
        return J(p)
    __rmul__=__mul__
    def inv(self):
        c=self.p.get(ZERO,0); need(c!=0,'JET_INVERSE')
        u=(self-c)*(1/c); term=J.const(1); out=term
        for _ in range(4): term=term*(-u); out+=term
        return out*(1/c)
    def __truediv__(self,v): return self*self.coerce(v).inv()
    def __rtruediv__(self,v): return self.coerce(v)*self.inv()
    def d(self,i):
        p={}
        for m,c in self.p.items():
            if m[i]:
                mm=list(m); mm[i]-=1; p[tuple(mm)]=c*m[i]
        return J(p)
    def value(self): return self.p.get(ZERO,F(0))

class S:
    """Affine-path epsilon series through epsilon squared (no factorials)."""
    __slots__=('v',)
    def __init__(self,a=0,b=0,c=0): self.v=(F(a),F(b),F(c))
    def coerce(self,v): return v if isinstance(v,S) else S(v)
    def __add__(self,v):
        v=self.coerce(v); return S(*(a+b for a,b in zip(self.v,v.v)))
    __radd__=__add__
    def __neg__(self): return S(*(-a for a in self.v))
    def __sub__(self,v): return self+-self.coerce(v)
    def __rsub__(self,v): return self.coerce(v)+-self
    def __mul__(self,v):
        v=self.coerce(v);a,b,c=self.v;d,e,f=v.v
        return S(a*d,a*e+b*d,a*f+b*e+c*d)
    __rmul__=__mul__
    def inv(self):
        a,b,c=self.v;return S(1/a,-b/a**2,b*b/a**3-c/a**2)
    def __truediv__(self,v):return self*self.coerce(v).inv()
    def sqrt(self,root):
        a,b,c=self.v;need(root*root==a,'EXACT_ROOT')
        return S(root,b/(2*root),c/(2*root)-b*b/(8*root**3))

def eps(a,b,c):
    if len({a,b,c})<3:return 0
    return 1 if (a,b,c) in ((0,1,2),(1,2,0),(2,0,1)) else -1

def matmul(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),S())
             for j in range(len(b[0]))] for i in range(len(a))]

def inv_affine(m):
    n=len(m); gi0=[[F(i==j)/m[i][i].v[0] for j in range(n)] for i in range(n)]
    # This oracle's background metric is diagonal; perturbed metric is full.
    need(all(m[i][j].v[0]==0 for i in range(n) for j in range(n) if i!=j),'DIAGONAL_BACKGROUND')
    q=[[m[i][j].v[1] for j in range(n)]for i in range(n)]
    return [[S(gi0[i][j],-gi0[i][i]*q[i][j]*gi0[j][j],
               sum(gi0[i][i]*q[i][k]*gi0[k][k]*q[k][j]*gi0[j][j] for k in range(n))
               -gi0[i][i]*m[i][j].v[2]*gi0[j][j])
             for j in range(n)] for i in range(n)]

def det3(m):
    return sum((eps(i,j,k)*m[0][i]*m[1][j]*m[2][k]
                for i,j,k in product(range(3),repeat=3)),S())

def serial_tensor(t): return [[str(v)for v in row]for row in t]

point=(F(0),F(0),F(1,5),F(2,7),F(-1,3)) # y, local t, stereographic x
coords=[]
for i,p in enumerate(point):
    m=list(ZERO);m[i]=1;coords.append(J({ZERO:p,tuple(m):1}))
y,t,*x=coords
rho=F(7,5)+F(2,5)*y-F(3,22)*y*y
c=F(5,4)+F(3,4)*t+F(5,8)*t*t
R=rho*c
den=1+sum((z*z for z in x),J.const(0)); omega=2/den
n=[(1-sum((z*z for z in x),J.const(0)))/den]+[2*z/den for z in x]
Y=n[0]*n[0]-n[1]*n[1]
v=[Y.d(i+2) for i in range(3)]
ksg={(e,a,b):(-2/den)*((x[b] if e==a else 0)+(x[a] if e==b else 0)-(x[e] if a==b else 0))
     for e,a,b in product(range(3),repeat=3)}
h=[[Y.d(a+2).d(b+2)-sum((ksg[e,a,b]*v[e] for e in range(3)),J.const(0))
    for b in range(3)]for a in range(3)]

MULTI=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
amplitudes={}
jets={}
for j,name in enumerate(('A','B','C','D','E','L','M')):
    vals=tuple(F((-1)**(j+k)*(3*j+2*k+1),11+j+2*k) for k in range(6))
    jets[name]=vals
    out=J.const(0)
    for (ny,nt),a in zip(MULTI,vals):
        mono=J.const(a/F(2 if ny==2 or nt==2 else 1))
        for _ in range(ny):mono*=y
        for _ in range(nt):mono*=t
        out+=mono
    amplitudes[name]=out
A,B,C,D,E,L,M=(amplitudes[k]for k in ('A','B','C','D','E','L','M'))
z=J.const(0)
bg=[[z for _ in range(5)]for _ in range(5)]
q=[[z for _ in range(5)]for _ in range(5)]
ksecond=[[z for _ in range(5)]for _ in range(5)]
bg[0][0]=J.const(1);bg[1][1]=-rho*rho
q[0][0]=A*Y;q[0][1]=q[1][0]=rho*B*Y;q[1][1]=rho*rho*D*Y
for a in range(3):
    bg[a+2][a+2]=R*R*omega*omega
    q[0][a+2]=q[a+2][0]=R*C*v[a]
    q[1][a+2]=q[a+2][1]=rho*R*E*v[a]
    for b in range(3):q[a+2][b+2]=R*R*(L*Y*omega*omega*(a==b)+M*h[a][b])

combined='--combined' in sys.argv
tt=[]
if combined:
    # The constants below reproduce the frozen ambient Weyl templates
    # independently; no angular module is imported.
    eps4={p:(-1)**sum(p[i]>p[j]for i in range(4)for j in range(i+1,4))
          for p in permutations(range(4))}
    weights={(0,1):F(-2,3),(2,3):F(-2,3),(0,2):F(1,3),
             (0,3):F(1,3),(1,2):F(1,3),(1,3):F(1,3)}
    ambient={k:F(0)for k in product(range(4),repeat=4)}
    for (i,j),w in weights.items():
        ambient[i,j,i,j]=ambient[j,i,j,i]=w
        ambient[i,j,j,i]=ambient[j,i,i,j]=-w
    star={k:sum((F(eps4.get((k[0],k[1],a,b),0))*ambient[a,b,k[2],k[3]]/2
                 for a,b in product(range(4),repeat=2)),F(0))for k in ambient}
    dn={(i,a):n[i].d(a+2)for i,a in product(range(4),range(3))}
    for sign in (1,-1):
        tpl=[[sum(((ambient[i,k,j,l]+sign*star[i,k,j,l])/2*n[k]*n[l]*dn[i,a]*dn[j,b]
                   for i,j,k,l in product(range(4),repeat=4)),J.const(0))
               for b in range(3)]for a in range(3)]
        tt.append(tpl)
    ut=F(7,17)
    u=F(5,13)+ut*t+F(11,19)*y+F(2,23)*y*t+F(3,29)*t*t
    for a,b in product(range(3),repeat=2):
        ksecond[a+2][b+2]=R*R*u*(tt[0][a][b]+tt[1][a][b])

g=[[S(bg[a][b].value(),q[a][b].value(),ksecond[a][b].value()/2)for b in range(5)]for a in range(5)]
dg={(i,a,b):S(bg[a][b].d(i).value(),q[a][b].d(i).value(),ksecond[a][b].d(i).value()/2)for i,a,b in product(range(5),repeat=3)}
ddg={(i,j,a,b):S(bg[a][b].d(i).d(j).value(),q[a][b].d(i).d(j).value(),ksecond[a][b].d(i).d(j).value()/2)for i,j,a,b in product(range(5),repeat=4)}
gi=inv_affine(g)
dgi={(i,a,b):-sum((gi[a][k]*dg[i,k,l]*gi[l][b]for k,l in product(range(5),repeat=2)),S())
     for i,a,b in product(range(5),repeat=3)}
Gamma={(a,b,d):sum((gi[a][e]*(dg[b,d,e]+dg[d,b,e]-dg[e,b,d])/2 for e in range(5)),S())
       for a,b,d in product(range(5),repeat=3)}
dGamma={(i,a,b,d):sum((dgi[i,a,e]*(dg[b,d,e]+dg[d,b,e]-dg[e,b,d])/2+
                           gi[a][e]*(ddg[i,b,d,e]+ddg[i,d,b,e]-ddg[i,e,b,d])/2
                          for e in range(5)),S())for i,a,b,d in product(range(5),repeat=4)}
Rup={(a,b,cc,d):dGamma[cc,a,d,b]-dGamma[d,a,cc,b]+
       sum((Gamma[a,cc,e]*Gamma[e,d,b]-Gamma[a,d,e]*Gamma[e,cc,b]for e in range(5)),S())
     for a,b,cc,d in product(range(5),repeat=4)}
Rlow={(a,b,cc,d):sum((g[a][e]*Rup[e,b,cc,d]for e in range(5)),S())for a,b,cc,d in product(range(5),repeat=4)}
Ric={(b,d):sum((Rup[a,b,a,d]for a in range(5)),S())for b,d in product(range(5),repeat=2)}
scal=sum((gi[a][b]*Ric[a,b]for a,b in product(range(5),repeat=2)),S())
P={(a,b):(Ric[a,b]-scal*g[a][b]/8)/3 for a,b in product(range(5),repeat=2)}
W={(a,b,cc,d):Rlow[a,b,cc,d]-g[a][cc]*P[b,d]+g[a][d]*P[b,cc]+g[b][cc]*P[a,d]-g[b][d]*P[a,cc]
   for a,b,cc,d in product(range(5),repeat=4)}
need(all(z.v[0]==0 for z in W.values()),'625_BACKGROUND_WEYL_COMPONENTS')

sp=[0,2,3,4];gam=[[g[i][j]for j in sp]for i in sp];gami=inv_affine(gam)
bt=[g[1][i]for i in sp]
beta=[sum((gami[i][j]*bt[j]for j in range(4)),S())for i in range(4)]
N2=-g[1][1]+sum((bt[i]*beta[i]for i in range(4)),S())
Ni=N2.sqrt(rho.value()).inv()
U=[S()for _ in range(5)];U[1]=Ni
for i,ii in enumerate(sp):U[ii]=-Ni*beta[i]
Utime=list(U)
# Registered outward spatial y-normal within t slice.
U=[S()for _ in range(5)]
radial_denom=gami[0][0].sqrt(F(1))
for i,ii in enumerate(sp):U[ii]=gami[i][0]/radial_denom
screen=[[g[a+2][b+2]for b in range(3)]for a in range(3)]
rootdet=det3(screen).sqrt((R.value()*omega.value())**3)
epsilon={(a,cc,d):sum((screen[a][e]*eps(e,cc,d)for e in range(3)),S())/rootdet
         for a,cc,d in product(range(3),repeat=3)}

def projected(curv, omit_u1=False,omit_eps1=False):
    out=[[F(0)for _ in range(3)]for _ in range(3)]
    first=[[F(0)for _ in range(3)]for _ in range(3)]
    for a,b in product(range(3),repeat=2):
        total=S()
        for aa,bb in ((a,b),(b,a)):
            for cc,d,mu in product(range(3),range(3),range(5)):
                e=epsilon[aa,cc,d];u=U[mu]
                if omit_u1:u=S(u.v[0],0,u.v[2])
                if omit_eps1:e=S(e.v[0],0,e.v[2])
                total+=e*curv[bb+2,mu,cc+2,d+2]*u/4
        first[a][b]=total.v[1]/(R.value()*omega.value())**2
        out[a][b]=2*total.v[2]/(R.value()*omega.value())**2
    return first,out

H1,H2=projected(W)
need(all(z==0 for row in H1 for z in row),'ENTIRE_SCALAR_H1_ZERO')
need(H2==projected(Rlow)[1],'EXACT_SCHOUTEN_SYMMETRIC_CANCELLATION')
values={k:val[0]for k,val in jets.items()}
A0,B0,C0,D0,E0,L0,M0=(values[k]for k in ('A','B','C','D','E','L','M'))
rho0=rho.value();R0=R.value();theta=F(3,5)
Ct=jets['C'][2];Ey=jets['E'][1];My=jets['M'][1];Mt=jets['M'][2];Lt=jets['L'][2]
eta=R0*(Ct-theta*C0)/(2*rho0)-R0*Ey/2-B0/2
Ly=jets['L'][1]
coef=((-A0+4*M0-L0)*My-2*M0*Ly)/4+C0*(A0-L0)/(2*R0)+(E0/R0-Mt/(2*rho0))*eta
axial=[[sum((F(eps(a,cc,d))*v[cc].value()*h[b][d].value()+
              F(eps(b,cc,d))*v[cc].value()*h[a][d].value()
              for cc,d in product(range(3),repeat=2)),F(0))/(2*omega.value()**3)
        for b in range(3)]for a in range(3)]
target=[[2*coef/R0*z for z in row]for row in axial]
checks=['625 background Weyl slots vanish','entire symmetric first H vanishes',
 'exact Schouten cancellation']
if combined:
    linear=[[3*F(11,19)/(2*R0)*(tt[0][a][b].value()-tt[1][a][b].value())/omega.value()**2
             for b in range(3)]for a in range(3)]
    wrong=[[target[a][b]+2*linear[a][b]for b in range(3)]for a in range(3)]
    target=[[target[a][b]+linear[a][b]for b in range(3)]for a in range(3)]
    need(H2!=wrong,'SECOND_DERIVATIVE_FACTOR_TWO_MUTATION_REJECTED')
    checks+=['direct frozen equal-TT second metric with epsilon squared k/2 normalization',
             'direct linear curl eigenvalues +3,-3 and radial coefficient uy/(2R)',
             'second-derivative factor-two mutation rejected']
need(H2==target,'FULL_DIRECT_WEYL_HESSIAN_EQUALS_COMPACT_TOTAL')
need(H2!=projected(W,omit_u1=True)[1],'OMITTED_OBSERVER_CORRECTION_REJECTED')
need(H2!=projected(W,omit_eps1=True)[1],'OMITTED_VOLUME_SCREEN_CORRECTION_REJECTED')
need(any(z!=0 for row in H2 for z in row),'GENERIC_SCALAR_QUADRATIC_NONZERO_WITNESS')
checks+=['complete direct five-dimensional curvature Hessian equals compact formula',
 'registered radial n1 omission rejected','screen epsilon1 omission rejected','generic scalar quadratic witness nonzero']
result={'status':'PASS_RADIAL_DIRECT_COMBINED_CURVATURE_WITNESS'if combined else 'PASS_RADIAL_DIRECT_CURVATURE_WITNESS',
 'arithmetic':'fractions.Fraction only; no project imports',
 'checks':checks,
 'coefficient_F':str(coef),'eta_Knormal':str(eta),
 'H2_orthonormal':serial_tensor(H2),
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'limits':['Finite generic curved nonnormal coordinate witness, not all-point proof.',
 'Registered radial and boosted comparison-frame calculations plus first relative Lorentz tilt.',
 'No actual first-wave values, physical observer transport, full nonlinear family or cosmology claim.']}
# First-Weyl scalar Hessian slots for the relative Lorentz-tilt calculation.
h01=h[0][1].value()/omega.value()**2
h12=h[1][2].value()/omega.value()**2
orth_ang=R0*omega.value()
gW=W[2,3,2,4].v[1]/orth_ang**4/h12
etW=W[2,1,3,1].v[1]/(orth_ang**2*rho0**2)/h01
eyW=W[2,0,3,0].v[1]/orth_ang**2/h01
kW=W[2,1,3,0].v[1]/(orth_ang**2*rho0)/h01
H=theta/rho0;a=rho.d(0).value()/rho0
m=Mt/(2*rho0)-E0/R0;j=My/2-C0/R0
mt=jets['M'][5]/(2*rho0)-(jets['E'][2]-theta*E0)/R0
my=(jets['M'][4]-a*Mt)/(2*rho0)-(jets['E'][1]-a*E0)/R0
jy=jets['M'][3]/2-(jets['C'][1]-a*C0)/R0
gRc=-L0/(2*R0**2)+H*m-a*j
etRc=-mt/rho0-2*H*m-D0/(2*R0**2)+a*j
eyRc=-jy-2*a*j-A0/(2*R0**2)
kRc=eta/R0**2-my-a*m-H*j
need(gW+etW==gRc+etRc,'FIRST_WEYL_G_PLUS_ET')
need(eyW-gW==eyRc-gRc,'FIRST_WEYL_EY_MINUS_G')
need(kW==kRc,'FIRST_WEYL_MIXED_HESSIAN')

# Transform the full first Weyl tensor by a relative orthonormal Lorentz
# generator, including every observer and screen slot. Angular v is measured
# in the unit-round-sphere orthonormal basis, matching the p,q definitions.
cboost,sboost=F(5,4),F(3,4)
p_tilt,q_tilt=F(2,7),F(-3,11)
frame_scales=[F(1),1/rho0]+[1/orth_ang]*3
Wframe={inds:W[inds].v[1]*frame_scales[inds[0]]*frame_scales[inds[1]]*frame_scales[inds[2]]*frame_scales[inds[3]] for inds in product(range(5),repeat=4)}
def vv_sum(*terms):
    return [sum((coef*vec[i] for coef,vec in terms),F(0))for i in range(5)]
def unit(i):return [F(i==j)for j in range(5)]
Ubar=[sboost,cboost,F(0),F(0),F(0)]
nbar=[cboost,sboost,F(0),F(0),F(0)]
Es=[unit(i+2)for i in range(3)]
vs=[v[i].value()/omega.value()for i in range(3)]
def curvature(a,b,c,d):
    vectors=(a,b,c,d)
    return sum((Wframe[ix]*vectors[0][ix[0]]*vectors[1][ix[1]]*vectors[2][ix[2]]*vectors[3][ix[3]]
                for ix in product(*(tuple(i for i,z in enumerate(vec)if z)for vec in vectors))),F(0))
def tilt_tensor(omit_screen=False, q_sign=1):
    du=vv_sum(*((p_tilt*vs[i],Es[i])for i in range(3)))
    de=[vv_sum((p_tilt*vs[i],Ubar),(q_sign*q_tilt*vs[i],nbar))for i in range(3)]
    out=[[F(0)for _ in range(3)]for _ in range(3)]
    for aa,bb in product(range(3),repeat=2):
        for a,b in ((aa,bb),(bb,aa)):
            for c,d in product(range(3),repeat=2):
                sign=eps(a,c,d)
                if not sign:continue
                value=curvature(Es[b],du,Es[c],Es[d])
                if not omit_screen:
                    value+=curvature(de[b],Ubar,Es[c],Es[d])
                    value+=curvature(Es[b],Ubar,de[c],Es[d])
                    value+=curvature(Es[b],Ubar,Es[c],de[d])
                out[aa][bb]+=F(sign,4)*value
    return out
G=cboost*cboost*(gW+etW)+sboost*sboost*(eyW-gW)+2*cboost*sboost*kW
K=cboost*sboost*(etW+eyW)+(cboost*cboost+sboost*sboost)*kW
tilt_expected=[[-(p_tilt*G+q_tilt*K)*z for z in row]for row in axial]
tilt_direct=tilt_tensor()
need(tilt_direct==tilt_expected,'FULL_RELATIVE_LORENTZ_GENERATOR_TILT')
need(tilt_tensor(omit_screen=True)!=tilt_expected,'OMITTED_SCREEN_GENERATOR_REJECTED')
need(tilt_tensor(q_sign=-1)!=tilt_expected,'WRONG_RADIAL_SCREEN_TILT_SIGN_REJECTED')

# Direct full boosted second tensor, including optional TT metric second jet.
Uradial=list(U)
radial_spacetime_norm=gi[0][0].sqrt(F(1))
U=[gi[i][0]/radial_spacetime_norm for i in range(5)]
need(projected(W)[1]==H2,'RADIAL_HYPERSURFACE_AND_REGISTERED_NORMAL_EQUIVALENCE')
U=[cboost*Utime[i]+sboost*Uradial[i]for i in range(5)]
boost_first,boost_direct=projected(W)
need(all(z==0 for row in boost_first for z in row),'BOOSTED_SCALAR_H1_ZERO')
Ft=((D0+4*M0-L0)*Mt-2*M0*Lt)/(4*rho0)-E0*(D0+L0)/(2*R0)+(C0/R0-My/2)*(eta+B0)
boost_factor=2*(cboost*Ft+sboost*coef)/R0
if combined:boost_factor+=3*(cboost*ut/rho0+sboost*F(11,19))/(4*R0)
boost_expected=[[boost_factor*z for z in row]for row in axial]
need(boost_direct==boost_expected,'FULL_BOOSTED_SECOND_TENSOR')
U=Uradial
result['checks']+=['first Weyl Hessian invariant g+et','first Weyl Hessian invariant ey-g','first Weyl mixed Hessian k',
 'full relative Lorentz generator with all four curvature slots','omitted screen generator mutation rejected',
 'wrong radial screen tilt sign mutation rejected','full direct boosted second tensor','radial hypersurface and registered normal equivalence']
result['first_Weyl']={k:str(v)for k,v in dict(g=gW,et=etW,ey=eyW,k=kW).items()}
result['candidate_Riemann']={k:str(v)for k,v in dict(g=gRc,et=etRc,ey=eyRc,k=kRc).items()}
result['relative_tilt']={k:str(v)for k,v in dict(c=cboost,s=sboost,p=p_tilt,q=q_tilt,G=G,K=K).items()}
result['tilt_first_operator_tensor']=serial_tensor(tilt_direct)
result['boosted_second_tensor']=serial_tensor(boost_direct)
if __name__=='__main__':print(json.dumps(result,sort_keys=True))
