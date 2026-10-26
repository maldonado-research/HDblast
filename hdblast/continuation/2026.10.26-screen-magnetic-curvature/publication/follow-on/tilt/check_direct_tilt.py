#!/usr/bin/env python3
"""Independent rational coordinate two-jet curvature fixtures for observer tilt.

Standalone stdlib; no accepted ADM coefficient code is imported. Synthetic
metric jets only, not actual fields or a new Einstein-matter solution.
"""
from fractions import Fraction as F
from itertools import product,permutations
import json

N=5;I=range(N);ZERO=(0,)*N;DEG=4

def add(*ps):
 out={}
 for p in ps:
  for e,c in p.items():out[e]=out.get(e,F(0))+c
 return {e:c for e,c in out.items() if c}
def scale(c,p):return {e:c*v for e,v in p.items() if c*v}
def mul(p,q):
 out={}
 for e,c in p.items():
  for f,d in q.items():
   g=tuple(x+y for x,y in zip(e,f))
   if sum(g)<=DEG:out[g]=out.get(g,F(0))+c*d
 return {e:c for e,c in out.items() if c}
def diff(p,k):
 out={}
 for e,c in p.items():
  if e[k]:f=list(e);f[k]-=1;out[tuple(f)]=c*e[k]
 return out
def value(p):return p.get(ZERO,F(0))
ONE={ZERO:F(1)}
X=[]
for i in I:
 e=[0]*N;e[i]=1;X.append({tuple(e):F(1)})

class D:
 __slots__=('a','b')
 def __init__(self,a=0,b=0):self.a=F(a);self.b=F(b)
 @staticmethod
 def cv(x):return x if isinstance(x,D) else D(x)
 def __add__(self,x):x=D.cv(x);return D(self.a+x.a,self.b+x.b)
 __radd__=__add__
 def __neg__(self):return D(-self.a,-self.b)
 def __sub__(self,x):return self+-D.cv(x)
 def __rsub__(self,x):return D.cv(x)+-self
 def __mul__(self,x):x=D.cv(x);return D(self.a*x.a,self.b*x.a+self.a*x.b)
 __rmul__=__mul__
 def __truediv__(self,x):
  x=D.cv(x)
  return self*D(1/x.a,-x.b/(x.a*x.a))

def eps(a,c,d):
 if len({a,c,d})<3:return 0
 return (-1)**sum(p>q for i,p in enumerate((a,c,d)) for q in (a,c,d)[i+1:])

def fixture(seed):
 rho=F(3,2);a=F(2,7);rhoy=a*rho;rhoyy=F(-4,9)
 c=F(5,4);theta=F(3,5);R=rho*c
 rj=add(scale(rho,ONE),scale(rhoy,X[0]),scale(rhoyy/2,mul(X[0],X[0])))
 cj=add(scale(c,ONE),scale(c*theta,X[1]),scale(c/2,mul(X[1],X[1])))
 Rj=mul(rj,cj)
 z2=add(*(mul(X[i],X[i]) for i in range(2,5)))
 den=add(ONE,scale(-1,z2),mul(z2,z2))
 normal=[F(3,5),0,F(4,5),0]
 triad=[[F(4,5),0,F(-3,5),0],[0,1,0,0],[0,0,0,1]]
 ambient=[mul(add(scale(normal[j],add(ONE,scale(-1,z2))),*(scale(2*triad[k][j],X[k+2]) for k in range(3))),den) for j in range(4)]
 Y=add(mul(ambient[0],ambient[0]),scale(-1,mul(ambient[1],ambient[1])))
 grad={i:diff(Y,i) for i in range(2,5)}
 loggrad={i:scale(-2,mul(X[i],den)) for i in range(2,5)}
 kappa={(i,j):scale(4,add(ONE,scale(-2,z2),scale(3,mul(z2,z2)))) if i==j else {} for i,j in product(range(2,5),repeat=2)}
 spheregamma={(e,i,j):add(loggrad[i] if e==j else {},loggrad[j] if e==i else {},scale(-1,loggrad[e]) if i==j else {}) for e,i,j in product(range(2,5),repeat=3)}
 hess={(i,j):add(diff(grad[j],i),scale(-1,add(*(mul(spheregamma[e,i,j],grad[e]) for e in range(2,5))))) for i,j in product(range(2,5),repeat=2)}
 jets={};amps={}
 for idx,key in enumerate('ABCDELM'):
  row=[F(((seed+2)*(idx+3)*(k+2)+3*k-7)%23-11,7+k) for k in range(6)]
  jets[key]=row
  amps[key]=add(scale(row[0],ONE),scale(row[1],X[0]),scale(row[2],X[1]),scale(row[3]/2,mul(X[0],X[0])),scale(row[4],mul(X[0],X[1])),scale(row[5]/2,mul(X[1],X[1])))
 bg={(i,j):{} for i,j in product(I,repeat=2)};q={k:{} for k in bg}
 bg[0,0]=ONE;bg[1,1]=scale(-1,mul(rj,rj))
 q[0,0]=mul(amps['A'],Y);q[0,1]=q[1,0]=mul(mul(rj,amps['B']),Y);q[1,1]=mul(mul(mul(rj,rj),amps['D']),Y)
 for i in range(2,5):
  q[0,i]=q[i,0]=mul(mul(Rj,amps['C']),grad[i])
  q[1,i]=q[i,1]=mul(mul(mul(rj,Rj),amps['E']),grad[i])
  for j in range(2,5):
   bg[i,j]=mul(mul(Rj,Rj),kappa[i,j])
   q[i,j]=mul(mul(Rj,Rj),add(mul(mul(amps['L'],Y),kappa[i,j]),mul(amps['M'],hess[i,j])))
 g={k:D(value(bg[k]),value(q[k])) for k in bg}
 dg={(k,i,j):D(value(diff(bg[i,j],k)),value(diff(q[i,j],k))) for k,i,j in product(I,repeat=3)}
 ddg={(k,l,i,j):D(value(diff(diff(bg[i,j],l),k)),value(diff(diff(q[i,j],l),k))) for k,l,i,j in product(I,repeat=4)}
 inv0={i:1/g[i,i].a for i in I}
 inv={(i,j):D(inv0[i] if i==j else 0,-inv0[i]*qv.b*inv0[j]) for (i,j),qv in g.items()}
 gamma={(e,i,j):sum((inv[e,k]*(dg[i,j,k]+dg[j,i,k]-dg[k,i,j])/2 for k in I),D()) for e,i,j in product(I,repeat=3)}
 riem={}
 for aa,bb,cc,dd in product(I,repeat=4):
  riem[aa,bb,cc,dd]=(ddg[bb,cc,aa,dd]+ddg[aa,dd,bb,cc]-ddg[bb,dd,aa,cc]-ddg[aa,cc,bb,dd])/2+sum((g[e,f]*(gamma[e,aa,dd]*gamma[f,bb,cc]-gamma[e,aa,cc]*gamma[f,bb,dd]) for e,f in product(I,repeat=2)),D())
 ric={(b,d):sum((inv[a,c0]*riem[a,b,c0,d] for a,c0 in product(I,repeat=2)),D()) for b,d in product(I,repeat=2)}
 scal=sum((inv[b,d]*ric[b,d] for b,d in product(I,repeat=2)),D())
 sch={(b,d):(ric[b,d]-scal*g[b,d]/8)/3 for b,d in product(I,repeat=2)}
 weyl={(a,b,c0,d):riem[a,b,c0,d]-g[a,c0]*sch[b,d]+g[a,d]*sch[b,c0]+g[b,c0]*sch[a,d]-g[b,d]*sch[a,c0] for a,b,c0,d in product(I,repeat=4)}
 if any(x.a for x in weyl.values()):raise ValueError('BACKGROUND_WEYL_NONZERO')
 frame=[F(1),1/rho,1/(2*R),1/(2*R),1/(2*R)]
 C={k:x.b*frame[k[0]]*frame[k[1]]*frame[k[2]]*frame[k[3]] for k,x in weyl.items()}
 v={i:value(grad[i])/2 for i in range(2,5)};h={(i,j):value(hess[i,j])/4 for i,j in product(range(2,5),repeat=2)}
 ax={(aa,bb):sum((F(eps(aa,cc,dd))*v[cc]*h[bb,dd] for cc,dd in product(range(2,5),repeat=2)),F(0)) for aa,bb in product(range(2,5),repeat=2)}
 ax={(aa,bb):(ax[aa,bb]+ax[bb,aa])/2 for aa,bb in product(range(2,5),repeat=2)}
 P=(jets['M'][5]+theta*jets['M'][2])/(2*rho*rho)-jets['E'][2]/(rho*R)+(jets['D'][0]+jets['L'][0])/(2*R*R)
 Q=(jets['M'][4]+theta*jets['M'][1])/(2*rho)-jets['E'][1]/(2*R)-(jets['C'][2]+theta*jets['C'][0])/(2*rho*R)+jets['B'][0]/(2*R*R)
 if not any(ax.values()):raise ValueError('DEGENERATE_ANGULAR_FIXTURE')
 if not P or not Q:raise ValueError('GENERIC_NEGATIVE_CONTROL_MUST_BE_NONZERO')
 tests=[]
 for b0,V,W,label in [(F(2,3),0,0,'radial_boost'),(0,F(3,7),0,'angular_observer'),(0,0,F(-2,5),'radial_screen'),(F(3,2),F(4,7),F(9,11),'combined')]:
  du={0:b0*value(Y),1:F(0),**{i:V*v[i] for i in range(2,5)}}
  de={i:{0:-W*v[i],1:V*v[i],2:F(0),3:F(0),4:F(0)} for i in range(2,5)}
  raw={}
  for aa,bb in product(range(2,5),repeat=2):
   raw[aa,bb]=sum((F(eps(aa,cc,dd))*sum((de[bb][e]*C[e,1,cc,dd]+du[e]*C[bb,e,cc,dd]+de[cc][e]*C[bb,1,e,dd]+de[dd][e]*C[bb,1,cc,e] for e in I),F(0)) for cc,dd in product(range(2,5),repeat=2)),F(0))
  for aa,bb in product(range(2,5),repeat=2):
   observed=(raw[aa,bb]+raw[bb,aa])/2
   expected=2*(V*P-W*Q)*ax[aa,bb]
   if observed!=expected:raise ValueError(('TILT_MISMATCH',seed,label,aa,bb,str(observed),str(expected)))
  tests.append(label)
 return {'seed':seed,'checks':tests,'P':str(P),'Q':str(Q),'all_screen_components_checked':9,'background_Weyl_zero':True}

if __name__=='__main__':
 rows=[fixture(s) for s in (1,4,9)]
 print(json.dumps({'status':'PASS_INDEPENDENT_DIRECT_COORDINATE_SCALAR_TILT_FIXTURES','arithmetic':'exact rational metric two-jets and parameter dual numbers','sphere_point':['3/5','0','4/5','0'],'rows':rows,'actual_field_evaluated':False,'field_equations_used':False},indent=2))
