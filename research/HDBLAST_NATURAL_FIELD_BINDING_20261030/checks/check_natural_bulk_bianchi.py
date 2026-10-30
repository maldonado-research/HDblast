#!/usr/bin/env python3
"""Private exact source-map Bianchi audit, coherent free ordinary jets."""
import argparse, ast, hashlib, importlib.util, json, sys
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
MAP='research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json'
MAP_PIN='34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'
CONTROLS=('reduced_instead_of_natural','forcing_instead_of_natural','half_second_variation','omit_deltaB_E1','omit_S1_deta','positive_time_contraction','missing_trace_half','missing_radial_connection','missing_FRW_connection','wrong_potential_chain','frozen_rho_radial','wrong_scalar_source_sign','premature_harmonic_scalar')
HELPER_PINS={'codazzi_laurent.py':'178b0a4e58db63eae6c24112ed9d91bf430f70dca9f6f8126c5e352e429f8fc2','angular_generic.py':'fd14d8cc3ac6f58c206b784c5c57d15430f63f38ca091ac03c7769f9bfc89df9'}
def need(ok,why):
 if not ok:raise ValueError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name):
 path=HERE/(name+'.py');spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def auth(root):
 out={}
 for name,pin in {MAP:MAP_PIN,**HELPER_PINS}.items():
  p=root/name if name==MAP else HERE/name
  need(p.is_file() and not p.is_symlink(),'SOURCE_TYPE:'+name)
  need(sha(p)==pin,'SOURCE_PIN:'+name);out[name]={'bytes':p.stat().st_size,'sha256':pin}
 for p in (Path(__file__),HERE/'codazzi_laurent.py',HERE/'angular_generic.py'):
  tree=ast.parse(p.read_text());need(not any(isinstance(x,ast.Assert) for x in ast.walk(tree)),'ASSERT_FORBIDDEN:'+p.name)
  permitted={'argparse','ast','hashlib','importlib','json','sys','fractions','itertools','math','pathlib','functools'}
  for node in ast.walk(tree):
   if isinstance(node,ast.Import):need(all(x.name.split('.')[0] in permitted for x in node.names),'IMPORT_CLOSURE:'+p.name)
   if isinstance(node,ast.ImportFrom):need(node.module and node.module.split('.')[0] in permitted,'IMPORT_CLOSURE:'+p.name)
   if isinstance(node,ast.Call) and isinstance(node.func,ast.Name):need(node.func.id not in {'eval','exec','compile','__import__'},'DYNAMIC_CODE_FORBIDDEN:'+p.name)
  out[p.name]={'bytes':p.stat().st_size,'sha256':sha(p)}
 return out

def calculate(root,mutation=None):
 before=auth(root);alg=load('codazzi_laurent');R=load('angular_generic');P=alg.P
 var=P.var;rho,a,p,c,th,Vp,Vpp,Vppp=[var(n) for n in ('rho','a','p','c','theta','Vp','Vpp','Vppp')]
 ay=-rho**-2-p*p/3;py=Vp-4*a*p;pot=6*rho**-2+p*p/2-6*a*a
 def deriv(x,d):
  x=P.cv(x);out=P()
  rules={'rho':a*rho,'a':ay,'p':py,'Vp':Vpp*p,'Vpp':Vppp*p} if d==0 else {'c':c*th,'theta':1-th*th}
  if mutation=='frozen_rho_radial' and d==0:rules['rho']=P()
  if mutation=='wrong_potential_chain' and d==0:rules['Vp']=Vpp
  for n in x.names():
   if n[0] in 'ABCDELMs' and (len(n)==1 or set(n[1:])<=set('yt')):
    suffix=n[1:]+('y' if d==0 else 't');nn=n[0]+'y'*suffix.count('y')+'t'*suffix.count('t')
    need(len(suffix)<= (2 if n[0]=='s' else 3),'JET_ORDER:'+nn);v=var(nn)
   else:v=rules.get(n,P())
   out+=x.diff(n)*v
  return out
 def quotient(x):
  out=P()
  for mon,coef in P.cv(x).t.items():
   z=dict(mon);n=z.pop('theta',0);need(n>=0,'NEGATIVE_THETA_POWER')
   for k in range(n//2+1):
    w=dict(z);w['theta']=n%2;w['c']=w.get('c',0)-2*k
    out+=P({tuple(sorted((key,e) for key,e in w.items() if e)):coef*comb(n//2,k)*(-1)**k})
  return out
 class AP:
  __slots__=('p',)
  def __init__(self,v=0):
   self.p=R.clean(v) if isinstance(v,dict) else ({R.ZERO:P.cv(v)} if v else {})
  @classmethod
  def cv(cls,v):return v if isinstance(v,cls) else cls(v)
  def __bool__(self):return bool(self.p)
  def __add__(self,x):return AP(R.add(self.p,self.cv(x).p))
  __radd__=__add__
  def __neg__(self):return AP(R.scale(-1,self.p))
  def __sub__(self,x):return self+-self.cv(x)
  def __rsub__(self,x):return self.cv(x)+-self
  def __mul__(self,x):return AP(R.mul(self.p,self.cv(x).p))
  __rmul__=__mul__
  def ordinary(self,d):return AP({e:deriv(v,d) for e,v in self.p.items()})
  def angular(self,d):return AP(R.diff(self.p,d))
 I=range(6);J=range(4);sig=(1,-1,1,1,1,1)
 X=[AP(v) for v in R.X];K={k:AP(v) for k,v in R.P.items()};Y=AP(R.Y)
 dY={k:AP(v) for k,v in R.grad(R.Y).items()};hY={k:AP(v) for k,v in R.hess(R.Y).items()}
 def cov(t,idx,d):
  val=t[idx]
  if d==0:return val.ordinary(0)
  if d==1:out=val.ordinary(1)*(1/rho)
  else:
   i=d-2;out=sum((K[i,k]*val.angular(k) for k in J),AP())
   for slot,v in enumerate(idx):
    if v>=2:
     repl=list(idx);repl[slot]=d;out+=X[v-2]*t[tuple(repl)]
   out=out*(1/(rho*c))
  for slot,v in enumerate(idx):
   repl=list(idx)
   if d==1:
    if v in (0,1):repl[slot]=1-v;out-=t[tuple(repl)]*a
   elif v==0:repl[slot]=d;out-=t[tuple(repl)]*a
   elif v==1:repl[slot]=d;out-=t[tuple(repl)]*(th/rho)
   else:
    repl[slot]=0;out+=K[d-2,v-2]*a*t[tuple(repl)]
    repl[slot]=1;out-=K[d-2,v-2]*(th/rho)*t[tuple(repl)]
  return out
 checks=[];angular_checks=[]
 def eq(name,left,right=0):
  left=P.cv(left);right=P.cv(right);r=quotient(left-right)
  need(not r,'NONZERO:'+name+':'+json.dumps(r.encode()[:3],separators=(',',':')))
  checks.append({'name':name,'left_collected_terms':len(left.t),'right_collected_terms':len(right.t),'residual_collected_terms':0})
 def eqAP(name,left,right=0):
  left=AP.cv(left);right=AP.cv(right);r=left-right
  for exp,x in r.p.items():
   reduced=quotient(x);need(not reduced,'NONZERO:'+name+':'+str(exp)+':'+json.dumps(reduced.encode()[:3],separators=(',',':')))
  angular_checks.append({'name':name,'left_angular_terms':len(left.p),'right_angular_terms':len(right.p),'residual_terms':0})
 eq('BACKGROUND_E0_yy',-4*(ay+a*a)-p*p-F(2,3)*pot)
 eq('BACKGROUND_E0_slice',3/rho**2-ay-4*a*a-F(2,3)*pot)
 eq('BACKGROUND_S0',py+4*a*p-Vp)
 # Compatibility is checked with the intended chain rule even in derivative controls.
 if mutation not in ('frozen_rho_radial','wrong_potential_chain'):eq('BACKGROUND_V_chain',deriv(pot,0),Vp*p)
 q={k:AP() for k in product(I,repeat=2)};g={k:AP() for k in q}
 q[0,0]=Y*var('A');q[0,1]=q[1,0]=Y*var('B');q[1,1]=Y*var('D')
 g[0,0]=AP(1);g[1,1]=AP(-1)
 for i in J:
  q[0,i+2]=q[i+2,0]=dY[i,]*var('C');q[1,i+2]=q[i+2,1]=dY[i,]*var('E')
  for j in J:q[i+2,j+2]=Y*var('L')*K[i,j]+hY[i,j]*var('M');g[i+2,j+2]=K[i,j]
 # P is the fully lower connection variation. Its ordinary derivatives are
 # generated by a derivation on free profile symbols, hence coherent through q3.
 dq={(d,i,j):cov(q,(i,j),d) for d,i,j in product(I,repeat=3)}
 conn={(k,i,j):(dq[i,j,k]+dq[j,i,k]-dq[k,i,j])*F(1,2) for k,i,j in product(I,repeat=3)}
 C1={k:sum((sig[i]*conn[k,i,i] for i in I),AP()) for k in I}
 ric={(i,j):sum((sig[k]*(cov(conn,(k,i,j),k)-cov(conn,(k,i,k),j)) for k in I),AP()) for i,j in product(I,repeat=2)}
 eta=Y*var('s');deta={(i,):cov({():eta},(),i) for i in I}
 E={(i,j):ric[i,j]-(deta[j,]*p if i==0 else AP())-(deta[i,]*p if j==0 else AP())-F(2,3)*(eta*Vp*g[i,j]+q[i,j]*pot) for i,j in product(I,repeat=2)}
 dE={(k,i,j):cov(E,(i,j),k) for k,i,j in product(I,repeat=3)}
 bgH={(i,j):g[i,j]*(a*p) for i,j in product(I,repeat=2)};bgH[0,0]=AP(py)
 # Full exact linear identity and independent first-gauge/scalar diagnostics.
 S=sum((sig[i]*cov(deta,(i,),i) for i in I),AP())-sum((sig[i]*sig[j]*q[i,j]*bgH[i,j] for i,j in product(I,repeat=2)),AP())-C1[0]*p-eta*Vpp
 for b in I:
  linear=sum((sig[k]*(dE[k,k,b]-F(1,2)*dE[b,k,k]) for k in I),AP())+(S*p if b==0 else AP())
  eqAP('LINEAR_BIANCHI:'+str(b),linear)
 first={n:var(n) for n in ('A','B','C','D','E','L','M')}
 A,B,C,D,EE,L,M=[first[n] for n in first]
 c1y=(var('Ay')+var('Dy')-3*var('Ly')+8*var('My'))/2-var('Bt')/rho+a*(4*A+D-3*L+8*M)-3*th*B/rho-8*C/(rho*c)
 c10=var('By')-(var('At')+var('Dt')+3*var('Lt')-8*var('Mt'))/(2*rho)+5*a*B-3*th*(D+L-F(8,3)*M)/rho-8*EE/(rho*c)
 c1a=var('Cy')-var('Et')/rho+5*a*C-4*th*EE/rho+((-A+D-L)/2-2*M)/(rho*c)
 eqAP('FIRST_C1_y',C1[0],Y*c1y);eqAP('FIRST_C1_0',C1[1],Y*c10)
 for i in J:eqAP('FIRST_C1_a:'+str(i),C1[i+2],dY[i,]*c1a)
 sdiag=var('syy')+4*a*var('sy')-(var('stt')+3*th*var('st')+8*var('s')/c**2)/rho**2-py*A-a*p*(-D+3*L-8*M)-p*c1y-Vpp*var('s')
 eqAP('FIRST_SCALAR_DIAGNOSTIC',S,Y*sdiag)
 for i,j in product(I,repeat=2):eqAP('E1_SYMMETRY:'+str((i,j)),E[i,j],E[j,i])
 correction={}
 for b in I:
  db=sum((sig[i]*sig[j]*q[i,j]*(-dE[i,j,b]+F(1,2)*dE[b,i,j]) for i,j in product(I,repeat=2)),AP())-sum((sig[k]*C1[k]*E[k,b] for k in I),AP())
  correction[b]=(AP() if mutation=='omit_deltaB_E1' else 2*db)+(AP() if mutation=='omit_S1_deta' else 2*S*deta[b,])
 def project(t,basis,label):
  coeff=[]
  for base in basis:
   numerator=sum((R.integrate(R.mul(t[k].p,v)) for k,v in base.items()),P())
   denominator=R.inner(base,base);coeff.append(numerator/denominator)
  for k,v in t.items():
   reconstructed=sum((AP(base[k])*cf for base,cf in zip(basis,coeff)),AP())
   eqAP(label+':'+str(k),v,reconstructed)
  return coeff
 corr={'y':project({():correction[0]},R.SCALARS,'CORRECTION_y_RECONSTRUCTION'),
       '0':project({():correction[1]},R.SCALARS,'CORRECTION_0_RECONSTRUCTION'),
       'a':project({(i,):correction[i+2] for i in J},R.GRADS,'CORRECTION_a_RECONSTRUCTION')}
 report=json.loads((root/MAP).read_text());need(report['status']=='PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION','MAP_STATUS');vals={n:var(n) for n in report['input_symbols']};vals.update(ay=ay,py=py)
 quantity={'reduced_instead_of_natural':'QEred','forcing_instead_of_natural':'metric_forcing'}.get(mutation,'QE')
 tensor={k:[alg.encoded(rows,vals) for rows in v] for k,v in report['maps'][quantity].items()}
 scalar_key='Qphi' if mutation=='premature_harmonic_scalar' else 'Qphi_general'
 scalar=[alg.encoded(rows,vals) for rows in report['maps'][scalar_key]]
 if mutation=='half_second_variation':tensor={k:[x/2 for x in v] for k,v in tensor.items()}
 B={'y':[],'0':[],'a':[]};H=th/rho
 for j,l in enumerate((0,2,4)):
  K0=tensor['yy'][j];N=tensor['00'][j];L0=tensor['y0'][j];J0=tensor['ab'][j]
  v=P() if j==0 else tensor['ya'][j-1];u=P() if j==0 else tensor['0a'][j-1]
  lam=l*(l+2);tf=F(1,2) if mutation!='missing_trace_half' else F(1)
  # Divergence minus trace-gradient written so trace-factor corruption is exact.
  trace=K0-N+3*J0
  by=deriv(K0,0)-tf*deriv(trace,0)-deriv(L0,1)/rho+a*(4*K0+N-3*J0)-3*H*L0-lam*v/(rho*c)
  bt=deriv(L0,0)-deriv(N,1)/rho-tf*deriv(trace,1)/rho+5*a*L0-3*H*(N+J0)-lam*u/(rho*c)
  if mutation=='positive_time_contraction':by+=2*deriv(L0,1)/rho;bt+=2*deriv(N,1)/rho
  if mutation=='missing_radial_connection':by-=a*(4*K0+N-3*J0);bt-=5*a*L0
  if mutation=='missing_FRW_connection':by+=3*H*L0;bt+=3*H*(N+J0)
  B['y'].append(by);B['0'].append(bt)
  if j:
   beta=2-F(2,3)*lam;h=tensor['ab'][j+2]
   ba=deriv(v,0)-deriv(u,1)/rho+5*a*v-4*H*u+(J0+beta*h-tf*trace)/(rho*c)
   if mutation=='positive_time_contraction':ba+=2*deriv(u,1)/rho
   if mutation=='missing_radial_connection':ba-=5*a*v
   if mutation=='missing_FRW_connection':ba+=4*H*u
   B['a'].append(ba)
 rows=[]
 for field in ('y','0','a'):
  for j,(bval,cval) in enumerate(zip(B[field],corr[field])):
   matter=(p*scalar[j] if field=='y' else P())*(-1 if mutation=='wrong_scalar_source_sign' else 1)
   eq('NATURAL_BIANCHI_'+field+':'+str(j),bval+matter+cval)
   rows.append({'row':field+':'+str(j),'B0_QE_terms':len(bval.t),'p_Qphi_terms':len(matter.t),'first_residual_correction_terms':len(cval.t),'residual_terms':0})
 need(auth(root)==before,'SOURCE_CHANGED_AFTER_USE')
 return {'status':'PASS_EXACT_RESIDUAL_CORRECTED_NATURAL_BULK_BIANCHI','source_authentication':before,'rows':rows,'checks':checks,'angular_coefficient_checks':angular_checks,
  'independent_first_gauge_amplitudes':{'y':c1y.encode(),'0':c10.encode(),'a':c1a.encode()},'independent_first_scalar_amplitude':sdiag.encode(),
  'correction_modal_coefficients':{k:[quotient(x).encode() for x in v] for k,v in corr.items()},
  'identity':'B0(QE)+p Qphi_general e^y+2 deltaB_q(E1)+2 S1 d eta=0',
  'independent_linear_residuals':'Ric1_ab=div P_ab-grad_b tr P_a; E1=Ric1-dphi deta-deta dphi-(2/3)(Vp eta g0+V q); S1=Box eta-q.Hess phi-C1.dphi-Vpp eta',
  'formal_metric_jets':70,'formal_scalar_jets':6,'background_quotient':'rho_y=a rho; a_y=-rho^-2-p²/3; p_y=Vp-4ap; V=6rho^-2+p²/2-6a²; c_t=c theta; theta_t=1-theta²; theta²=1-c^-2',
  'homogeneous_specialization':'Requires E1=0 and nabla E1=0, S1=0. Replacing Qphi_general by Qphi additionally requires C1=0. These are conditional premises, not supplied by this audit.',
  'scope':'Exact finite registered natural source/proof join; no actual master reconstruction membership, new evolution theorem, nonlinear family, accuracy or cosmological result.'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--mutation',choices=CONTROLS);args=ap.parse_args()
 need(not args.output.exists(),'REFUSE_OVERWRITE')
 result=calculate(args.repo_root.resolve(),args.mutation)
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
 print(json.dumps({'status':result['status'],'rows':len(result['rows']),'checks':len(result['checks'])}))
if __name__=='__main__':main()
