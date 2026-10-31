#!/usr/bin/env python3
"""Exact reduced source/operator binding on coherent free ordinary jets."""
import argparse, ast, hashlib, importlib.util, json, sys
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
MAP='research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json'
MAP_PIN='34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'
OLD_CAPSULE='research/HDBLAST_NATURAL_FIELD_BINDING_20261030'
OLD_MANIFEST_PIN='5121e4304fad7ab5856977b38c6ebd9fdc205900b543bc1c3d97384d246f378f'
CONTROLS=('premature_QC_background','q_missing_two','q_positive_time','PC_missing_two','PC_positive_time','drop_PC1_role','wrong_gauge_gradient_sign','unnormalized_symgrad','omit_TT_PC1','drop_BPC1','QE_instead_general','QEred_instead_general','forcing_instead_general','premature_Qphi','omit_deltaB_E1','omit_S1_deta','wrong_wave_sign','missing_half_wave','frozen_rho_radial','wrong_potential_chain')
HELPER_PINS={'covector_wave_rows.py':'cb025648deda6eec4cfd77b3bb96c6cf6e0713f2b55987b3fdfff22dbb066bbf','codazzi_laurent.py':'178b0a4e58db63eae6c24112ed9d91bf430f70dca9f6f8126c5e352e429f8fc2','angular_generic.py':'fd14d8cc3ac6f58c206b784c5c57d15430f63f38ca091ac03c7769f9bfc89df9'}
def need(ok,why):
 if not ok:raise ValueError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name):
 path=HERE/(name+'.py');spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def auth(root):
 out={}
 cap=root/OLD_CAPSULE;manifest=cap/'MANIFEST.json'
 need(manifest.is_file() and not manifest.is_symlink() and sha(manifest)==OLD_MANIFEST_PIN,'OLD_MANIFEST_PIN')
 data=json.loads(manifest.read_text());need(len(data['files'])==data['file_count_excluding_this_manifest'],'OLD_MANIFEST_COUNT')
 for row in data['files']:
  rel=Path(row['path']);need(not rel.is_absolute() and '..' not in rel.parts,'OLD_MANIFEST_PATH')
  f=cap/rel;need(f.is_file() and not f.is_symlink() and f.stat().st_size==row['bytes'] and sha(f)==row['sha256'],'OLD_CAPSULE_SOURCE:'+str(rel))
 out[OLD_CAPSULE+'/MANIFEST.json']={'sha256':OLD_MANIFEST_PIN,'file_count':len(data['files'])}

 for name,pin in {MAP:MAP_PIN,**HELPER_PINS}.items():
  p=root/name if name==MAP else HERE/name
  need(p.is_file() and not p.is_symlink(),'SOURCE_TYPE:'+name)
  need(sha(p)==pin,'SOURCE_PIN:'+name);out[name]={'bytes':p.stat().st_size,'sha256':pin}
 for p in (Path(__file__),*(HERE/name for name in HELPER_PINS)):
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

 def project(t,basis,label):
  coeff=[]
  for base in basis:
   numerator=sum((R.integrate(R.mul(t[k].p,v)) for k,v in base.items()),P())
   denominator=R.inner(base,base);coeff.append(numerator/denominator)
  for k,v in t.items():
   reconstructed=sum((AP(base[k])*cf for base,cf in zip(basis,coeff)),AP())
   eqAP(label+':'+str(k),v,reconstructed)
  return coeff
 def project_vector(t,label):
  return {'y':project({():t[0]},R.SCALARS,label+'_y'),'0':project({():t[1]},R.SCALARS,label+'_0'),'a':project({(i,):t[i+2] for i in J},R.GRADS,label+'_a')}
 def project_tensor(t,label):
  return {'yy':project({():t[0,0]},R.SCALARS,label+'_yy'),'y0':project({():t[0,1]},R.SCALARS,label+'_y0'),'00':project({():t[1,1]},R.SCALARS,label+'_00'),'ya':project({(i,):t[0,i+2] for i in J},R.GRADS,label+'_ya'),'0a':project({(i,):t[1,i+2] for i in J},R.GRADS,label+'_0a'),'ab':project({(i,j):t[i+2,j+2] for i,j in product(J,repeat=2)},R.TENSORS,label+'_ab')}
 def unproject_vector(t):
  return {0:sum((AP(base[()])*v for base,v in zip(R.SCALARS,t['y'])),AP()),1:sum((AP(base[()])*v for base,v in zip(R.SCALARS,t['0'])),AP()),**{i+2:sum((AP(base[i,])*v for base,v in zip(R.GRADS,t['a'])),AP()) for i in J}}
 def unproject_tensor(t):
  out={k:AP() for k in product(I,repeat=2)}
  for key,ij in [('yy',(0,0)),('y0',(0,1)),('00',(1,1))]:out[ij]=sum((AP(base[()])*v for base,v in zip(R.SCALARS,t[key])),AP())
  out[1,0]=out[0,1]
  for key,k in [('ya',0),('0a',1)]:
   for i in J:out[k,i+2]=out[i+2,k]=sum((AP(base[i,])*v for base,v in zip(R.GRADS,t[key])),AP())
  for i,j in product(J,repeat=2):out[i+2,j+2]=sum((AP(base[i,j])*v for base,v in zip(R.TENSORS,t['ab'])),AP())
  return out
 def modes_eq(label,left,right):
  for field in left:
   need(len(left[field])==len(right[field]),'MODE_COUNT:'+label+':'+field)
   for j,(x,y) in enumerate(zip(left[field],right[field])):eq(label+':'+field+':'+str(j),x,y)
 def modes_add(*ts):return {field:[sum((t[field][j] for t in ts),P()) for j in range(len(ts[0][field]))] for field in ts[0]}
 def modes_scale(s,t):return {field:[s*x for x in values] for field,values in t.items()}
 def B0(tensor):
  B={'y':[],'0':[],'a':[]};H=th/rho
  for j,l in enumerate((0,2,4)):
   K0=tensor['yy'][j];N=tensor['00'][j];L0=tensor['y0'][j];J0=tensor['ab'][j]
   v=P() if j==0 else tensor['ya'][j-1];u=P() if j==0 else tensor['0a'][j-1]
   lam=l*(l+2);trace=K0-N+3*J0
   B['y'].append(deriv(K0,0)-deriv(trace,0)/2-deriv(L0,1)/rho+a*(4*K0+N-3*J0)-3*H*L0-lam*v/(rho*c))
   B['0'].append(deriv(L0,0)-deriv(N,1)/rho-deriv(trace,1)/(2*rho)+5*a*L0-3*H*(N+J0)-lam*u/(rho*c))
   if j:
    beta=2-F(2,3)*lam;h=tensor['ab'][j+2]
    B['a'].append(deriv(v,0)-deriv(u,1)/rho+5*a*v-4*H*u+(J0+beta*h-trace/2)/(rho*c))
  return B
 raw=json.loads((root/MAP).read_text());need(raw['status']=='PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION','MAP_STATUS')
 vals={n:var(n) for n in raw['input_symbols']};vals.update(ay=ay,py=py)
 maps={}
 for name in ('QC','QC_actual','QC_background','PC1','gauge_gradient','QE','QEred','QEred_general','metric_forcing'):
  maps[name]={k:[alg.encoded(rows,vals) for rows in v] for k,v in raw['maps'][name].items()}
 modes_eq('QC_ALIAS',maps['QC'],maps['QC_actual'])
 q_contract={b:sum(((1 if mutation=='q_positive_time' else sig[d])*q[b,d]*C1[d] for d in I),AP())*(1 if mutation=='q_missing_two' else 2) for b in I}
 q_modes=project_vector(q_contract,'qC1_RECONSTRUCTION')
 qc=maps['QC_background'] if mutation=='premature_QC_background' else maps['QC_actual']
 modes_eq('ACTUAL_LOWERING',qc,modes_add(maps['QC_background'],q_modes))
 QC=unproject_vector(maps['QC_actual']);QCT={(i,):x for i,x in QC.items()}
 symgrad={(i,j):(cov(QCT,(j,),i)+cov(QCT,(i,),j))*(1 if mutation=='unnormalized_symgrad' else F(1,2)) for i,j in product(I,repeat=2)}
 gradient_modes=project_tensor(symgrad,'SYMNABLA_QC_RECONSTRUCTION')
 modes_eq('GAUGE_GRADIENT',maps['gauge_gradient'],gradient_modes)
 pc={(i,j):sum(((1 if mutation=='PC_positive_time' else sig[d])*conn[d,i,j]*C1[d] for d in I),AP())*(1 if mutation=='PC_missing_two' else 2) for i,j in product(I,repeat=2)}
 pc_modes=project_tensor(pc,'PC1_RECONSTRUCTION')
 if mutation=='omit_TT_PC1':pc_modes['ab'][5]=P();pc_modes['ab'][6]=P()
 modes_eq('PC1_CONTRACTION',maps['PC1'],pc_modes)
 rhs=modes_add(maps['QE'],modes_scale(1 if mutation=='wrong_gauge_gradient_sign' else -1,maps['gauge_gradient']),modes_scale(0 if mutation=='drop_PC1_role' else 1,maps['PC1']))
 modes_eq('REDUCED_GENERAL_ROLE',maps['QEred_general'],rhs)
 modes_eq('REDUCED_ROLE',maps['QEred'],modes_add(maps['QE'],modes_scale(-1,maps['gauge_gradient'])))
 # Reconstruct all full tensor components, including the two TT output modes.
 for name,rhs in [('REDUCED_GENERAL_ROLE',modes_add(maps['QE'],modes_scale(-1,maps['gauge_gradient']),maps['PC1'])),('REDUCED_ROLE',modes_add(maps['QE'],modes_scale(-1,maps['gauge_gradient'])))]:
  left=unproject_tensor(maps['QEred_general' if name=='REDUCED_GENERAL_ROLE' else 'QEred']);right=unproject_tensor(rhs)
  for k in left:eqAP(name+'_FULL:'+str(k),left[k],right[k])
 ric={(i,j):sum((sig[k]*(cov(conn,(k,i,j),k)-cov(conn,(k,i,k),j)) for k in I),AP()) for i,j in product(I,repeat=2)}
 eta=Y*var('s');deta={(i,):cov({():eta},(),i) for i in I}
 scalar_gap=q_contract[0]*p-2*sum((sig[k]*C1[k]*deta[k,] for k in I),AP())
 scalar_gap_modes=project({():scalar_gap},R.SCALARS,'SCALAR_HARMONIC_DIFFERENCE_RECONSTRUCTION')
 for j,gap in enumerate(scalar_gap_modes):eq('SCALAR_HARMONIC_ROLE:'+str(j),alg.encoded(raw['maps']['Qphi_general'][j],vals)-alg.encoded(raw['maps']['Qphi'][j],vals),gap)
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

 corr=project_vector(correction,'FIRST_RESIDUAL_CORRECTION_RECONSTRUCTION')
 wave=load('covector_wave_rows')
 LQC=wave.wave_rows(maps['QC_actual'],deriv,rho,a,p,c,th)
 Bgradient=B0(maps['gauge_gradient'])
 modes_eq('COMPOSED_COVECTOR_OPERATOR',Bgradient,modes_scale(F(1,2),LQC))
 reconstructed_wave=unproject_vector(LQC)
 for b in I:
  direct_Bgradient=sum((sig[k]*(cov(symgrad,(k,b),k)-F(1,2)*cov(symgrad,(k,k),b)) for k in I),AP())
  eqAP('COMPOSED_COVECTOR_OPERATOR_FULL:'+str(b),direct_Bgradient,reconstructed_wave[b]*F(1,2))
 # This is a collected map substitution in the independently derived wave.
 # Neither LQC nor E1/S1 is defined by rearranging the reduced target identity.
 tensor_name={'QE_instead_general':'QE','QEred_instead_general':'QEred','forcing_instead_general':'metric_forcing'}.get(mutation,'QEred_general')
 Bred=B0(maps[tensor_name]);Bpc=B0(maps['PC1']);Bnatural=B0(maps['QE'])
 scalar_name='Qphi' if mutation=='premature_Qphi' else 'Qphi_general'
 scalar=[alg.encoded(rows,vals) for rows in raw['maps'][scalar_name]]
 rows=[]
 for field in ('y','0','a'):
  for j,bval in enumerate(Bred[field]):
   matter=p*scalar[j] if field=='y' else P()
   natural=Bnatural[field][j]+matter+corr[field][j]
   if mutation not in ('premature_Qphi','omit_deltaB_E1','omit_S1_deta'):eq('REPLAY_NATURAL_BIANCHI:'+field+':'+str(j),natural)
   w=LQC[field][j]*(-F(1,2) if mutation=='wrong_wave_sign' else 1 if mutation=='missing_half_wave' else F(1,2))
   pc=P() if mutation=='drop_BPC1' else Bpc[field][j]
   total=bval+matter+corr[field][j]+w-pc
   eq('REDUCED_GENERAL_BIANCHI:'+field+':'+str(j),total)
   rows.append({'row':field+':'+str(j),'B0_QEred_general_terms':len(bval.t),'p_Qphi_general_terms':len(matter.t),'first_residual_correction_terms':len(corr[field][j].t),'half_L0_QC_actual_terms':len(w.t),'B0_PC1_terms':len(pc.t),'residual_terms':0})
 need(auth(root)==before,'SOURCE_CHANGED_AFTER_USE')
 return {'status':'PASS_EXACT_REDUCED_SOURCE_BINDING','source_authentication':before,'rows':rows,'checks':checks,'angular_coefficient_checks':angular_checks,
  'role_counts':{'QC_alias':8,'actual_lowering':8,'gauge_gradient':20,'PC1':20,'QEred_general':20,'QEred':20,'Qphi_general_minus_Qphi':3,'composed_covector_operator':8},
  'independent_first_gauge_amplitudes':{'y':c1y.encode(),'0':c10.encode(),'a':c1a.encode()},
  'wave_composition':'Authenticated independent direct rough-wave contraction; L0=Box0+(2/3)V+p^2 radial projector; QC_actual and coherent derivatives substituted after arbitrary eight-mode operator proof.',
  'identity':'B0(QEred_general)+p Qphi_general e^y+2 deltaB_q(E1)+2 S1 d eta+(1/2)L0 QC_actual-B0(PC1)=0',
  'independent_linear_residuals':'Ric1_ab=div P_ab-grad_b tr P_a; E1=Ric1-dphi deta-deta dphi-(2/3)(Vp eta g0+V q); S1=Box eta-q.Hess phi-C1.dphi-Vpp eta',
  'formal_metric_jets':70,'formal_scalar_jets':6,'background_quotient':'rho_y=a rho; a_y=-rho^-2-p²/3; p_y=Vp-4ap; V=6rho^-2+p²/2-6a²; c_t=c theta; theta_t=1-theta²; theta²=1-c^-2',
  'homogeneous_specialization':{'identity':'B0(QEred)=-p Qphi e^y-(1/2)L0 QC_actual','first_natural_premises':['E1=0','nabla E1=0','S1=0'],'harmonic_premises':['C1=0','nabla C1=0'],'notes':'First-natural master PDE as a neighbourhood identity and its prolongation supply first-natural premises; harmonic-neighbourhood identity supplies C1 prolongations. QC_actual is retained; exchanging to QC_background inside L0 additionally requires nabla^2 C1=0.'},
  'scope':'Exact finite registered reduced source/operator/proof join. No actual numerical q3 adapter, subsidiary uniqueness proof, new evolution theorem, nonlinear family, accuracy or cosmological result.'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--mutation',choices=CONTROLS);args=ap.parse_args()
 need(not args.output.exists(),'REFUSE_OVERWRITE')
 result=calculate(args.repo_root.resolve(),args.mutation)
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
 print(json.dumps({'status':result['status'],'rows':len(result['rows']),'checks':len(result['checks'])}))
if __name__=='__main__':main()
