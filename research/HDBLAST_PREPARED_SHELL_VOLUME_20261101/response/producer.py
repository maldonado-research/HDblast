#!/usr/bin/env python3
"""New-family five-mean resolvent and weighted-history enclosure.

The pilot flag distinguishes exploration from the frozen final campaign.
No old-family initial/total-response bounds are used.
"""
import argparse, hashlib, importlib.util, json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
CAP=HERE.parent
SOURCE=None
RAW=None
CONDITIONED=None
ORACLE_CALLS=0
INJECT='none'
class Failure(Exception):pass
ELL=F(1,10**6); T=F(5,4); DELTA=T-1
WEIGHTS=(1,2,1,3,1)
FH=F(17387513156573); GH=F(3851400)
Q=T*(F(2,10000)+DELTA*F(1,10**6))+DELTA*F(1,2000)

class I:
 def __init__(self,a=0,b=None):
  self.lo=F(a); self.hi=self.lo if b is None else F(b)
  need(self.lo<=self.hi,'interval order')
 @staticmethod
 def cv(x):return x if isinstance(x,I) else I(x)
 def __add__(self,x):x=I.cv(x);return I(self.lo+x.lo,self.hi+x.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,x):return self+-I.cv(x)
 def __rsub__(self,x):return I.cv(x)+-self
 def __mul__(self,x):
  x=I.cv(x); y=[a*b for a in (self.lo,self.hi) for b in (x.lo,x.hi)]
  return I(min(y),max(y))
 __rmul__=__mul__
 def __truediv__(self,x):
  x=I.cv(x);need(x.lo>0 or x.hi<0,'zero denominator')
  return self*I(1/x.hi,1/x.lo)
 def __rtruediv__(self,x):return I.cv(x)/self
 def __pow__(self,n):
  if n<0:return 1/self**(-n)
  if n==0:return I(1)
  if n%2==0:return I(0 if self.lo<=0<=self.hi else min(self.lo**n,self.hi**n),max(self.lo**n,self.hi**n))
  return self*self**(n-1)
 def mag(self):return max(abs(self.lo),abs(self.hi))
 def row(self):return [str(self.lo),str(self.hi)]

def sqrtup(x):
 grid=1<<80
 n=(x*grid*grid).__ceil__()
 return F(math.isqrt(n)+1,grid)

def hnorm(xs):return sqrtup(sum(w*x*x for w,x in zip(WEIGHTS,xs)))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def auth():
 result=guard.authenticate(REPO,CAP)
 supplied={}
 for role,path,reference in [('initial',RAW,CAP/'initial/references/RAW_BOXES.json'),('conditioned',CONDITIONED,CAP/'initial/references/CONDITIONED_RAW_BOXES.json')]:
  need(guard.regular(path),'SUPPLIED_REPORT_REGULAR')
  need(path.read_bytes()==reference.read_bytes(),'SUPPLIED_REPORT_CUSTODY')
  supplied[role]={'bytes':path.stat().st_size,'sha256':sha(path)}
 result['supplied_report_bytes']=supplied
 return result

def need(ok,message):
 if not ok:raise Failure(message)

def interval(x):return I(*x['exact'])

def mean_operator_check():
 path=SOURCE/'projection/inputs/linear_operator.py'
 spec=importlib.util.spec_from_file_location('authenticated_mean_operator',path)
 op=importlib.util.module_from_spec(spec);spec.loader.exec_module(op)
 indices=[op.IDX[x+':0']for x in ('K','L','N','J','w')]
 matrices=op.matrices()
 for key,matrix in matrices.items():
  need(all(matrix[i][j]==matrix[j][i]==0 for i in indices for j in range(op.N)if j not in indices),'mean block closure:'+key)
 need([op.gram()[i]for i in indices]==list(WEIGHTS),'exact mean Gram')
 return {'indices':indices,'primitive_matrices_checked':len(matrices),'exact_invariant_block':True}

def adapted(profile,py):
 d={k:interval(v) for k,v in profile.items()}
 rho,p,a,u,ut,uy,qm,f,fy,xi,xiy=[d[k]for k in ('rho','p','a','u','ut','uy','QM','f','fy','xi','xiy')]
 fr=-ELL*rho*fy
 U=[f,I(0),-u*u/3,f,-p*xi]
 Ur=[fr,I(0),2*ELL*rho*u*uy/3,fr,ELL*rho*(py*xi+p*xiy)]
 V=[ELL*(2*rho*qm/(3*a)+F(32,3)*u*ut),
    ELL*rho*(-fy+a*f-F(13,3)*u*uy-a*u*u/3),
    ELL*(-2*rho*qm/(3*a)+F(4,3)*u*ut),
    -F(4,3)*ELL*u*ut,I(0)]
 if INJECT=='drop_conditioned_scalar':U[4]=I(0)
 if INJECT=='reverse_conditioned_scalar':U[4]=-U[4]
 return U,Ur,V

def first_bounds():
 global ORACLE_CALLS
 # Authenticate complete source closure before importing its actual jet oracle.
 sys.path.insert(0,str(SOURCE/'numeric'))
 import actual_first_jets as jets
 ORACLE_CALLS+=1
 m,s,b=jets.actual_jet(jets.I(1,T),jets.I(0),mode='actual')
 m={k:I(v.lo,v.hi)for k,v in m.items()};s={k:I(v.lo,v.hi)for k,v in s.items()}
 rho=I(b['rho'].lo,b['rho'].hi); A=rho*I(1,1+F(1,10**11))
 H=I(0,ELL*T)/rho; Hp=I(0,1)/rho**2
 X=(-m['D']/2).mag();Z=(m['L']/2).mag();d=(m['M']/2).mag();V=m['E'].mag()
 C=rho.hi*ELL*DELTA*X
 E=Z+H.mag()*C
 Ep=(m['Lt']/(2*rho)).mag()+Hp.mag()*C+H.mag()*X
 Dp=(m['Mt']/(2*rho)).mag()+V/A.lo+C/A.lo**2
 psi=s['s'].mag();psip=(s['st']/rho).mag()
 correction=C*psip/3+(E+F(8,3)*d)*psi
 shell=json.loads((SOURCE/'numeric/inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_text())['boxes']
 pick=lambda k:interval(shell[k])
 sigma,sp,spp,Vp,Vpp=[pick(k)for k in ('sigma_b','sigma_prime_b','sigma_second_b','U_prime_b','U_second_b')]
 Wp=2*Vp-sp*spp/2+F(2,3)*sigma*sp
 if INJECT=='reverse_volume_source':Wp=-Wp
 # sigma'''=4 for the retained cubic tension.
 Wpp=2*Vpp-(spp*spp+4*sp)/2+F(2,3)*(sp*sp+sigma*spp)
 qfirst=(Wpp.mag()*psi**2/12+psip**2/6+F(4,3)*psi**2/A.lo**2
       +F(5,3)*E**2/A.lo**2+Ep**2+F(16,3)*Ep*Dp+F(8,3)*Dp**2)
 return dict(rho=rho,A=A,Wprime=Wp,Wsecond=Wpp,psi=psi,psi_s=psip,
             clock=C,E=E,E_s=Ep,D_s=Dp,M2_correction=correction,Q1=qfirst,
             raw_actual_shell_2jet={'x':['1','5/4'],'r':['0','0'],'m':m,'s':s,
              'background':{k:I(v.lo,v.hi) if hasattr(v,'lo') else v for k,v in b.items()}})

def run(pilot=True,freeze=None):
 before=auth();raw=json.loads(RAW.read_text())
 if not pilot:
  need(freeze is not None,'final campaign needs freeze')
  frozen=json.loads(freeze.read_text())
  need(frozen['inputs']==before,'FROZEN_SOURCE_CODE_BINDING')
  need(frozen['T']=='5/4' and frozen['N']==1024 and frozen['delta']=='1/4','frozen target and grid')
  need(frozen['observables']==['M2','V2_over_V0','Theta2'],'frozen observables')
 mean_operator=mean_operator_check()
 bg=json.loads((SOURCE/'numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json').read_text())
 py=interval(bg['y_jets']['p'][1])
 cond=json.loads(CONDITIONED.read_text())
 need(len(raw['profile_cells'])==len(cond['profile_cells']),'conditioned coverage')
 inclusion_checks=0
 for old,new in zip(raw['profile_cells'],cond['profile_cells']):
  need(old['r']==new['r'],'conditioned same cells')
  for name,vector in zip(('U0','U0_r','V0'),adapted(old['profiles'],py)):
   for key,v in zip(('K0','L0','N0','J0','w0'),vector):
    rounded=interval(new['conditioned'][name][key])
    need(rounded.lo<=v.lo<=v.hi<=rounded.hi,'CONDITIONED_ADAPTER_INCLUSION')
    inclusion_checks+=1
 cells=[x for x in raw['profile_cells'] if F(x['r'][1])<=T]
 maxima=[[F(0)]*5 for _ in range(3)];history=I(0);history_abs=F(0);history0=I(0);history0_abs=F(0)
 for cell in cells:
  l,r=map(F,cell['r']); vectors=adapted(cell['profiles'],py)
  for a,v in zip(maxima,vectors):
   for j,z in enumerate(v):a[j]=max(a[j],z.mag())
  weight=(T*(r-l)-(r*r-l*l)/2)
  history+=weight*vectors[0][4];history_abs+=weight*vectors[0][4].mag()
  history0+=(r-l)*vectors[0][4];history0_abs+=(r-l)*vectors[0][4].mag()
 d0,dr,dv=map(hnorm,maxima)
 free=max((d0+T*dv)/3,dr+dv)
 if INJECT=='old_total_cap_substitution':free=F(103*10**31)
 S=T*FH+GH
 if INJECT=='drop_actual_source_history':S=F(0)
 error=(Q*free+S)/(1-Q)
 target=next(v for v in raw['targets']if F(v['r'])==T)
 u0=adapted(target['profiles'],py)[0][4]
 fb=first_bounds();correction=fb['M2_correction']
 if INJECT=='omit_proper_mean_correction':correction=F(0)
 M2=u0+I(-DELTA*error-correction,DELTA*error+correction)
 history_error=DELTA**3*error/6+DELTA**2*correction/2
 mhistory=history+I(-history_error,history_error)
 # Exact de Sitter Green: k(T,x)=sinh(2ell(T-x))/(2ell)*cosh²(ell x)/cosh²(ell T).
 # On this tiny interval, |k/(T-x)-1|<1e-11 by cosh Taylor bounds.
 kernel_error=F(1,10**11)*(history_abs+history_error)
 ghistory=mhistory+I(-kernel_error,kernel_error)
 rho=fb['rho'];vol=(rho*ELL)**2*(fb['Wprime']*ghistory/2+I(-fb['Q1'],fb['Q1'])*DELTA**2/2)
 history0_error=DELTA**2*error/2+DELTA*correction
 mhistory0=history0+I(-history0_error,history0_error)
 derivative_kernel_error=F(1,10**11)*(history0_abs+history0_error)
 ghistory0=mhistory0+I(-derivative_kernel_error,derivative_kernel_error)
 theta=rho*ELL*(fb['Wprime']*ghistory0/2+I(-fb['Q1'],fb['Q1'])*DELTA*(1+F(1,10**11)))
 if not pilot:need(M2.lo>0 and vol.hi<0 and theta.hi<0,'REGISTERED_SIGN_GATES')
 after=auth();need(before==after,'source or code changed during evaluation')
 def enc(x):
  if isinstance(x,I):return {'exact':x.row(),'display_only':[float(x.lo),float(x.hi)]}
  if isinstance(x,F):return {'exact':str(x),'display_only':float(x)}
  if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
  if isinstance(x,list):return [enc(v)for v in x]
  return x
 return enc({'status':'EXPLORATORY_BOUNDED_GAUGE_CONDITIONED_MEAN_PILOT' if pilot else 'PASS_CONDITIONAL_TOTAL_MEAN_SCALAR_VOLUME_EXPANSION_ENCLOSURES',
  'source_authentication':before,'actual_first_oracle_calls':ORACLE_CALLS,'cone':{'T':T,'delta':DELTA,'zero_wedge':'x+r<=1','domain':'0<=x<=5/4,0<=r<=5/4-x'},
  'mean_operator':mean_operator,
  'conditioned_adapter_independent_inclusion_checks':inclusion_checks,
  'q':Q,'source_history_X_cap':S,'initial_component_caps':maxima,'initial_norm_caps':{'U0':d0,'U0_r':dr,'V0':dv},
  'free_X_cap':free,'full_resolvent_minus_free_X_cap':error,'first':fb,
  'free_scalar_target':u0,'proper_M2_target':M2,'linear_weighted_free_history':history,
  'linear_weighted_history_error':history_error,'linear_weighted_M2_history':mhistory,
  'Green_kernel_comparison_bound':kernel_error,'Green_weighted_M2_history':ghistory,'proper_v2':vol,
  'free_unweighted_scalar_history':history0,'unweighted_history_error':history0_error,
  'Green_derivative_weighted_M2_history':ghistory0,'proper_Theta2':theta,
  'M2_strict_positive':M2.lo>0,'v2_strict_positive':vol.lo>0,'v2_strict_negative':vol.hi<0,
  'Theta2_strict_negative':theta.hi<0,
  'limits':['Retained existence, operator, source and gauge-completion premises remain conditional.',
   'No old-family Ucap or initial error is reused.','All actual bulk and shell source histories are charged once via S.',
   'Formal parameter derivatives of the chosen exact preparation; no effective nonzero parameter radius.',
   'Negative second changes do not establish negative total volume or contraction.',
   'No physical radiation or cosmological-origin conclusion.',
   'This is an exploratory pilot.' if pilot else 'Conditional enclosure; external review and nonlinear quantitative continuation are not supplied.']})

def main():
 global REPO,SOURCE,RAW,CONDITIONED,guard,INJECT
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--repo-root',type=Path,required=True)
 p.add_argument('--initial-report',type=Path,default=CAP/'initial/references/RAW_BOXES.json')
 p.add_argument('--conditioned-report',type=Path,default=CAP/'initial/references/CONDITIONED_RAW_BOXES.json')
 p.add_argument('--output-dir',type=Path,required=True)
 p.add_argument('--freeze',type=Path,default=HERE/'CAMPAIGN_FREEZE.json')
 p.add_argument('--inject',choices=('none','drop_conditioned_scalar','reverse_conditioned_scalar','reverse_volume_source','omit_proper_mean_correction','drop_actual_source_history','old_total_cap_substitution'),default='none')
 a=p.parse_args();REPO=a.repo_root.resolve(strict=True);SOURCE=REPO/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
 RAW=a.initial_report.resolve(strict=True);CONDITIONED=a.conditioned_report.resolve(strict=True);INJECT=a.inject
 spec=importlib.util.spec_from_file_location('response_authentication',HERE/'authenticate.py')
 guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
 guard.output_guard(a.output_dir,[REPO,CAP]);auth()
 need(RAW.read_bytes()==(CAP/'initial/references/RAW_BOXES.json').read_bytes(),'INITIAL_REFERENCE_BINDING')
 need(CONDITIONED.read_bytes()==(CAP/'initial/references/CONDITIONED_RAW_BOXES.json').read_bytes(),'CONDITIONED_REFERENCE_BINDING')
 result=run(pilot=False,freeze=a.freeze)
 reference=json.loads((HERE/'FROZEN_SCIENTIFIC_REFERENCE.json').read_text())
 science={k:v for k,v in result.items()if k not in ('source_authentication','actual_first_oracle_calls')}
 need(science==reference,'FROZEN_SCIENTIFIC_REFERENCE')
 result['unchanged_frozen_scientific_values']=True
 a.output_dir.mkdir(parents=True,exist_ok=False)
 (a.output_dir/'REPORT.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'status':'PASS','witness':'NEW_TOTAL_MEAN_ENCLOSURES','actual_first_oracle_calls':ORACLE_CALLS},sort_keys=True))

if __name__=='__main__':
 try:main()
 except (Failure,ValueError,OSError,KeyError) as exc:
  witness=str(exc) if isinstance(exc,Failure) else getattr(exc,'witness','INPUT_OR_AUTHENTICATION_FAILURE')
  print(json.dumps({'status':'FAIL','witness':witness,'actual_first_oracle_calls':ORACLE_CALLS},sort_keys=True));sys.exit(2)
