#!/usr/bin/env python3
"""Bounded differences for the frozen original or explicitly recentered candidate.

Each primitive carries its actual/nominal jets and a three-channel exact
telescoping difference: shell Taylor coefficients, rational packet/cutoff,
and validated ODE Taylor remainder. Never subtract full profile ranges.
"""
import argparse,hashlib,importlib.util,json,math,sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
OLD=Path('/workspace/hdblast-research-work/total-mean-response-20261023')
BRANCH_PIN='024acc3720696322e9fa2f864de4c75b17c88d81ec1a6b80c020e928cefd9cd0'
def loadmod(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def strict_tree(root,declared):
 for p in root.rglob('*'):
  if p.is_symlink():raise ValueError('unexpected symbolic path '+str(p))
  if p.is_file() and str(p.relative_to(root)) not in declared|{'MANIFEST.json'}:raise ValueError('undeclared regular file '+str(p))
def authenticate(dependency_map=None):
 # The current authenticator is read only after independently checking its
 # enclosing complete 49-file capsule against an external published pin.
 root=Path(dependency_map['T']).resolve() if dependency_map else OLD/'capsule-source'; pin='cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72'
 if sha(root/'MANIFEST.json')!=pin:raise ValueError('current manifest')
 if dependency_map is None and json.loads((OLD/'SOURCE_FREEZE.json').read_text())['manifest']['sha256']!=pin:raise ValueError('external current pin')
 seen=set()
 for e in json.loads((root/'MANIFEST.json').read_text())['files']:
  p=Path(e['path'])
  if p.is_absolute() or '..' in p.parts or str(p) in seen:raise ValueError('manifest path')
  seen.add(str(p)); q=root/p
  if q.is_symlink() or q.stat().st_size!=e['bytes'] or sha(q)!=e['sha256']:raise ValueError('current file '+str(p))
 strict_tree(root,seen)
 audit=loadmod('audit',root/'initial/authenticate_inputs.py'); r=audit.run({k:dependency_map[k] for k in ('P','N','G','C')} if dependency_map else None)
 # The older authenticator permits cache paths. Tighten the imported closure
 # before any numerical import so a stale unhashed pyc can never be consumed.
 for v in r['closures'].values():strict_tree(Path(v['root']),{e['path'] for e in v['entries']})
 r['current']={'manifest_sha256':pin,'files_including_manifest':len(seen)+1,'external_pin':str(OLD/'SOURCE_FREEZE.json') if dependency_map is None else 'embedded independent manifest hash in direct_difference_pilot.py'}
 return r

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--case',choices=['old','conditional','branch'],default='old');ap.add_argument('--beta-center',type=F);ap.add_argument('--beta-radius',type=F,default=F(1,10**11));ap.add_argument('--beta-certificate',type=Path);ap.add_argument('--recenter',action='store_true');ap.add_argument('--dependency-map',type=Path);ap.add_argument('--start',type=int,default=128);ap.add_argument('--stop',type=int,default=224);args=ap.parse_args()
 dependency_map=json.loads(args.dependency_map.read_text()) if args.dependency_map else None
 capsule_root=Path(dependency_map['T']).resolve() if dependency_map else OLD/'capsule-source'
 if args.beta_certificate is None and dependency_map and 'B' in dependency_map:args.beta_certificate=Path(dependency_map['B']).resolve()
 before=authenticate(dependency_map) # before every numerical import
 init=loadmod('initial_numeric',Path(before['closures']['N']['root'])/'numeric/evaluate_initial_profiles.py')
 I,J,ell=init.I,init.J,init.ELL;data,pins=init.load();power=1<<40
 report0=json.loads((capsule_root/'initial/INITIAL_GATE_NORMAL.json').read_text())
 nominal_cs={k:list(map(F,v)) for k,v in report0['representation']['background_coefficients_ascending'].items()}
 def const(v,n=4):return J([v]+[0]*n)
 def jp(cs,z):
  out=const(0,z.n)
  for c in reversed(cs):out=out*z+c
  return out
 def rhs(state):
  rho,e,p,a,alpha=[state[k] for k in ('rho','eta','p','a','alpha')];h=init.poly(init.UP[1],e)/p
  return {'rho':-rho*rho*a,'eta':-rho*p,'p':-rho*(init.poly(init.UP[1],e)-4*a*p),'a':1/rho+rho*p*p/3,'alpha':rho*(F(5,2)*a-h)*alpha}
 def expand(vals,n):
  state={k:J([v]+[0]*n) for k,v in vals.items()}
  for j in range(n):
   rr=rhs(state)
   for k in state:state[k].c[j+1]=rr[k].c[j]/(j+1)
  return state
 def sqrti(x):
  if x.lo<=0:raise ValueError('sqrt positive')
  g=1<<160;l=math.isqrt((x.lo*g*g).__floor__());h=math.isqrt((x.hi*g*g).__ceil__())+1
  return I(F(l,g),F(h,g))
 beta_old=I(*data['shell']['boxes']['eta_b']['exact']);branch_pin=None
 bc=args.beta_center if args.beta_center is not None else nominal_cs['eta'][0]
 if args.case=='branch':
  if args.beta_certificate is None:raise ValueError('branch certificate required')
  branch_pin=sha(args.beta_certificate)
  if branch_pin!=BRANCH_PIN:raise ValueError('independently frozen branch report pin')
  branch_certificate=json.loads(args.beta_certificate.read_text())
  if branch_certificate['status']!='PASS_RETAINED_BRANCH_ANALYTIC_COMPARISON_ENCLOSURE':raise ValueError('branch status')
  beta=I(*branch_certificate['beta'])
 else:beta=beta_old if args.case=='old' else I(bc-args.beta_radius,bc+args.beta_radius)
 if beta.lo<beta_old.lo or beta.hi>beta_old.hi:raise ValueError('conditional beta outside retained theorem')
 model_c=F(data['shell']['parameters_exact']['c']);delta=F(1,1000)
 def shell_at(bb):
  W=F(1,3)+bb**2+bb**3/3;Wp=2*bb+bb**2;ff=1+model_c+model_c*bb
  H2=delta*(W*ff/9-Wp*model_c/12)+delta**2*(ff**2/36-model_c**2/48)
  return {'rho':1/sqrti(H2),'eta':bb,'p':-Wp-delta*model_c/2,'a':(2*W+delta*ff)/6,'alpha':I(1)}
 shell=shell_at(beta)
 if args.recenter:
  if args.case!='branch':raise ValueError('repair requires branch certificate')
  center=(F(branch_certificate['beta'][0])+F(branch_certificate['beta'][1]))/2
  centered=shell_at(I(center));mid=lambda v:(v.lo+v.hi)/2
  nominal={k:const(mid(v)) for k,v in centered.items()}
  for j in range(4):
   rr=rhs(nominal)
   for k in nominal:nominal[k].c[j+1]=I(mid(rr[k].c[j]/(j+1)))
  nominal_cs={k:[mid(x) for x in v.c] for k,v in nominal.items()}
 for k in ('rho','eta','p','a'):
  old=I(*data['shell']['boxes'][k+'_b']['exact'])
  if shell[k].lo<old.lo or shell[k].hi>old.hi:
   # The correlated actual shell is in both exact consequences. Intersect.
   shell[k]=I(max(shell[k].lo,old.lo),min(shell[k].hi,old.hi))
 shell_jets=expand(shell,4)
 # The retained first-exit enclosure covers the whole local interval. It
 # bounds the fifth derivative of each true state component, via the same
 # exact analytic ODE recurrence, independently of the conditional beta.
 coarse=init.background(I(0,F(7,4)),data,'none')
 # Recover rho,eta,p,a,alpha value bounds directly (alpha is C*rho^2).
 coarse_vals={k:v.val() for k,v in zip(('rho','eta','p','a'),coarse[:4])}
 kg=I(*data['coefficients']['alpha_global']['log_derivative']['exact'])
 coarse_vals['alpha']=init.expbox(-kg*I(0,ell*F(7,4)))
 fifth={k:v.c[5] for k,v in expand(coarse_vals,5).items()}
 class P:
  def __init__(self,a,b,d=None):self.a=a;self.b=b;self.d=d if d is not None else [const(0,a.n) for _ in range(3)]
  @property
  def n(self):return self.a.n
  def cv(self,b):return b if isinstance(b,P) else P(const(b,self.n),const(b,self.n))
  def __add__(self,b):b=self.cv(b);return P(self.a+b.a,self.b+b.b,[x+y for x,y in zip(self.d,b.d)])
  __radd__=__add__
  def __neg__(self):return P(-self.a,-self.b,[-x for x in self.d])
  def __sub__(self,b):return self+-self.cv(b)
  def __rsub__(self,b):return self.cv(b)+-self
  def __mul__(self,b):
   b=self.cv(b);return P(self.a*b.a,self.b*b.b,[x*b.b+self.a*y for x,y in zip(self.d,b.d)])
  __rmul__=__mul__
  def inv(self):
   ia,ib=self.a.inv(),self.b.inv();return P(ia,ib,[-ia*x*ib for x in self.d])
  def __truediv__(self,b):return self*self.cv(b).inv()
  def __rtruediv__(self,b):return self.cv(b)*self.inv()
  def __pow__(self,n):
   if n<0:return self.inv()**(-n)
   out=self.cv(1)
   for _ in range(n):out=out*self
   return out
  def dz(self):return P(self.a.dz(),self.b.dz(),[x.dz() for x in self.d])
 def intersect(a,b):return J([I(max(x.lo,y.lo),min(x.hi,y.hi)) for x,y in zip(a.c,b.c)])
 def bg_pair(r):
  z=ell*J([r,1/ell,0,0,0]);out={}
  # Differentiated Taylor theorem, normalized z-jet coefficient n.
  for k in shell:
   nom=jp(nominal_cs[k],z);br=[];rem=[]
   for n in range(5):
    br.append(sum((math.comb(j,n)*(shell_jets[k].c[j]-nominal_cs[k][j])*z.c[0]**(j-n) for j in range(n,5)),I(0)))
    rem.append(math.comb(5,n)*fifth[k]*z.c[0]**(5-n))
   ds=[J(br),const(0),J(rem)];act=nom+ds[0]+ds[2]
   out[k]=P(act,nom,ds)
  rho,e,p,a,alpha=[out[k] for k in ('rho','eta','p','a','alpha')]
  h=init.poly(init.UP[1],e)/p;k=rho*(F(5,2)*a-h)
  return rho,e,p,a,h,alpha/rho**2,3*alpha/(rho**3*p),k,init.poly(init.UP[2],e)
 def square(j):
  for _ in range(40):j=j*j
  return j
 def kernel(t):
  # Same rational R_n(t) as frozen candidate; scalar differences composed
  # by normalized Taylor coefficients. No exp/R whole-range subtraction.
  b=1+t/power;rat=square(1/b);t0=t.val();b0=1+t0/power;R=rat.val()
  dif=[];A=F(1)
  for k in range(5):
   if k:A*=1+F(k-1,power)
   # b^k-A is expanded exactly before interval evaluation.
   coeff=[F(math.comb(k,j),power**j) for j in range(k+1)];coeff[0]-=A
   poly_diff=init.poly(coeff,t0)
   D0=I(-(R*t0**2/(2*power)).hi,0)
   dif.append((-1)**k*(D0+R*poly_diff/(b0**k)))
  dt=J([I(0)]+t.c[1:]);dj=const(0);actual=const(0);power_dt=const(1)
  E=init.expbox(t0)
  for k in range(5):
   dj=dj+power_dt*dif[k]/math.factorial(k)
   actual=actual+power_dt*((-1)**k*E)/math.factorial(k)
   power_dt=power_dt*dt
  return P(actual,rat,[const(0),dj,const(0)])
 def packet_pair(r,tail):
  rr=J([r,1/ell,0,0,0])
  if tail:
   actual=init.bump(r,data,True);return P(actual,const(0),[const(0),actual,const(0)])
  # Centered denominator prevents interval d>1/4 around r=3/2.
  dr=rr-F(3,2);d=F(1,4)-dr*dr
  # Replace constant square range by the exact square operation; convolution
  # remains the exact derivative formula on the entire closed cell.
  d.c[0]=F(1,4)-r.__sub__(F(3,2))**2
  t=1/d-4;q=kernel(t)
  if r.lo<F(66,64):
   tt=64*rr-65;S=126*tt**5-420*tt**6+540*tt**7-315*tt**8+70*tt**9
   cut=P(const(1),S,[const(0),const(1)-S,const(0)])
   q=q*cut
  # The actual packet derivative polynomial oracle is an independent
  # enclosure. Intersection tightens ranges but not difference bounds.
  q.a=intersect(q.a,init.bump(r,data,False));return q
 def registered(bg,E,q):
  rho,e,p,a,h,C,B,k,U2=bg;ir=1/rho;ip=1/p;dy=lambda j:-ir*j.dz()
  u=C*q;v=B*(q.dz()+k*q);vy=dy(v);ut=E*u;d=3*ir*ip;vt=E*(v-d*u)
  A=-2*vy**2-8*p*u*vy+8*p*(h+2*a)*u*v+(16*p*p-192*ir**2)*u*u-(2*U2+F(8,3)*p*p)*v*v-ir**2*(12*ut*ut+2*vt*vt)
  D=ir**2*(36*u*u-2*v*v);f=rho**2*(A/6+F(4,3)*D)/6;fy=dy(f)
  K=E*ir/3*(18*u*u-v*v+d*u*v)
  Pi=E*ir*ip*(d*u*vy+(d*(h-ir)+10*p)*u*v+(30*a-6*ir)*u*u+(2*ir-a)*v*v)
  zero=P(const(0),const(0))
  U=[zero,zero,-u*u/3,f,-F(3,2)*ip*fy]
  V=[ell*(2*rho*K-8*u*ut)/6,ell*rho*((-26*u*dy(u)-2*a*u*u)/6-F(3,2)*fy-3*a*f),ell*(-8*rho*K+44*u*ut)/6,ell*(2*rho*K+4*u*ut)/6,ell*(rho*Pi+2*u*vt)/6]
  return [U,[ell*x.dz() for x in U],V]
 def sqrtup(q):
  g=1<<112;return F(math.isqrt((q*g*g).__ceil__())+1,g)
 def norm(caps):return sqrtup(sum(w*c*c for w,c in zip((1,2,1,3,1),caps)))
 caps=[[[F(0)]*5 for _ in range(3)] for _ in range(3)];rows=[]
 if not 128<=args.start<args.stop<=224:raise ValueError('96-cell pilot cap')
 for j in range(args.start,args.stop):
  r=I(F(j,128),F(j+1,128));tail=r.hi<=F(65,64)
  profiles=registered(bg_pair(r),kernel(ell*J([r,1/ell,0,0,0])),packet_pair(r,tail))
  err=[[[x.d[c].val().mag() for x in row] for row in profiles] for c in range(3)]
  for c in range(3):
   for kind in range(3):caps[c][kind]=[max(a,b) for a,b in zip(caps[c][kind],err[c][kind])]
  rows.append({'index':j,'r':r.row(),'channel_component_error_caps':[[[str(x) for x in row] for row in chan] for chan in err]})
  if (j-args.start)%8==0:print('cell',j,flush=True)
 hs=[[norm(row) for row in chan] for chan in caps];tot=[sum(hs[c][k] for c in range(3)) for k in range(3)];free=lambda h:max((h[0]+F(3,2)*h[2])/3,h[1]+h[2])
 after=authenticate(dependency_map)
 if before!=after:raise ValueError('dependencies changed')
 if branch_pin is not None and sha(args.beta_certificate)!=branch_pin:raise ValueError('branch certificate changed')
 report={'status':'CONDITIONAL_BETA_DESIGN_ONLY' if args.case=='conditional' else 'CERTIFIED_RETAINED_THEOREM_CONDITIONAL_INITIAL_ERROR_BOUND','case':args.case,'beta_interval':beta.row(),'beta_is_newly_validated_root':False,'scope':'unchanged rational five-mean C1 initial/free-wave candidate; no B1 or bulk source work','old_candidate_sha256':sha(capsule_root/'initial/INITIAL_GATE_NORMAL.json'),'candidate_changed':False,'means':['K0','L0','N0','J0','w0'],'gram':[1,2,1,3,1],'ell':str(ell),'channels':['inherited_shell_Taylor_coefficient_difference','rational_packet_cutoff_and_E_approximation','validated_fifth_ODE_derivative_Taylor_remainder_outward_arithmetic'],'channel_H_caps':[[str(x) for x in row] for row in hs],'channel_H_caps_display':[[float(x) for x in row] for row in hs],'channel_free_caps_display':[float(free(row)) for row in hs],'total_H_caps':[str(x) for x in tot],'total_H_caps_display':[float(x) for x in tot],'free_wave_X_cap':str(free(tot)),'free_wave_X_cap_display':float(free(tot)),'passes_1e28':free(tot)<=10**28,'actual_error_lower_bound_proved':False,'whole_cells':len(rows),'partition_indices':[args.start,args.stop],'full_footprint_covered':args.start==128 and args.stop==224,'shell_values':{k:v.row() for k,v in shell.items()},'shell_Taylor_coefficients':{k:[x.row() for x in v.c] for k,v in shell_jets.items()},'uniform_fifth_normalized_ODE_derivative':{k:v.row() for k,v in fifth.items()},'rows':rows,'authentication_summary':{k:{q:v[q] for q in ('manifest_sha256','files_verified','declared_bytes_verified')} for k,v in before['closures'].items()},'current_authentication':before['current'],'reauthenticated_after':True,'arithmetic':'Inherited 112-bit outward rational intervals; all algebraic endpoint operations exact Fraction before outward rounding. Interval rounding is included once in each reported channel; no additional uncharged floating error. Third channel is the independent ODE Taylor remainder evaluated with the same arithmetic, not a separately isolated roundoff-only estimate.','limits':['Conditional beta report does not assert root inclusion.','Paired algebra keeps exact small differences but may still overestimate due to interval product dependence.','Failed upper cap is not an actual error lower bound.','No evolved U or invariant observable has been evaluated.']}
 if args.case=='branch':
  report.update(status='INITIAL_ERROR_BOUND_JOINED_TO_RETAINED_BRANCH_COMPARISON_CERTIFICATE',beta_is_newly_validated_root=False,branch_input_certificate={'path':str(args.beta_certificate),'sha256':branch_pin,'route':branch_certificate['route'],'scope':branch_certificate['scope']})
 if args.recenter:
  report.update(candidate_changed=True,scope='one explicitly registered recentered rational five-mean C1 initial/free-wave candidate; original candidate and report frozen',representation_repair='replace only degree4 rational nominal background coefficients, selected once from one shared certified beta midpoint and exact shell laws; packet/cutoff/E unchanged')
 report['representation_background_coefficients_ascending']={k:[str(x) for x in v] for k,v in nominal_cs.items()}
 report['channel_component_error_caps']=[[[str(x) for x in row] for row in chan] for chan in caps]
 report['representation_background_coefficients_sha256']=hashlib.sha256(json.dumps(report['representation_background_coefficients_ascending'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
 report['beta_shell_dependence']='All shell intervals derive from the same beta interval and exact constraint/junction laws. They are then projected to rectangular intervals for the jet recurrence; this valid enclosure loses beta correlations and is not a full beta Taylor model.'
 report['numeric_cap_below_allocation']=report['passes_1e28']
 report['passes_1e28']=bool(report['full_footprint_covered'] and args.case=='branch' and report['numeric_cap_below_allocation'])
 args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ('status','case','whole_cells','channel_H_caps_display','free_wave_X_cap_display','passes_1e28')}))
if __name__=='__main__':main()
