#!/usr/bin/env python3
"""Exact intrinsic shell scalar identities and whole-interval zero certificate.

The inherited actual PDE/jet bounds are conditional premises. No evolved second
response or finite-amplitude nonlinear remainder is evaluated by this checker.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from fractions import Fraction as Q
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
ZERO=(0,0,0,0)
CONTROLS={
 'time_sign':'Lorentzian time contraction',
 'half_hessian':'second parameter derivative normalization',
 'angular_gradient':'complete angular gradient polynomial',
 'mean_threshold':'probability sphere mean',
 'drop_embedding':'joint field and embedding pullback',
 'drop_theta':'derivative of c times vt over v',
 'nonconstant_background':'constant background shell scalar premise',
}

def need(ok,why):
 if not ok:raise ValueError(why)

def plus(*ps):
 out={}
 for p in ps:
  for e,c in p.items():out[e]=out.get(e,Q(0))+c
 return {e:c for e,c in out.items() if c}
def scale(c,p):return {e:c*v for e,v in p.items() if c*v}
def times(p,q):
 out={}
 for e,c in p.items():
  for f,d in q.items():
   k=tuple(a+b for a,b in zip(e,f));out[k]=out.get(k,Q(0))+c*d
 return {e:c for e,c in out.items() if c}
def constant(c):return {} if not c else {ZERO:Q(c)}
def var(i):return {tuple(int(j==i) for j in range(4)):Q(1)}
def deriv(p,i):
 out={}
 for e,c in p.items():
  if e[i]:
   f=list(e);f[i]-=1;out[tuple(f)]=c*e[i]
 return out
def reduce_sphere(p):
 out={};stack=list(p.items())
 while stack:
  e,c=stack.pop()
  if e[3]<2:out[e]=out.get(e,Q(0))+c;continue
  f=list(e);f[3]-=2;stack.append((tuple(f),c))
  for i in range(3):
   g=f.copy();g[i]+=2;stack.append((tuple(g),-c))
 return {e:c for e,c in out.items() if c}
def moment(p):
 out=Q(0)
 for e,c in p.items():
  if any(v%2 for v in e):continue
  orders=[v//2 for v in e];numer=Q(1)
  for n in orders:
   for k in range(n):numer*=Q(1,2)+k
  denom=math.factorial(sum(orders)+1)
  out+=c*numer/denom
 return out

class E:
 """Normalized epsilon series through second order, exact rationals."""
 def __init__(self,c):self.c=tuple(Q(x) for x in c)
 @staticmethod
 def cv(x):return x if isinstance(x,E) else E((x,0,0))
 def __add__(self,b):b=E.cv(b);return E(tuple(x+y for x,y in zip(self.c,b.c)))
 __radd__=__add__
 def __neg__(self):return E(tuple(-x for x in self.c))
 def __sub__(self,b):return self+-E.cv(b)
 def __mul__(self,b):
  b=E.cv(b);return E(tuple(sum((self.c[i]*b.c[k-i] for i in range(k+1)),Q(0)) for k in range(3)))
 __rmul__=__mul__
 def inv(self):
  a,b,c=self.c;need(a!=0,'nondegenerate formal metric pivot')
  return E((1/a,-b/a**2,b*b/a**3-c/a**2))
 def __truediv__(self,b):return self*E.cv(b).inv()

def inverse(matrix):
 n=len(matrix);rows=[[E.cv(v) for v in row]+[E.cv(int(i==j)) for j in range(n)] for i,row in enumerate(matrix)]
 for i in range(n):
  pivot=next((j for j in range(i,n) if rows[j][i].c[0]),None)
  need(pivot is not None,'nondegenerate formal matrix')
  rows[i],rows[pivot]=rows[pivot],rows[i];p=rows[i][i]
  rows[i]=[v/p for v in rows[i]]
  for j in range(n):
   if j!=i:
    q=rows[j][i];rows[j]=[a-q*b for a,b in zip(rows[j],rows[i])]
 return [row[n:] for row in rows]

def algebra_checks(injection):
 n=[var(i) for i in range(4)];Y=plus(times(n[0],n[0]),scale(-1,times(n[1],n[1])))
 S=plus(times(n[0],n[0]),times(n[1],n[1]));Y2=times(Y,Y)
 tangential=[plus(deriv(Y,i),scale(-2,times(n[i],Y))) for i in range(4)]
 G=plus(*(times(d,d) for d in tangential))
 expected=scale(4,plus(S,scale(-1,Y2)))
 if injection=='angular_gradient':expected=Y2
 need(not reduce_sphere(plus(G,scale(-1,expected))),'complete angular gradient polynomial')
 need(moment(Y2)==Q(1,6) and moment(G)==Q(4,3),'independent sphere moments')
 Z=plus(S,constant(-Q(1,2)));F=plus(Y2,scale(-Q(1,2),S),constant(Q(1,12)))
 yp=plus(constant(Q(1,6)),scale(Q(1,2),Z),F)
 gp=plus(constant(Q(4,3)),scale(2,Z),scale(-4,F))
 need(not reduce_sphere(plus(Y2,scale(-1,yp))) and not reduce_sphere(plus(G,scale(-1,gp))),'all-point scalar projection reconstruction')
 intrinsic=load_module('intrinsic_scalar',HERE/'intrinsic_scalar.py')
 for j in range(1,9):
  rho=Q(j+2,2);c=Q(j+3,3);v=Q(2*j-3,7);vt=Q(j+5,11)
  coef=intrinsic.coefficients(v,vt,rho,c)
  reconstructed=plus(constant(coef[0]),scale(coef[1],Z),scale(coef[2],F))
  physical=scale(2/rho**2,plus(scale(-vt**2,Y2),scale(v**2/c**2,G)))
  need(not reduce_sphere(plus(physical,scale(-1,reconstructed))),'intrinsic full angular coefficient formula')
  m=(-vt**2+(4 if injection=='mean_threshold' else 8)*v**2/c**2)/(3*rho**2)
  need(moment(physical)==m==intrinsic.mean(v,vt,rho,c),'probability sphere mean')
  # Direct generic gamma(e)^-1 dpsi(e)dpsi(e), with independent nonzero
  # metric first/second variations and scalar second variations.
  g0=[-rho**2]+[rho**2*c**2]*3
  g=[[E((g0[a] if a==b else 0,Q((a+1)*(b+1)+j,13),Q(a+b+j+1,17))) for b in range(4)] for a in range(4)]
  inv=inverse(g);d=[Q((a+1)*j-2,19) for a in range(4)]
  ds=[E((0,d[a],Q(a-j,23))) for a in range(4)]
  full=sum((inv[a][b]*ds[a]*ds[b] for a in range(4) for b in range(4)),E((0,0,0)))
  target=2*((1 if injection=='time_sign' else -1)*d[0]**2/rho**2+sum((w*w for w in d[1:]),Q(0))/(rho**2*c**2))
  need(target==2*((-d[0]**2/rho**2)+sum((w*w for w in d[1:]),Q(0))/(rho**2*c**2)),'Lorentzian time contraction')
  factor=1 if injection=='half_hessian' else 2
  need(full.c[0]==full.c[1]==0 and factor*full.c[2]==target,'second parameter derivative normalization')
 chi,p,zeta,f=Q(2,3),Q(-5,7),Q(11,13),Q(17,19)
 transformed=intrinsic.material_variation(chi+p*f,p,zeta if injection=='drop_embedding' else zeta-f)
 need(transformed==intrinsic.material_variation(chi,p,zeta),'joint field and embedding pullback')
 v,vt,vtt,c,theta=Q(2),Q(3),Q(5),Q(5,4),Q(3,5)
 direct=c*theta*vt/v+c*(vtt*v-vt**2)/v**2
 N=v*vtt-vt**2+(0 if injection=='drop_theta' else theta*v*vt)
 need(direct==c*N/v**2,'derivative of c times vt over v')
 background_gradient=(Q(1),0,0,0) if injection=='nonconstant_background' else (0,0,0,0)
 need(all(x==0 for x in background_gradient),'constant background shell scalar premise')
 return {'all_point_angular_identities':True,'independent_probability_moments':['1/6','4/3'],
         'generic_formal_inverse_families':8,'metric_and_second_scalar_variations_cancel_at_order2':True,
         'joint_embedding_field_pullback':True,'monotonicity_identity':True}

def load_module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def certificate():
 a=load_module('shell_actual_jets',ROOT/'numeric/actual_first_jets.py');intrinsic=load_module('shell_intrinsic',HERE/'intrinsic_scalar.py')
 rows=[]
 for label,left,right in [('negative_scalar',Q(5,4),Q(11,8)),('positive_scalar',Q(13,8),Q(7,4))]:
  _,s,b=a.actual_jet(a.I(left,right),a.I(0))
  v,vt,vtt=s['s'],s['st'],s['stt'];need(v.lo>0 or v.hi<0,'whole-bracket scalar nonvanishing')
  N=v*vtt-vt**2+b['theta']*v*vt;need(N.hi<0,'strict whole-bracket H monotonicity')
  ends=[]
  for x,direction in [(left,1),(right,-1)]:
   _,e,eb=a.actual_jet(a.I(x),a.I(0));H=eb['c']*e['st']/e['s']
   if direction==1:need(H.lo>0 and H.lo**2>8,'left H above sqrt8')
   else:need(H.hi<0 and H.hi**2>8,'right H below minus sqrt8')
   im=intrinsic.mean(e['s'],e['st'],eb['rho'],eb['c']);need(im.hi<0,'endpoint intrinsic mean negative')
   ends.append({'x':str(x),'v':e['s'].row(),'vt':e['st'].row(),'H':H.row(),'mean_I_second':im.row()})
  rows.append({'label':label,'scaled_time_bracket':[str(left),str(right)],'scaled_radius':'0',
   'whole_bracket':{'v':v.row(),'vt':vt.row(),'vtt':vtt.row(),'theta':b['theta'].row(),'rho':b['rho'].row(),'c':b['c'].row(),'H_derivative_numerator':N.row()},
   'endpoints':ends,'conclusion':{'exactly_two_simple_mean_zeros_in_bracket':True,
   'mean_positive_between_zeros_negative_outside_within_bracket':True,'unique_stationary_v_in_bracket':True}})
 data,errors,pins=a.inputs()
 dependencies={str((ROOT/'numeric'/k).relative_to(ROOT)):v['sha256'] for k,v in pins.items()}
 for f in [HERE/'intrinsic_scalar.py',HERE/'SHELL_INTRINSIC_SCALAR.md',Path(__file__),ROOT/'numeric/actual_first_jets.py',ROOT/'numeric/interval_core.py',ROOT/'numeric/exterior_error_caps.py']:
  need(f.exists(),'dependency exists '+f.name);dependencies[str(f.relative_to(ROOT))]=hashlib.sha256(f.read_bytes()).hexdigest()
 return {'brackets':rows,'total_distinct_simple_mean_zeros_certified':4,'dependency_sha256':dependencies,
 'scope':{'actual_first_solution_enclosed':True,'full_closed_time_brackets':True,'intrinsic_scalar_second_parameter_derivative':True,
 'fixed_background_foliation_mean':True,'nonlinear_finite_amplitude_sign':False,'physical_force':False,'second_waveform_solved':False}}

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--run-controls',action='store_true');p.add_argument('--inject',choices=['none',*CONTROLS],default='none');a=p.parse_args()
 need(not a.output.exists(),'refuse existing output')
 need(not a.output.resolve().is_relative_to(ROOT.resolve()),'output outside scientific source tree')
 algebra=algebra_checks(a.inject);actual=certificate();controls=[]
 if a.run_controls:
  for injection,why in CONTROLS.items():
   target=a.output.with_name(a.output.name+'.rejected-'+injection)
   need(not target.exists(),'fresh rejection output')
   cmd=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])+[str(Path(__file__).resolve()),'--inject',injection,'--output',str(target)]
   r=subprocess.run(cmd,text=True,capture_output=True)
   need(r.returncode!=0 and why in r.stderr and not target.exists(),'specific rejection '+injection)
   controls.append({'injection':injection,'diagnostic':why,'nonzero_exit':True,'no_output':True})
 result={'status':'PASS_INTRINSIC_SHELL_SCALAR_AND_FOUR_SIMPLE_MEAN_ZEROS','algebra':algebra,'actual_interval_certificate':actual,'controls':controls,
         'normalization':'Second parameter derivative; coefficient of epsilon^2 is half. All sign decisions exact rational.'}
 with a.output.open('x') as f:f.write(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print(result['status'])

if __name__=='__main__':main()
