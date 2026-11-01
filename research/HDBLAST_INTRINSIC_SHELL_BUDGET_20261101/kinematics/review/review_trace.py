#!/usr/bin/env python3
"""Independent exact Fraction interval audit; no producer/oracle imports."""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import argparse,hashlib,json
class AuditError(ValueError):pass
def need(ok,msg):
 if not ok:raise AuditError(msg)
class I:
 def __init__(self,a,b=None):
  if isinstance(a,I):self.lo,self.hi=a.lo,a.hi;return
  if b is None and isinstance(a,(list,tuple)):a,b=a
  self.lo,self.hi=F(a),F(a if b is None else b);need(self.lo<=self.hi,'ORDERED_BOX')
 def __add__(self,b):b=I(b);return I(self.lo+b.lo,self.hi+b.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,b):return self+-I(b)
 def __rsub__(self,b):return I(b)+-self
 def __mul__(self,b):
  b=I(b);v=[a*c for a in (self.lo,self.hi) for c in (b.lo,b.hi)];return I(min(v),max(v))
 __rmul__=__mul__
 def __truediv__(self,b):
  b=I(b);need(not b.lo<=0<=b.hi,'DIVISOR_NONZERO');return self*I(1/b.hi,1/b.lo)
 def __rtruediv__(self,b):return I(b)/self
 def sq(self):return I(0 if self.lo<=0<=self.hi else min(self.lo**2,self.hi**2),max(self.lo**2,self.hi**2))
 def mag(self):return max(abs(self.lo),abs(self.hi))
 def meet(self,b):b=I(b);return I(max(self.lo,b.lo),min(self.hi,b.hi))
def encode(x):
 if isinstance(x,I):return [str(x.lo),str(x.hi)]
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [encode(v) for v in x]
 return x
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def projection(r):
 return {'error_input_pins':{k:{'bytes':v['bytes'],'sha256':v['sha256']} for k,v in r['actual_error_input_pins'].items()},'shell':{k:v['exact'] for k,v in r['raw_actual_shell_coefficient_boxes'].items()},
  'history':[{'x':v['scaled_x'],'r':v['scaled_r'],'scalar':v['actual_scalar'],'background':v['actual_background']} for v in r['history_rows']],
  'targets':[{'x':v['scaled_x'],'r':v['scaled_r'],'scalar':v['actual_scalar'],'background':v['actual_background'],'master':v['actual_master']} for v in r['target_rows']]}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def regular(p):return p.is_file() and not p.is_symlink() and all(not q.is_symlink() for q in p.absolute().parents)
def guard(r,pin):
 need(digest(projection(r))==pin,'RAW_BOX_PROJECTION_PIN')
 need(r['actual_oracle_calls']==73 and len(r['history_rows'])==64 and len(r['target_rows'])==9,'FROZEN_64_8_1')
 for key in ['actual_scalar_source_with_complete_inherited_errors','whole_history_from_zero_including_precontact_errors']:
  need(r['scope'][key] is True,'ACTUAL_SOURCE_ERROR_GUARD:'+key)
 need(r['scope']['finite_propagation_zero_before_one_used'] is False,'NO_PRECONTACT_ZERO_SHORTCUT')
 need(r['all_sources_authenticated_before_import_and_after'] is True,'SOURCE_AUTHENTICATION')
 need(r['frozen_campaign']['no_adaptive_refinement'] is True,'NO_ADAPTATION')
def contained(raw,own,label):
 c=I(raw);need(c.lo<=own.lo<=own.hi<=c.hi,'ENCLOSURE:'+label)
def upper(raw,own,label):need(F(raw)>=own,'UPPER:'+label)
def lower(raw,own,label):need(F(raw)<=own,'LOWER:'+label)
def weight(x,a,b,k):
 # Binomial expansion integrated in source coordinate; independently differs from producer antiderivative.
 end=min(x,b)
 if end<=a:return F(0)
 return sum((F(comb(k,j))*(-1)**j*x**(k-j)*(end**(j+1)-a**(j+1))/(j+1) for j in range(k+1)),F(0))
def decimal(i,d=42):
 n=10**d;v=[(i.lo*n).__floor__(),(i.hi*n).__ceil__()]
 return [('' if a>=0 else '-')+str(abs(a)//n)+'.'+str(abs(a)%n).zfill(d) for a in v]
def reconstruct(r):
 ell=F(1,10**6);T=F(63,50);boxes=r['raw_actual_shell_coefficient_boxes']
 p=I(boxes['p_b']['exact']);bc=I(boxes['B_b']['exact']);a=I(boxes['a_b']['exact']);vp=I(boxes['U_prime_b']['exact']);sig2=I(boxes['sigma_second_b']['exact'])
 direct=2*vp-8*a*p+sig2*p;fp=(2*p*bc).meet(direct)
 need(p.hi<0 and bc.lo>0 and I(boxes['mu']['exact']).lo>0 and fp.hi<0,'TRACE_FORCE_NONZERO')
 contained(r['F_prime'],fp,'Fprime');contained(r['F_prime_direct'],direct,'Fprime_direct')
 histories=[];acap=bcap=force_mass=F(0);end=F(0)
 for i,row in enumerate(r['history_rows']):
  lo,hi=map(F,row['scaled_x']);need(lo==end and lo<hi,'HISTORY_COVER');end=hi
  need(row['scaled_r']==['0','0'],'PHYSICAL_SHELL')
  bg=row['actual_background'];rho=I(bg['rho']);f=ell**2*rho.sq()*fp*I(row['actual_scalar']['s'])/6
  aa=3*ell*I(bg['theta']);bb=ell**2*(8/I(bg['c']).sq()-4)
  for key,val in [('trace_forcing',f),('trace_damping_coefficient',aa),('trace_mass_coefficient',bb)]:contained(row[key],val,f'history{i}:{key}')
  acap=max(acap,aa.mag());bcap=max(bcap,bb.mag());force_mass+=(hi-lo)*f.mag();histories.append((lo,hi,f))
 need(end==T,'HISTORY_ENDPOINT')
 den=1-acap*T-bcap*T*T/2;need(den>0,'FEEDBACK_DENOMINATOR')
 V=force_mass/den;U=T*V;G=acap*V+bcap*U
 scalars={'damping_abs_cap':acap,'mass_abs_cap':bcap,'source_abs_integral':force_mass,'u_x_abs_cap':V,'u_abs_cap':U,'feedback_acceleration_error_cap':G}
 for key,val in scalars.items():upper(r[key],val,key)
 lower(r['global_denominator'],den,'denominator');scalars['global_denominator']=den
 target=[];weight_count=0
 for i,row in enumerate(r['target_rows']):
  xl,xh=map(F,row['scaled_x']);bg=row['actual_background'];rho=I(bg['rho']);forces=[]
  for k,termkey in [(0,'forcing_ux_terms'),(1,'forcing_u_terms'),(2,'forcing_integrated_u_terms')]:
   total=I(0);terms=[]
   for lo,hi,f in histories:
    if lo>=xh:break
    wl,wh=weight(xl,lo,hi,k),weight(xh,lo,hi,k);need(0<=wl<=wh,'POSITIVE_MONOTONE_WEIGHT');w=I(wl,wh);term=w*f;total+=term;terms.append((lo,hi,w,term))
   need(len(terms)==len(row[termkey]),'WEIGHT_TERM_COUNT');weight_count+=len(terms)
   for (lo,hi,w,term),stored in zip(terms,row[termkey]):
    need([lo,hi]==list(map(F,stored['source_cell'])),'SOURCE_LABEL');contained(stored['exact_weight'],w,f'weight{i}:{k}');contained(stored['term'],term,f'term{i}:{k}')
   forces.append(total)
  force_ux,force_u,force_int=forces[0],forces[1],forces[2]/2
  ux=force_ux+I(-xh*G,xh*G);u=force_u+I(-xh*xh*G/2,xh*xh*G/2);integral=force_int+I(-xh**3*G/6,xh**3*G/6)
  P=rho*ell*integral;A=rho*I(bg['c']);H=I(bg['theta'])/rho;Hp=1/A.sq();up=ux/(rho*ell);E=u-H*P;Ep=up-H*u-Hp*P;L=P/A.sq();theta=3*Ep+8*L
  vals={'forcing_ux_integral':force_ux,'forcing_u_integral':force_u,'forcing_integrated_u':force_int,'recovered_u':u,'recovered_u_x':ux,'history_integral_u_dx':integral,'recovered_u_t':ux/ell,'P':P,'A':A,'H':H,'H_prime':Hp,'u_prime':up,'E':E,'E_prime':Ep,'L':L,'theta_first_amplitude':theta,'curvature_defect_second':F(10,3)*E.sq()/A.sq(),'mean_squared_shear_second':F(80,9)*L.sq(),'expansion_variance_second':theta.sq()/3,'buchert_Q_second':F(2,9)*(theta.sq()-40*L.sq())}
  for key,val in vals.items():contained(row[key],val,f'target{i}:{key}')
  for key,val in [('ux_feedback_error_cap',xh*G),('u_feedback_error_cap',xh*xh*G/2),('integrated_u_feedback_error_cap',xh**3*G/6)]:upper(row[key],val,f'target{i}:{key}')
  old=I(row['actual_master']['u']);oldux=ell*I(row['actual_master']['u_t'])
  need(not(u.hi<old.lo or u.lo>old.hi),'GLOBAL_MASTER_COMPATIBILITY');need(not(ux.hi<oldux.lo or ux.lo>oldux.hi),'GLOBAL_DERIVATIVE_COMPATIBILITY')
  target.append({'scaled_x':[xl,xh],**vals})
 need(target[0]['scaled_x'][0]==F(31,25) and target[7]['scaled_x'][1]==T,'TARGET_ENDPOINTS')
 need(all(target[i]['scaled_x'][1]==target[i+1]['scaled_x'][0] for i in range(7)),'TARGET_COVER')
 need(target[-1]['scaled_x']==[F(5,4),F(5,4)],'SELECTED_POINT')
 keys=['recovered_u','recovered_u_x','P','E','E_prime','theta_first_amplitude','curvature_defect_second','mean_squared_shear_second','expansion_variance_second','buchert_Q_second']
 hulls={k:I(min(row[k].lo for row in target[:8]),max(row[k].hi for row in target[:8])) for k in keys}
 for key,val in hulls.items():contained(r['continuous_interval_hulls'][key],val,'hull:'+key)
 flags={k:hulls[k].lo>0 for k in ['recovered_u','E_prime','curvature_defect_second','expansion_variance_second','buchert_Q_second','mean_squared_shear_second']}
 need(all(flags[k] for k in ['recovered_u','E_prime','curvature_defect_second','expansion_variance_second','buchert_Q_second']),'STRICT_MAIN_SIGNS')
 return {'raw_box_projection_sha256':digest(projection(r)),'F_prime':fp,'F_prime_direct':direct,'feedback':scalars,'exact_kernel_weights_reviewed':weight_count,'continuous_hulls':hulls,'strict_positive_flags':flags,'targets':target,'continuous_outward_decimal':{k:decimal(v) for k,v in hulls.items()},'selected_target_outward_decimal':{k:decimal(target[-1][k]) for k in keys}}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--raw-report',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 pins=json.loads((Path(__file__).resolve().parent/'SOURCE_PINS.json').read_text());need(regular(args.raw_report),'REGULAR_RAW_REPORT');before=sha(args.raw_report)
 raw=json.loads(args.raw_report.read_text());guard(raw,pins['raw_box_projection_sha256']);science=reconstruct(raw);need(before==sha(args.raw_report),'RAW_POSTCHECK')
 need(not args.output.exists() and not args.output.is_symlink(),'FRESH_OUTPUT');need(all(not p.is_symlink() for p in args.output.absolute().parents),'OUTPUT_ANCESTORS')
 result={'status':'PASS_INDEPENDENT_EXACT_FRACTION_TRACE_AGGREGATION_REVIEW','producer_or_oracle_imports':False,'third_party_imports':False,'arithmetic':'Unrounded exact Fraction intervals reconstructed inside producer outward dyadic enclosures','shared_premise':'Authenticated raw full actual scalar/background boxes with inherited continuous errors from x=0; no new oracle campaign or analytic PDE certificate','science':science}
 args.output.write_text(json.dumps(encode(result),sort_keys=True,indent=2)+'\n');print(json.dumps({'status':result['status'],'report_sha256':sha(args.output),'strict_positive_flags':science['strict_positive_flags']},sort_keys=True))
if __name__=='__main__':main()
