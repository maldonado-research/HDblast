#!/usr/bin/env python3
"""Independent H1-to-shell bound from raw retained inputs, without producer imports.

Exact rational consequences only. Continuous operator and PDE premises are
explicit shared conditional assumptions, not re-proved by this program.
"""
import argparse
from fractions import Fraction as F
import hashlib,json
from math import factorial
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
PINS={
 'premises':('research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json','436305e19f375ddead2958dd44ea48dd080302852ef056631efd3ae2e6b433c1'),
 'static':('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json','3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22'),
 'base':('research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json','1a9ebffb6bee9332bc0fd0df0f5d1d9f7f2a796825c2a563e4c644974ba4d8a3'),
 'witness':('research/HDBLAST_INTRINSIC_PREPARATION_ROBUSTNESS_20261104/independent/references/CONTROLS.artifacts/WITNESS_NORMAL.json','ba3973210fb364f40cafaada84d3dd04abb29b385c594d70603903817ad31e5c')}
PRIMARY_KEYS=set('authenticated_sources class continuous_operator_caps_are_retained_conditional_assertions design_sha256 energy initial_L2_caps_per_R limits mean_Gram new_domain observable_difference_caps_per_R retained_domain source_pin_file_sha256 status verification whole_class_absolute_difference_caps whole_class_exact_intervals whole_class_outward_intervals'.split())
ENERGY_KEYS=set('boundary_alpha boundary_growth_exact_upper boundary_scalar_factor boundary_value_square_guard bulk_growth_exact_upper exponential_Taylor_order exponential_Taylor_tail_ratio exponential_argument exponential_guard exponential_rational_upper growth_guard initial_energy_guard initial_energy_per_R_squared tilt_k'.split())
class Rejected(Exception):pass
class Parser(argparse.ArgumentParser):
 def error(self,message):raise Rejected('ARGUMENTS')
def need(x,w):
 if not x:raise Rejected(w)
def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+'\n').encode()
def digest(raw):return hashlib.sha256(raw).hexdigest()
def regular(p):return p.is_file() and not p.is_symlink() and all(not q.is_symlink() for q in p.absolute().parents)
def safe_path(p):return not any(q.is_symlink() for q in [p,*p.absolute().parents])
def interval(n):return tuple(F(x) for x in n['exact'])
def exp_cap(x,n):return sum((x**j/F(factorial(j)) for j in range(n+1)),F())+(x**(n+1)/factorial(n+1))/(1-x/F(n+2))
def main():
 ap=Parser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--candidate',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
 need(a.repo_root.is_absolute() and safe_path(a.repo_root),'REPO_PATH');repo=a.repo_root.resolve(strict=True)
 need(a.output.is_absolute() and safe_path(a.output),'OUTPUT_PATH');out=a.output.resolve();need(not out.exists(),'OUTPUT_EXISTS');need(out.parent.is_dir(),'OUTPUT_PARENT');need(all(base!=out and base not in out.parents for base in [repo,HERE]),'OUTPUT_LOCATION')
 need(a.candidate.is_absolute() and regular(a.candidate),'CANDIDATE_REGULAR')
 docs={};auth={}
 for alias,(name,pin) in PINS.items():
  p=repo/name;need(regular(p),'INPUT_REGULAR:'+alias);raw=p.read_bytes();need(digest(raw)==pin,'INPUT_PIN:'+alias);docs[alias]=json.loads(raw);auth[alias]={'path':name,'sha256':pin,'bytes':len(raw)}
 # All retained bytes are authenticated before numerical reconstruction.
 checks=[]
 def gate(x,n):need(x,'MATH:'+n);checks.append(n)
 prem,st,base,witness=(docs[k] for k in ['premises','static','base','witness'])
 gate(prem['status']=='EXPLICIT_CONDITIONAL_PREMISES_NOT_REPROVED_BY_PORTABLE_CHECKER','explicit_conditional_operator_status')
 gate(prem['mean_gram']==[1,2,1,3,1] and prem['coefficient_norm_upper']=={'A':'1/10000','B':'1/10000','D':'1/1000000','M':'1/2000'},'same_Gram_and_normalized_operator_caps')
 gate(prem['domain']=='0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x','retained_operator_domain')
 ell,T,Delta=F(1,10**6),F(5,4),F(1,4)
 gate(T<=F(3,2) and T==F(5,4) and T<=F(7,4),'backward_cone_is_subset')
 gate(F(st['domain']['t'])==0 and F(st['domain']['ell'])==ell and F(st['domain']['z_min'])<=ell and F(st['domain']['z_max'])>=ell*T and st['domain']['coverage']=='Every point of the original packet annulus; no samples.','raw_continuous_initial_domain')
 gate(st['normalization']=='coefficient j = D_z^j(value)/j! or D_y^j(value)/j!, as named','raw_first_y_derivative_normalization')
 rho,pM,aL,aM,pyM,ayM=F(131),F(1,7000),F(111,1000),F(14,125),F(1,4000),F(1,10000)
 for key,j,lo,hi in [('rho',0,F(129),rho),('p',0,-pM,-F(3,25000)),('a',0,aL,aM),('p',1,-pyM,pyM),('a',1,-ayM,ayM)]:
  box=interval(st['y_jets'][key][j]);gate(lo<=box[0]<=box[1]<=hi,'raw_static_'+key+'_'+str(j))
 gate(F(base['cone']['T']['exact'])==T and F(base['cone']['delta']['exact'])==Delta and base['cone']['domain']=='0<=x<=5/4,0<=r<=5/4-x' and base['cone']['zero_wedge']=='x+r<=1','same_base_target_and_zero_wedge')
 for name,lo,hi in [('rho',F(129),rho),('Wprime',-F(1,1000),-F(1,1200))]:
  box=interval(base['first'][name]);gate(lo<=box[0]<=box[1]<=hi,'raw_shell_'+name)
 # ODE cap, using only D rather than D_r in the forcing.
 b=ell*rho*aM;g=ell*rho*pM*pM/(3*aL);mu=F(20001,20000)
 gate(max(1+3*b,3*ell*ell+b+g)<=mu and exp_cap(mu*Delta,6)<F(4,3),'continuous_ODE_exponential_majorant')
 H=Delta*F(4,3)*F(2,3)*g
 # Deliberately different from the producer: supD<=3R, direct L2 forcing,
 # and actual Gram-square sums instead of squaring a Gram triangle cap.
 gate(F(17,4)<9,'independent_H1_to_C0_constant_3')
 h=F(3,2)*H
 hr=F(3,2)*(3*ell*ell+b+g)*H+F(2,3)*g
 kappa=pM/(2*aL);kr=ell*rho*(pyM/(2*aL)+pM*ayM/(2*aL*aL))
 w=kappa*(h+F(2,3));wr=kr*(h+F(2,3))+kappa*(hr+F(2,3));vel=hr+b*h
 E0=(4*h*h+w*w+4*hr*hr+wr*wr+2*vel*vel)/2
 energy_guard=F(1,5000000);gate(0<E0<energy_guard,'independent_initial_energy_guard')
 # Tilted bulk energy and nonsymmetric Robin estimate.
 k=alpha=F(1,2);A=B=F(1,10000);Dm=F(1,1000000);M=F(1,2000)
 cb=max(5*A+2*B+2*Dm+2,4*B+A+Dm,F(3,2)*Dm+1)
 boundary=(2/alpha)*((alpha+(1-k)*M)**2/k+(1-k)*M*M+alpha*M)
 C=F(201,100);gate(max(cb,boundary)<C,'independent_bulk_and_Robin_growth')
 exponential=exp_cap(C*T,16);gate(C*T==F(201,80) and exponential<13,'independent_order16_exponential_guard')
 scalar_square=(2/alpha)*13*energy_guard;gate(scalar_square<F(1,90000),'boundary_scalar_factor_1_over_300')
 radius=F(2*10**31);Cv=(rho*ell)**2*F(1,1000)*Delta*Delta/F(1200);Ct=rho*ell*F(1,1000)*Delta/F(600)
 baseboxes={'v2':interval(base['proper_v2']),'Theta2':interval(base['proper_Theta2'])};operator={'v2':Cv,'Theta2':Ct}
 actual={key:(box[0]-operator[key]*radius,box[1]+operator[key]*radius) for key,box in baseboxes.items()}
 outward={'v2':(F('-7.2400e16'),F('-8.6400e15')),'Theta2':(F('-1.0113e22'),F('-3.3540e21'))}
 for key,box in actual.items():gate(outward[key][0]<=box[0]<=box[1]<=outward[key][1]<0,'closed_H1_class_'+key)
 # Reuse prior exact witness output; do not repeat its analytic bump evaluation.
 gate(witness['schema']=='HDBLAST_INDEPENDENT_GEOMETRIC_CLASS_WITNESS_EXCLUSION_V1' and F(witness['lambda'])==-10**21,'retained_scaled_witness')
 oldD=interval(witness['direction_D']);sobolev=F(21,10);witness_lower=abs(F(witness['lambda']))*oldD[0]/sobolev
 gate(oldD[0]>F('3.7621e12') and witness_lower>F('1.7914e33')>radius,'old_witness_outside_H1_class_via_endpoint_trace')
 expected={
 'class_radius':str(radius),'class_norm':'sqrt(integral_[1,T](D^2+D_r^2)dr) of the complete preparation increment',
 'mean_Gram':[1,2,1,3,1],'new_domain':'0<=x<=T,0<=r<=T-x,T=5/4','retained_domain':prem['domain'],
 'conditional_operator_assertions':True,'energy':{'tilt_k':str(k),'boundary_alpha':str(alpha),'bulk_growth_exact_upper':str(cb),'boundary_growth_exact_upper':str(boundary),'growth_guard':str(C),'exponential_argument':str(C*T),'exponential_guard':'13','initial_energy_guard':str(energy_guard),'boundary_scalar_factor':'1/300','boundary_value_square_guard':str(scalar_square)},
 'observable_difference_caps_per_R':{key:str(val) for key,val in operator.items()},
 'whole_class_absolute_difference_caps':{key:str(val*radius) for key,val in operator.items()},
 'whole_class_exact_intervals':{key:[str(x) for x in box] for key,box in actual.items()},
 'whole_class_outward_intervals':{key:[str(x) for x in box] for key,box in outward.items()}}
 # Compare all declared shared consequences, retaining distinct initial majorants.
 primary=json.loads(a.candidate.read_bytes());need(set(primary)==PRIMARY_KEYS,'PRIMARY_SCHEMA');need(set(primary['energy'])==ENERGY_KEYS,'PRIMARY_SCHEMA')
 observed={'class_radius':primary['class']['radius_closed'],'class_norm':primary['class']['norm'],'mean_Gram':primary['mean_Gram'],'new_domain':primary['new_domain'],'retained_domain':primary['retained_domain'],'conditional_operator_assertions':primary['continuous_operator_caps_are_retained_conditional_assertions'],'energy':{key:primary['energy'][key] for key in expected['energy']},'observable_difference_caps_per_R':primary['observable_difference_caps_per_R'],'whole_class_absolute_difference_caps':primary['whole_class_absolute_difference_caps'],'whole_class_exact_intervals':primary['whole_class_exact_intervals'],'whole_class_outward_intervals':primary['whole_class_outward_intervals']}
 need(primary['status']=='PASS_CONDITIONAL_H1_PREPARATION_RESPONSE_BOUNDS' and primary['class']['parameter_normalization_fixed'] is True and observed==expected,'PRIMARY_COMPARISON')
 need(0<F(primary['energy']['initial_energy_per_R_squared'])<energy_guard,'PRIMARY_INITIAL_ENERGY_GUARD')
 degree=primary['energy']['exponential_Taylor_order'];need(type(degree)is int and degree>=3 and degree<=100,'PRIMARY_EXPONENTIAL')
 need(primary['energy']['exponential_Taylor_tail_ratio']==str(C*T/F(degree+2)) and primary['energy']['exponential_rational_upper']==str(exp_cap(C*T,degree)) and F(primary['energy']['exponential_rational_upper'])<13,'PRIMARY_EXPONENTIAL')
 gate(True,'shared_primary_projection_matches_independent_math')
 report={'schema':'HDBLAST_INDEPENDENT_H1_PREPARATION_RESPONSE_V1','status':'PASS_INDEPENDENT_H1_PREPARATION_RESPONSE','authenticated_inputs':auth,'shared_projection':expected,'independent_initial_L2_caps_per_R':{key:str(val) for key,val in {'h':h,'h_r':hr,'w':w,'w_r':wr,'V_L':vel}.items()},'independent_initial_energy_per_R_squared':str(E0),'independent_exponential':{'degree':16,'rational_upper':str(exponential)},'old_witness_H1_radius_lower':str(witness_lower),'old_witness_lower_display_guard':'1791400000000000000000000000000000','primary_comparison':{'scope':'Shared mathematical consequences and reported guards; intermediate initial majorants are intentionally different. Root source authentication and canonical byte readback separately protect complete producer reports.','shared_projection_matches':True,'primary_initial_energy_below_shared_guard':True,'primary_exponential_tail_reconstructed':True,'unknown_top_level_fields_rejected':True},'verification':{'finite_checks':len(checks),'check_names':checks},'activity':{'PDE_solves':0,'background_or_evolved_field_queries':0,'bump_queries':0,'analytic_target_bump_replays':0},'limits':['Continuous operator/background/development premises remain conditional. Byte pins do not prove them.','Tilted energy is an analytic estimating device, not an identified positive physical energy.','No physical selection of the class or radius, no uniform nonlinear amplitude, and no Big Bang-origin evidence.','H1 trace theorem bypasses H1-to-C1 embedding failure; it asserts no containment of the two fixed-radius balls.']}
 raw=canonical(report);out.write_bytes(raw);print(json.dumps({'status':'PASS','witness':'INDEPENDENT_H1_RECONSTRUCTION','report_sha256':digest(raw)},sort_keys=True))
if __name__=='__main__':
 try:main()
 except Rejected as e:print(json.dumps({'status':'FAIL','witness':str(e)},sort_keys=True));sys.exit(2)
 except (OSError,ValueError,KeyError,TypeError,ArithmeticError) as e:print(json.dumps({'status':'FAIL','witness':'INVALID_INPUT_STRUCTURE'},sort_keys=True));sys.exit(2)
