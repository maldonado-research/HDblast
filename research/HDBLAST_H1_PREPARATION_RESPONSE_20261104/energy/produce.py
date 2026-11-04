#!/usr/bin/env python3
"""Exact rational bound for a conditional tilted five-mean H1 response estimate."""
import argparse,hashlib,json,sys
from fractions import Fraction as F
from pathlib import Path,PurePosixPath
HERE=Path(__file__).resolve().parent
PINS_SHA='ad8796d5fa3a18642819bee3ca0baf7118581d40df6e7e02f917aae24c9e3b7d'
DESIGN_SHA='991bf44eb5b3f92865d072c5893b7347a89fdd9965e9b9617df02f0749f2d158'
class Failure(Exception):pass
def need(x,w):
 if not x:raise Failure(w)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def regular(p):return p.is_file() and not p.is_symlink() and all(not q.is_symlink() for q in p.absolute().parents)
def read(p):return json.loads(p.read_bytes())
def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+'\n').encode()
def scalar(v):return F(v['exact'])
def interval(v):return [F(x) for x in v['exact']]
def encode(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {k:encode(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [encode(x)for x in v]
 return v
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--design',type=Path,default=HERE/'FROZEN_DESIGN.json');a=ap.parse_args()
 repo=a.repo_root.absolute();out=a.output.absolute();design=a.design.absolute()
 need(all(not q.is_symlink()for q in [repo,*repo.parents,out,*out.parents,design,*design.parents]),'PATH_SYMLINK')
 repo=repo.resolve(strict=True);out=out.resolve();design=design.resolve(strict=True)
 need(not out.exists(),'FRESH_OUTPUT_FILE');need(out.parent.is_dir(),'OUTPUT_PARENT_EXISTS');need(all(base!=out and base not in out.parents for base in [repo,HERE]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
 need(regular(HERE/'SOURCE_PINS.json') and sha(HERE/'SOURCE_PINS.json')==PINS_SHA,'SOURCE_PIN_FILE_BYTES');need(regular(design) and sha(design)==DESIGN_SHA,'FROZEN_DESIGN_BYTES')
 pins=read(HERE/'SOURCE_PINS.json');need(pins['schema']=='HDBLAST_H1_ENERGY_SELECTED_SOURCE_PINS_V1','SOURCE_PIN_SCHEMA')
 for name,row in pins['files'].items():
  rel=PurePosixPath(name);need(not rel.is_absolute() and '..'not in rel.parts and rel.as_posix()==name,'SOURCE_RELATIVE_PATH');p=repo/name
  need(regular(p) and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'AUTHENTICATED_SOURCE_BYTES')
 d=read(design);par=d['parameters'];checks=[]
 def check(name,x):need(x,name);checks.append(name)
 T,delta,ell,k,alpha,C,expo,Eguard,mguard,radius=[F(par[n])for n in ['T','Delta','ell','tilt_k','boundary_alpha','growth_guard','exponential_guard','initial_energy_guard','boundary_scalar_factor','closed_H1_radius']]
 check('predeclared_norm_radius_and_target',(T,delta,ell,k,alpha,C,expo,Eguard,mguard,radius)==(F(5,4),F(1,4),F(1,10**6),F(1,2),F(1,2),F(201,100),F(13),F(1,5*10**6),F(1,300),F(2*10**31)))
 premises=read(repo/'research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json')
 check('explicit_conditional_continuous_operator_premise',premises['status']=='EXPLICIT_CONDITIONAL_PREMISES_NOT_REPROVED_BY_PORTABLE_CHECKER' and premises['mean_gram']==[1,2,1,3,1] and premises['coefficient_norm_upper']=={'A':'1/10000','B':'1/10000','D':'1/1000000','M':'1/2000'})
 check('retained_domain_declaration',premises['domain']=='0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x')
 check('new_cone_contained_in_retained_domain',T<=F(3,2) and T<=F(5,4) and T<=F(7,4))
 amax,bmax,dmax,mmax=[F(premises['coefficient_norm_upper'][n])for n in ['A','B','D','M']]
 check('predeclared_operator_caps_match',all(F(par[n])==v for n,v in [('bulk_A_cap',amax),('bulk_B_cap',bmax),('bulk_D_cap',dmax),('Robin_M_cap',mmax)]))
 vcoef=amax+(bmax+k*amax)/2+dmax/2+F(1,2)
 rcoef=k*bmax+(bmax+k*amax)/2+k*dmax/2
 ucoef=dmax/2+k*dmax/2+F(1,2)
 bulkC=max(2*vcoef/(1-k),2*rcoef/(1-k),2*ucoef)
 beta=alpha+(1-k)*mmax;boundaryC=2*(beta*beta/k+(1-k)*mmax*mmax+alpha*mmax)/alpha
 check('positive_tilt_and_boundary_state',0<k<1 and alpha>0)
 check('bulk_energy_growth_guard',bulkC<C);check('nonsymmetric_Robin_growth_guard',boundaryC<C)
 x=C*T;n=24;term=F(1);partial=term
 for j in range(1,n+1):term=term*x/j;partial+=term
 tailfirst=term*x/(n+1);ratio=x/(n+2);expupper=partial+tailfirst/(1-ratio)
 check('exponential_series_tail_ratio',0<ratio<1);check('strict_rational_exponential_guard',expupper<expo)
 prior=read(repo/'research/HDBLAST_INTRINSIC_PREPARATION_ROBUSTNESS_20261104/bounds/references/NORMAL.json');const=prior['constants'];ode=prior['ode_caps'];component=prior['conditioned_component_caps_per_R']
 check('retained_exact_profile_bound_status',prior['status']=='PASS_CONDITIONAL_PHYSICAL_INITIAL_PROFILE_CLASS_ROBUSTNESS')
 rhoM=interval(const['rho'])[1];pM=interval(const['minus_p'])[1];aL=interval(const['a'])[0];aM=interval(const['a'])[1];pyM=scalar(const['abs_p_y']);ayM=scalar(const['abs_a_y']);Fprime=scalar(const['abs_Fprime'])
 check('retained_scales_and_Gram',scalar(const['T'])==T and scalar(const['delta'])==delta and scalar(const['ell'])==ell)
 bM=ell*rhoM*aM;gM=ell*rhoM*pM*pM/(3*aL);H=delta*F(4,3)*F(2,3)*gM;Hr=3*ell*ell*H+(bM+gM)*H+F(2,3)*gM;V=Hr+bM*H;kM=pM/(2*aL);krM=ell*rhoM*(pyM/(2*aL)+pM*ayM/(2*aL*aL))
 check('reconstructed_old_continuous_ODE_caps',H==scalar(ode['H_per_R']) and Hr==scalar(ode['Hr_per_R']) and bM==scalar(ode['bmax']) and gM==scalar(ode['gmax']))
 check('reconstructed_scalar_pullback_and_velocity_caps',kM==scalar(component['k']) and krM==scalar(component['abs_k_r']) and V==scalar(component['V_L']))
 sob=F(21,10);check('H1_to_C0_exact_square_guard',F(17,4)<sob*sob)
 # sqrt(length)=1/2; h and h_r use C0(D), while w_r retains the L2 D_r term.
 factor=sob/2;h2=factor*H;hr2=factor*Hr;vel2=factor*V;w2=kM*(h2+F(2,3));wr2=krM*(h2+F(2,3))+kM*(hr2+F(2,3))
 D0=2*h2+w2;Dr=2*hr2+wr2;Dv=2*vel2;E0=(D0*D0+Dr*Dr+Dv*Dv)/2
 velocity_slots=[F(0),vel2,F(0),F(0),F(0)];radial_slots=[hr2,F(0),F(0),hr2,wr2];check('initial_cross_term_slots_disjoint',sum(weight*v*r for weight,v,r in zip([1,2,1,3,1],velocity_slots,radial_slots))==0)
 check('initial_energy_fixed_guard',E0<Eguard)
 trace_square=2*expo*Eguard/alpha;check('boundary_scalar_fixed_factor_square',trace_square<mguard*mguard)
 base=read(repo/'research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json');bv=interval(base['proper_v2']);bt=interval(base['proper_Theta2'])
 check('same_absolute_base_intervals',bv==interval(prior['base_exact_intervals']['v2']) and bt==interval(prior['base_exact_intervals']['Theta2']))
 check('base_coefficients_strictly_negative',bv[1]<0 and bt[1]<0)
 Cv=(rhoM*ell)**2*Fprime*delta*delta*mguard/4;Ct=rhoM*ell*Fprime*delta*mguard/2;dv=Cv*radius;dt=Ct*radius;ov=[bv[0]-dv,bv[1]+dv];ot=[bt[0]-dt,bt[1]+dt]
 check('whole_H1_class_volume_strict_negative',ov[1]<0);check('whole_H1_class_expansion_strict_negative',ot[1]<0)
 # Fixed outward displays chosen before final report freeze after predeclared radius.
 displayv=[F('-7.2400e16'),F('-8.6400e15')];displayt=[F('-1.0113e22'),F('-3.3540e21')]
 check('whole_H1_class_volume_outward_display',displayv[0]<=ov[0]<=ov[1]<=displayv[1]);check('whole_H1_class_expansion_outward_display',displayt[0]<=ot[0]<=ot[1]<=displayt[1])
 need(all(regular(repo/name) and (repo/name).stat().st_size==row['bytes'] and sha(repo/name)==row['sha256'] for name,row in pins['files'].items()),'AUTHENTICATED_SOURCE_POSTCHECK')
 result=encode({'status':'PASS_CONDITIONAL_H1_PREPARATION_RESPONSE_BOUNDS','class':{'norm':'sqrt(integral_[1,T](D^2+D_r^2)dr) of the complete preparation increment','radius_closed':radius,'parameter_normalization_fixed':True},'design_sha256':DESIGN_SHA,'source_pin_file_sha256':PINS_SHA,'authenticated_sources':pins['files'],'continuous_operator_caps_are_retained_conditional_assertions':True,'retained_domain':premises['domain'],'new_domain':'0<=x<=T,0<=r<=T-x,T=5/4','mean_Gram':[1,2,1,3,1],'energy':{'tilt_k':k,'boundary_alpha':alpha,'bulk_growth_exact_upper':bulkC,'boundary_growth_exact_upper':boundaryC,'growth_guard':C,'exponential_argument':x,'exponential_Taylor_order':n,'exponential_Taylor_tail_ratio':ratio,'exponential_rational_upper':expupper,'exponential_guard':expo,'initial_energy_per_R_squared':E0,'initial_energy_guard':Eguard,'boundary_value_square_guard':trace_square,'boundary_scalar_factor':mguard},'initial_L2_caps_per_R':{'h':h2,'h_r':hr2,'V_L':vel2,'w':w2,'w_r':wr2,'Gram_U0':D0,'Gram_U0_r':Dr,'Gram_U0_x':Dv,'initial_tilt_cross_term':F(0),'initial_boundary_value':F(0)},'observable_difference_caps_per_R':{'v2':Cv,'Theta2':Ct},'whole_class_absolute_difference_caps':{'v2':dv,'Theta2':dt},'whole_class_exact_intervals':{'v2':ov,'Theta2':ot},'whole_class_outward_intervals':{'v2':displayv,'Theta2':displayt},'verification':{'finite_checks':len(checks),'check_names':checks,'new_PDE_solves':0,'new_background_or_evolved_field_queries':0,'new_bump_queries':0},'limits':['Analytic tilted energy is an estimating device, not an identified positive physical energy.','Continuous operator/background/linear uniqueness and exact-family premises remain conditional and are not reproved.','Class radius and preparation are not physically selected.','Coefficient-level control does not give a uniform effective nonlinear parameter window.','No claim of set containment between the former C1 ball and this H1 ball.','No physical force, radiation, SI scale, observational result or higher-dimensional Big Bang origin.']})
 with out.open('xb')as f:f.write(canonical(result))
 print(json.dumps({'status':'PASS','witness':'EXACT_H1_ENERGY_RESPONSE_BOUNDS','finite_checks':len(checks),'output_sha256':sha(out)},sort_keys=True))
if __name__=='__main__':
 try:main()
 except (Failure,ValueError,OSError,KeyError,TypeError) as e:
  print(json.dumps({'status':'FAIL','witness':str(e)},sort_keys=True));sys.exit(2)
