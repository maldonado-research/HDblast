#!/usr/bin/env python3
"""Independent direct-ODE audit of inner initial-curvature smoothing diagnostics."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode=True
PINS={
 'producer':'2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1',
 'initial':'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad',
 'bump':'add927ba57117a967da7d6725821d95e894ef954405d503731b3055fc64377bd',
}
JOIN_PIN='e1d0d6e3d8650cbf402ccb4f5a87fd20c20db651e3bd323c68f0bdea0219dd35'
def need(ok,message):
 if not ok: raise RuntimeError(message)
def pinned(path,digest):
 raw=path.read_bytes();need(hashlib.sha256(raw).hexdigest()==digest,'INPUT_HASH:'+str(path));return raw
def encode(z):
 if isinstance(z,Q):return str(z)
 if isinstance(z,dict):return {k:encode(v) for k,v in z.items()}
 if isinstance(z,(list,tuple)):return [encode(v) for v in z]
 return z
def pder(p):return [k*p[k] for k in range(1,len(p))] or [Q(0)]
def evalp(p,x):return sum((v*x**i for i,v in enumerate(p)),Q(0))
def ptrim(p):
 p=list(p)
 while len(p)>1 and p[-1]==0:p.pop()
 return p
def padd(p,q):return ptrim([(p[i] if i<len(p) else 0)+(q[i] if i<len(q) else 0) for i in range(max(len(p),len(q)))])
def pmul(p,q):
 out=[Q(0)]*(len(p)+len(q)-1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):out[i+j]+=a*b
 return ptrim(out)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path)
 ap.add_argument('--producer-report',required=True,type=Path)
 args=ap.parse_args();root=args.repo_root.resolve()
 producer_join_raw=pinned(args.producer_report,JOIN_PIN);producer_join=json.loads(producer_join_raw)
 paths={'producer':root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/full23_initial.py',
 'initial':root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results/INITIAL_REPORT.json',
 'bump':root/'research/HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016/numeric/inputs/BUMP_DERIVATIVE_BOUNDS.json'}
 raw={k:pinned(p,PINS[k]) for k,p in paths.items()};report=json.loads(raw['initial']);bump=json.loads(raw['bump'])
 spec=importlib.util.spec_from_file_location('independent_join_registered_initial',paths['producer']);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 c=m.Context(root);I,J,ell=c.I,c.J,c.ell;delta=Q(1,2**30);order=5
 groups=[];controls=[]
 def ck(name,ok):need(ok,name);groups.append(name)
 ck('continuous_initial_report',report['status']=='PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER')
 ck('coordinate_scaling',ell==Q(1,10**6))
 ck('exact_positive_Gram',len(c.weights)==23 and min(c.weights)>0)
 polys=[list(map(Q,row)) for row in bump['polynomial_coefficients_ascending']]
 ck('actual_fifth_packet_jet_available',bump['maximum_derivative_order']>=5 and len(polys)>=6)
 ck('packet_polynomial_base',polys[0]==[1])
 for n in range(5):
  diff=padd(polys[n],[-z for z in pder(polys[n])])
  want=pmul([0,0,1] if n%2==0 else [0,-4,1],diff)
  if n%2:want=padd(want,[-2*z for z in polys[n]])
  ck('actual_packet_derivative_recurrence_'+str(n+1),ptrim(polys[n+1])==want)
 test=J([I(0),I(1/ell),I(0),I(0)])
 test2,test3=test**2,test**3
 ck('ordinary_scaled_second_derivative_factorial',2*ell**2*test2.c[2].lo==2*ell**2*test2.c[2].hi==2)
 ck('ordinary_scaled_third_derivative_factorial',6*ell**3*test3.c[3].lo==6*ell**3*test3.c[3].hi==6)
 need(ell**3*test3.c[3].lo!=6,'omitted_third_factorial');controls.append('omitted_third_factorial')
 # Polynomial exact integral of normalized even kernel on positive half.
 kp=[Q(35,32)*v for v in [1,0,-3,0,3,0,-1]]
 moment=lambda j:sum((v/Q(n+j+1) for n,v in enumerate(kp)),Q(0))
 ck('kernel_positive_half_mass',moment(0)==Q(1,2));ck('kernel_positive_half_first',moment(1)==Q(35,256))
 ck('kernel_vanishes_at_boundary',evalp(kp,Q(1))==0)
 for wrong,name in ((Q(1),'unit_instead_of_half_mass'),(Q(35,128),'full_instead_of_half_first')):
  need(wrong!=(moment(0) if wrong==1 else moment(1)),name);controls.append(name)
 S=list(map(Q,[0,0,0,0,0,126,-420,540,-315,70]));D=S
 for j in range(5):
  ck('packet_cutoff_C4_match_'+str(j),evalp(D,Q(0))==0 and evalp(D,Q(1))==(1 if j==0 else 0));D=pder(D)
 ck('packet_fifth_join_jump',evalp(D,Q(0))==evalp(D,Q(1))==15120)
 need(evalp(D,Q(1))!=0,'wrong_C5_cutoff');controls.append('assumed_C5_existing_cutoff')
 # Nominal expression evaluation with forced exact branch; no numerical side stencil.
 def nominal(r,branch):
  if branch=='zero':return [c.constant(0,6)]*23
  bg,rr,z=c.nominal_background(r,6);d=rr-Q(3,2);den=Q(1,4)-d*d;den.c[0]=Q(1,4)-(r-Q(3,2))**2
  packet=c.square40(1/(1+(1/den-4)/c.power))
  if branch=='cutoff':
   t=64*rr-65;packet*=126*t**5-420*t**6+540*t**7-315*t**8+70*t**9
  elif branch!='plain':raise RuntimeError('unknown_branch')
  E=c.square40(1/(1+z/c.power));return c._registered(bg,E,packet)[0]
 # Independent actual jet by direct autonomous ODE recurrence, not paired Taylor differences.
 def actual(rc):
  vals=c.actual_background(I(rc));state={k:J([v]+[0]*order) for k,v in vals.items()}
  for n in range(order):
   rho,eta,p,a,alpha=[state[k] for k in ('rho','eta','p','a','alpha')];V1=c.init.poly(c.init.UP[1],eta)
   rhs={'rho':-rho*rho*a,'eta':-rho*p,'p':-rho*(V1-4*a*p),'a':1/rho+rho*p*p/3,'alpha':rho*(Q(5,2)*a-V1/p)*alpha}
   for k in state:state[k].c[n+1]=rhs[k].c[n]/(n+1)
  s=1/(Q(1,4)-(rc-Q(3,2))**2);v=3-2*rc;ep=c.init.expbox(I(s-4))
  packet=J([c.init.poly(polys[n],I(s))*ep*(v if n%2 else 1)/(ell**n*math.factorial(n)) for n in range(order+1)])
  e=c.init.expbox(I(ell*rc));E=J([(-1)**n*e/math.factorial(n) for n in range(order+1)])
  out=c._registered(state,E,packet)[0];ck('actual_order5_output_order_'+str(rc),all(z.n>=2 for z in out))
  return [2*ell**2*z.c[2] for z in out],vals
 joins=[]
 for rc,left,right in ((Q(65,64),'zero','cutoff'),(Q(66,64),'cutoff','plain')):
  ck('isolated_join_'+str(rc),2*delta<Q(1,64) and 1<rc-delta<rc+delta<Q(7,4))
  au,state=actual(rc);nl,nr=nominal(I(rc),left),nominal(I(rc),right)
  sl,sr=nominal(I(rc-delta,rc),left),nominal(I(rc,rc+delta),right)
  rows=[];lower2=Q(0);nonzero=[]
  for i,label in enumerate(c.op.LABELS):
   L=2*ell**2*nl[i].c[2];R=2*ell**2*nr[i].c[2]
   Mminus=(6*ell**3*sl[i].c[3]).mag();Mplus=(6*ell**3*sr[i].c[3]).mag();rem=moment(1)*delta*(Mminus+Mplus)
   average=(L+R)/2;smooth=average+I(-rem,rem);err=au[i]-smooth
   lo=Q(0) if err.lo<=0<=err.hi else min(abs(err.lo),abs(err.hi))
   lower2+=c.weights[i]*lo**2
   if lo:nonzero.append(label)
   rows.append({'label':label,'actual_u_rr_direct_ODE':au[i].row(),'nominal_left_u_rr':L.row(),'nominal_right_u_rr':R.row(),'nominal_u_rr_jump':(R-L).row(),
    'left_local_M3':Mminus,'right_local_M3':Mplus,'half_kernel_remainder_abs_upper':rem,'smoothed_nominal_u_rr_enclosure':smooth.row(),'actual_minus_smoothed_u_rr_direct_enclosure':err.row(),'component_abs_lower':lo})
   ck('valid_error_interval_'+str(rc)+'_'+label,err.lo<=err.hi and rem>=0)
  joins.append({'r':rc,'left_branch':left,'right_branch':right,'rows':rows,'actual_error_sign_excluding_labels':nonzero,'H_error_squared_lower':lower2,'actual_background_values':{k:v.row() for k,v in state.items()}})
 ck('first_join_nominal_curvature_identically_zero',all(all(Q(v)==0 for v in row[key]) for row in joins[0]['rows'] for key in ('nominal_left_u_rr','nominal_right_u_rr')))
 scalar=[z for z in joins[1]['rows'] if z['label'] in ('w:0','w:2','w:4')]
 ck('second_join_nominal_scalar_curvature_jump',all(Q(z['nominal_u_rr_jump'][0])*Q(z['nominal_u_rr_jump'][1])>0 for z in scalar))
 ck('producer_report_successful',producer_join['status']=='PASS_INNER_JOIN_ACTUAL_INITIAL_CURVATURE_ERROR')
 ck('same_fixed_smoothing_width',Q(producer_join['delta_u'])==delta)
 ck('same_unmodified_target_and_candidate',producer_join['actual_target_unchanged'] and producer_join['initial_candidate_unchanged'] and producer_join['original_bulk_and_shell_sources_unchanged'])
 ck('same_complete_Gram',list(map(Q,producer_join['positive_Gram']))==c.weights)
 ck('same_join_count',len(producer_join['joins'])==len(joins)==2)
 for produced,independent in zip(producer_join['joins'],joins):
  rc=independent['r'];ck('same_join_and_branches_'+str(rc),Q(produced['r'])==rc and all(produced[k]==independent[k] for k in ('left_branch','right_branch')))
  ck('same_sign_excluding_labels_'+str(rc),produced['fixed_width_nonzero_component_labels']==independent['actual_error_sign_excluding_labels'])
  producer_lower2=Q(0)
  for i,(p,z) in enumerate(zip(produced['rows'],independent['rows'])):
   tag=str(rc)+'_'+z['label'];ck('same_label_'+tag,p['label']==z['label'])
   ck('nominal_sides_and_remainders_exact_match_'+tag,all(p[k]==encode(z[k]) for k in ('nominal_left_u_rr','nominal_right_u_rr','left_local_M3','right_local_M3','half_kernel_remainder_abs_upper')))
   pa=list(map(Q,p['actual_u_rr']));ia=list(map(Q,z['actual_u_rr_direct_ODE']))
   ck('independent_actual_routes_overlap_'+tag,max(pa[0],ia[0])<=min(pa[1],ia[1]))
   pe=list(map(Q,p['actual_minus_smoothed_u_rr']));ie=list(map(Q,z['actual_minus_smoothed_u_rr_direct_enclosure']))
   ck('independent_error_routes_overlap_'+tag,max(pe[0],ie[0])<=min(pe[1],ie[1]))
   minimum=Q(0) if pe[0]<=0<=pe[1] else min(abs(pe[0]),abs(pe[1]))
   ck('producer_component_lower_recomputed_'+tag,Q(p['fixed_width_error_component_abs_lower'])==minimum)
   if minimum:
    plim=list(map(Q,p['zero_width_limit_actual_minus_trace_average']))
    ck('nonzero_width_limit_same_sign_'+tag,(pe[0]>0 and plim[0]>0) or (pe[1]<0 and plim[1]<0))
   producer_lower2+=c.weights[i]*minimum**2
  ck('producer_complete_H_lower_recomputed_'+str(rc),Q(produced['fixed_width_error_H_norm_squared_lower'])==producer_lower2 and 0<Q(produced['fixed_width_error_H_norm_strict_lower'])**2<producer_lower2)
 # All scientific upstream roots are read back after evaluation.
 c.reauthenticate()
 for k,p in paths.items():pinned(p,PINS[k])
 pinned(args.producer_report,JOIN_PIN)
 print(json.dumps(encode({'status':'PASS','candidate':'unchanged W1_quadratic initial displacement','delta_u':delta,'actual_route':'independent autonomous ODE Taylor recurrence at actual point state; direct exact bump/exponential jets through5; no producer paired actual derivatives imported',
 'joins':joins,'groups':groups,'group_count':len(groups),'wrong_formula_controls':controls,'control_count':len(controls),'input_sha256':PINS,'producer_join_report_sha256':JOIN_PIN,
 'scope':'coordinate initial radial curvature error only; isolated fixed-width lower bounds, no propagated spacetime or physical-relative accuracy claim'}),sort_keys=True))
if __name__=='__main__':main()
