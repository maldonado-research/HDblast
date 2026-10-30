#!/usr/bin/env python3
"""Independent exact reconstruction of the original-source C2 error gate."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


PINS={
 'initial':'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207',
 'repair':'4a6e4dcc02c8e0ece0d7452b45b4b2f1e37e8221303d20b6fbe8ec6496004002',
 'transfer':'4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3',
 'time':'6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d',
 'refined':'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84',
 'data':'44d05eadf32770e152ba15e8482f2d290a8baafdcf59f20ac281841960b76c1f',
 'raw_initial':'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad',
}
checks=[];controls=[];digests={}
def need(ok,label):
 if not ok:raise RuntimeError(label)
 checks.append(label)
def load(path,role,digest=None):
 raw=Path(path).read_bytes();sha=hashlib.sha256(raw).hexdigest()
 if digest or role in PINS:need(sha==(digest or PINS[role]),'authenticate '+role)
 digests[role]=sha
 return json.loads(raw)
def reject(label,bad,good):
 need(bad!=good,'reject '+label);controls.append(label)
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v) for v in x]
 return x

ap=argparse.ArgumentParser()
ap.add_argument('--repo-root',type=Path,required=True)
ap.add_argument('--higher-capsule',type=Path)
ap.add_argument('--causal-capsule',type=Path)
ap.add_argument('--refined-response',type=Path)
ap.add_argument('--gate-report',type=Path,required=True)
ap.add_argument('--fx-report',type=Path)
ap.add_argument('--fx-sha256')
a=ap.parse_args()
hi=a.higher_capsule or a.repo_root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030'
ca=a.causal_capsule or a.repo_root/'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030'
i=load(hi/'reference/INITIAL.json','initial')
rp=load(hi/'reference/REPAIR.json','repair')
t=load(hi/'reference/TRANSFER.json','transfer')
tm=load(hi/'reference/TIME.json','time')
r=load(a.refined_response or a.repo_root/'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json','refined')
d=load(ca/'reference/CANDIDATE.json','data')
raw=load(a.repo_root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results/INITIAL_REPORT.json','raw_initial')
g=load(a.gate_report,'gate')
need(g['candidate'].startswith('W1_quadratic') and g['target'].startswith('same original actual U'),'original candidate and target identity')
need(g['boundary_error']=={'g':'0','g_x':'0','M_x':'0'},'same-source homogeneous differentiated boundary')
need(t['labels']==tm['labels']==i['labels']==rp['labels'] and
 len(t['labels'])==23 and min(map(F,t['Gram']))>0,'consistent complete positive full23 norm')

act={k:F(v) for k,v in r['actual_coefficient_H_norm_caps'].items()}
frozen={k:F(v) for k,v in t['frozen_H_operator_caps'].items()}
da={k:F(v) for k,v in r['actual_coefficient_error_H_norm_caps'].items()}
der={k:F(v) for k,v in tm['global']['H_operator_caps'].items()}
need(all(act[k]==frozen[k]+da[k] for k in act),'same actual versus frozen operator caps')
need(der['B_x']==der['M_x']==0,'actual time-independent B and M')
span=F(3,2);b=act['M']/4;inv=1/(1-b)
q=F(3409,8000000)
coupled=act['A']+act['B']+F(3,4)*act['D']
comm=der['A_x']+der['B_x']+F(3,4)*der['D_x']
need(b==F(t['b_actual_boundary_X_operator_cap']) and inv==F(t['P_boundary_inverse_X_cap']) and
 span*coupled+b<=q==F(r['retained_coupled_X_operator_cap'])<1,'actual Robin inverse and full feedback contraction')
need(comm==F(tm['coefficient_commutator_X_to_H_upper']),'full actual time commutator')

sv={k:F(v) for k,v in d['initial_absolute_H_caps'].items()}
bs=F(t['direct_reproof']['new_smooth_Bs_X_cap'])
feedback=F(t['direct_reproof']['residual_majorants']['W_minus_B0_whole_cone_X_cap'])
correction=F(t['direct_reproof']['new_actual_bulk_correction_X_cap'])
targeterr=F(t['best_reproved_global_C1_error_X_cap'])
f0=F(t['direct_reproof']['new_reduced_forcing_defect_H_cap'])
need(f0==coupled*correction and bs+feedback==F(d['W0_X_cap']),'same original finite-candidate bounds')

free_second=sv['u_rr']+sv['v_r']
frozen_initial=sum(frozen[k]*sv[field] for k,field in [('A','v'),('B','u_r'),('D','u')])
frozen_time=(frozen['A']+frozen['B'])*free_second+frozen['D']*bs
feedback_time=inv*(frozen_initial+span*frozen_time+act['M']*bs+F(d['source_caps']['G_x_H']))
source_derivative_constant=(act['A']+act['B'])*feedback_time+(act['D']+comm)*feedback+(
 da['A']+da['B'])*free_second+(da['D']+comm)*bs
source_initial=sum(da[k]*sv[field] for k,field in [('A','v'),('B','u_r'),('D','u')])+F(d['source_caps']['F_H'])
correction_time0=inv*(source_initial+span*source_derivative_constant)
correction_time1=span*inv
defect_time0=(act['A']+act['B'])*correction_time0+(act['D']+comm)*correction
defect_time1=(act['A']+act['B'])*correction_time1

# Independently reassemble each initial C1 error from interior nominal,
# P1 interpolation and unchanged PR100 convolution allowances.
nominal=list(map(F,r['refined_actual_minus_nominal_initial_H_caps']))
sm=t['direct_reproof']
errors=[nominal[0]+F(sm['displacement_smoothing_value_H_cap']),
 nominal[1]+F(sm['displacement_smoothing_radial_derivative_H_cap']),
 nominal[2]+F(raw['P1_H_error_caps'][2])+F(sm['velocity_smoothing_H_cap'])]
need(errors==list(map(F,r['initial_smoothed_data_H_error_caps'])),'independent three initial C1 errors')
eu,eur,ev=errors
err2=F(rp['new_uniform_initial_u_rr_H_upper'])
errvr=F(i['continuous']['actual_minus_smoothed_P1_v_r_H_upper'])
accel=err2+act['A']*ev+act['B']*eur+act['D']*eu
initial_free=max((ev+span*accel)/3,errvr+accel)
z0=(initial_free+span*(defect_time0+comm*targeterr))/(1-q)
z1=span*defect_time1/(1-q)
rr0=z0+coupled*targeterr+f0
need(z1==span*span*(act['A']+act['B'])/((1-b)*(1-q)),'independent affine original Fx sensitivity')

for name,value in [('W0_minus_Bs_time_derivative_X_cap',feedback_time),
 ('r0_initial_H_cap',source_initial),('initial_differentiated_acceleration_H_cap',accel),
 ('differentiated_initial_free_wave_X_cap',initial_free)]:
 need(F(g[name])==value,'independent '+name)
affine={
 'r0_time_derivative_H_affine_in_Fx':(source_derivative_constant,F(1)),
 'bulk_correction_time_derivative_X_affine_in_Fx':(correction_time0,correction_time1),
 'reduced_defect_time_derivative_H_affine_in_Fx':(defect_time0,defect_time1),
 'error_xx_xr_H_affine_in_Fx':(z0,z1),
 'error_rr_H_affine_in_Fx':(rr0,z1),
}
for name,(v0,v1) in affine.items():
 need(F(g[name]['constant'])==v0 and F(g[name]['Fx_coefficient'])==v1,'independent '+name)
need(g['initial_error_H_caps']=={k:str(v) for k,v in zip(('e0','e0_r','v0','e0_rr','v0_r'),(eu,eur,ev,err2,errvr))},'all five initial error inputs agree')

# Polynomial forcing r1=A*C_x+B*C_r+D*C with nonconstant coefficients.
# This witness independently checks the complete product-rule structure.
x,r0=F(2,7),F(3,11)
A=2+x+x*x;Ax=1+2*x;B=3-2*x;Bx=-2;D=1+3*x;Dx=3
C=x*x*r0+3*x*r0*r0+5*x**3
Cx=2*x*r0+3*r0*r0+15*x*x;Cr=x*x+6*x*r0
Cxx=2*r0+30*x;Cxr=2*x+6*r0
parts=[A*Cxx,B*Cxr,D*Cx,Ax*Cx,Bx*Cr,Dx*C]
product_derivative=sum(parts)
# Expand r1 as a polynomial in x before differentiating.
def mult(p,q):
 z=[F(0)]*(len(p)+len(q)-1)
 for i,a0 in enumerate(p):
  for j,b0 in enumerate(q):z[i+j]+=a0*b0
 return z
poly=[F(0)]*6
for a0,b0 in [([F(2),F(1),F(1)],[3*r0*r0,2*r0,F(15)]),
 ([F(3),F(-2)],[F(0),6*r0,F(1)]),
 ([F(1),F(3)],[F(0),3*r0*r0,r0,F(5)])]:
 for j,c0 in enumerate(mult(a0,b0)):poly[j]+=c0
expanded=sum(j*coef*x**(j-1) for j,coef in enumerate(poly) if j)
need(expanded==product_derivative,'independent complete product differentiation fixture')
for j,label in [(2,'omit D C_x'),(3,'omit A_x C_x'),(4,'omit B_x C_r'),(5,'omit D_x C')]:
 reject(label,product_derivative-parts[j],expanded)
reject('omit free initial forcing acceleration',F(0),F(1))
reject('omit differentiated boundary source',feedback_time-inv*(act['M']*bs+F(d['source_caps']['G_x_H'])),feedback_time)
reject('use boundary inverse in place of full error inverse',(initial_free+span*(defect_time0+comm*targeterr))*inv,z0)

numeric=None
if g['new_Fx_input_supplied']:
 need(a.fx_report is not None and a.fx_sha256 is not None,'actual full Fx evidence required for numerical C2 audit')
 fxd=load(a.fx_report,'actual_Fx',a.fx_sha256)
 need(fxd['status'].startswith('PASS') and fxd['scope']['actual_original_full23_F_x_bound'] is True and fxd['source_coordinate']=='x=t/ell','actual original full scaled Fx scope')
 fx=F(fxd['genuine_scaled_F_x_H_cap'])
 numeric={'xx':z0+z1*fx,'xr':z0+z1*fx,'rr':rr0+z1*fx}
 need(g['actual_original_bulk_Fx_report_sha256']==a.fx_sha256 and
 F(g['actual_original_bulk_Fx_H_cap'])==fx and fx>=0,'same authenticated actual Fx bound')
 need(g['original_target_second_coordinate_error_H_caps']=={k:str(v) for k,v in numeric.items()},'all three independent original-target C2 numeric caps')
 need(F(g['original_reduced_defect_time_derivative_H_cap'])==defect_time0+defect_time1*fx,'independent actual differentiated reduced defect')
 need(max(numeric.values())<F('6.222e33'),'all three strict readable original-target C2 error caps')
 share=initial_free/(initial_free+span*(defect_time0+defect_time1*fx+comm*targeterr))
 need(F(g['initial_free_wave_fraction_of_differentiated_majorant_numerator'])==share and share>F('0.999999999'),'independent initial-majorant share, not a true-error lower bound')
 need(F(g['actual_Fx_contribution_to_second_error_upper'])==z1*fx,'independent original Fx contribution to second-error majorant')
else:
 need(g['status']=='PASS_AFFINE_ORIGINAL_TARGET_C2_GATE_FX_INPUT_OPEN','affine gate explicitly open without actual Fx')

out={'status':'PASS_INDEPENDENT_ORIGINAL_TARGET_C2_GATE_AUDIT',
 'input_sha256':digests,'same_original_candidate_and_sources':True,
 'affine_error_xx_xr_H':{'constant':z0,'Fx_coefficient':z1},
 'affine_error_rr_H':{'constant':rr0,'Fx_coefficient':z1},
 'actual_Fx_numeric_caps_audited':numeric,
 'named_checks':checks,'wrong_formula_controls_rejected':controls,
 'scope':'Independent internal audit of conditional original-target physical-side C2 coordinate error bound; no across-shell C2, physical relative accuracy, natural constraints, nonlinear family, or cosmological evidence. Numerical error gate stays open if actual full Fx input is absent.'}
print(json.dumps(enc(out),indent=2,sort_keys=True))
