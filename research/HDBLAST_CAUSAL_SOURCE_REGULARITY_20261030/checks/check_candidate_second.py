#!/usr/bin/env python3
"""Absolute second coordinate jets of a distinct time-filtered candidate.

All arithmetic is rational. This is not a C2 accuracy estimate to the unsmoothed
target. The analytic differentiated principal-wave construction is in the note.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse, hashlib, json

CHECKS=[]
def need(q, reason):
    if not q:
        raise RuntimeError(reason)
    CHECKS.append(reason)

PINS={
 'initial_absolute':'3f6b0a62c0cd21762655353a99e68a6b68b3dc774686d670388e003a257e59a4',
 'higher_transfer':'4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3',
 'higher_time':'6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d',
 'actual_response':'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84',
 'actual_bulk_source':'45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f',
 'G_x':'e9487c9bece0b5df69547a08f74ace87948f890234afb48c4cf785ed60c574c0',
}
def load(path,key):
    b=path.read_bytes()
    need(hashlib.sha256(b).hexdigest()==PINS[key], 'authenticated '+key)
    return json.loads(b)

ap=argparse.ArgumentParser()
ap.add_argument('--repo-root',type=Path,default=Path('.'))
ap.add_argument('--higher-capsule',type=Path)
ap.add_argument('--initial-absolute-report',type=Path)
ap.add_argument('--response-report',type=Path)
ap.add_argument('--gx-report',type=Path)
ap.add_argument('--epsilon',type=Q,default=Q(1,1024))
a=ap.parse_args()
this=Path(__file__).resolve().parent.parent
hi=a.higher_capsule or a.repo_root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030'
i=load(a.initial_absolute_report or this/'reference/INITIAL_ABSOLUTE.json','initial_absolute')
t=load(hi/'reference/TRANSFER.json','higher_transfer')
dt=load(hi/'reference/TIME.json','higher_time')
r=load(a.response_report or a.repo_root/'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json','actual_response')
f=load(a.repo_root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_SOURCES_NORMAL.json','actual_bulk_source')
gx=load(a.gx_report or this/'reference/GX.json','G_x')
need(i['status']=='PASS_ABSOLUTE_QUADRATIC_INITIAL_JETS','absolute smoothed quadratic initial data')
need(i['labels']==t['labels']==dt['labels'] and i['Gram']==t['Gram']==dt['Gram'],'same full23 positive norm')
need(t['direct_reproof']['extension_component_curvature_caps_unchanged'],'same nominal curvature caps for quadratic continuation')
need(gx['status']=='PASS' and gx['rows']==23 and gx['source']=='actual registered full 23 coefficient G, not an error','full23 actual shell derivative')
Gx=Q(gx['H_norm_rational_upper'])
need(0<Gx<=58807538,'actual G_x H cap')
G=Q(r['actual_affine_boundary_G_H_cap'])
FH=Q(f['global_cover']['summary']['positive_bundle_bulk_forcing_H_bound']['F_H_upper'])
need(FH==17387513156573,'genuine forcing magnitude, not response units')
eps=a.epsilon
need(0<eps<=Q(1,4),'positive retained smoothing width')
T=Q(3,2)
op={k:Q(v) for k,v in r['actual_coefficient_H_norm_caps'].items()}
frozen={k:Q(v) for k,v in t['frozen_H_operator_caps'].items()}
der={k:Q(v) for k,v in dt['global']['H_operator_caps'].items()}
need(der['B_x']==0 and der['M_x']==0,'time independent B and M')
need(all(z>=0 for z in list(op.values())+list(frozen.values())+list(der.values())),'nonnegative actual/frozen operator caps')
b=Q(t['b_actual_boundary_X_operator_cap'])
need(b==op['M']/4 and b<1,'actual boundary inverse contraction')
Pinv=1/(1-b)
need(Pinv==Q(t['P_boundary_inverse_X_cap']),'authenticated P inverse')
c=op['A']+op['B']+Q(3,4)*op['D']
c0=frozen['A']+frozen['B']+Q(3,4)*frozen['D']
need(T*c==Q(t['k_actual_bulk_X_operator_cap']),'actual coupled bulk bound')
need(T*c0==Q(t['k0_frozen_bulk_X_operator_cap']),'frozen coupled bulk bound')
data={k:Q(v) for k,v in i['absolute_initial_H_caps'].items()}
Bs=Q(t['direct_reproof']['new_smooth_Bs_X_cap'])
need(max(data['u_r'],data['v'])<=Bs and data['u']<=Q(3,4)*Bs,'initial data compatible with free-wave bound')
free_second=data['u_rr']+data['v_r']
# f0=L0 Bs; constant frozen matrices.  Differentiate before norm bounds.
f0=c0*Bs
f0_x=(frozen['A']+frozen['B'])*free_second+frozen['D']*Bs
f00=frozen['A']*data['v']+frozen['B']*data['u_r']+frozen['D']*data['u']
need(f00<=f0,'initial frozen forcing covered by uniform cap')
W0=Bs+Q(t['direct_reproof']['residual_majorants']['W_minus_B0_whole_cone_X_cap'])
def differentiated_bound(initial_f, forcing_x):
    accel=data['u_rr']+initial_f
    free=max((data['v']+T*accel)/3,data['v_r']+accel)
    return accel,free,Pinv*(free+T*forcing_x+Gx)
acc0,J0,V0=differentiated_bound(f00,f0_x)
# This is d/dx of the COMPLETE coupled candidate forcing before reflection.
comm=der['A_x']+der['B_x']+Q(3,4)*der['D_x']
need(comm==Q(dt['coefficient_commutator_X_to_H_upper']),'actual coefficient derivative commutator')
filter_x=Q(15,4)*FH/eps
f1=c*W0+FH
f1_x=(op['A']+op['B'])*V0+(comm+op['D'])*W0+filter_x
f10=op['A']*data['v']+op['B']*data['u_r']+op['D']*data['u']
# F_eps(0)=0; fixed actual F(0) is not set to zero in the target.
need(f10<=c*Bs<=f1,'filtered candidate initial forcing, no F(0) cancellation')
acc1,J1,V1=differentiated_bound(f10,f1_x)
Z=Pinv*(Bs+T*f1+G)
charge=3*FH*Pinv
exact_error=Q(t['best_reproved_global_C1_error_X_cap'])+charge
need(exact_error<Q('3.541e25'),'explicit source approximation fits old readable C1 threshold')
# The boundary-inverted initial/frozen forcing construction is finite; this
# bound covers the resulting candidate, not a fixed point with La Z in forcing.
caps={'value':Q(3,4)*Z,'x':Z,'r':Z,'xx':V1,'xr':V1,'rr':V1+f1}
readable={'value':Q('9.37e32'),'x':Q('1.250e33'),'r':Q('1.250e33'),
          'xx':Q('1.315e35'),'xr':Q('1.315e35'),'rr':Q('1.315e35')}
for k,v in caps.items():
    need(v<readable[k],'candidate readable '+k)
controls=[]
for name,wrong,correct in [
 ('missing initial f0 acceleration',data['u_rr'],acc0),
 ('missing initial actual lower-order acceleration',data['u_rr'],acc1),
 ('half filter total variation',Q(15,8)*FH/eps,filter_x),
 ('zero actual source magnitude',Q(0),FH),
 ('omit source approximation charge',Q(t['best_reproved_global_C1_error_X_cap']),exact_error),
 ('omit actual coefficient derivative term',(op['A']+op['B'])*V0+op['D']*W0+filter_x,f1_x),
]:
    need(wrong!=correct,'reject '+name)
    controls.append(name)
report={
 'status':'PASS_FILTERED_CANDIDATE_ABSOLUTE_SECOND_COORDINATE_JETS',
 'candidate':'Z_epsilon=P[Bs_quadratic+K(La W0_quadratic+F_epsilon)+LG]',
 'fixed_target_source_unchanged':True,'epsilon':str(eps),'input_sha256':PINS,
 'coordinate_domain':'D+:0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x',
 'norm':'positive full23 H; X=max(sup norm(V)/3,sup norm(V_x),sup norm(V_r)); zero wedge implies norm(V)<=3X/4',
 'initial_absolute_H_caps':{k:str(v) for k,v in data.items()},
 'source_caps':{'F_H':str(FH),'G_H':str(G),'G_x_H':str(Gx),'filtered_F_x_H':str(filter_x)},
 'boundary_inverse_X_cap':str(Pinv),
 'free_wave_X_cap':str(Bs),'free_wave_xx_xr_H_cap':str(free_second),
 'W0_X_cap':str(W0),'W0_initial_forcing_H_cap':str(f00),'W0_initial_acceleration_H_cap':str(acc0),
 'frozen_forcing_H_cap':str(f0),'frozen_forcing_x_H_cap':str(f0_x),
 'differentiated_W0_free_wave_X_cap':str(J0),'W0_x_X_cap':str(V0),
 'W0_rr_H_cap':str(V0+f0),
 'filtered_candidate_forcing_H_cap':str(f1),'filtered_candidate_forcing_x_H_cap':str(f1_x),
 'filtered_candidate_initial_forcing_H_cap':str(f10),'filtered_candidate_initial_acceleration_H_cap':str(acc1),
 'differentiated_candidate_free_wave_X_cap':str(J1),
 'candidate_X_cap':str(Z),'candidate_absolute_coordinate_jet_H_caps':{k:str(v) for k,v in caps.items()},
 'candidate_absolute_coordinate_jet_H_readable_strict_uppers':{k:str(v) for k,v in readable.items()},
 'candidate_C1_change_from_W1_quadratic_X_cap':str(charge),
 'conditional_target_C1_error_X_cap':str(exact_error),'conditional_target_C1_readable_strict_upper':'3.541e25',
 'initial_error_forcing':'r_epsilon(0,r)=F(0,r), generally nonzero; not discarded in any differentiated target-error estimate',
 'regularity_premises':'same qualitative smooth physical-side data/source and flat corner premises as registered construction; quantitative source time derivative supplied by causal convolution, not actual F_x; derivatives at r=0 are one-sided physical derivatives',
 'scope':'absolute coordinate jets of distinct auxiliary numerical candidate on the physical side only; no C2 across doubled shell, no C2 accuracy to unsmoothed target, no completed natural residual or constraint propagation theorem, no nonlinear family or physical evidence',
 'named_checks':CHECKS,'wrong_formula_controls_rejected':controls,
}
print(json.dumps(report,sort_keys=True,indent=2))
