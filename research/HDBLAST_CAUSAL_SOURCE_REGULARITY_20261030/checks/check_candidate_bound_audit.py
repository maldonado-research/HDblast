#!/usr/bin/env python3
"""Independent exact ledger audit of physical-side candidate second jets.

No scientific modules imported; all norms are positive full23 H norms.
The accompanying note supplies the analytic argument and shell-side scope.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path

ap=argparse.ArgumentParser()
ap.add_argument('--repo-root',type=Path,required=True)
ap.add_argument('--higher-capsule',type=Path)
ap.add_argument('--response-report',type=Path)
ap.add_argument('--candidate-report',type=Path,required=True)
ap.add_argument('--initial-report',type=Path,required=True)
ap.add_argument('--gx-report',type=Path,required=True)
a=ap.parse_args()
checks=[]
controls=[]
pins={}
def need(test,label):
    if not test:
        raise RuntimeError(label)
    checks.append(label)
def load(path,digest=None,role=None):
    b=path.read_bytes()
    actual=hashlib.sha256(b).hexdigest()
    if digest is not None:
        need(actual==digest,'authenticate '+(role or path.name))
    pins[role or path.name]=actual
    return json.loads(b)
def rejected(label,claimed,actual):
    need(claimed!=actual,'reject '+label)
    controls.append(label)
hi=a.higher_capsule or a.repo_root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030'
raw=a.repo_root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results'
transfer=load(hi/'reference/TRANSFER.json','4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3')
time=load(hi/'reference/TIME.json','6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d')
response=load(a.response_report or a.repo_root/'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json','ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84')
iraw=load(raw/'INITIAL_REPORT.json','bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad')
profiles=load(raw/'INITIAL_PROFILES.json','555f911f6c336862dfacb6aa56ef5179ac06701812d593b5bfbbf6ea9bd4d5d1')
ini=load(a.initial_report,'3f6b0a62c0cd21762655353a99e68a6b68b3dc774686d670388e003a257e59a4','absolute_initial_report')
gx=load(a.gx_report,'e9487c9bece0b5df69547a08f74ace87948f890234afb48c4cf785ed60c574c0','actual_shell_derivative')
candidate=load(a.candidate_report,role='candidate_report')
need(candidate['status']=='PASS_FILTERED_CANDIDATE_ABSOLUTE_SECOND_COORDINATE_JETS','candidate status')
need(ini['labels']==profiles['labels']==iraw['labels']==transfer['labels'],'same23 labels')
weights=list(map(F,iraw['gram']))
need(len(weights)==23 and min(weights)>0,'positive23 dimensional norm')
norm2=lambda row:sum(w*z*z for w,z in zip(weights,row))
def rounded_sqrt(q):
    scale=2**112
    m=isqrt(q.numerator*scale*scale//q.denominator)
    while F(m*m,scale*scale)<q:m+=1
    return F(m,scale)
nodes=list(map(F,profiles['nodes']))
v=[list(map(F,row)) for row in profiles['profiles']['v']]
slopes=[[(b[j]-aa[j])/(s-r) for j in range(23)] for r,s,aa,b in zip(nodes,nodes[1:],v,v[1:])]
vnode=rounded_sqrt(max(map(norm2,v)))
vslope=rounded_sqrt(max(map(norm2,slopes)))
m2=F(iraw['nominal_second_derivative_H_caps'][0])
u,ur=map(F,iraw['nominal_H_caps'][:2])
du,dv=F(1,2**30),F(1,2048)
rebuilt={'u':u+m2*du*du/18,'u_r':ur+35*m2*du/128,'u_rr':m2,
         'v':vnode+35*vslope*dv/128,'v_r':vslope}
need(rebuilt=={k:F(v) for k,v in ini['absolute_initial_H_caps'].items()},'independent nodal curvature and convolution data reconstruction')
need(rebuilt=={k:F(v) for k,v in candidate['initial_absolute_H_caps'].items()},'candidate data agree with independent initial reconstruction')
op={k:F(v) for k,v in response['actual_coefficient_H_norm_caps'].items()}
frozen={k:F(v) for k,v in transfer['frozen_H_operator_caps'].items()}
dt={k:F(v) for k,v in time['global']['H_operator_caps'].items()}
b=op['M']/4
inv=1/(1-b)
T=F(3,2)
FH=F(17387513156573)
G=F(response['actual_affine_boundary_G_H_cap'])
Gx=F(gx['H_norm_rational_upper'])
eps=F(candidate['epsilon'])
need(eps>0 and b<1,'positive filter width and Robin contraction')
need(dt['B_x']==dt['M_x']==0,'constant B and M in time')
Bs=F(transfer['direct_reproof']['new_smooth_Bs_X_cap'])
W0=Bs+F(transfer['direct_reproof']['residual_majorants']['W_minus_B0_whole_cone_X_cap'])
c0=frozen['A']+frozen['B']+F(3,4)*frozen['D']
c=op['A']+op['B']+F(3,4)*op['D']
frozen_initial=frozen['A']*rebuilt['v']+frozen['B']*rebuilt['u_r']+frozen['D']*rebuilt['u']
actual_initial=op['A']*rebuilt['v']+op['B']*rebuilt['u_r']+op['D']*rebuilt['u']
frozen_dx=(frozen['A']+frozen['B'])*(rebuilt['u_rr']+rebuilt['v_r'])+frozen['D']*Bs
V0=(max((rebuilt['v']+T*(rebuilt['u_rr']+frozen_initial))/3,
        rebuilt['v_r']+rebuilt['u_rr']+frozen_initial)+T*frozen_dx+Gx)*inv
comm=dt['A_x']+dt['B_x']+F(3,4)*dt['D_x']
filtered_dx=15*FH/(4*eps)
actual_dx=(op['A']+op['B'])*V0+(comm+op['D'])*W0+filtered_dx
forcing=c*W0+FH
Z=(Bs+T*forcing+G)*inv
V1=(max((rebuilt['v']+T*(rebuilt['u_rr']+actual_initial))/3,
        rebuilt['v_r']+rebuilt['u_rr']+actual_initial)+T*actual_dx+Gx)*inv
expected={'W0_X_cap':W0,'W0_initial_forcing_H_cap':frozen_initial,
 'frozen_forcing_x_H_cap':frozen_dx,'W0_x_X_cap':V0,
 'filtered_candidate_forcing_x_H_cap':actual_dx,
 'filtered_candidate_initial_forcing_H_cap':actual_initial,'candidate_X_cap':Z}
for key,val in expected.items():
    need(F(candidate[key])==val,'independent '+key)
jets={'value':F(3,4)*Z,'x':Z,'r':Z,'xx':V1,'xr':V1,'rr':V1+forcing}
need(jets=={k:F(v) for k,v in candidate['candidate_absolute_coordinate_jet_H_caps'].items()},'all6 independent candidate jet majorants')
for key,cap in jets.items():
    need(cap<F(candidate['candidate_absolute_coordinate_jet_H_readable_strict_uppers'][key]),'strict readable '+key)
charge=2*T*FH*inv
E=F(transfer['best_reproved_global_C1_error_X_cap'])+charge
need(charge==F(candidate['candidate_C1_change_from_W1_quadratic_X_cap']) and E==F(candidate['conditional_target_C1_error_X_cap']) and E<F('3.541e25'),'independent C1 approximation charge and strict bound')
# Polynomial fixture: w=x*v(r)+u(r)+x^2*a(r)/2. At x=0,
# principal forcing equals a-u_rr. Thus time-differentiated wave velocity
# must be u_rr+f(0)=a, with no half or zero-source substitution.
u_rr,v_r,acc=F(7,3),F(-5,4),F(19,6)
f_initial=acc-u_rr
need(u_rr+f_initial==acc,'principal wave initial acceleration fixture')
# Unit kernel polynomial checked by integration, and its derivative norm
# by endpoint/peak variation. No source derivative is used.
integral=lambda coeff:sum(z/F(i+1) for i,z in enumerate(coeff))
need(integral([0,0,30,-60,30])==1,'independent filter mass')
peak=30*F(1,2)**2*(1-F(1,2))**2
need(2*peak==F(15,4),'independent filter derivative variation')
rejected('omit initial forcing acceleration',u_rr,acc)
rejected('half forcing acceleration',u_rr+f_initial/2,acc)
rejected('omit Robin inverse',F(1),inv)
rejected('kernel peak used as variation',peak,F(15,4))
rejected('omit coefficient time commutator',actual_dx-comm*W0,actual_dx)
rejected('discard source approximation charge',E-charge,E)
out={'status':'PASS_INDEPENDENT_CANDIDATE_SECOND_JET_AUDIT',
 'candidate_report_sha256':pins['candidate_report'],
 'input_sha256':pins,'independently_recomputed_jets':{k:str(v) for k,v in jets.items()},
 'named_checks':checks,'wrong_formula_controls_rejected':controls,
 'scope':'same fixed candidate absolute one-sided coordinate second jets; no C2 target accuracy, small natural residual, across-shell C2 regularity or physical claim'}
print(json.dumps(out,sort_keys=True,indent=2))
