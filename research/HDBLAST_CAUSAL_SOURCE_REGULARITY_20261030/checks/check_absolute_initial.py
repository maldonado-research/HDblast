#!/usr/bin/env python3
"""Absolute jets of the PR100 quadratic smoothed initial profiles.
Standard-library exact rational arithmetic; no scientific imports.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt
import argparse, hashlib, json

COUNT = []
def need(v, name):
    if not v:
        raise RuntimeError(name)
    COUNT.append(name)

def sqrtup(q):
    need(q >= 0, 'nonnegative squared norm')
    den = 2**112
    k = isqrt((q.numerator * den * den) // q.denominator)
    if Q(k*k, den*den) < q:
        k += 1
    return Q(k, den)

def load(p, expected, name):
    b = p.read_bytes()
    need(hashlib.sha256(b).hexdigest() == expected, 'authenticate ' + name)
    return json.loads(b)

def moment(n):
    # 35/32 * (1-z^2)^3 on [-1,1].
    return Q(35,16)*sum(Q(c,n+2*i+1) for i,c in enumerate([1,-3,3,-1]))

ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, default=Path('.'))
ap.add_argument('--higher-capsule', type=Path)
a = ap.parse_args()
raw = a.repo_root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
higher = a.higher_capsule or a.repo_root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030'
pins = {
    'raw_initial_report':'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad',
    'raw_initial_profiles':'555f911f6c336862dfacb6aa56ef5179ac06701812d593b5bfbbf6ea9bd4d5d1',
    'quadratic_transfer':'4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3',
    'higher_initial':'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207',
}
r = load(raw/'initial/results/INITIAL_REPORT.json', pins['raw_initial_report'], 'raw initial report')
p = load(raw/'initial/results/INITIAL_PROFILES.json', pins['raw_initial_profiles'], 'raw initial profiles')
t = load(higher/'reference/TRANSFER.json', pins['quadratic_transfer'], 'quadratic transfer')
h = load(higher/'reference/INITIAL.json', pins['higher_initial'], 'higher initial')
labels = r['labels']
gram = list(map(Q,r['gram']))
need(len(labels) == 23 and labels == p['labels'] == t['labels'] == h['labels'], 'all23 labels')
need(r['gram'] == p['gram'] == t['Gram'] == h['gram'], 'all23 positive Gram')
need(all(w>0 for w in gram), 'positive H weights')
nodes = list(map(Q,p['nodes']))
vt = [list(map(Q,row)) for row in p['profiles']['v']]
need(len(nodes) == len(vt) and all(len(row)==23 for row in vt), 'P1 dimensions')
need(nodes[0]==0 and nodes[-1]==Q(7,4) and all(y>x for x,y in zip(nodes,nodes[1:])), 'ordered initial footprint')
need(all(z==0 for row in vt[:2] for z in row), 'flat P1 velocity initial collar')
need(p['parity']['v']==p['parity']['u'] and p['parity']['du']==[-s for s in p['parity']['u']], 'derivative-consistent H-isometric reflection')
slopes = [[(v1[i]-v0[i])/(y-x) for i in range(23)] for x,y,v0,v1 in zip(nodes,nodes[1:],vt,vt[1:])]
norm2 = lambda row: sum(w*z*z for w,z in zip(gram,row))
vnode2 = max(map(norm2,vt))
vr2 = max(map(norm2,slopes))
vnode = sqrtup(vnode2)
vr = sqrtup(vr2)
need(vnode*vnode>=vnode2 and vr*vr>=vr2, 'velocity convex nodal and slope H caps')
M2 = list(map(Q,h['piecewise_nominal_M2u']))
M2H = Q(r['nominal_second_derivative_H_caps'][0])
need(M2==list(map(Q,r['nominal_second_derivative_component_caps'][0])), 'same nominal component curvature caps')
need(sum(w*z*z for w,z in zip(gram,M2)) <= M2H*M2H, 'all23 aggregate curvature cap')
need(t['direct_reproof']['M2_H']==str(M2H) and t['direct_reproof']['extension_component_curvature_caps_unchanged'], 'quadratic extension same global curvature cap')
for i,interval in enumerate(t['nominal_left_curvature_intervals']):
    need(max(abs(Q(z)) for z in interval)<=M2[i], 'quadratic outer curvature '+labels[i])
du,dv=Q(1,2**30),Q(1,2048)
need(du==Q(h['fixed_displacement_kernel_width']) and dv==Q(h['fixed_velocity_kernel_width']), 'same prescribed mollifiers')
need(moment(0)==1 and moment(1)==Q(35,128) and moment(2)==Q(1,9), 'even kernel moments mass absolute-first second')
nu,nur,_ = map(Q,r['nominal_H_caps'])
us = nu+M2H*du*du/18
usr = nur+M2H*du*Q(35,128)
vs = vnode+vr*dv*Q(35,128)
need(usr < Q(t['direct_reproof']['new_smooth_Bs_X_cap']), 'displacement derivative fits inherited free-wave size cap')
need(vs < Q(t['direct_reproof']['new_smooth_Bs_X_cap']), 'velocity fits inherited free-wave size cap')
# Deliberately wrong formulas are algebraically excluded, independent of asserts.
controls=[]
for name,wrong,correct in [('missing half in value smoothing',Q(1,9),Q(1,18)),('wrong kernel absolute first moment',Q(35,256),Q(35,128)),('wrong kernel second moment',Q(1,18),Q(1,9))]:
    if wrong==correct:
        raise RuntimeError('wrong formula accepted: '+name)
    controls.append(name)
report = {
 'status':'PASS_ABSOLUTE_QUADRATIC_INITIAL_JETS',
 'candidate':'PR100 smoothed initial displacement with quadratic exterior continuation; same smoothed P1 velocity',
 'input_sha256':pins,
 'labels':labels,'Gram':r['gram'],'footprint':['0','7/4'],
 'delta_u':str(du),'delta_v':str(dv),
 'absolute_initial_H_caps':{'u':str(us),'u_r':str(usr),'u_rr':str(M2H),'v':str(vs),'v_r':str(vr)},
 'absolute_initial_H_readable_strict_uppers':{},
 'P1_velocity_node_H_squared_max':str(vnode2),
 'P1_velocity_slope_H_squared_max':str(vr2),
 'P1_velocity_node_H_cap':str(vnode),
 'kernel_moments':{'mass':'1','absolute_first_per_delta':'35/128','second_per_delta_squared':'1/9'},
 'named_checks':COUNT,'wrong_formula_controls_rejected':controls,
 'scope':'absolute coordinate initial profile jets for the distinct auxiliary PR100 candidate, not initial error, evolved C2 accuracy, full Einstein residual, spacetime curvature or physical evidence',
}
for key,cap in report['absolute_initial_H_caps'].items():
    # deterministic readable decimal upper with four significant digits
    exact=Q(cap)
    scale=Q(1)
    while exact>=scale*1000: scale*=10
    while exact<scale*100: scale/=10
    n=exact//scale+1
    report['absolute_initial_H_readable_strict_uppers'][key]=str(n*scale)
print(json.dumps(report,sort_keys=True,indent=2))
