#!/usr/bin/env python3
"""Direct physical-frame harmonic-constraint value intervals at one shell point.

The same fixed actual quadratic gauge source is retained using its genuine
uniform component caps. No source point value, nonlinear spacetime, complete
Einstein residual, or quantitative second-jet error is inferred.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode = True


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:'+str(path))
    return raw


def box(z):
    a,b = map(Q,z)
    need(a <= b,'ordered interval')
    return a,b


def enc(z):
    return [str(z.lo),str(z.hi)]


def sq_upper(z):
    return max(abs(z.lo),abs(z.hi))**2


def ceil_sqrt(q):
    need(q >= 0,'nonnegative norm square')
    n=math.isqrt(q.numerator//q.denominator)
    if Q(n*n)<q:
        n+=1
    need(Q(n*n)>=q,'outward integer norm square root')
    return n


PINS={
 'driver':'efa6a6c2a9c7db596fe024246684f760b807e727fd2ddac8e8126182ba941a73',
 'point':'ab80240cedec3807d40ac960457b384e68857acafa78ee9792730eb75a3dda99',
 'source':'45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f',
 'source_map':'34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615',
 'profiles':'00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc',
 'shell_report':'0b4a7fd639deec8bfbe65f7b2427e6b084f18195ef954d0add39471528b6f6fe',
 'shell_proof':'518926555ba97b62e5b4219d098ad1bed68e8bb4dec72e877c8d4fcdeacae979',
}
ap=argparse.ArgumentParser()
ap.add_argument('--repo-root',type=Path,required=True)
ap.add_argument('--quadrature-capsule',type=Path)
args=ap.parse_args()
root=args.repo_root.resolve()
base=root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
cap=args.quadrature_capsule or root/'research/HDBLAST_CAUSAL_BULK_QUADRATURE_20261030'
paths={'driver':base/'compute/full23_response.py',
 'point':cap/'reference/TRANSFER.json',
 'source':root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_SOURCES_NORMAL.json',
 'source_map':root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json',
 'profiles':base/'compute/response-normal/RESPONSE_PROFILES.json',
 'shell_report':base/'shell/verified/FULL23_SHELL_NORMAL.json',
 'shell_proof':base/'shell/FULL23_SHELL_AND_SUBSIDIARY.md'}
raw={k:pin(p,PINS[k]) for k,p in paths.items()}
point,source,maps,profiles,shell=[json.loads(raw[k]) for k in ('point','source','source_map','profiles','shell_report')]
need(point['status']=='PASS' and point['point']==['9/8','0'],'registered successful W1 point')
need(source['status']=='PASS_CONTINUOUS_ACTUAL_EVOLVED_GEOMETRIC_SOURCE_CELL_ENCLOSURES','actual source certificate')
need(source['source_map']=={'bytes':len(raw['source_map']),'sha256':PINS['source_map'],'status':'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION'},'complete source-map premise')
need(maps['gauge_convention']=='QC=second variation of full metric-lowered gauge covector; QC_background also retained','full lowered quadratic constraint convention')
need(maps['variation_convention']=='second parameter derivative d2/deps2 along affine g0+eps*q,phi0+eps*eta','second variation convention')
need(maps['input_units']['ell']=='t=ell*x and dy=-rho*ell*dr' and maps['input_units']['derivatives']=='ordinary physical y,t derivatives; no factorial normalization','coordinate and derivative units')
need(source['domain']['causal_triangle']=='0<=x<=3,r>=0,x+r<=8','uniform source domain')
need('QC is the actual-metric-lowered quadratic gauge Hessian.' in source['normalization'],'QC has physical frame units already')
need(shell['status']=='PASS' and all('NORMAL_GAUGE_TRACE_JOIN:'+str(l) in shell['check_names'] for l in (0,2,4)),'exact normal gauge source identity at all degrees')

spec=importlib.util.spec_from_file_location('direct_constraint_authenticated_driver',paths['driver'])
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
response=mod.Response(root) # authenticates seven complete scientific roots before scientific imports
I=response.I
need(response.labels==profiles['labels']==[row['label'] for row in point['full23_point_table']],'full registered mode order')
need(response.weights==list(map(Q,profiles['gram'])),'same positive frame input norm')
need(len(response.labels)==23 and response.labels[-2:]==['T:+','T:-'],'full modes retained')
W={row['label']:{k:I(*box(v)) for k,v in row['enclosed_actual_bulk_candidate_W1'].items()} for row in point['full23_point_table']}
need(all(set(z)=={'value','dx','dr'} for z in W.values()),'complete candidate first jets')
need(all(W[label][k].lo==W[label][k].hi==0 for label in response.labels if label.startswith(('L:','v:')) for k in ('value','dx')),'Dirichlet value/time identities only')
# Dirichlet normal derivatives are retained. They are not zero just because
# their boundary value and tangential derivative vanish.
need(any(W[label]['dr'].lo!=0 or W[label]['dr'].hi!=0 for label in response.labels if label.startswith(('L:','v:'))),'nontrivial Dirichlet normal jets retained')
bg=response.context.actual_background(I(0))
rho,warp=bg['rho'],bg['a']
ell=response.ell
need(ell==Q(1,1000000) and rho.lo>0,'positive actual physical normalization')
x=Q(9,8);t=ell*x
ht=I(t-t**3/3,t);invc=I(1-t*t/2,1)
alpha=ell*rho*warp
eta=ell*ht
zeta=ell*invc
physical=1/(ell*rho)
need(physical.lo>0,'physical covector conversion')
QC=source['global_cover']['summary']['absolute_coefficient_caps']
scalar_gram={0:Q(1),2:Q(1,12),4:Q(19,720)}
rows=[]
raw_square=refined_square=Q(0)
nonzero=[]
for l in (0,2,4):
    def get(ch,jet='value'):
        key=ch+':'+str(l)
        return W[key][jet] if key in W else I(0)
    K,N,J,L,v,u,H=[get(ch) for ch in ('K','N','J','L','v','u','H')]
    lam=l*(l+2);beta=Q(2)-Q(2,3)*lam
    Cy=-(get('K','dr')+get('N','dr')-3*get('J','dr'))/2-get('L','dx')+alpha*(4*K+N-3*J)-3*eta*L-lam*zeta*v
    Ct=-get('L','dr')-(get('K','dx')+get('N','dx')+3*get('J','dx'))/2+5*alpha*L-3*eta*(N+J)-lam*zeta*u
    Ca=-get('v','dr')-get('u','dx')+5*alpha*v-4*eta*u+zeta*((-K+N-J)/2+beta*H)
    for component,linear in [('y',Cy),('0',Ct)]+([('a',Ca)] if l else []):
        slot=(0,2,4).index(l) if component!='a' else (2,4).index(l)
        cap=Q(QC['QC/'+component+'/'+str(slot)])
        need(cap>=0,'genuine nonnegative QC component cap')
        gauge=I(-cap,cap)
        direct=physical*linear+gauge
        need(direct.lo<=0<=direct.hi,'direct interval includes zero; no nonzero constraint sign established')
        refined=I(0) if component=='y' else direct
        weight=scalar_gram[l]*(lam if component=='a' else 1)
        need(weight>0,'positive five-covector sphere Gram')
        raw_square+=weight*sq_upper(direct)
        refined_square+=weight*sq_upper(refined)
        rows.append({'label':{'y':'R','0':'T','a':'A'}[component]+':'+str(l),
          'component':component,'degree':l,'output_probability_sphere_Gram_weight':str(weight),
          'scaled_linear_DC':enc(linear),'physical_linear_DC':enc(physical*linear),
          'genuine_uniform_actual_QC_interval':enc(gauge),
          'raw_direct_second_constraint':enc(direct),
          'boundary_refined_second_constraint':enc(refined),
          'raw_interval':enc(direct),'refined_interval':enc(refined),
          'refinement':'exact common-source GN=-2tr(QR)-2QC_y identity and actual shell law' if component=='y' else None})
        if not(direct.lo<=0<=direct.hi):nonzero.append([component,l])
need(len(rows)==8 and len(nonzero)==0,'all three scalar pairs and both gradients')
raw_cap=ceil_sqrt(raw_square);refined_cap=ceil_sqrt(refined_square)
readable=Q('1.319e24')
need(raw_square<readable**2 and refined_square<readable**2,'strict readable complete covector point norm')
# A wrong Dirichlet-normal omission changes the displayed complete expression.
# This check must not infer normal zero from the trace/time boundary data.
l=2
proper=-W['L:2']['dr']-(W['K:2']['dx']+W['N:2']['dx']+3*W['J:2']['dx'])/2-3*eta*(W['N:2']['value']+W['J:2']['value'])-8*zeta*W['u:2']['value']
wrong=proper+W['L:2']['dr']
need((proper.lo,proper.hi)!=(wrong.lo,wrong.hi),'dropped Dirichlet normal corruption rejected')
response.context.reauthenticate()
need(mod.authenticate_local()==response.local_before,'local scientific input readback')
for k,path in paths.items():pin(path,PINS[k])
print(json.dumps({'status':'PASS','point':['9/8','0'],
 'constraint_definition':'mathfrak_C^(2)[W1]=DC[W1]+QC[q,q], same fixed actual first pair and full metric-lowered second-parameter QC',
 'method':'direct full23 candidate first-jet interval assembly; genuine uniform physical QC component caps; exact normal boundary-source refinement',
 'output_norm':'positive orthonormal five-covector frame norm integrated against unit S3 probability measure; scalar Gram g_l, gradient Gram lambda_l*g_l; not Lorentzian signed norm',
 'background':{'ell':str(ell),'rho':enc(rho),'alpha=ell*rho*a':enc(alpha),'eta=ell*tanh(ell*x)':enc(eta),'zeta=ell*sech(ell*x)':enc(zeta),'physical_conversion_1_over_ell_rho':enc(physical)},
 'full8_modal_constraint_table':rows,'rows':rows,
 'raw_direct_covector_norm_squared_upper':str(raw_square),
 'boundary_refined_covector_norm_squared_upper':str(refined_square),
 'raw_direct_covector_norm_integer_upper':str(raw_cap),
 'boundary_refined_covector_norm_integer_upper':str(refined_cap),
 'positive_L2_H_cap':str(refined_cap),
 'covector_point_norm_readable_strict_upper':str(readable),
 'nonzero_constraint_intervals':nonzero,
 'normal_components_exact_zero_by_boundary_identity':[0,2,4],
 'TT_and_scalar_w_contribution':'zero identically in linear DC by transverse-tracelessness and metric-only gauge definition; both TT slots and all scalar field slots remain in authenticated candidate input',
 'Dirichlet_normal_derivatives_retained':True,
 'input_sha256':PINS,
 'scope':'direct conditional harmonic-constraint value enclosure of one finite reduced-system candidate at one shell event; does not require formal target U to have zero constraint; zero is included and no nonzero defect sign is proved; not a constraint derivative bound, full Einstein/Hamiltonian/momentum residual, C2 error, nonlinear family, observable, or Big Bang origin evidence'},sort_keys=True))
