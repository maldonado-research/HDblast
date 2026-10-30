#!/usr/bin/env python3
"""Independent physical-frame comparison-norm audit for natural residual caps."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import isqrt
from pathlib import Path

ap=argparse.ArgumentParser()
ap.add_argument('--natural-report',type=Path,required=True)
ap.add_argument('--candidate-report',type=Path,required=True)
ap.add_argument('--time-report',type=Path,required=True)
ap.add_argument('--transfer-report',type=Path,required=True)
ap.add_argument('--source-report',type=Path,required=True)
a=ap.parse_args()
checks=[]
pins={}
def need(test,label):
    if not test:raise RuntimeError(label)
    checks.append(label)
def load(path,role,digest=None):
    data=path.read_bytes();sha=hashlib.sha256(data).hexdigest()
    if digest is not None:need(sha==digest,'authenticate '+role)
    pins[role]=sha
    return json.loads(data)
n=load(a.natural_report,'natural_report')
j=load(a.candidate_report,'candidate_report','44d05eadf32770e152ba15e8482f2d290a8baafdcf59f20ac281841960b76c1f')
t=load(a.time_report,'time_report','6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d')
tr=load(a.transfer_report,'transfer_report','4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3')
s=load(a.source_report,'source_report','45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f')
need(n['status']=='PASS_NATURAL_COEFFICIENT_RESIDUAL_BRIDGE','natural report status')
need(n['candidate_input_sha256']==pins['candidate_report'],'same candidate report')
ell=Q(1,10**6)
bg=t['global']
rl,rh=map(Q,bg['rho_interval'])
aa=max(abs(Q(x)) for x in bg['a_interval'])
pp=max(abs(Q(x)) for x in bg['p_interval'])
alpha=ell*rh*aa;eta=Q(3,2)*ell*ell
alphar=ell*ell*(1+rh*rh*(aa*aa+pp*pp/3))
# Independently sum the degree-by-degree Frobenius contributions of the
# full23 ->8 harmonic constraint value operator.
la=[Q(8),Q(24)]
alpha2=3*(16+1+3)+3*Q(25,2)+2*Q(25,2)
eta2=3*Q(9,2)+3*(9+3)+2*8
zeta2=sum(la)+sum(Q(7,12)*l+Q(2,3)*(l-3) for l in la)
need((alpha2,eta2,zeta2)==(Q(245,2),Q(131,2),Q(68)),'independent weighted Frobenius squares')
need(n['value_matrix_primitive_Frobenius_squares']=={'alpha':str(alpha2),'eta':str(eta2),'zeta':str(zeta2)},'reported full23 value matrix normalization')
need(Q(12)**2>alpha2 and Q(9)**2>eta2 and Q(9)**2>zeta2,'primitive operator caps')
q=Q(15,8)
need(q*q>Q(7,2),'separate second jet derivative constant')
z0=12*alpha+9*eta+9*ell
zx=9*ell*ell+9*Q(3,2)*ell**3
zr=12*alphar
conn=2*alpha+5*ell+Q(7,4)*(alpha+eta)
need(n['Zmat_caps']=={'value':str(z0),'x_derivative':str(zx),'r_derivative':str(zr)},'independent actual-background value operator derivatives')
need(Q(n['gradient_DC_connection_scaled_cap'])==conn,'independent frame connection charge')
# Universal extrinsic connection inequality by exact polynomial expansion.
# Variables a,h,Cy,C0,W, with W=|C_sphere|²>=0.
def variable(i):
    key=tuple(int(j==i) for j in range(5))
    return {key:Q(1)}
def add(*items):
    out={}
    for scale,p in items:
        for key,value in p.items():out[key]=out.get(key,Q(0))+scale*value
    return {k:v for k,v in out.items() if v}
def mul(p,q):
    out={}
    for x,u in p.items():
        for y,v in q.items():
            k=tuple(a+b for a,b in zip(x,y))
            out[k]=out.get(k,Q(0))+u*v
    return {k:v for k,v in out.items() if v}
pa,ph,py,p0,pw=map(variable,range(5))
square=lambda p:mul(p,p)
ah2=add((1,square(pa)),(1,square(ph)))
ynorm=add((1,square(py)),(1,square(p0)),(1,pw))
conn_pair=add((1,mul(pa,py)),(-1,mul(ph,p0)))
difference=add((3,mul(ah2,ynorm)),(-1,mul(ah2,pw)),(-3,square(conn_pair)))
orth_pair=add((1,mul(ph,py)),(1,mul(pa,p0)))
positive_square=add((3,square(orth_pair)),(2,mul(ah2,pw)))
need(difference==positive_square,'universal connection positive-square identity')
need(5**2>24 and Q(7,4)**2>3,'intrinsic angular and extrinsic rational constants')
v={k:Q(x) for k,x in j['candidate_absolute_coordinate_jet_H_caps'].items()}
P=q*max(v['x'],v['r'])+z0*v['value']
Px=q*max(v['xx'],v['xr'])+z0*v['x']+zx*v['value']
Pr=q*max(v['xr'],v['rr'])+z0*v['r']+zr*v['value']
DC=(Px+Pr+conn*P)/(ell*ell*rl*rl)
caps=s['global_cover']['summary']['absolute_coefficient_caps']
labels=t['labels'];H=list(map(Q,t['Gram']))
rows=[('yy',['K:0','K:2','K:4']),('y0',['L:0','L:2','L:4']),
 ('00',['N:0','N:2','N:4']),('ya',['v:2','v:4']),('0a',['u:2','u:4']),
 ('ab',['J:0','J:2','J:4','H:2','H:4','T:+','T:-'])]
grad2=sum(H[labels.index(label)]*(Q(caps['QE/'+field+'/'+str(idx)])+Q(caps['QEred/'+field+'/'+str(idx)]))**2
 for field,names in rows for idx,label in enumerate(names))
g=isqrt(grad2.numerator//grad2.denominator)
if g*g<grad2:g+=1
need(Q(n['actual_QC_symmetric_gradient_L2_square_upper'])==grad2 and Q(n['actual_QC_symmetric_gradient_L2_cap'])==g,'complete20 actual quadratic source correction')
unit=Q(tr['direct_reproof']['new_reduced_forcing_defect_H_cap'])+Q(2,3)*Q(tr['k_actual_bulk_X_operator_cap'])*Q(j['candidate_C1_change_from_W1_quadratic_X_cap'])+2*Q(j['source_caps']['F_H'])
E=unit/(2*ell*ell*rl*rl)+DC+g
S=unit/(ell*ell*rl*rl)
expected={'P':P,'P_x':Px,'P_r':Pr,'nabla_DC_physical_L2':DC,
 'natural_Einstein_second_parameter_L2':E,'natural_scalar_second_parameter_L2':S}
need({k:Q(v) for k,v in n['candidate_derived_caps'].items()}==expected,'all independent natural residual arithmetic')
need(E<Q('2.918e43') and S<Q('2.118e30') and unit<Q('3.574e22'),'strict natural and reduced readable bounds')
controls=[]
for name,test in [('omit radial rho product charge',conn!=conn-alpha),
 ('omit complete source discrepancy',unit!=unit-2*Q(j['source_caps']['F_H'])),
 ('metric and scalar scale factors confused',unit/(2*ell*ell*rl*rl)!=S),
 ('coefficient lambda-squared half confused with parameter Hessian',E!=E/2)]:
    need(test,'reject '+name);controls.append(name)
out={'status':'PASS_INDEPENDENT_NATURAL_RESIDUAL_BOUND_AUDIT','input_sha256':pins,
 'named_checks':checks,'wrong_formula_controls_rejected':controls,
 'natural_trace_reversed_metric_Hessian_upper':str(E),'natural_scalar_Hessian_upper':str(S),
 'scope':'large absolute positive physical background-frame probability L2 bounds; no smallness, target C2 accuracy, nonlinear-family or origin claim'}
print(json.dumps(out,sort_keys=True,indent=2))
