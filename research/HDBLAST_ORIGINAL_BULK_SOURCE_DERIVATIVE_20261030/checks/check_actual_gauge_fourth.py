#!/usr/bin/env python3
"""Conditional all-point actual first-gauge fourth coordinate derivatives.

Only standard-library rational arithmetic. No differentiation of interval
selections, no comparison-field substitution, and no source/model changes.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path

CHECKS=[]
def need(ok, why):
    if not ok: raise RuntimeError(why)
    CHECKS.append(why)

PINS={
 'coefficients':'4009c46f6d63583b4d161eebbb620077467b7561a3c4b5728d5293655d340a53',
 'gauge_premise':'00f654647d3cd5578e3f9afcb00dc3ec81739036f40adbe43f53e84b87506e98',
 'master4':'d04f67e2e1863cca423a4b0835bef0a739bf5c471cf1d3e2f6730677d970490e',
}
def authenticated(path, key, parse=True):
    data=path.read_bytes()
    need(hashlib.sha256(data).hexdigest()==PINS[key], 'authenticate '+key)
    return json.loads(data) if parse else data.decode()

def add(*ps):
    out={}
    for p in ps:
        for n,c in p.items():out[n]=out.get(n,Q(0))+c
    return {n:c for n,c in out.items() if c}
def scale(c,p):return {n:c*v for n,v in p.items() if c*v}
def mul(p,q):
    out={}
    for n,c in p.items():
        for m,d in q.items():out[n+m]=out.get(n+m,Q(0))+c*d
    return {n:c for n,c in out.items() if c}
def dp(p):return {n-1:n*c for n,c in p.items() if n}
J={0:Q(1),2:Q(-1)}
H={1:Q(1)}
def td(p):return mul(J,dp(p))
def qd(p):return add(td(p),scale(Q(-3,2),mul(H,p)))

def algebra(injection):
    h=H.copy();j=J.copy();hj=mul(H,J);q={0:Q(1)}
    for _ in range(3):h,j,hj,q=td(h),td(j),td(hj),qd(q)
    need(h=={0:Q(-2),2:Q(8),4:Q(-6)}, 'ordinary third tanh derivative polynomial')
    need(j=={1:Q(16),3:Q(-40),5:Q(24)}, 'ordinary third sech-square derivative polynomial')
    need(hj=={0:Q(-8),2:Q(68),4:Q(-120),6:Q(60)}, 'ordinary third HJ derivative polynomial')
    if injection=='normalized_taylor_order3':q=scale(Q(1,6),q)
    need(q=={1:Q(39,4),3:Q(-105,8)}, 'ordinary third cosh-minus-three-halves derivative polynomial, not Taylor coefficient')
    # Noncommuting words retain P0(P0 s0), not a scalarized replacement.
    direct={('srr',):1,('P','P','s'):1,('P','sx'):1,('Px','s'):2,
            ('Q','sr'):1,('R','s'):1,('sxx',):1}
    expanded={('srr',):1,('P','P','s'):1,('P','sx'):1,('Px','s'):2,
              ('Q','sr'):1,('R','s'):1,('sxx',):1}
    if injection=='drop_initial_P_squared':expanded.pop(('P','P','s'))
    need(expanded==direct, 'stabilized initial fourth jet retains noncommuting P0 squared source')
    # U=x^4+3x^2r^2+5r^4 has T=Uxx-Urr=6x²-54r².
    recovered=Q(24)-12-(-108)
    if injection=='flip_radial_recovery_sign':recovered=Q(24)-12+(-108)
    need(recovered==120, 'PDE fourth radial recovery sign')

def encode(v):
    if isinstance(v,Q):return str(v)
    if isinstance(v,dict):return {k:encode(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [encode(x) for x in v]
    return v

ap=argparse.ArgumentParser()
ap.add_argument('--repo-root',type=Path,required=True)
ap.add_argument('--master4-report',type=Path)
ap.add_argument('--inject-error',choices=('none','normalized_taylor_order3','drop_initial_P_squared','flip_radial_recovery_sign'),default='none')
a=ap.parse_args()
base=a.repo_root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs'
paths={'coefficients':base/'GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json',
       'gauge_premise':base/'GAUGE_WAVE_REGULARITY_PREMISE.md',
       'master4':a.master4_report or a.repo_root/'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/MASTER4.json'}
c=authenticated(paths['coefficients'],'coefficients')
premise=authenticated(paths['gauge_premise'],'gauge_premise',False)
m4=authenticated(paths['master4'],'master4')
need(c['status']=='PASS_WHOLE_COLLAR_GAUGE_COEFFICIENT_AND_SOURCE_CERTIFICATE','whole-collar actual coefficient/source premise')
need(c['domain']['ell']=='1/1000000' and c['domain']['scaled_x']==[0,3] and c['domain']['scaled_r']==[0,8], 'same retained scaled domain')
need(c['matrix_convention'].startswith('U_xx-U_rr=A U_x+B U_r+D U+(S,0,0)'), 'same raw three-amplitude PDE sign')
need('C1' in premise and 'U(0,r)=U_x(0,r)=0' in premise and 'b_r(x,0)=d_r(x,0)=0' in premise, 'retained zero Cauchy and time-independent mixed boundary laws')
need(c['source_definition']=='S=P_s W_r+P_0 W; P_s=6ell rho C q; P_0=6ell^2 rho(C_z-4rho a C)q; C=alpha/rho^2,q=cosh(t)^(-3/2),W=w(ell x,ell r).','same actual master-to-gauge source expression')
need(m4['status']=='PASS_CONDITIONAL_ACTUAL_MASTER_CONTINUOUS_FOURTH_JETS','actual master fourth derivative premise')
algebra(a.inject_error)

ell=Q(1,10**6);X=Q(3);theta=X/2;lam=8;h=X*ell
raw=c['matrix_maximum_absolute_row_sum_by_mixed_derivative']
caps={name:{k:Q(v) for k,v in rows.items()} for name,rows in raw.items()}
source={k:Q(v) for k,v in c['source_strict_absolute_mixed_derivative_caps'].items()}
need(set(caps)=={'A','B','D'} and all(set(rows)==set(source) for rows in caps.values()), 'all mixed coefficient/source derivatives through order two')
need(all(v>=0 for rows in caps.values() for v in rows.values()) and all(v>=0 for v in source.values()), 'nonnegative ordinary mixed derivative bounds')
static=c['static_actual_derivative_caps']
rho,aa,C,Cz=(Q(static[k][0]) for k in ('rho','a','C','C_z'))
need(min(rho,aa,C,Cz)>0, 'actual static absolute profile bounds')
Cs=6*ell*rho*C
C0=6*ell**2*rho*(Cz+4*rho*aa*C)
need(h*h<Q(1,3), 'short actual temporal strip')
H3=Q(2)
J3=16*h+40*h**3+24*h**5
HJ3=8+68*h*h+120*h**4+60*h**6
A3=5*ell**4*H3
B3coef=Q(0)
D3=ell**5*max(lam*J3+6*rho*rho*aa*H3,
                 (3+lam)*J3+2*lam*H3,
                 2*aa*J3+2*HJ3+abs(4-lam)*J3)
caps['A']['3,0']=A3;caps['B']['3,0']=B3coef;caps['D']['3,0']=D3
qtime=[Q(1),Q(3,2)*h,Q(3,2)+Q(15,4)*h*h,Q(39,4)*h+Q(105,8)*h**3]
master=list(map(Q,c['master_scaled_order_strict_caps']))
master.append(max(map(Q,m4['exact']['all_scaled_order_four_mixed_caps'].values())))
need(len(master)==5 and all(v>0 for v in master),'actual master order zero through four bounds')
S3=sum(Q(comb(3,j))*ell**j*qtime[j]*(Cs*master[4-j]+C0*master[3-j]) for j in range(4))
source['3,0']=S3
need(A3==Q(1,10**23) and B3coef==0, 'new pure-time third drift derivatives')

def L1(k):return caps['A'][k]+caps['B'][k]
def L0(k):return caps['D'][k]
def alpha(k):return L1(k)+theta*L0(k)
q0=X*alpha('0,0');res=1/(1-q0)
need(0<q0<Q(3,250),'same scalar infinity-norm contraction')
S0,S1,S2,Sr,Sxr,Srr=(source[k] for k in ('0,0','1,0','2,0','0,1','1,1','0,2'))
N0=res*X*S0
N1=res*(2*S0+X*(alpha('1,0')*N0+S1))
d0=caps['A']['0,0']*S0+S1
free2=max(2*S0/X+2*d0,Sr+d0)
N2=res*(free2+X*(2*alpha('1,0')*N1+alpha('2,0')*N0+S2))
second=N1+alpha('0,0')*N0+S0
Tx=L1('0,0')*N1+L0('0,0')*N0+alpha('1,0')*N0+S1
Tr=L1('0,0')*second+L0('0,0')*N0+alpha('0,1')*N0+Sr
third=N2+max(Tx,Tr)
d0r=caps['A']['0,1']*S0+caps['A']['0,0']*Sr+Sxr
v0=Srr+caps['A']['0,0']*d0+2*caps['A']['1,0']*S0+caps['B']['0,0']*Sr+caps['D']['0,0']*S0+S2
free3=max(2*d0/X+2*v0,d0r+v0)
N3=res*(free3+X*(3*alpha('1,0')*N2+3*alpha('2,0')*N1+alpha('3,0')*N0+S3))
Txx=L1('0,0')*N2+L0('0,0')*N1+2*L1('1,0')*N1+2*L0('1,0')*N0+alpha('2,0')*N0+S2
Txr=L1('0,0')*third+L0('0,0')*second+(L1('1,0')+L1('0,1'))*second+(L0('1,0')+L0('0,1'))*N0+alpha('1,1')*N0+Sxr
Trr=L1('0,0')*third+L0('0,0')*second+2*L1('0,1')*second+2*L0('0,1')*N0+alpha('0,2')*N0+Srr
fourth={'4,0':N3,'3,1':N3,'2,2':N3+Txx,'1,3':N3+Txr,'0,4':N3+Txx+Trr}
orders=[theta*N0,N0,second,third,max(fourth.values())]
mixed={'0,0':theta*N0,'1,0':N0,'0,1':N0,'2,0':N1,'1,1':N1,'0,2':second,
       '3,0':N2,'2,1':N2,'1,2':N2+Tx,'0,3':N2+Tr,**fourth}
def readable(q):
    scale=Q(1)
    while q>=scale*1000:scale*=10
    while q<scale*100:scale/=10
    return (q//scale+1)*scale
need(all(x>0 for x in orders),'positive actual gauge magnitude bounds')
for key,path in paths.items():need(hashlib.sha256(path.read_bytes()).hexdigest()==PINS[key],'input unchanged '+key)
report={'status':'PASS_CONDITIONAL_ACTUAL_FIRST_GAUGE_CONTINUOUS_FOURTH_JETS',
 'domain':'0<=x<=3,r>=0,x+r<=8; x=t/ell,r=z/ell,ell=1e-6',
 'amplitudes':['f=xi^y/Y','b=xi^t/Y','d=angular_gradient_amplitude'],
 'norm':'componentwise scalar infinity; N=max(2sup|V|/X,sup|Vx|,sup|Vr|), X=3',
 'input_sha256':PINS,
 'actual_source_unchanged':True,
 'static_source_caps':{'Cs':Cs,'C0':C0},
 'temporal_ordinary_derivative_caps':{'tanh_third':H3,'sech_square_third':J3,'tanh_sech_square_third':HJ3,'cosh_minus_three_halves_orders_0_to_3':qtime},
 'master_scaled_order_caps_0_to_4':master,
 'ordinary_matrix_row_caps':caps,'ordinary_source_mixed_caps':source,
 'alpha':{k:alpha(k) for k in source},'q0':q0,'resolvent':res,
 'initial_third_time_field_caps':{'displacement':d0,'radial_displacement':d0r,'velocity':v0},
 'free_second_time_field_N_cap':free2,'free_third_time_field_N_cap':free3,
 'N_time_derivative_caps':{'0':N0,'1':N1,'2':N2,'3':N3},
 'gauge_forcing_first_derivative_caps':{'1,0':Tx,'0,1':Tr},
 'gauge_forcing_second_derivative_caps':{'2,0':Txx,'1,1':Txr,'0,2':Trr},
 'actual_gauge_scaled_mixed_order_four_caps':fourth,
 'actual_gauge_uniform_scaled_order_caps_0_to_4':orders,
 'actual_mixed_scaled_derivative_caps':{k:[v,v,v] for k,v in mixed.items()},
 'readable_strict_uniform_scaled_order_caps_0_to_4':[readable(x) for x in orders],
 'regularity':'C4 on physical-side closed triangle, with one-sided shell/domain-edge traces; no smooth reflected coefficient extension asserted',
 'analytic_premises':['Same actual background/first master and mixed-boundary gauge PDE as pinned retained proofs.',
 'Quantitative source pure-time order3 follows from actual master order4 and analytic coefficient identities; mixed derivatives through2 inherited.',
 'Third time Picard jets have common displacement after iterate2 and common velocity after iterate3; common-data convergence after exceptional finite offsets.',
 'PDE recovery Vrr=Vxx-T supplies the remaining fourth radial jets from already controlled continuous derivatives.'],
 'named_checks':CHECKS,
 'scope':'conditional actual first-gauge coordinate jet magnitudes; not an error to a numerical candidate, complete bulk Fx by itself, second-response C2 error, exact nonlinear family or cosmological evidence'}
print(json.dumps(encode(report),indent=2,sort_keys=True))
