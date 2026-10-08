#!/usr/bin/env python3
"""Exact symbolic/rational checks only: no physical callbacks or saved arrays.

This verifies a finite representation theorem, not actual incoming-state error.
The exact Taylor coefficients specified by the theorem are NOT evaluated here.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_CEILING
import hashlib
import json
from pathlib import Path
import sys
import sympy as S

checks = []
controls = []

def require(value, name):
    if not value:
        raise RuntimeError(name)
    checks.append(name)

def eq(actual, expected, name):
    require(S.simplify(actual - expected) == 0, name)

def reject(actual, expected, name):
    if S.simplify(actual - expected) == 0:
        raise RuntimeError('Mutation accepted: ' + name)
    controls.append(name)

def encode(value):
    if isinstance(value, F):
        return {'numerator': str(value.numerator), 'denominator': str(value.denominator)}
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, list):
        return [encode(v) for v in value]
    return value

def upper_decimal(value):
    with localcontext() as c:
        c.prec = 12
        c.rounding = ROUND_CEILING
        return str(Decimal(value.numerator) / Decimal(value.denominator))

s, k, h, hp, hpp, E, Phi = S.symbols('s k h hp hpp E Phi', nonzero=True)
I = S.I
L = -1/s

def ds(expr):
    return (S.diff(expr, s) + S.diff(expr, h)*hp + S.diff(expr, hp)*hpp
            - 2*I*k*E*S.diff(expr, E) - E*S.diff(expr, Phi))

g = 4*L**2*h - 2*L*hp - hpp
BW = -E*(hp + 2*L*h + 2*I*k*h)
IW = E*(4*k*k - 4*I*k*L + 6*L*L)*h
eq(ds(BW)+IW, E*g, 'W integration-by-parts identity')
BU = -Phi*(hp+2*L*h)-E*h
IU = (6*Phi*L*L-2*(L+I*k)*E)*h
eq(ds(BU)+IU, Phi*g, 'U integration-by-parts identity')
reject(ds(BW)+E*(4*k*k-4*I*k*L+4*L*L)*h, E*g,
       'reject 6L2 to 4L2 cap mutation')
reject(ds(BW)+E*(4*k*k+4*I*k*L+6*L*L)*h, E*g,
       'reject imaginary sign mutation')
reject(ds(BU)+(6*Phi*L*L-2*(L-I*k)*E)*h, Phi*g,
       'reject U imaginary sign mutation')

z = S.symbols('z', real=True)
d = 1-z*z
B = S.exp(1-1/d)
bp = -2*z/d**2
bpp = -2/d**2-8*z*z/d**3+4*z*z/d**4
eq(S.diff(B,z)/B, bp, 'B derivative rational factor')
eq(S.diff(B,z,2)/B, bpp, 'B second derivative rational factor')
eq(bpp, (6*z**4-2)/d**4, 'B second derivative reduced numerator')
GB = 4/s**2 + 2*bp/s - bpp
GZ = z*GB+2/s-2*bp
eq((4*z*B/s**2+2*S.diff(z*B,z)/s-S.diff(z*B,z,2))/B,
   GZ, 'signed source rational factor')
reject(z*GB+2/s+2*bp, GZ, 'reject signed source cross derivative')

tau = S.symbols('tau', real=True)
phi = (S.exp(2*I*k*tau)-1)/(2*I*k)
eq(S.limit(phi,k,0),tau,'Phi removable limit')
# Real symbols are intentional for the state and direct-operator checks.
kr, tr, lr, p, q, c = S.symbols('kr tr lr p q c', real=True, nonzero=True)
real_phi = S.sin(2*kr*tr)/(2*kr)
eq(S.limit(real_phi,kr,0),tr,'real Phi removable limit')
eq(-real_phi-(-S.sin(2*kr*tr))/(2*kr),0,
   'real source kernel has zero canonical c')
R = (kr+3*lr*lr/(2*kr))*c+3*lr*lr*p/(2*kr)+lr*q
P = (kr/3-lr*lr/(2*kr))*c-(2*kr/3+lr*lr/(2*kr))*p+lr*q
normP2 = (2*kr/3+lr*lr/(2*kr))**2+lr*lr
eq((2*kr/3+5*lr*lr/(4*kr))**2-normP2,
   21*lr**4/(16*kr**2), 'pressure positive norm majorant')
eq(S.diff(R,lr)*lr**2+S.diff(R,p)*(-2*kr*q)+S.diff(R,q)*(2*kr*p),
   lr*(R-3*P),'unchanged direct Ward relation')

delta = F(1,128)
D = delta*(2-delta)
H = F(3,8)**63
lam = 1/(5-delta)
pdelta = 2*(1-delta)/D**2
require(D == F(255,16384),'exact cap denominator')
require(1-1/D == -F(16129,255),'exact cap exponent')
require(1-1/D < -63,'cap exponent below minus63')
require(sum(F(1,int(S.factorial(j))) for j in range(5)) > F(8,3),
        'positive exp series proves e greater than8/3')
# For n>=1, n!>=2**(n-1), strict for n>=3. Thus the full
# positive exponential series is strictly below 1+sum_j>=0 2**(-j)=3.
require(1+1/(1-F(1,2)) == 3 and F(3) < 4,
        'factorial geometric majorant proves e less than4')

degree = 112
panels = []
left = delta
while left < F(1,2):
    right = min(F(3,2)*left,F(1,2))
    center = (left+right)/2
    halfwidth = (right-left)/2
    radius = left/2
    Q = 1-center+radius
    T = 5-center-radius
    Ddisk = 1-Q*Q
    require(0 < halfwidth <= radius/2,'panel ratio '+str(len(panels)))
    require(Q < 1 and T > 0 and Ddisk > 0,'panel holomorphic '+str(len(panels)))
    MB = 2*(4/T**2+4*Q/(T*Ddisk**2)+2/Ddisk**2
            +8*Q*Q/Ddisk**3+4*Q*Q/Ddisk**4)
    MZ = Q*MB+2*(2/T+4*Q/Ddisk**2)
    panels.append(dict(left=left,right=right,center_eta=-5+center,
                       halfwidth=halfwidth,radius=radius,Q=Q,T=T,
                       positive_bound=MB,signed_bound=MZ))
    left = right
require(len(panels)==11,'eleven panels')
require(sum(p['right']-p['left'] for p in panels)==F(1,2)-delta,
        'complete interior cover')
interior = {
    source: sum((p['right']-p['left'])*p[key] for p in panels)/2**degree
    for source,key in [('positive_B','positive_bound'),('signed_uB','signed_bound')]
}

rows=[]
for source,rho in [('positive_B',0),('signed_uB',1)]:
    for K in (64,128,256):
        cap = H*(pdelta+rho+2*lam+2*K
                  +delta*(4*K*K+4*K*lam+6*lam*lam))
        total=cap+interior[source]
        rcoef=F(K*K,252)+F(K,294)
        pcoef=F(K**3,162)+F(5*K,1764)
        require(cap > 0 and interior[source] > 0,'positive bounds '+source+str(K))
        require(total*pcoef < F(1,10**17),'pressure analytic truncation below1e-17 '+source+str(K))
        rows.append(dict(source=source,K=K,cap_W_upper=cap,
                         interior_W_upper=interior[source],total_W_upper=total,
                         density_upper=total*rcoef,pressure_upper=total*pcoef,
                         density_upper_decimal=upper_decimal(total*rcoef),
                         pressure_upper_decimal=upper_decimal(total*pcoef),
                         cap_W_upper_decimal=upper_decimal(cap),
                         interior_W_upper_decimal=upper_decimal(interior[source])))

# Independently integrate the two norm majorants in k before substituting
# the all-time rational bounds; no source values enter these checks.
Ksym,Csym,Pisym = S.symbols('Ksym Csym Pisym',positive=True)
rawR = S.integrate(Csym*(lr*kr+3*lr*lr/2)/(4*Pisym**2),(kr,0,Ksym))
rawP = S.integrate(Csym*(2*kr*kr/3+5*lr*lr/4)/(4*Pisym**2),(kr,0,Ksym))
eq(rawR.subs({lr:S.Rational(2,7),Pisym:3})/Csym,
   Ksym*Ksym/252+Ksym/294,'integrated density coefficient')
eq(rawP.subs({lr:S.Rational(2,7),Pisym:3})/Csym,
   Ksym**3/162+5*Ksym/1764,'integrated pressure coefficient')
reject(Ksym**3/324+5*Ksym/1764,
       rawP.subs({lr:S.Rational(2,7),Pisym:3})/Csym,
       'reject halved pressure coefficient')

output = dict(schema='incoming-representation-exact-checks-v1',
              status='CONDITIONAL_ANALYTIC_REPRESENTATION_TRUNCATION_PROVED',
              actual_incoming_state='NOT_ENCLOSED',
              physical_callbacks=0,saved_array_decodes=0,
              taylor_coefficients_evaluated=0,
              full_twelve_case_certificate='UNRESOLVED',metric_calibration='FAIL',
              big_bang_cause='NOT_ESTABLISHED',
              external_novelty='NOT_ASSESSED',
              polynomial_degree=degree,panel_count=len(panels),
              exact_coefficients_per_source=len(panels)*(degree+1),
              checks=checks,checks_passed=len(checks),
              rejected_controls=controls,rejected_controls_count=len(controls),
              panels=panels,bounds=rows,
              source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
dest = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('PRIMARY_EXACT_RECEIPT.json')
dest.write_text(json.dumps(encode(output),indent=2,sort_keys=True)+'\n')
print(json.dumps({k:output[k] for k in ('checks_passed','rejected_controls_count',
                                     'panel_count','polynomial_degree','actual_incoming_state')}))
for row in rows:
    print(row['source'], row['K'], row['density_upper_decimal'],row['pressure_upper_decimal'])
