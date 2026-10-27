#!/usr/bin/env python3
"""Exact full23 shell G_x bound from two authenticated inherited reports.
This does not bound G-minus-an-arbitrary-approximation or solve the field equations.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, math


def need(q,reason):
    if not q: raise RuntimeError(reason)


def read(path, digest):
    b=path.read_bytes()
    need(hashlib.sha256(b).hexdigest()==digest,'authenticated report hash mismatch: '+path.name)
    return json.loads(b)

# Exactly differentiate the eight primitive profiles via independent dual-number
# product rules at a generic rational fixture; validate a distinct formula for
# all eight, including source adapters and physical-to-scaled time conversion.
def plus(a,b):return (a[0]+b[0],a[1]+b[1])
def mul(a,b):return (a[0]*b[0],a[1]*b[0]+a[0]*b[1])
def scale(a,k):return (k*a[0],k*a[1])
def power2(a):return mul(a,a)
def comb(*args):
    out=(F(0),F(0))
    for x in args:out=plus(out,x)
    return out

# Dual second entry here is d/dt. Pair sigma1 and sigma2 constants; no shell
# background derivative is introduced because the registered shell is fixed.
h, v, X, L, V, d, tau = [(F(i),F(i+7)) for i in range(1,8)]
a,s1,s2,p,s3,ell,rho=map(F,(2,3,5,7,11,13,17))
beta=comb(scale(h,a),scale(v,s1/3))
alpha=comb(scale(power2(h),a/4),scale(mul(h,v),-s1/6),scale(power2(v),-s2/6))
pi=comb(scale(power2(h),p/4),scale(mul(h,v),s2/2),scale(power2(v),s3/2))
g0=comb(scale(power2(v),4*s2/3),scale(mul(h,v),4*s1),scale(power2(h),-10*a),scale(mul(h,tau),8*a),scale(power2(X),-16*a),scale(power2(L),-48*a),scale(mul(L,d),256*a))
g1=scale(power2(V),8*a)
g2=scale(power2(d),-16*a)
profiles=[scale(g0,-ell*rho),scale(g1,-ell*rho),scale(g2,-ell*rho),scale(comb(scale(alpha,-1),scale(mul(beta,X),2)),2*ell*rho),scale(pi,ell*rho),scale(mul(beta,V),-2*ell*rho),scale(comb(alpha,scale(mul(beta,L),-2)),2*ell*rho),scale(mul(beta,d),-4*ell*rho)]
h0,ht=h;v0,vt=v;x,xt=X;l,lt=L;z,zt=V;d0,dt=d;t,tt=tau
b=a*h0+s1*v0/3;bt=a*ht+s1*vt/3
alphat=a*h0*ht/2-s1*(ht*v0+h0*vt)/6-s2*v0*vt/3
pit=p*h0*ht/2+s2*(ht*v0+h0*vt)/2+s3*v0*vt
expected=[-ell**2*rho*((8*s2*v0*vt/3)+4*s1*(ht*v0+h0*vt)-20*a*h0*ht+8*a*(ht*t+h0*tt)-4*a*(8*x*xt+24*l*lt-64*(lt*d0+l*dt))),-ell**2*rho*16*a*z*zt,ell**2*rho*32*a*d0*dt,2*ell**2*rho*(-alphat+2*bt*x+2*b*xt),ell**2*rho*pit,-2*ell**2*rho*(bt*z+b*zt),2*ell**2*rho*(alphat-2*bt*l-2*b*lt),-4*ell**2*rho*(bt*d0+b*dt)]
for k in range(8):need(ell*profiles[k][1]==expected[k],'dual product-rule fixture profile '+str(k))
need(ell*profiles[7][1]!=-4*ell**2*rho*bt*d0,'reject dropped d_t')
need(ell*profiles[5][1]!=-2*ell**2*rho*bt*z,'reject dropped V_t')
need(ell*profiles[6][1]!=2*ell**2*rho*(alphat-2*bt*l),'reject dropped L_t')

parser=argparse.ArgumentParser();parser.add_argument('--repo-root',type=Path,required=True);args=parser.parse_args()
base=args.repo_root/'research'
shell=read(base/'HDBLAST_COMMON_GAUGE_FORCING_20261016/shell-producer/COMMON_HARMONIC_SHELL_NORMAL.json','559f7587feca262cc68990f2ee79a3aeea22df6955e5aafb698f3326fdacd94a')
full=read(base/'HDBLAST_FULL_ANGULAR_RESPONSE_20261026/shell/verified/FULL23_SHELL_NORMAL.json','0b4a7fd639deec8bfbe65f7b2427e6b084f18195ef954d0add39471528b6f6fe')
need(shell['status']=='PASS_CONTINUOUS_FULL_COMMON_FIRST_HARMONIC_PHYSICAL_SHELL_FORCING' and shell['domain']['scaled_shell']=='0<=x<=3,r=0; entire unit S3','shell scope/domain')
need(full['status']=='PASS' and full['robin_dimension']==18 and full['universal_boundary_source_rank']==8,'full23 packing scope')
B=[[F(v) for v in row] for row in full['source_basis_B_rows_in_LABELS_order']]
S=[[F(v) for v in row] for row in full['source_gram_Btranspose_H_B']]
labels=full['labels'];need(len(labels)==23 and len(B)==23 and all(len(row)==8 for row in B) and len(S)==8 and all(len(row)==8 for row in S),'packing shape')
expected_labels=[f'{ch}:{n}' for ch in ('K','L','N','J','w') for n in (0,2,4)] + [f'{ch}:{n}' for ch in ('v','u','H') for n in (2,4)] + ['T:+','T:-']
need(labels==expected_labels,'exact 23 label ordering')
# Rebuild the positive section-space diagonal in the frozen normalization,
# then check every entry of the reported source Gram independently of S.
def section_weight(label):
    if label.startswith('T:'):return F(1,6)
    ch,l=label.split(':'); l=int(l)
    w={0:F(1),2:F(1,12),4:F(19,720)}[l]
    if ch=='L':w*=2
    if ch=='J':w*=3
    if ch in ('v','u'):w*=2*l*(l+2)
    if ch=='H':w*=F(2,3)*l*(l+2)*(l*(l+2)-3)
    return w
weights=list(map(section_weight,labels))
for i in range(8):
    for j in range(8):
        need(S[i][j]==sum(weights[k]*B[k][i]*B[k][j] for k in range(23)), 'Gram B^T H B '+str((i,j)))
need(S[7][7]!=sum((weights[k] if labels[k] != 'T:+' else 3*weights[k])*B[k][7]*B[k][7] for k in range(23)), 'wrong TT gram weight control')
A=shell['bounds']['actual_field_amplitude_caps']; P=shell['bounds']['parameters']
get=lambda k:F(A[k]['exact_upper']); par=lambda k:F(P[k])
ell=par('ell');rho=max(map(F,P['rho_b_interval']));a=par('a_abs');s1=par('sigma_prime_abs');s2=par('sigma_second_abs');s3=par('sigma_third');p=par('p_abs')
need(ell==F(1,10**6) and min(rho,a,s1,s2,s3,p)>=0,'background parameters')
h,ht,v,vt,X,Xt,L,Lt,V,Vt,d,dt,tau,taut=map(get,['J','J_t','v','v_t','X','X_t','Z','Z_t','V_gauge','V_gauge_t','d','d_t','tau','tau_t'])
need(min(h,ht,v,vt,X,Xt,L,Lt,V,Vt,d,dt,tau,taut)>=0,'nonnegative amplitude bounds')
beta=a*h+s1*v/3; betat=a*ht+s1*vt/3
alphat=a*h*ht/2+s1*(ht*v+h*vt)/6+s2*v*vt/3
pit=p*h*ht/2+s2*(ht*v+h*vt)/2+s3*v*vt
g0t=F(8,3)*s2*v*vt+4*s1*(ht*v+h*vt)+20*a*h*ht+8*a*(ht*tau+h*taut)+4*a*(8*X*Xt+24*L*Lt+64*(Lt*d+L*dt))
g1t=16*a*V*Vt
g2t=32*a*d*dt
profile_bounds=[ell**2*rho*g0t,ell**2*rho*g1t,ell**2*rho*g2t,2*ell**2*rho*(alphat+2*betat*X+2*beta*Xt),ell**2*rho*pit,2*ell**2*rho*(betat*V+beta*Vt),2*ell**2*rho*(alphat+2*betat*L+2*beta*Lt),4*ell**2*rho*(betat*d+beta*dt)]
rows=[sum(abs(B[i][j])*profile_bounds[j] for j in range(8)) for i in range(23)]
# Absolute-value gram majorant, valid for all independent signs and correlations.
norm2=sum(abs(S[i][j])*profile_bounds[i]*profile_bounds[j] for i in range(8) for j in range(8))
need(norm2>=0,'positive H norm square upper bound')
# Choose an integer ceiling from exact rationals, not floating point.
rootceil=math.isqrt(norm2.numerator//norm2.denominator)
if rootceil**2*norm2.denominator < norm2.numerator:rootceil+=1
need(F(rootceil)**2 > norm2,'strict exact square root cap')
need(rootceil <= 58807538 and norm2 < F(58807538)**2, 'documented full norm bound')
need(max(rows) < F(182093220), 'documented modal row bound')
need(max(rows) < math.ceil(max(rows)), 'strict modal integer ceiling')
need(labels[21:]==['T:+','T:-'] and rows[21]==rows[22] and rows[21]<F(301,10**21),'TT coefficient normalization')
need(all(rows[i]==0 for i in (3,4,5,15,16)),'Dirichlet rows zero')
# Explicitly preserve terms in signed matrix before taking absolute values.
need(any(B[0][i]*B[2][i]<0 for i in range(8)) and any(S[i][j]<0 for i in range(8) for j in range(8)),'signed matrix/Gram controls')
print(json.dumps({'status':'PASS','uniform_scaled_shell_domain':'0<=x<=3,r=0','source':'actual registered full 23 coefficient G, not an error','profile_count':8,'rows':len(rows),'dirichlet_rows_exactly_zero':5,'max_modal_row_exact':str(max(rows)),'max_modal_row_ceiling':str(math.ceil(max(rows))),'H_norm_exact_square_upper':str(norm2),'H_norm_rational_upper':str(rootceil),'TT_row_exact':str(rows[21]),'TT_row_strict_upper':'301/10^21','dual_number_checks':8,'independent_gram_entries_checked':64,'wrong_formula_controls':['dropped d_t','dropped V_t','dropped L_t','wrong TT Gram weight'],'scope':'conditional on inherited master/gauge shell caps and fixed rational source factorization; not an independent PDE replay, source approximation error, or C2 solution certificate'},sort_keys=True))
