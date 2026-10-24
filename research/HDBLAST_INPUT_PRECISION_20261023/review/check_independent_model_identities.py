#!/usr/bin/env python3
"""Independent exact polynomial/model checks; no scientific capsule imports."""
import argparse, hashlib, json
from fractions import Fraction as F
from pathlib import Path

def add(p,q):
 r=[F(0)]*max(len(p),len(q))
 for i,x in enumerate(p):r[i]+=x
 for i,x in enumerate(q):r[i]+=x
 return trim(r)
def trim(p):
 while len(p)>1 and p[-1]==0:p.pop()
 return p

def mul(p,q):
 r=[F(0)]*(len(p)+len(q)-1)
 for i,x in enumerate(p):
  for j,y in enumerate(q):r[i+j]+=x*y
 return trim(r)
def scale(p,c):return trim([x*c for x in p])
def diff(p):return trim([i*p[i] for i in range(1,len(p))] or [F(0)])
def symmetric_range(p,r):
 e=sum(abs(p[i])*r**i for i in range(1,len(p)))
 return p[0]-e,p[0]+e

def req(ok,label):
 if not ok:raise ValueError(label)
 checks.append(label)

checks=[]
c=F(96714321425692,161839258930731);delta=F(1,1000);k=F(1,9);X=F(9,100000);h0=F(1,100);amin=k-F(1,1000000);amax=F(3,20)
req(c==2/F('1.0357712571566784')-F(4,3),'exact_decimal_model_identity')
req(0<delta<=F(1,500),'actual_delta_within_proved_uniform_interval')
W=[F(1,3),F(0),F(1),F(1,3)];f=[1+c,c]
U=add(scale(mul(diff(W),diff(W)),F(1,2)),scale(mul(W,W),-F(2,3)))
expected=[-F(2,27),0,F(14,9),F(50,27),-F(1,6),-F(4,9),-F(2,27)]
req(U==expected,'bulk_potential_polynomial_identity')
A=add(scale(mul(W,f),F(1,9)),scale(mul(diff(W),diff(f)),-F(1,12)))
B=add(scale(mul(f,f),F(1,36)),scale(mul(diff(f),diff(f)),-F(1,48)))
sigma=add(scale(W,2),scale(f,delta));p=scale(diff(sigma),-F(1,2));a=scale(sigma,F(1,6))
Hconstraint=add(add(mul(a,a),scale(mul(p,p),-F(1,12))),scale(U,F(1,6)))
H=add(scale(A,delta),scale(B,delta**2))
req(H==Hconstraint,'exact_correlated_shell_constraint_identity')
beta_center=-F(9,64)*c*delta;beta_radius=F(799,250)*delta**2
req(beta_center+beta_radius<0,'retained_beta_box_strictly_negative')
req(abs(beta_center)+beta_radius<X,'retained_absolute_beta_less_than_9e_minus_5')
M=diff(U)[1:];Mlo,Mhi=symmetric_range(M,X)
req(0<Mlo and Mhi<4,'riccati_divided_potential_between_zero_and_four')
# U-U0=x² times this positive polynomial on the entire comparison domain.
V=U[2:];Vlo,Vhi=symmetric_range(V,X)
req(Vlo>0,'bulk_potential_above_U0_on_comparison_domain')
radicand_min=k*k-X*X*Vhi/F(6)
radicand_max=h0+X*X/F(3)+k*k
req(0<amin and amin*amin<radicand_min,'positive_a_lower_bound_for_all_z')
req(radicand_max<amax*amax,'actual_a_upper_bound_for_z_in_0_2')
Hlo,Hhi=symmetric_range(H,X)
req(0<Hlo<=Hhi<h0,'correlated_shell_curvature_inside_comparison_domain')
req(F(37,20)-12*amax==F(1,20),'degree_six_memory_margin_given_kappa')
req(F(37,20)-10*amax==F(7,20),'degree_five_residual_margin_given_kappa')
report={'status':'PASS_INDEPENDENT_EXACT_MODEL_AND_SHELL_IDENTITIES','arithmetic':'fractions.Fraction only; no scientific module imports','checks':checks,'conditional_comparison_checks':'The kappa>=37/20 graph damping premise is not proved here; only its proposed margins are checked. z in [0,2] follows by regular-cone first crossing after importing retained theorem existence/sign/positive a premises.','exact_constants':{'c':str(c),'delta':str(delta),'old_beta_center':str(beta_center),'old_beta_radius':str(beta_radius),'comparison_eta_abs':str(X),'h0':str(h0),'amin':str(amin),'amax':str(amax),'M_range':[str(Mlo),str(Mhi)],'shell_H_range_symmetric_beta_domain':[str(Hlo),str(Hhi)]},'polynomials_ascending':{'U':list(map(str,U)),'A':list(map(str,A)),'B0':list(map(str,B)),'shell_H':list(map(str,H))},'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limitations':['No endpoint graph Z or error barrier checked yet.','No numerical branch enclosure or initial gate supplied.','No global uniqueness or physical observable claim.']}
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();args.output.write_text(json.dumps(report,indent=2)+'\n');print(report['status'],len(checks))
