#!/usr/bin/env python3
"""Independent exact first-only curvature-defect and Buchert algebra."""
import json
from pathlib import Path
import sympy as s

def need(x,what):
    if not x: raise RuntimeError(what)
e,E,A,V0,T,S,V=s.symbols('e E A V0 T S V',positive=True)
# Round criticality of C=S/V-6/A^2*(V0/V)^(2/3).
C=S/V-6/A**2*(V0/V)**s.Rational(2,3)
dC=s.diff(C,S)*T/A**2+s.diff(C,V)*T/2
need(s.simplify(dC.subs({S:6*V0/A**2,V:V0}))==0,'round criticality for arbitrary trace integral')
# Conformal scalar-curvature law, integrated by parts on closed S3:
# integral R[e^(2f)k] dV[e^(2f)k] = integral e^f(6+2|Df|^2).
Y2=s.Rational(1,6); DY2=s.Rational(4,3)
D=1+s.Rational(9,2)*e**2*E**2*Y2
N=6+e**2*E**2*(3*Y2+2*DY2)
meanR=(N/D/A**2).series(e,0,3).removeO()
ref=(6/A**2*D**(-s.Rational(2,3))).series(e,0,3).removeO()
meanR2=s.expand(meanR).coeff(e,2)*2
ref2=s.expand(ref).coeff(e,2)*2
C2=s.simplify(meanR2-ref2)
need(C2==10*E**2/(3*A**2),'curvature defect coefficient')
# Moving-volume average commutator and Lorentzian Gauss identity.
M,Mp,T2,shear2,meanR4,meanR3=s.symbols('M Mp T2 shear2 meanR4 meanR3')
# M=<theta>; <theta'>=M'-<theta^2>+M^2.
Rdiff=2*(Mp-T2+M**2)+s.Rational(4,3)*T2+shear2
QB=s.Rational(2,3)*(T2-M**2)-shear2
need(s.expand(Rdiff+QB-(2*Mp+s.Rational(4,3)*M**2))==0,'Buchert identity')
# Gaussian first metric B=2(E Y I-J HessY).
Ep,j,W=s.symbols('Ep j W',real=True)
# <HessY:HessY>=8, DeltaY=-8Y, tracefree Hess norm=8-(64/3)<Y²>=40/9.
Hess2=s.Rational(8)
shear_first_norm=s.simplify(Hess2-s.Rational(64,3)*Y2)
need(shear_first_norm==s.Rational(40,9),'tracefree Hess norm')
theta1coef=3*Ep+8*j
QB2=s.expand(2*(s.Rational(2,3)*theta1coef**2*Y2-j**2*shear_first_norm))
expect=2*Ep**2+s.Rational(32,3)*Ep*j+s.Rational(16,3)*j**2
need(s.expand(QB2-expect)==0,'Buchert first-only coefficient')
QB2W=s.expand(QB2.subs(Ep,W-j))
need(s.expand(QB2W-(2*W**2+s.Rational(20,3)*W*j-s.Rational(10,3)*j**2))==0,'Buchert u,P coefficient')
# Background de Sitter check of curvature sign convention.
H,Hp=s.symbols('H Hp')
background=6/A**2+6*Hp+12*H**2
need(s.simplify(background.subs(Hp,1/A**2)-12*(H**2+1/A**2))==0,'de Sitter background curvature sign')
print(json.dumps({'status':'PASS','imports_source_or_producer':False,'C2':str(C2),'mean_R3_second_conformal_path':str(meanR2),'round_reference_second_conformal_path':str(ref2),'tracefree_Hessian_mean_norm':str(shear_first_norm),'Buchert_Q2_Eprime_j':str(QB2),'Buchert_Q2_W_j':str(QB2W),'Buchert_convention':'Q_B=(2/3)Var(theta)-<|shear|^2>; norm has no factor1/2','scope':'standard integrated scalar-curvature first/second variation and averaged Gauss identity applied to retained first Gaussian shell metric; no actual u/P evaluation'},indent=2,sort_keys=True))
