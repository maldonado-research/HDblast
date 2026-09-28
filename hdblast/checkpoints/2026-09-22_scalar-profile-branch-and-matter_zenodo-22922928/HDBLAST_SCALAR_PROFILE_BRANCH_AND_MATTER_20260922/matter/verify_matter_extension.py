#!/usr/bin/env python3
"""Exact algebra checks for the conditional brane-matter extension.

Requires SymPy. These checks do not solve a quantum state or a bulk evolution.
Run with python -B verify_matter_extension.py; explicit guards survive -O.
"""
import json
from pathlib import Path
import sympy as s

out = []
def check(name, lhs, rhs=0):
    residual = s.simplify(s.expand(lhs-rhs))
    if residual != 0:
        raise RuntimeError(f"{name}: residual {residual}")
    out.append({"name": name, "pass": True, "residual": str(residual)})

# k denotes kappa_5^2 (dimension mass^-3), not an AdS curvature.
k, sig, sp, rho, p, j, v, w, U, Up, H, E = s.symbols(
    "k sigma sigma_prime rho p j v normal_phi U U_prime H Weyl")
q0, qs = (sig-k*(2*rho+3*p))/6, (sig+k*rho)/6
traceS = -4*sig/k-rho+3*p
Sdiag = [-sig/k-rho]+[-sig/k+p]*3
Kdiag = [q0]+[qs]*3
for i in range(4):
    check(f"Israel mixed component {i}", Kdiag[i]-sum(Kdiag), k*Sdiag[i]/2)
check("pure tension spatial", qs.subs(rho,0), sig/6)
check("pure tension temporal", q0.subs({rho:0,p:0}), sig/6)
wn = -(sp+k*j)/2
check("two bulk scalar boundary variation", -2*wn/k-sp/k-j)
check("pure tension scalar", wn.subs(j,0), -sp/2)
check("Codazzi signed exchange", -sp/k-j, 2*wn/k)
check("total brane energy flux", (sp/k+j)*v, -2*wn*v/k)

# Direct homogeneous chi energy identity, including mass variation.
chi, chidot, mass2, mass2p = s.symbols("chi chidot mass_squared mass_squared_prime")
chidotdot = -3*H*chidot-mass2*chi
rhochi = (chidot**2+mass2*chi**2)/2
pchi = (chidot**2-mass2*chi**2)/2
drhochi = chidot*chidotdot+mass2*chi*chidot+mass2p*v*chi**2/2
check("chi signed energy exchange", drhochi+3*H*(rhochi+pchi), mass2p*chi**2*v/2)
J, Q = s.symbols("J Q")
check("decay energy transfers cancel", (J*v-Q)+Q, J*v)

# Project bulk stress directly using Gauss, in orthonormal shell coordinates.
X = -v**2+w**2
Ttrace = -s.Rational(3,2)*X-5*U  # k*T^A_A
Tnn = w**2-X/2-U               # k*T_nn
T00 = v**2+X/2+U
Tii = -X/2-U
F00 = s.Rational(2,3)*(T00-(Tnn-Ttrace/4))
Fii = s.Rational(2,3)*(Tii+(Tnn-Ttrace/4))
check("Gauss scalar stress temporal", F00, v**2/4-w**2/4+U/2)
check("Gauss scalar stress spatial", Fii, 5*v**2/12+w**2/4-U/2)
Ktrace = sum(Kdiag)
Ksqtrace = sum(x*x for x in Kdiag)
Qmixed = [Ktrace*x-x*x-(Ktrace*Ktrace-Ksqtrace)/2 for x in Kdiag]
check("Gauss extrinsic temporal", -Qmixed[0], 3*qs**2)
check("Gauss extrinsic spatial", Qmixed[1], -qs**2-2*qs*q0)
Fried = qs**2+(v**2-w**2)/12+U/6+E
Hdot = qs*(q0-qs)-v**2/3-2*E
check("spatial Einstein to Hdot", -(2*Hdot+3*Fried), Qmixed[1]+Fii+E)
check("Hdot matter form", Hdot, -k*(sig+k*rho)*(rho+p)/12-v**2/3-2*E)
check("acceleration", Fried+Hdot, qs*q0-v**2/4-w**2/12+U/6-E)
check("Gauss scalar trace", 6*(Hdot+2*Fried),
      -v**2-w**2+2*U+Ktrace**2-Ksqtrace)
check("scalar substituted Friedmann", Fried.subs(w,wn),
      (sig+k*rho)**2/36+v**2/12-(sp+k*j)**2/48+U/6+E)

# Differentiate Friedmann, imposing the signed matter exchange and junction.
vd, wd = s.symbols("v_dot normal_phi_dot")
rhod = j*v-3*H*(rho+p)
qsd = (sp*v+k*rhod)/6
qsd_normal = -w*v/3-k*H*(rho+p)/2
check("junction time derivative", qsd.subs(sp,-2*w-k*j), qsd_normal)
Edot = (4*qs*w*v-4*H*v**2-v*vd+w*wd-Up*v)/6-4*H*E
check("Weyl balance from differentiated Friedmann",
      2*H*Hdot-(2*qs*qsd_normal+(v*vd-w*wd)/6+Up*v/6+Edot))

# Radiation thresholds under the stated frozen scalar, zero Weyl assumptions.
hv, ss, kk = s.symbols("Hvac2 sigma_positive k_positive", positive=True)
R = s.symbols("rho_radiation", nonnegative=True)
Frad = kk*ss*R/18+kk**2*R**2/36
Rcrit = (s.sqrt(ss**2+36*hv)-ss)/kk
Rdec = (s.sqrt(ss**2+108*hv)-ss)/(3*kk)
check("radiation dominates Friedmann threshold", Frad.subs(R,Rcrit), hv)
check("radiation deceleration threshold", (kk*ss*R/18+kk**2*R**2/12).subs(R,Rdec), hv)
check("radiation acceleration expression", (Fried+Hdot).subs({p:rho/3,v:0,E:0}),
      (sig**2/36-w**2/12+U/6)-k*sig*rho/18-k**2*rho**2/12)

# The occupation below is the after-crossing asymptotic particle occupation.
# Weighting it by massless momentum gives a formal energy proxy, NOT the
# renormalized stress during the nonadiabatic event or a validated decay budget.
momentum, q, m0 = s.symbols("momentum q m0", positive=True)
occupation = s.exp(-s.pi*(momentum**2+m0**2)/q)
number = s.integrate(momentum**2*occupation/(2*s.pi**2), (momentum,0,s.oo))
energy0 = s.integrate(momentum**3*occupation.subs(m0,0)/(2*s.pi**2), (momentum,0,s.oo))
check("single real scalar number integral", number, q**s.Rational(3,2)*s.exp(-s.pi*m0**2/q)/(8*s.pi**3))
check("formal massless weighting of asymptotic occupation", energy0, q*q/(4*s.pi**4))

# Registered potential and the constant-scalar endpoint obstruction.
phi, delta, c = s.symbols("phi delta c")
W = 1-phi+phi**3/3
pot = s.diff(W,phi)**2/2-s.Rational(2,3)*W**2
tension = 2*W+delta*(1+c*phi)
check("AdS plus bulk stationary point", s.diff(pot,phi).subs(phi,1))
check("AdS plus bulk curvature", pot.subs(phi,1), -s.Rational(2,27))
check("AdS plus mass", s.diff(pot,phi,2).subs(phi,1), s.Rational(28,9))
check("constant phi plus junction obstruction", s.diff(tension,phi).subs(phi,1), delta*c)
eta = -s.Rational(9,64)*delta*c
check("leading radial matching", s.Rational(14,9)*eta, -2*eta-delta*c/2)
Hstatic = tension**2/36-s.diff(tension,phi)**2/48+pot/6
Hseries = s.series(Hstatic.subs(phi,1+eta),delta,0,3).removeO()
Hseries_target = delta*(1+c)/27+delta**2*((1+c)**2/36-c**2/384)
check("near plus vacuum H squared through second order", Hseries,Hseries_target)
x = s.symbols("x")
gegen = s.gegenbauer(14,2,x)
f = gegen/gegen.subs(x,1)
check("regular radial polynomial normalized", f.subs(x,1), 1)
check("regular radial Gegenbauer equation", (x*x-1)*s.diff(f,x,2)+5*x*s.diff(f,x)-252*f)
check("regular radial asymptotic exponent", s.limit(x*s.diff(f,x)/f,x,s.oo), 14)

# Exact mass dimensions (exponents). Coordinate x has dimension -1.
d = {"kappa5":s.Rational(-3,2),"kappa5_squared":-3,"phi":0,"Phi":s.Rational(3,2),
     "U":2,"sigma":1,"lambda":4,"chi":1,"gb":s.Rational(-1,2),"gbar":1,
     "mass2":2,"j":4,"rho":4,"H":1,"Weyl":2,"q":2,"fermion":s.Rational(3,2),"yukawa":0}
checks_dim = {
 "bulk canonical field":(d['phi']-d['kappa5'],d['Phi']),
 "physical brane tension":(d['sigma']-d['kappa5_squared'],d['lambda']),
 "brane mass interaction":(2*d['gb']+2*d['Phi']+2*d['chi'],4),
 "geometric coupling":(d['gb']-d['kappa5'],d['gbar']),
 "production q":(d['gb']+d['Phi']+1,d['q']),
 "matter junction dimension":(d['kappa5_squared']+d['rho'],d['sigma']),
 "scalar source dimension":(d['kappa5_squared']+d['j'],1),
 "Yukawa interaction":(d['yukawa']+d['chi']+2*d['fermion'],4),
 "energy exchange":(d['j']+d['phi']+1,d['rho']+1),
 "Gaussian production energy":(2*d['q'],d['rho'])}
for name,(lhs,rhs) in checks_dim.items():
    check("dimensions: "+name,s.sympify(lhs),s.sympify(rhs))

# Negative controls: the verifier must reject common sign/factor mistakes.
negative = {}
for name,expr in {
    "missing factor two scalar junction": -2*(-sp-k*j)/k-sp/k-j,
    "wrong matter exchange sign": -j*v-j*v,
    "constant phi=1 with registered tension": delta*c,
    "missing quadratic matter term": (sig+k*rho)**2/36-(sig**2/36+k*sig*rho/18),
}.items():
    if s.simplify(expr)==0:
        raise RuntimeError("Negative control unexpectedly vanishes: "+name)
    negative[name] = str(s.factor(expr))

# Numerical seed estimates ONLY; no shooting/BVP or endpoint simulation here.
cn = s.Rational(2)/s.Float("1.0357712571566784",50)-s.Rational(4,3)
dn = s.Rational(1,1000)
Hrs2 = dn*(1+cn)/27+dn**2*(1+cn)**2/36
curv = s.Rational(1,9)
rhostar = 1/s.sqrt(Hrs2)
xshell = s.sqrt(1+(curv*rhostar)**2)
gain = f.subs(x,xshell)
gdn = curv*s.sqrt(xshell**2-1)*s.diff(f,x).subs(x,xshell)/gain
etab = -dn*cn/(2*(gdn+2))
numerical = {"delta":"0.001","c":str(s.N(cn,30)),
    "comparison_rho_shell":str(s.N(rhostar,30)),
    "comparison_y_shell":str(s.N(s.asinh(curv*rhostar)/curv,30)),
    "linear_regular_radial_gain":str(s.N(gain,30)),
    "linear_finite_curvature_eta_shell":str(s.N(etab,30)),
    "linear_finite_curvature_eta_cone":str(s.N(etab/gain,30)),
    "H_squared_second_order":str(s.N(Hseries_target.subs({delta:dn,c:cn}),30)),
    "constant_scalar_RS_H_squared":str(s.N(Hrs2,30)),
    "status":"seed estimates, not solutions of the full static BVP; no physical length scale assigned"}
receipt = {"status":"PASS","sympy_version":s.__version__,"exact_checks":len(out),
    "checks":out,"negative_controls":negative,"dimensions":{k:str(v) for k,v in d.items()},
    "near_AdS_plus_seed_estimates":numerical,
    "scope":"conditional classical algebra and Gaussian integrals. The massless-weighted asymptotic occupation integral is only a formal energy proxy, not the instantaneous renormalized stress or a validated prompt-decay budget; no full evolution, quantum renormalization, branch existence, heating, or observation is established"}
Path(__file__).with_name("MATTER_EXTENSION_EXACT_CHECKS.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"status":"PASS","exact_checks":len(out),"negative_controls":len(negative),"near_AdS_plus_seed_estimates":numerical},indent=2))
