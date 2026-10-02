#!/usr/bin/env python3
"""Independent exact variational/boundary checks; no physical evolution.

Derive FRW local sources from a lapse-retaining minisuperspace action, rather
than defining sources by the desired Ward identity. Then check the inherited
doubled-bulk / single-shell junction algebra and adverse source mutations.
All failures raise explicitly, including under python -O.
"""
import hashlib
import json
from pathlib import Path

import sympy as s

HERE = Path(__file__).resolve().parent
checks = []
controls = []


def reduce(expr):
    return s.factor(s.cancel(s.together(s.expand(expr))))


def check(name, expression, note=""):
    residual = reduce(expression)
    if residual != 0:
        raise AssertionError(f"{name}: {residual}")
    checks.append({"name": name, "passed": True, "residual": "0", "note": note})


def control(name, expression, note=""):
    residual = reduce(expression)
    if residual == 0:
        raise AssertionError(f"Control escaped detection: {name}")
    controls.append({"name": name, "rejected": True, "residual": str(residual), "note": note})


t = s.symbols("tau", real=True)
a = s.Function("a", positive=True)(t)
N = s.Function("N", positive=True)(t)
phi = s.Function("phi", real=True)(t)
V = s.Function("V")(phi)
F = s.Function("F")(phi)
alpha, kap2 = s.symbols("alpha kappa5_squared", nonzero=True)


def D(expr, order=1):
    return s.diff(expr, t, order)


def euler(lagrangian, field, highest):
    return sum((-1)**q * D(s.diff(lagrangian, D(field, q)), q)
               for q in range(highest + 1))


R_lapse = 6 * (D(a, 2)/(a*N**2) + D(a)**2/(a**2*N**2)
               - D(a)*D(N)/(a*N**3))
lagrangian = N*a**3 * (-V + F*R_lapse + alpha*R_lapse**2)
gauge_fix = {N: 1, **{D(N, q): 0 for q in range(1, 6)}}
rho_variation = (-euler(lagrangian, N, 1)/a**3).subs(gauge_fix)
p_variation = (euler(lagrangian, a, 2)/(3*N*a**2)).subs(gauge_fix)
J_variation = (-euler(lagrangian, phi, 0)/(N*a**3)).subs(gauge_fix)

H = D(a)/a
v = D(phi)
R = 6*(D(H) + 2*H**2)
rho = V - 6*F*H**2 - 6*H*D(F) + alpha*(36*D(H)**2 - 216*H**2*D(H) - 72*H*D(H, 2))
p = -V + 2*F*(2*D(H) + 3*H**2) + 2*D(F, 2) + 4*H*D(F) + alpha*(108*D(H)**2 + 216*H**2*D(H) + 144*H*D(H, 2) + 24*D(H, 3))
J = s.diff(V, phi) - s.diff(F, phi)*R
check("lapse variation independently yields local energy density", rho_variation-rho)
check("scale-factor variation independently yields local pressure", p_variation-p)
check("scalar variation independently yields V_phi-F_phi R", J_variation-J)
check("variational local Ward identity", D(rho)+3*H*(rho+p)-J*v)

# Covariant R^2 tensor's FRW reduction provides another expression.
R00 = -3*(D(H)+H**2)
Rii = D(H)+3*H**2
E00 = 2*R*R00 + R**2/2 + 6*H*D(R)
Eii = 2*R*Rii - R**2/2 - 2*D(R, 2) - 4*H*D(R)
check("covariant R_squared 00 variation agrees with lapse variation", -2*alpha*E00 - alpha*(36*D(H)**2-216*H**2*D(H)-72*H*D(H, 2)))
check("covariant R_squared spatial variation agrees with scale variation", -2*alpha*Eii - alpha*(108*D(H)**2+216*H**2*D(H)+144*H*D(H, 2)+24*D(H, 3)))

sigma = s.Function("sigma")(phi)
rt = s.Function("rho_total")(t)
pt = s.Function("p_total")(t)
jt = s.Function("J_total")(t)
ks = (sigma + kap2*rt)/6
k0 = (sigma - kap2*(2*rt+3*pt))/6
w = -(s.diff(sigma, phi) + kap2*jt)/2
ward = D(rt)+3*H*(rt+pt)-jt*v
momentum = -3*D(ks)-3*H*(ks-k0)-v*w
check("shell momentum residual equals minus kappa_squared Ward residual /2", momentum + kap2*ward/2,
      "Coordinate conformal-gauge momentum residual has an extra e^(2 B_b).")
check("doubled bulk signed flux budget", D(sigma/kap2+rt)+3*H*(rt+pt)+2*w*v/kap2-ward)
check("radiation temporal junction reduces to inherited A1 path", k0.subs(pt, rt/3)-(sigma-3*kap2*rt)/6)
check("vacuum junctions recover identical metric normal derivatives", (ks-k0).subs({rt:0, pt:0}))

# Quantum chain rule; x and Q here are arbitrary state/geometry-dependent data.
x = s.Function("x")(phi)
Q = s.Function("Q")(t)
check("mass-current chain rule J_phi=x_phi Q/2", D(x)*Q/2-s.diff(x, phi)*Q*v/2)
G, phistar, m0 = s.symbols("G phi_star m0", real=True)
mass = m0**2 + G**2*(phi-phistar)**2
check("quadratic mass current", s.diff(mass, phi)*Q/2-G**2*(phi-phistar)*Q)

# Negative controls act on fully derived sources, with no desired identity used
# to manufacture the mutated pressure/current.
control("reverse scalar-source sign", D(rho)+3*H*(rho+p)+J*v)
control("omit F_phi R current but retain F metric variation", D(rho)+3*H*(rho+p)-s.diff(V,phi)*v)
control("omit F_ddot spatial variation", D(rho)+3*H*(rho+p-2*D(F, 2))-J*v)
control("omit R_squared H_triple_dot pressure term", D(rho)+3*H*(rho+p-24*alpha*D(H, 3))-J*v)
control("double local metric terms while leaving current single", D(2*rho)+3*H*(2*rho+2*p)-J*v)
control("drop scalar matter source from junction", momentum.subs(jt, 0) + kap2*ward/2,
        "For a physical Ward source jt, removing it leaves -kappa_squared jt v/2.")

# Detect double-counting by original-action junction coefficients, because
# doubling stress AND current consistently still passes Ward.
check("paired double counting can evade the Ward test", D(2*rho)+3*H*(2*rho+2*p)-2*J*v,
      "PASS is intentional: coefficient bookkeeping is required to detect model change.")
control("source duplication changes spatial junction by kappa_squared rho_local/6", (sigma+kap2*(rt+2*rho))/6-(sigma+kap2*(rt+rho))/6)
control("source duplication changes temporal junction", (sigma-kap2*(2*(rt+2*rho)+3*(pt+2*p)))/6-(sigma-kap2*(2*(rt+rho)+3*(pt+p)))/6)

# Conservation alone cannot detect a pressure error at a static instant.
hjet, delta_p = s.symbols("H_at_instant delta_pressure", real=True)
check("static instant hides a pressure-only error in Ward", (3*hjet*delta_p).subs(hjet,0),
      "Use direct spatial junction and a separately evaluated trace check at H=0.")
control("temporal junction detects static pressure-only error", -kap2*delta_p/2)

# Determine differential order without relying on terminology.
check("energy contains scale-factor third derivative", s.diff(rho, D(a,3)) + 72*alpha*D(a)/a**2)
check("pressure contains scale-factor fourth derivative", s.diff(p, D(a,4)) - 24*alpha/a)
check("F pressure contains scalar acceleration", s.diff(p, D(phi,2))-2*s.diff(F,phi))
check("F scalar current contains scale-factor acceleration", s.diff(J, D(a,2))+6*s.diff(F,phi)/a)

# The positive-reference recipe itself requires high jets even when the
# additional matched finite alpha is set to zero. These are direct integrand
# coefficients, not a physical-field heavy-mass approximation.
eta = s.symbols("eta", real=True)
ae = s.Function("a_conformal", positive=True)(eta)
xe = s.Function("x_conformal", real=True)(eta)
ke, r = s.symbols("k reference_mass_squared", positive=True)
de = lambda z, q=1: s.diff(z, eta, q)
wr = s.sqrt(ke**2+ae**2*r)
Le = de(ae)/ae
De = ae**2*(xe-r)
u2 = (De-de(ae,2)/ae)/(2*wr)-de(wr,2)/(4*wr**2)+3*de(wr)**2/(8*wr**3)
u4 = -u2**2/(2*wr)-de(u2,2)/(4*wr**2)+de(wr,2)*u2/(4*wr**3)+3*de(wr)*de(u2)/(4*wr**3)-3*de(wr)**2*u2/(4*wr**4)
be = -ke**2/3-ae**2*r
ce = 1-be/wr**2
J4 = Le*de(u2)/wr**2-2*Le*de(wr)*u2/wr**3+de(wr)*de(u2)/(2*wr**3)-3*de(wr)**2*u2/(4*wr**4)
p4 = (ce*u4+be*u2**2/wr**3-(Le**2-De)*u2/wr**2+J4)/(4*ae**4)
Q2 = -u2/(2*ae**2*wr**2)
check("reference pressure subtraction has fourth conformal scale derivative", s.diff(p4,de(ae,4))-ce*(1/(32*ae**5*wr**3)+r/(64*ae**3*wr**5)))
check("reference pressure subtraction has second conformal mass derivative", s.diff(p4,de(xe,2))+ce/(32*ae**2*wr**3))
check("reference Q subtraction has second conformal scale derivative", s.diff(Q2,de(ae,2))-(1/(4*ae**3*wr**3)+r/(8*ae*wr**5)))

out = {"scope":"exact local variation and boundary compatibility only; no physical evolution or certified bulk closure",
       "checks_passed":len(checks), "negative_controls_rejected":len(controls),
       "checks":checks, "negative_controls":controls,
       "sympy_version":s.__version__,
       "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/"INDEPENDENT_JUNCTION_CHECKS.json").write_text(json.dumps(out, indent=2)+"\n")
print(f"{len(checks)} exact checks passed; {len(controls)} negative controls rejected")
