#!/usr/bin/env python3
"""Exact independent checks; no numerical sampling or main-draft mutation."""
import hashlib
import json
from pathlib import Path
import platform

import sympy as s


checks = {}


def zero(name, expr):
    residual = s.factor(s.cancel(expr))
    checks[name] = {"passed": residual == 0, "residual": str(residual)}


def positive(name, expr):
    expr = s.factor(expr)
    checks[name] = {"passed": expr.is_positive is True, "expression": str(expr)}


# Formal nonlinear radial vector field. H2 denotes rho**(-2).
a, p, h, q, v, b, H2, U1, U2 = s.symbols("a p h q v b H2 U1 U2")
flow = {
    a: -H2 - p**2 / 3,
    p: U1 - 4*a*p,
    h: q,
    q: U2*h - 4*a*q - 4*p*b,
    v: b,
    b: 2*H2*v - s.Rational(2, 3)*p*q,
    H2: -2*a*H2,
    U1: U2*p,
}


def radial(expr):
    return s.expand(sum(s.diff(expr, variable)*rate for variable, rate in flow.items()))


constraint = 2*a*b + 2*H2*v - (p*q - U1*h)/6
j = b + p*h/3
D = h*flow[p] - p*q  # Wronskian; d below is reserved for a0-k.
zero("j_identity_mod_linearized_constraint", radial(j) + 4*a*j + 2*H2*v - 2*constraint)
zero("auxiliary_constraint_identity", p*q-flow[p]*h-12*a*j-12*H2*v+6*constraint)
zero("wronskian_D_identity", radial(D)+4*a*D-4*p*(H2*h+p*j))

# Differentiating a'=-a**2-p**2/4-U/6 gives this damping equation.
bprime_damped = -2*a*b-p*q/2-U1*h/6
zero("metric_damping_mod_constraint", flow[b]-bprime_damped-constraint)
for epsilon in (-1, 1):
    zmetric = -epsilon*b
    abs_p, abs_U1 = epsilon*p, epsilon*U1
    zero(f"metric_positive_system_sign_{epsilon}",
         -epsilon*bprime_damped+2*a*zmetric-abs_p*q/2-abs_U1*h/6)
    zero(f"scalar_positive_system_sign_{epsilon}",
         flow[q]+4*a*q-U2*h-4*abs_p*zmetric)

beta_T, D_junction = s.symbols("beta_T D_junction", nonzero=True)
Q = j/h
A = flow[p]-p*q/h
R_T = flow[a]+b*(beta_T-p)/h+p*beta_T/3  # sigma'=-2p
zero("endpoint_residual_identity", R_T+H2-Q*(beta_T-p))
zero("fixed_beta_scalar_identity", A-D/h)
zero("fixed_beta_metric_identity", flow[a]-p*b/h+H2+p*Q)
zero("junction_residual_decomposition",
     R_T.subs(beta_T, -A/D_junction)+H2+p*Q+Q*A/D_junction)

# Center coefficients are solved from the scalar and metric equations.
y = s.symbols("y", positive=True)
U0, U3, r3, e2 = s.symbols("U0 U3 r3 e2")
rho_trial, eta_minus_c = y+r3*y**3, e2*y**2
p_trial = s.diff(eta_minus_c, y)
a_trial = s.diff(rho_trial, y)/rho_trial
scalar_constant = s.limit(s.diff(p_trial, y)+4*a_trial*p_trial-U1, y, 0)
metric_linear = s.limit((s.diff(rho_trial, y, 2)+rho_trial*(p_trial**2/4+U0/6))/y, y, 0)
center_solution = s.solve([scalar_constant, metric_linear], [e2, r3], dict=True)[0]
zero("center_eta_quadratic", center_solution[e2]-U1/10)
zero("center_rho_cubic", center_solution[r3]+U0/36)
rho = rho_trial.subs(center_solution)
eta_minus_c = eta_minus_c.subs(center_solution)
p_center = s.diff(eta_minus_c, y)
a_center = s.diff(rho, y)/rho


def central(expr):
    return s.diff(expr, U0)*U1+s.diff(expr, U1)*U2+s.diff(expr, U2)*U3


h_center = 1+central(eta_minus_c)
q_center = s.diff(h_center, y)
v_center = central(rho)/rho
b_center = central(a_center)
j_center = b_center+p_center*h_center/3
D_center = h_center*s.diff(p_center, y)-p_center*q_center
zero("center_h_quadratic", s.limit((h_center-1)/y**2, y, 0)-U2/10)
zero("center_q_linear", s.limit(q_center/y, y, 0)-U2/5)
zero("center_v_quadratic", s.limit(v_center/y**2, y, 0)+U1/36)
zero("center_b_linear", s.limit(b_center/y, y, 0)+U1/18)
zero("center_j_linear", s.limit(j_center/y, y, 0)-U1/90)
zero("center_D_constant", s.limit(D_center, y, 0)-U1/5)
zero("center_rho4_j_boundary", s.limit(rho**4*j_center, y, 0))
zero("center_rho4_D_boundary", s.limit(rho**4*D_center, y, 0))

# Linear Riccati variables; parametrization w=lambda+4k enforces w>4k.
k, lam, d = s.symbols("k lambda d", positive=True)
z, E = s.symbols("z E", real=True)
w = lam+4*k
m2 = lam*w
dprime = -d*(d+2*k)
zprime = m2-4*(k+d)*z-z**2
e = lam-z
zero("linear_e_identity", -zprime+(w+z)*e-4*d*z)
Eprime = (-zprime*d-e*dprime)/d**2
Eflow = 4*lam-2*(w-3*k)*E+d*E*(E-3)
zero("linear_E_identity", Eprime.subs(z, lam-d*E)-Eflow)
C = 2*lam/(w-3*k)
barrier = s.factor(Eflow.subs(E, C))
zero("barrier_vector_field", barrier-d*C*(C-3))
positive("barrier_C_positive", C)
positive("barrier_2_minus_C_positive", 2-C)
positive("barrier_inward_margin_positive", -barrier)
CH = lam/(k*(w-3*k))
H02 = d*(d+2*k)
zero("linear_curvature_identity", H02-((k+d)**2-k**2))
positive("linear_H2_bound_margin", CH*H02-C*d)

F_center = 1+m2*y**2/10
z_center = s.diff(F_center, y)/F_center
d_center = k*s.coth(k*y)-k
E_center = (lam-z_center)/d_center
zero("linear_center_z_linear", s.limit(z_center/y, y, 0)-m2/5)
zero("linear_center_E_zero", s.limit(E_center, y, 0))
zero("linear_center_E_linear", s.limit(E_center/y, y, 0)-lam)

# Contraction inequalities give M3*B*b <= gap/2; check its strict margin.
nu = lam/2
mu_nu = nu*(nu+4*k)
gap = m2-mu_nu
positive("green_gap_positive", gap)
positive("convexity_lower_bound_positive", m2-gap/2)
zero("convexity_lower_bound_identity", m2-gap/2-(m2+mu_nu)/2)

source = Path(__file__).resolve()
report = {
    "all_passed": all(check["passed"] for check in checks.values()),
    "check_count": len(checks),
    "python_version": platform.python_version(),
    "sympy_version": s.__version__,
    "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "assumptions": [
        "Nonlinear j identity uses the linearized constraint; Wronskian identity does not.",
        "Center checks solve and differentiate the displayed leading regular series.",
        "Riccati checks use k>0, lambda>0, d>0, w=lambda+4k.",
        "First-crossing barrier argument is analytic; symbolic checks verify its exact identities and strict margins.",
        "Transfer from H0^2 to H^2 additionally uses the analytic conclusion 0<R=rho/r0<=1.",
        "Convexity additionally uses CG*LN<=1/2, CG>=1/gap, and LN>=M3*B*b.",
    ],
    "derived_constants": {
        "barrier_C": str(C),
        "barrier_vector_field": str(barrier),
        "linear_curvature_constant": str(CH),
        "convexity_lower_bound": str(s.factor(m2-gap/2)),
    },
    "checks": checks,
}
destination = source.with_suffix(".json")
destination.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
print(json.dumps({"all_passed": report["all_passed"], "check_count": len(checks), "report": str(destination)}))
raise SystemExit(0 if report["all_passed"] else 1)
