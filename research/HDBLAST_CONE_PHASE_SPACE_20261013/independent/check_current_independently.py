"""Exact audit from a nonlinear metric-family presymplectic potential.

This checks algebraic steps, not the global physical phase-space hypotheses.
The companion note derives the potential from its connection contraction.
No filesystem input besides this script and SymPy; explicit failures survive -O.
"""
import argparse
import json
from pathlib import Path

import sympy as s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    passed, rejected = [], []

    def equal(name, left, right=0):
        residual = s.factor(s.simplify(left - right))
        if residual != 0:
            raise RuntimeError(f"{name}: {residual}")
        passed.append(name)

    def unequal(name, wrong, right):
        residual = s.factor(s.simplify(wrong - right))
        if residual == 0:
            raise RuntimeError(f"Control failed to distinguish {name}")
        rejected.append({"name": name, "nonzero_residual": str(residual)})

    e1, e2 = s.symbols("e1 e2")
    rho, a, p, py = s.symbols("rho a p py", nonzero=True)
    q1, q2, q1y, q2y, x1, x2, x1y, x2y = s.symbols(
        "Psi1 Psi2 Psi1_y Psi2_y X1 X2 X1_y X2_y"
    )
    q1m, q2m, x1m, x2m = s.symbols("Psi1_m Psi2_m X1_m X2_m")
    psi = e1*q1 + e2*q2
    psi_y = e1*q1y + e2*q2y
    scalar_y = p + e1*x1y + e2*x2y
    psi_m = e1*q1m + e2*q2m
    scalar_m = e1*x1m + e2*x2m
    # sqrt(-gamma) is factored out throughout. n=-2 psi exactly on this
    # commuting nonlinear family; antisymmetric second variations cancel.
    def theta_y(dq, dqy, dx):
        return -rho**4*s.exp(6*psi)*(4*dqy + 12*(a+psi_y)*dq + scalar_y*dx)

    def theta_m(dq, dqm, dx):
        return -rho**2*(dqm + 6*dq*psi_m + scalar_m*dx)

    def curl(t2, t1):
        return s.expand((s.diff(t2,e1)-s.diff(t1,e2)).subs({e1:0,e2:0}))

    jy = curl(theta_y(q2,q2y,x2), theta_y(q1,q1y,x1))
    jm = curl(theta_m(q2,q2m,x2), theta_m(q1,q1m,x1))
    wq = q1*q2y-q2*q1y
    wx = x1*x2y-x2*x1y
    cross = q1*x2-q2*x1
    raw = rho**4*(wx-12*wq-6*p*cross)
    equal("nonlinear_potential_radial_curl", jy, raw)
    equal("nonlinear_potential_tangent_curl", jm,
          rho**2*(x1*x2m-x2*x1m+6*(q1*q2m-q2*q1m)))
    constraint = {q1y:-2*a*q1-p*x1/3, q2y:-2*a*q2-p*x2/3}
    current = rho**4*(wx-2*p*cross)
    equal("radial_current_after_momentum_constraint", jy.subs(constraint), current)
    unequal("drop_density_scalar_cross_variation", rho**4*(wx-12*wq), raw)
    unequal("drop_nonlinear_metric_connection_variation", rho**4*(wx-24*wq-6*p*cross), raw)
    unequal("omit_constraint_cross_term", rho**4*wx, current)

    l1,l2 = s.symbols("Lambda1 Lambda2")
    equations = {
        q1:-2*a*q1-p*x1/3, q2:-2*a*q2-p*x2/3,
        x1:py*x1/p+(3*l1/(rho**2*p)-2*p)*q1,
        x2:py*x2/p+(3*l2/(rho**2*p)-2*p)*q2,
        rho:a*rho, p:py,
    }
    reduced = s.factor(current.subs({x1y:equations[x1],x2y:equations[x2]}))
    equal("on_shell_radial_current_first_order", reduced,
          3*rho**2/p*(l2*q2*x1-l1*q1*x2))
    radial_div = sum(s.diff(reduced,var)*deriv for var,deriv in equations.items())
    tangential_div = rho**2*(l2-l1)*(x1*x2+6*q1*q2)
    equal("mixed_eigenvalue_bulk_conservation", radial_div+tangential_div)

    P1,P2,P1y,P2y = s.symbols("P1 P2 P1_y P2_y")
    toP = {q1:P1/rho**2,q2:P2/rho**2,x1:-3*P1y/(rho**2*p),x2:-3*P2y/(rho**2*p)}
    alpha = 1/(rho**2*p**2)
    master = 9*alpha*(l1*P1*P2y-l2*P2*P1y)
    equal("mixed_eigenvalue_master_current", reduced.subs(toP), master)
    equal("equal_eigenvalue_green_current", master.subs(l2,l1),
          9*l1*alpha*(P1*P2y-P2*P1y))
    unequal("use_single_eigenvalue_for_packet", 9*l1*alpha*(P1*P2y-P2*P1y),master)
    beta, co = s.symbols("beta c_o", nonzero=True)
    equal("even_shell_flux_zero", master.subs({P1y:l1*P1/(alpha*beta),P2y:l2*P2/(alpha*beta)}))
    odd = (2*master).subs({P1y:-co*P1/alpha,P2y:-co*P2/alpha})
    equal("odd_shell_current_balance", odd,18*co*(l2-l1)*P1*P2)
    unequal("odd_shell_balance_wrong_orientation", odd,-18*co*(l2-l1)*P1*P2)

    u,v = s.symbols("u v", positive=True)
    du,dv = s.symbols("du dv")
    A = s.symbols("A")
    # Leading flat-cone metric: dy² + y² ds²_dS = du dv + (u+v)² dΩ³/4.
    dy = (v*du+u*dv)/(2*s.sqrt(u*v))
    equal("smooth_leading_metric_null_covectors", -6*u*v*A*dy**2,
          -s.Rational(3,2)*A*(v*du+u*dv)**2)
    # Radial and tangential densities near the cone, per unit solid angle,
    # for X_i=F_i(v), Psi=O(uv), with Weyl-gauge gravitational term subleading.
    W = s.symbols("F1_F2v_minus_F2_F1v")
    jy0 = (u+v)**3*v*W/8
    jt0 = -(u+v)**3*v*W/(8*s.sqrt(u*v))
    Ju = (s.sqrt(u/v)*jy0-u*jt0)/(2*s.sqrt(u*v))
    equal("characteristic_density_pullback", Ju,(u+v)**3*W/8)
    equal("characteristic_density_limit", s.limit(Ju,u,0,dir="+"),v**3*W/8)
    unequal("omit_tangential_current_in_pullback",s.limit(s.sqrt(u/v)*jy0/(2*s.sqrt(u*v)),u,0,dir="+"),v**3*W/8)

    args.output.parent.mkdir(parents=True,exist_ok=True)
    result = {
        "status":"PASS_EXACT_LOCAL_CURRENT_AND_NULL_LIMIT_ALGEBRA",
        "scope":"Algebra from the independently derived nonlinear-family Lee-Wald potential; conditional cone expansions; no global solution existence, gauge-charge or completeness proof.",
        "positive_checks":passed,"wrong_formula_controls":rejected,
        "positive_check_count":len(passed),"wrong_formula_control_count":len(rejected),
        "sympy":s.__version__,
    }
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:result[k] for k in ("status","positive_check_count","wrong_formula_control_count")},sort_keys=True))


if __name__ == "__main__":
    main()
