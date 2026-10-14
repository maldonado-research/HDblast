#!/usr/bin/env python3
"""Independent exact symbolic checks, not a PDE existence certificate."""
import argparse
import json
from pathlib import Path
import sympy as S


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--inject-failure", action="store_true")
    args = parser.parse_args()
    t, chi, x = S.symbols("tau chi x", real=True)
    b, c = S.symbols("b c", positive=True)
    d, lam = S.symbols("d lambda", real=True)
    checks = []

    def simp(expr):
        return S.trigsimp(S.simplify(S.expand_trig(expr)))

    def check(name, expr):
        residual = simp(expr)
        ok = residual == 0
        checks.append({"name": name, "passed": bool(ok), "residual": str(residual)})

    def reject(name, expr):
        residual = simp(expr)
        # A regular-point witness proves nonidentity; an expression merely
        # left unsimplified would not be an adequate negative control.
        witness = residual.xreplace({jet: S.Integer(1) for jet in residual.atoms(S.Derivative)})
        witness = witness.xreplace({fun: S.Integer(1) for fun in witness.atoms(S.Function)
                                    if fun.func not in (S.sin, S.cos, S.tan, S.cot, S.sec, S.csc)})
        witness = simp(witness.subs({t:S.pi/6, chi:S.pi/4, x:S.pi/6, b:1, c:1, d:1, lam:1}))
        ok = witness.is_zero is False
        checks.append({"name": name, "passed": bool(ok), "wrong_residual": str(residual),
                       "nonzero_witness": str(witness)})

    q = b / (S.cos(t) + S.cos(chi))
    TT, RR = q * S.sin(t), q * S.sin(chi)
    check("metric_tau_tau", -S.diff(TT, t)**2 + S.diff(RR, t)**2 + q**2)
    check("metric_chi_chi", -S.diff(TT, chi)**2 + S.diff(RR, chi)**2 - q**2)
    check("metric_cross", -S.diff(TT, t)*S.diff(TT, chi) + S.diff(RR, t)*S.diff(RR, chi))
    check("angular_metric", RR**2 - q**2*S.sin(chi)**2)
    ss = b**2*(S.cos(t)-S.cos(chi))/(S.cos(t)+S.cos(chi))
    check("actual_s", RR**2-TT**2-ss)
    check("s_at_center", ss.subs(chi, 0)+b**2*S.sin(t)**2/(1+S.cos(t))**2)
    check("s_at_shell", ss.subs(chi, S.pi/2)-b**2)

    # Direct conformal scalar wave conjugation in 1+4 dimensions.
    lt, lc = S.diff(q, t)/q, S.diff(q, chi)/q
    k = S.Rational(3, 2)
    zeroth = k*(3-k)*(lt**2-lc**2) + k*(S.diff(lt,t)-S.diff(lc,chi)-3*S.cot(chi)*lc)
    check("bulk_conformal_potential_9_over_4", zeroth+S.Rational(9,4))
    check("bulk_time_first_derivative_cancels", (2*k-3)*lt)
    check("bulk_radial_first_derivative_cancels", (3-2*k)*lc)
    reject("negative_control_bulk_curvature_2", zeroth+2)

    qb = b/S.cos(t)
    kb = S.diff(qb,t)/qb
    Qb = S.Rational(3,2)*S.diff(kb,t)+S.Rational(3,4)*kb**2
    mu = c/qb
    beta = (d-S.Rational(3,2)/b)*qb-mu*Qb
    check("normal_log_weight", lc.subs(chi,S.pi/2)-qb/b)
    check("boundary_mass_derivative", S.diff(mu,t)+mu*S.tan(t))
    z = S.Function("z")(t)
    weight = qb**(-k)
    # Boundary dimension 4: -dtt - 2(log q)' dt - lambda.
    boundary_direct = (-S.diff(weight*z,t,2)-2*kb*S.diff(weight*z,t)-lam*weight*z)/weight
    boundary_expected = -S.diff(z,t,2)+kb*S.diff(z,t)-lam*z+Qb*z
    check("boundary_wave_conjugation", boundary_direct-boundary_expected)
    check("boundary_beta_conversion", (d*qb-k*qb/b-mu*Qb)-beta)
    reject("negative_control_omit_boundary_friction", boundary_direct-(-S.diff(z,t,2)-lam*z+Qb*z))
    reject("negative_control_flip_boundary_mass", mu+c/qb)

    v = S.Function("v")(x)
    radial_weight = S.cos(x)**(-k)
    radial_direct = (-S.diff(radial_weight*v,x,2)+3*S.tan(x)*S.diff(radial_weight*v,x)
                     +(lam/S.cos(x)**2+S.Rational(9,4))*radial_weight*v)/radial_weight
    radial_expected = -S.diff(v,x,2)+(lam+S.Rational(3,4))/S.cos(x)**2*v
    check("radial_mode_conjugation", radial_direct-radial_expected)
    reject("negative_control_omit_radial_three_quarters", radial_direct-(-S.diff(v,x,2)+lam/S.cos(x)**2*v))

    p, w, zp, zpp = S.symbols("p w zprime zsecond")
    vx = p-w
    check("even_boundary_ode", (-vx+mu*zpp+S.diff(mu,t)*w+(mu*lam+beta)*z)
          -(mu*zpp+(S.diff(mu,t)+1)*w+(mu*lam+beta)*z-p))
    alpha = qb*(d-k/b)
    check("odd_boundary_ode", (-vx+alpha*z) - (w+alpha*z-p))

    a = S.symbols("a", real=True)
    # Rational tangent representation avoids branch-dependent simplification.
    y = S.symbols("tan_half_a", positive=True)
    R0 = b*(1-y)/(1+y)
    contact = 2*b*y/(1-y**2)
    check("ambient_contact_time", (b**2-R0**2)/(2*R0)-contact)
    check("outgoing_null_contact", (R0+contact)**2-(b**2+contact**2))

    ambient_t = S.symbols("T", real=True)
    A = S.sqrt(b**2+ambient_t**2)
    ur, ut = S.symbols("u_r u_T")
    h = ambient_t/A**2
    direct_normal = (ambient_t*(ut-h*ur)+ur)/b
    check("fixed_ball_normal", direct_normal-(ambient_t*ut/b+b*ur/A**2))
    check("fixed_ball_radial_coercivity", A**2-ambient_t**2-b**2)

    if args.inject_failure:
        check("injected_wrong_bulk_curvature", zeroth+S.Rational(5,2))

    report = {
        "kind": "exact_symbolic_identity_checks_not_PDE_certificate",
        "sympy_version": S.__version__,
        "passed": all(row["passed"] for row in checks),
        "checks": checks,
    }
    payload = json.dumps(report, indent=2)
    if args.output:
        args.output.write_text(payload+"\n")
    print(payload)
    raise SystemExit(0 if report["passed"] else 1)


if __name__ == "__main__":
    main()
