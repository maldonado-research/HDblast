#!/usr/bin/env python3
"""Exact algebra checks; not a PDE existence or numerical-response certificate."""
import argparse
import json
from pathlib import Path
import sympy as s


def run(corrupt):
    checks = []

    def identity(name, lhs, rhs):
        residual = s.simplify(s.expand(lhs - rhs))
        checks.append({"name": name, "passed": residual == 0,
                       "residual": str(residual)})

    rho, a, p, py, B, sigma2 = s.symbols("rho a p py B sigma2", nonzero=True)
    V, Vt, Vtt, Vz, Vzz, Py, Pyy = s.symbols("V Vt Vtt Vz Vzz Py Pyy")
    t = s.symbols("t", real=True)
    lam = s.Integer(8)
    A = 3*s.tanh(t)
    L4P = -Vtt-A*Vt+(4-lam/s.cosh(t)**2)*V
    y1 = -Vz/rho
    y2 = Vzz/rho**2+a*Vz/rho
    raw_master = y2-(2*a+2*py/p)*y1+L4P/rho**2-s.Rational(2, 3)*p**2*V
    C = rho*(3*a+2*py/p)
    claimed_C = -C if corrupt == "reverse_spatial_drift" else C
    D = lam/s.cosh(t)**2-4+s.Rational(2, 3)*rho**2*p**2
    claimed_master = Vtt-Vzz+A*Vt-claimed_C*Vz+D*V
    identity("registered_master_in_inward_depth", -rho**2*raw_master, claimed_master)
    raw_boundary = B*y1-L4P/rho**2
    claimed_boundary = Vtt+A*Vt+(lam/s.cosh(t)**2-4)*V-rho*B*Vz
    identity("full_even_junction_to_dynamic_boundary", rho**2*raw_boundary, claimed_boundary)
    claimed_scalar = (3 if corrupt != "reverse_observation" else -3)*Vz/(rho**3*p)
    identity("physical_scalar_observation", -3*y1/(rho**2*p), claimed_scalar)

    Psi = V/rho**2
    Psiy = Py/rho**2-2*a*Psi
    chi = -3*Py/(rho**2*p)
    chiy = -3*(Pyy-(2*a+py/p)*Py)/(rho**2*p)
    identity("all_even_tensor_junction_diagonal_rows", Psiy+2*a*Psi+p*chi/3, 0)
    lapse_term = 0 if corrupt == "drop_normal_lapse" else 2*p*Psi
    scalar_jump = 2*(chiy+lapse_term+sigma2*chi/2)
    onshell_Pyy = (2*a+2*py/p)*Py-L4P/rho**2+s.Rational(2, 3)*p**2*V
    claimed_jump = -6*((py/p+sigma2/2)*Py-L4P/rho**2)/(rho**2*p)
    identity("complete_scalar_junction_including_normal_variation",
             scalar_jump.subs(Pyy, onshell_Pyy), claimed_jump)

    c, cp = s.symbols("c cp")
    w, wt, wtt, wz, wzz = s.symbols("w wt wtt wz wzz")
    lt = -A/2
    lz = -c/2
    transformed = (wtt+2*lt*wt+(s.diff(lt,t)+lt**2)*w
                   -wzz-2*lz*wz-(-cp/2+lz**2)*w
                   +A*(wt+lt*w)-c*(wz+lz*w)+D*w)
    potential = s.Rational(35,4)/s.cosh(t)**2-s.Rational(25,4)+cp/2+c*c/4+s.Rational(2,3)*rho*rho*p*p
    identity("time_and_radial_conjugation_potential", transformed,
             wtt-wzz+potential*w)
    mu = s.symbols("mu", nonzero=True)
    boundary_conjugated = (mu*(wtt+2*lt*wt+(s.diff(lt,t)+lt**2)*w
                             +A*(wt+lt*w)+(lam/s.cosh(t)**2-4)*w)
                          -(wz+lz*w))
    beta = mu*(s.Rational(35,4)/s.cosh(t)**2-s.Rational(25,4))+c/2
    identity("positive_mass_boundary_ODE_after_conjugation",
             boundary_conjugated, -wz+mu*wtt+beta*w)

    xi, zeta, ch = s.symbols("xi zeta ch")
    embedding_shift = zeta if corrupt == "drop_embedding_shift" else zeta-xi
    identity("physical_pullback_gauge_cancellation", ch+p*xi+p*embedding_shift, ch+p*zeta)

    r, q, qvel = s.symbols("r q qvel", positive=True)
    identity("ambient_to_deSitter_zero_slice_displacement", r**4*(r**2*q), r**6*q)
    identity("ambient_to_deSitter_zero_slice_velocity", r**4*r*(r**2*qvel), r**7*qvel)

    x = s.symbols("x0:4")
    H = x[0]**2-x[1]**2
    identity("quadrupole_is_harmonic_homogeneous_polynomial", sum(s.diff(H, xx, 2) for xx in x), 0)
    identity("quadrupole_homogeneity", sum(xx*s.diff(H,xx) for xx in x), 2*H)
    identity("sphere_eigenvalue", s.Integer(2)*(2+4-2), lam)
    area = 2*s.pi**2
    E_n1four = s.Rational(3, 4*6)
    E_n1sq_n2sq = s.Rational(1, 4*6)
    normY = area*(2*E_n1four-2*E_n1sq_n2sq)
    identity("sphere_quadrupole_L2_norm", normY, s.pi**2/3)

    eps, s1, s2, m1, m2, m0, clock1 = s.symbols("eps s1 s2 m1 m2 m0 clock1")
    residual = eps*s1+eps**2*s2/2
    measure = m0*(1+eps*m1+eps**2*m2/2)
    qsecond = s.diff(residual**2*measure, eps, 2).subs(eps, 0)
    claimed_second = 2*m0*s1**2
    if corrupt == "insert_measure_variation":
        claimed_second += m0*m1*s1**2
    identity("exact_observable_second_variation", qsecond, claimed_second)
    phi_b = s.symbols("phi_b")
    identity("intrinsic_clock_shift_of_constant_background", clock1*s.diff(phi_b,t), 0)
    scalar_amp = 3*Vz/(rho**3*p)
    observable = 2*rho**3*s.cosh(t)**3*normY*scalar_amp**2
    expected = 6*s.pi**2*s.cosh(t)**3*Vz**2/(rho**3*p**2)
    if corrupt == "double_count_common_shell":
        expected *= 2
    identity("one_common_shell_observable_normalization", observable, expected)
    return {"status": "PASS" if all(c["passed"] for c in checks) else "FAIL",
            "corruption": corrupt, "checks": checks,
            "scope": "Exact symbolic algebra only; analytic PDE/ODE arguments remain in the note."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--corrupt", default="none", choices=[
        "none", "reverse_spatial_drift", "reverse_observation", "drop_normal_lapse",
        "drop_embedding_shift", "insert_measure_variation", "double_count_common_shell"])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.corrupt)
    encoded = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(encoded)
    print(json.dumps({"status": result["status"], "corruption": args.corrupt,
                      "checks": len(result["checks"]),
                      "failures": [c["name"] for c in result["checks"] if not c["passed"]]}))
    raise SystemExit(0 if result["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
