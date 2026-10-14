#!/usr/bin/env python3
"""Independent algebra checks for the finite-angular contact argument.

This checks identities, not the analytic existence theorem. No inherited
scientific checker or serialized formula is imported; no assert is used.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys

import sympy as sp


class VerificationError(RuntimeError):
    pass


def run(corrupt=None):
    checks = []

    def equal(name, expression):
        residual = sp.simplify(sp.expand(expression))
        if residual != 0:
            raise VerificationError(f"{name}: nonzero residual {residual}")
        checks.append(name)

    # Direct coordinate Christoffel contraction at a generic Lorentzian
    # tangential metric; mixed normal derivatives and tangential lapse jets
    # are independent and are not silently set to zero.
    N, Ny, a, sig = sp.symbols("N Ny a sigma", nonzero=True)
    h = sp.Matrix([[-4, 1, 0, 1], [1, 3, 1, 0],
                   [0, 1, 4, 1], [1, 0, 1, 5]])
    h0 = sp.diag(-1, 1, 1, 1)
    g = sp.zeros(5)
    g[0, 0] = N**2
    g[1:5, 1:5] = h
    gi = g.inv()
    dg = [[[sp.Integer(0) for _ in range(5)] for _ in range(5)]
          for _ in range(5)]  # coordinate, lower, lower
    dg[0][0][0] = 2*N*Ny
    for i in range(4):
        dg[0][0][i+1] = dg[0][i+1][0] = sp.Symbol(f"normal_mixed_{i}")
        dg[i+1][0][0] = sp.Symbol(f"tangential_lapse_{i}")
        for j in range(4):
            dg[0][i+1][j+1] = sig*N*h[i, j]/3
    for k in range(4):
        for i in range(4):
            for j in range(i, 4):
                value = sp.Symbol(f"tangential_h_{k}_{i}_{j}")
                dg[k+1][i+1][j+1] = dg[k+1][j+1][i+1] = value

    def christoffel(A, B, C):
        return sp.simplify(sum(gi[A, L] *
            (dg[B][C][L]+dg[C][B][L]-dg[L][B][C])/2
            for L in range(5)))

    Gamma = [[[christoffel(A, B, C) for C in range(5)]
              for B in range(5)] for A in range(5)]
    reference_y = sp.zeros(5)
    reference_y[1:5, 1:5] = -a*h0
    Cy_upper = sum(gi[B, C]*(Gamma[0][B][C]-reference_y[B, C])
                   for B in range(5) for C in range(5))
    trace = sp.trace(h.inv()*h0)
    equal("direct_relative_normal_gauge_contraction",
          Cy_upper-(Ny/N**3-2*sig/(3*N)+a*trace))
    lapse_coefficient = sp.Rational(1, 2) if corrupt == "lapse" else sp.Rational(2, 3)
    Ny_candidate = lapse_coefficient*sig*N**2-a*N**3*trace
    equal("physical_rows_plus_Cy_give_lapse_Robin",
          Cy_upper.subs(Ny, Ny_candidate))
    equal("relative_gauge_accepts_background",
          (sp.Rational(2, 3)*sig*N**2-a*N**3*sp.Symbol("trace"))
          .subs({sig:6*a, N:1, sp.Symbol("trace"):4}))

    # Direct second differentiation of the three exact boundary maps. The
    # noncommuting tangential perturbation matrices exercise inverse-metric
    # trace terms; scalar/lapse jets and tension derivatives remain symbolic.
    eps = sp.Symbol("eps")
    nu1, nu2, chi, scalar2, sig1, sig2 = sp.symbols("nu1 nu2 chi scalar2 sig1 sig2")
    H1 = sp.Matrix([[2, 1, 0, -1], [1, 3, 2, 0],
                    [0, 2, 4, 1], [-1, 0, 1, 2]])
    K2 = sp.Matrix([[1, 2, -1, 0], [2, -1, 0, 1],
                    [-1, 0, 2, 3], [0, 1, 3, -2]])
    he = h0+eps*H1+eps**2*K2/2
    Ne = 1+eps*nu1+eps**2*nu2/2
    scalar_increment = eps*chi+eps**2*scalar2/2
    tension = 6*a+sig1*scalar_increment+sig2*scalar_increment**2/2+4*scalar_increment**3/6
    tension_prime = sig1+sig2*scalar_increment+2*scalar_increment**2
    metric_exact = tension*Ne*he/3
    metric_R = (2*sig1*chi/3+4*a*nu1)*H1+(sig2*chi**2/3+2*sig1*chi*nu1/3)*h0
    metric_candidate = 2*a*K2+(sig1*scalar2/3+2*a*nu2)*h0+metric_R
    for i in range(4):
        for j in range(i, 4):
            equal(f"second_boundary_metric_source_{i}_{j}",
                  sp.diff(metric_exact[i,j],eps,2).subs(eps,0)-metric_candidate[i,j])
    scalar_R = -2*chi**2-sig2*chi*nu1
    equal("second_boundary_scalar_source",
          sp.diff(-tension_prime*Ne/2,eps,2).subs(eps,0)
          -(-sig2*scalar2/2-sig1*nu2/2+scalar_R))
    A1 = h0.inv()*H1
    quadratic_lapse = 12 if corrupt == "second-source" else 16
    lapse_R = (2*sig2*chi**2/3+8*sig1*chi*nu1/3
        -quadratic_lapse*a*nu1**2-2*a*sp.trace(A1*A1)+6*a*nu1*sp.trace(A1))
    lapse_exact = 2*tension*Ne**2/3-a*Ne**3*sp.trace(he.inv()*h0)
    equal("second_boundary_relative_lapse_source",
          sp.diff(lapse_exact,eps,2).subs(eps,0)
          -(-4*a*nu2+a*sp.trace(h0.inv()*K2)+2*sig1*scalar2/3+lapse_R))
    equal("second_metric_lapse_conversion",
          sp.diff(Ne**2,eps,2).subs(eps,0)-(2*nu2+2*nu1**2))

    # Compute the covariant mixed subsidiary expression explicitly, with
    # C_y=0 and its tangential derivatives zero as boundary traces.
    Ct = sp.symbols("C0:4")
    Cty = sp.symbols("Cy0:4")
    factor = 1 if corrupt == "subsidiary" else 2
    for i in range(4):
        covariant_sum = Cty[i] - sum(
            (Gamma[j+1][0][i+1]+Gamma[j+1][i+1][0])*Ct[j]
            for j in range(4))
        equal(f"subsidiary_coordinate_Robin_{i}",
              covariant_sum-(Cty[i]-factor*N*sig*Ct[i]/6))

    # Independent first longitudinal gauge defect: the covariant divergence
    # is evaluated directly in a tangential normal frame, with warp jets
    # retained. Angular and time derivatives of u are independent.
    eta = sp.diag(1, -1, 1, 1, 1)
    u, uy = sp.symbols("u uy")
    ua = sp.symbols("u0:4")
    perturb = sp.diag(-4*u, -2*u, 2*u, 2*u, 2*u)
    conn = [[[sp.Integer(0) for _ in range(5)] for _ in range(5)]
            for _ in range(5)]
    for i in range(1, 5):
        conn[0][i][i] = -a*eta[i, i]
        conn[i][0][i] = conn[i][i][0] = a
    dh = [[[sp.Integer(0) for _ in range(5)] for _ in range(5)]
          for _ in range(5)]
    dh[0][0][0] = -4*uy
    for i in range(1, 5):
        dh[0][i][i] = (2*uy+4*a*u)*eta[i, i]
    for k in range(1, 5):
        dh[k][0][0] = -4*ua[k-1]
        for i in range(1, 5):
            dh[k][i][i] = 2*ua[k-1]*eta[i, i]

    def nabla(k, i, j):
        return dh[k][i][j] - sum(conn[L][k][i]*perturb[L, j]
                + conn[L][k][j]*perturb[i, L] for L in range(5))

    first_C = [sum(eta[B, B]*nabla(B, A, B)
                 - eta[B, B]*nabla(A, B, B)/2 for B in range(5))
               for A in range(5)]
    expected_y = -6*(uy+(3 if corrupt == "first-gauge" else 4)*a*u)
    equal("longitudinal_first_normal_gauge_defect", first_C[0]-expected_y)
    for i in range(1, 5):
        equal(f"longitudinal_first_tangential_gauge_defect_{i}", first_C[i])

    # Actual mixed Lie derivative for a boundary-tangent covector.
    xiy = sp.symbols("xiy0:4")
    xia = sp.symbols("xia0:4")
    for i in range(1, 5):
        lie_mixed = xiy[i-1] - sum((conn[j][0][i]+conn[j][i][0])*xia[j-1]
                                   for j in range(1, 5))
        equal(f"boundary_preserving_first_gauge_{i}",
              lie_mixed-(xiy[i-1]-2*a*xia[i-1]))

    # Temporal gauge-jet Jacobian from the defining DC formula, not copied
    # from the proposed solved expressions. Signature here is time-first.
    etat = sp.diag(-1, 1, 1, 1, 1)
    jets = {}
    for d in range(5):
        for i in range(5):
            for j in range(i, 5):
                jets[d, i, j] = sp.Symbol(f"dk_{d}_{i}_{j}")
    def jh(d, i, j):
        return jets[d, min(i, j), max(i, j)]
    DC = sp.Matrix([sum(etat[B, B]*jh(B, A, B)
                    - etat[B, B]*jh(A, B, B)/2 for B in range(5))
                    for A in range(5)])
    unknown = sp.Matrix([jh(0, 0, i) for i in range(5)])
    Jacobian = DC.jacobian(unknown)
    equal("initial_lapse_shift_gauge_Jacobian",
          Jacobian.det()+sp.Rational(1, 2))
    for i in range(5):
        equal(f"initial_gauge_Jacobian_diagonal_{i}",
              Jacobian[i, i]+(sp.Rational(1, 2) if i == 0 else 1))
    Cdot = sp.symbols("Cdot0:5")
    E = sp.zeros(5)
    E[0, 0] = Cdot[0]
    for i in range(1, 5):
        E[0, i] = E[i, 0] = Cdot[i]/2
    traceE = sp.trace(etat*E)
    H = E-etat*traceE/2
    initial_factor = 1 if corrupt == "initial-constraint" else sp.Rational(1, 2)
    for i in range(5):
        equal(f"physical_initial_constraint_gauge_derivative_{i}",
              H[0, i]-initial_factor*Cdot[i])

    # Character check for the actual complete symmetric square. Laurent
    # variables parametrize the two maximal tori in SU(2)xSU(2).
    z, w = sp.symbols("z w", nonzero=True)
    def char(j, v):
        return sum(v**(2*m) for m in range(-j, j+1))
    vchar = char(1, z)*char(1, w)
    symmetric_square = (vchar**2+vchar.subs({z:z**2, w:w**2}, simultaneous=True))/2
    decomposition = char(2, z)*char(2, w)+char(2, z)+char(2, w)+vchar+1
    if corrupt == "angular":
        decomposition -= 1
    equal("SO4_complete_symmetric_square_character", symmetric_square-decomposition)
    equal("SO4_complete_symmetric_square_dimension", symmetric_square.subs({z:1,w:1})-45)

    # Solve a genuinely nonsymmetric two-component Robin trace ODE by its
    # fundamental matrix. Check both the ODE and the incoming reconstruction.
    t, s = sp.symbols("t s", real=True)
    M = sp.Matrix([[1, 3], [0, 1]])
    fundamental = sp.exp(t)*sp.Matrix([[1, 3*t], [0, 1]])
    for i in range(2):
        for j in range(2):
            equal(f"nonsymmetric_Robin_fundamental_{i}_{j}",
                  (sp.diff(fundamental, t)-M*fundamental)[i, j])
    actual_z = sp.Matrix([t**2+1, 2*t-1])
    outgoing = sp.Matrix([t, 1-t**2])
    boundary_source = sp.diff(actual_z, t)-M*actual_z-outgoing
    integrand = fundamental.inv()*(outgoing+boundary_source)
    integral = sp.Matrix([sp.integrate(v.subs(t, s), (s, 0, t)) for v in integrand])
    recovered_z = sp.simplify(fundamental*(actual_z.subs(t, 0)+integral))
    reflection_factor = 1 if corrupt == "Robin" else 2
    incoming = outgoing+reflection_factor*(M*recovered_z+boundary_source)
    for i in range(2):
        equal(f"nonsymmetric_Robin_trace_ODE_solution_{i}", recovered_z[i]-actual_z[i])
        equal(f"nonsymmetric_Robin_incoming_trace_{i}",
              (incoming[i]+outgoing[i])/2-sp.diff(actual_z[i], t))
        equal(f"nonsymmetric_Robin_normal_trace_{i}",
              (incoming[i]-outgoing[i])/2-(M*recovered_z+boundary_source)[i])

    # Candidate principal normal map: 4 mixed metric values and independent
    # normal derivatives of 10 tangent metric, normal lapse, and scalar.
    principal = sp.diag(*([1]*4+[sp.Rational(1, 2)]*10+[1, 1]))
    equal("sixteen_row_principal_boundary_map_determinant",
          principal.det()-sp.Rational(1, 1024))
    return checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output")
    parser.add_argument("--inject-corruption", choices=["lapse", "subsidiary",
        "first-gauge", "initial-constraint", "angular", "Robin", "second-source"])
    args = parser.parse_args()
    checks = run(args.inject_corruption)
    result = {
        "status": "pass", "check_count": len(checks), "checks": checks,
        "python": platform.python_version(), "sympy": sp.__version__,
        "optimized": bool(sys.flags.optimize),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "Independent exact algebra identities and manufactured Robin ODE; not an existence proof.",
        "corruption": args.inject_corruption,
    }
    rendered = json.dumps(result, indent=2)+"\n"
    if args.output:
        Path(args.output).write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
