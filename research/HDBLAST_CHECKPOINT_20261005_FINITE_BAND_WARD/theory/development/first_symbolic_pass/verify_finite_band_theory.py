#!/usr/bin/env python3
"""Exact symbolic Ward/kernel checks. No source callback or saved array input."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

import sympy as s


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), "A fresh receipt path is required")
    checks, mutations = [], []

    def zero(name, expression):
        require(s.cancel(s.expand(expression)) == 0, name)
        checks.append(name)

    def reject(name, expression):
        require(s.cancel(s.expand(expression)) != 0, "Mutation survived: " + name)
        mutations.append(name)

    L, k, K, Pi, lag = s.symbols("L k K Pi lag", positive=True)
    X, Y, Z, T, f = s.symbols("X Y Z T f", real=True)
    ruR, ruI, rwR, rwI = s.symbols("ru_real ru_imag rw_real rw_imag", real=True)
    h = s.symbols("h0:6", real=True)
    g = 4 * L**2 * h[0] - 2 * L * h[1] - h[2]

    def D(expr):
        return s.expand(
            L**2 * s.diff(expr, L)
            + (Z + ruR) * s.diff(expr, X)
            + (T + ruI) * s.diff(expr, Y)
            + (-2 * k * T - f + rwR) * s.diff(expr, Z)
            + (2 * k * Z + rwI) * s.diff(expr, T)
            + sum(h[j + 1] * s.diff(expr, h[j]) for j in range(5))
        )

    Rm = ((2 * k**2 + 3 * L**2) * X - k * T - L * Z) / (2 * k)
    Pm = ((2 * k**2 / 3 - L**2) * X - k * T - L * Z) / (2 * k)
    G = (3 * L**2 * X - L * Z) / (2 * k)
    c = X - T / (2 * k)
    R0 = (2 * k**2 + 3 * L**2) / (4 * k)
    P0 = (2 * k**2 / 3 - L**2) / (4 * k)
    CRb = L * h[1] / (2 * k) - 2 * k * h[0] - 2 * L**2 * h[0] / k
    CPb = L * h[1] / (2 * k) - 2 * k * h[0] / 3
    Fm = L * (Rm - 3 * Pm)
    FCb = L * (CRb - 3 * CPb) - 3 * h[1] * (R0 + P0)
    projector = ((2 * k**2 + 3 * L**2) * ruR - k * rwI - L * rwR) / (2 * k)
    source_defect = L * (f - g) / (2 * k)
    zero("direct mode density equals canonical invariant plus G", Rm - k * c - G)
    zero("canonical drift includes both relevant trajectory residuals", D(c) - ruR + rwI / (2 * k))
    zero("G derivative includes its own projected residual", D(G) - Fm - L * f / (2 * k) - (3 * L**2 * ruR - L * rwR) / (2 * k))
    zero("full direct mode derivative includes invariant drift", D(Rm) - Fm - L * f / (2 * k) - projector)
    zero("bare independently defined action-contact derivative", D(CRb) - FCb + L * g / (2 * k))
    zero("complete bare off-shell Ward identity", D(Rm + CRb) - Fm - FCb - projector - source_defect)
    zero("on-shell bare Ward with consistent forcing", (D(Rm + CRb) - Fm - FCb).subs({ruR: 0, ruI: 0, rwR: 0, rwI: 0, f: g}))
    zero("direct baseline Ward", D(R0) - L * (R0 - 3 * P0))

    CR, CP, CRprime, B0 = s.symbols("CR CP CRprime B0", real=True)
    zeta = CRprime - L * (CR - 3 * CP) + 3 * h[1] * B0 + L * g / (2 * k)
    Ffull = L * (Rm + CR - 3 * (Pm + CP)) - 3 * h[1] * B0
    zero("full arbitrary contact-defect Ward identity", D(Rm) + CRprime - Ffull - projector - source_defect - zeta)

    # Conditional subtraction transfer: rederive the algebra after assuming
    # the inherited separately checked full inventory Ward relations.
    SR, SP, dSR, dSP = s.symbols("SR SP dSR dSP", real=True)
    SRprime = L * (SR - 3 * SP)
    dSRprime = L * (dSR - 3 * dSP) + h[1] * (SR - 3 * SP)
    dCR = -dSR + 4 * h[0] * SR
    dCP = -dSP + 4 * h[0] * SP
    dCRprime = -dSRprime + 4 * h[1] * SR + 4 * h[0] * SRprime
    zero("conditional complete fixed-reference subtraction transfer", dCRprime - L * (dCR - 3 * dCP) - 3 * h[1] * (SR + SP))

    Mnu, MK = s.symbols("M_nu M_K", positive=True)
    zero("mixed-measure source defect decomposition", Mnu * L * f - MK * L * g - L * Mnu * (f - g) - L * (Mnu - MK) * g)
    pressure_shift = (MK - Mnu) * g / 3
    zero("declared mixed-pressure target shifts direct integrand", -3 * L * pressure_shift - (Mnu - MK) * L * g)

    reject("omit full-R canonical drift when differencing only G", k * D(c))
    reject("omit forcing mismatch", source_defect)
    reject("reverse forcing mismatch", 2 * source_defect)
    reject("omit baseline work", 3 * h[1] * B0)
    reject("omit independently measured contact defect", zeta)
    reject("omit imaginary w residual", -rwI / 2)
    reject("omit real u residual", (2 * k**2 + 3 * L**2) * ruR / (2 * k))
    reject("omit mixed-measure source moment", (Mnu - MK) * L * g)

    phase = s.exp(2 * s.I * k * lag)
    phi = (phase - 1) / (2 * s.I * k)
    zero("exact phase ODE", s.diff(phase, lag) - 2 * s.I * k * phase)
    zero("entire Duhamel drift derivative", s.diff(phi, lag) - phase)
    zero("Duhamel drift zero lag", phi.subs(lag, 0))
    zero("Duhamel drift continuous zero momentum", s.limit(phi, k, 0) - lag)

    H0 = (1 - s.cos(2 * K * lag)) / (2 * lag)
    H1 = K * s.sin(2 * K * lag) / (2 * lag) + (s.cos(2 * K * lag) - 1) / (4 * lag**2)
    H2 = -K**2 * s.cos(2 * K * lag) / (2 * lag) + K * s.sin(2 * K * lag) / (2 * lag**2) + (s.cos(2 * K * lag) - 1) / (4 * lag**3)
    zero("finite-band sine primitive differentiated in K", s.diff(H0, K) - s.sin(2 * K * lag))
    zero("finite-band k cosine primitive differentiated in K", s.diff(H1, K) - K * s.cos(2 * K * lag))
    zero("finite-band k2 sine primitive differentiated in K", s.diff(H2, K) - K**2 * s.sin(2 * K * lag))
    for name, expression in (("H0", H0), ("H1", H1), ("H2", H2)):
        zero(name + " zero lower momentum endpoint", expression.subs(K, 0))
    zero("H0 zero-lag removable limit", s.limit(H0, lag, 0))
    zero("H1 zero-lag removable limit", s.limit(H1, lag, 0) - K**2 / 2)
    zero("H2 zero-lag removable limit", s.limit(H2, lag, 0))
    zero("sine-kernel lag derivative", s.diff(H0, lag) - 2 * H1)
    zero("k-cosine-kernel lag derivative", s.diff(H1, lag) + 2 * H2)

    impulse_subs = {X: -s.sin(2 * k * lag) / (2 * k), Z: -s.cos(2 * k * lag), T: -s.sin(2 * k * lag)}
    Rimpulse = L * s.cos(2 * k * lag) / (2 * k) - 3 * L**2 * s.sin(2 * k * lag) / (4 * k**2)
    Pimpulse = s.sin(2 * k * lag) / 3 + L**2 * s.sin(2 * k * lag) / (4 * k**2) + L * s.cos(2 * k * lag) / (2 * k)
    zero("density impulse derived directly from density operator", Rm.subs(impulse_subs) - Rimpulse)
    zero("pressure impulse derived directly from pressure operator", Pm.subs(impulse_subs) - Pimpulse)
    PsiR = L * H1 / (4 * Pi**2) - 3 * L**2 * H0 / (8 * Pi**2)
    PsiP = H2 / (6 * Pi**2) + L * H1 / (4 * Pi**2) + L**2 * H0 / (8 * Pi**2)
    zero("density band kernel fundamental theorem in K", s.diff(PsiR, K) - (k**2 * Rimpulse / (2 * Pi**2)).subs(k, K))
    zero("pressure band kernel fundamental theorem in K", s.diff(PsiP, K) - (k**2 * Pimpulse / (2 * Pi**2)).subs(k, K))
    zero("density entire band kernel zero K", PsiR.subs(K, 0))
    zero("pressure entire band kernel zero K", PsiP.subs(K, 0))
    zero("kernel Ward control with independently derived pressure", L**2 * s.diff(PsiR, L) + s.diff(PsiR, lag) - L * (PsiR - 3 * PsiP))
    zero("density zero-lag kernel reproduces source moment", s.limit(PsiR, lag, 0) - L * K**2 / (8 * Pi**2))
    reject("omit direct pressure k2 sine term", L * H2 / (2 * Pi**2))
    reject("reverse density sine term", 3 * L**2 * H0 / (4 * Pi**2))
    reject("replace canonical frequency 2k by k", s.sin(2 * K * lag) - s.sin(K * lag))

    # Explicit integrated projector bounds: derivative and lower endpoint fix
    # each primitive, with the radical term kept as one complex norm bound.
    au = (K**4 + 3 * L**2 * K**2) / (8 * Pi**2)
    aw = ((K**2 + L**2)**s.Rational(3, 2) - L**3) / (12 * Pi**2)
    ap = K**4 / (24 * Pi**2) + L**2 * K**2 / (8 * Pi**2)
    for name, expression, integrand in (
        ("density u-error primitive", au, (k + 3 * L**2 / (2 * k)) * k**2 / (2 * Pi**2)),
        ("complex w-error primitive", aw, s.sqrt(k**2 + L**2) * k / (4 * Pi**2)),
        ("pressure u-error primitive", ap, (k / 3 + L**2 / (2 * k)) * k**2 / (2 * Pi**2)),
    ):
        zero(name + " derivative", s.diff(expression, K) - integrand.subs(k, K))
        zero(name + " lower endpoint", expression.subs(K, 0))
    zero("weighted density impulse finite infrared limit", s.limit(k**2 * Rimpulse, k, 0))
    zero("weighted pressure impulse finite infrared limit", s.limit(k**2 * Pimpulse, k, 0))

    A = s.symbols("A", real=True)
    zero("wrong incoming homogeneous state still conserves direct stress", (D(Rm) - Fm).subs({X: A, Z: 0, T: 0, ruR: 0, ruI: 0, rwR: 0, rwI: 0, f: 0}))
    reject("wrong incoming state can change density arbitrarily", Rm.subs({X: A, Z: 0, T: 0}))

    # Exact Cauchy-tail and geometric constants only; no physical source
    # value, mode, arithmetic coefficient or measured response is supplied.
    tail_uniform = F(1, 15 * 2**90)
    alpha = tail_uniform / 26
    beta = alpha / 2
    require(alpha == F(1, 390 * 2**90), "Integrated inherited Cauchy-tail constant")
    require(beta == F(1, 780 * 2**90), "Symmetric time moment constant")
    checks.extend(["inherited integrated exact tail constant", "exact symmetric time moment constant"])
    Lmax = F(2, 7)
    Pi2_min = F(9)
    rows = []
    for cutoff in (64, 128, 256):
        q = F(cutoff)
        moment = q**2 / (8 * Pi2_min)
        rho_time = moment * (Lmax * alpha + 3 * Lmax**2 * beta)
        rho_band = (Lmax * moment + 3 * Lmax**2 * q / (8 * Pi2_min)) * alpha
        rho = min(rho_time, rho_band)
        pressure = Lmax * moment * alpha + min(q**3 * alpha / (18 * Pi2_min), q**4 * beta / (12 * Pi2_min)) + Lmax**2 * min(q * alpha, q**2 * beta) / (8 * Pi2_min)
        ward = Lmax * moment * alpha
        require(rho == min(q**2 * F(20, 72 * 49), q**2 / 252 + q / 294) * alpha, "Rational density formula")
        require(pressure == (q**3 / 162 + q**2 / 252 + q / 882) * alpha, "Rational pressure formula")
        require(ward == q**2 * alpha / 252, "Rational source-Ward formula")
        require(rho < F(6, 10**28), "Uniform density truncation ceiling")
        require(pressure < F(22, 10**26), "Uniform pressure truncation ceiling")
        require(ward < F(6, 10**28), "Integrated source-Ward truncation ceiling")
        checks.append("K=" + str(cutoff) + " exact rational truncation-only bounds")
        rows.append({"K": cutoff, "density_upper_rational": str(rho), "density_time_moment_upper_rational": str(rho_time), "density_band_upper_rational": str(rho_band), "pressure_upper_rational": str(pressure), "integrated_source_ward_upper_rational": str(ward)})
    require(F(14488038916154245685, 4611686018427387904) > 3, "Represented producer Pi meets common bound")
    checks.append("represented producer Pi exceeds 3")

    receipt = {
        "status": "PASS_PURE_FINITE_BAND_SYMBOLIC_THEORY",
        "identity_and_rational_check_count": len(checks),
        "checks": checks,
        "mutation_count": len(mutations),
        "mutations_rejected": mutations,
        "source_evaluations": 0,
        "saved_array_values_loaded": 0,
        "remote_writes": 0,
        "numeric_source_coefficients_loaded": 0,
        "python_optimization": sys.flags.optimize,
        "sympy_version": s.__version__,
        "python_version": sys.version.split()[0],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "derivation_sha256": hashlib.sha256(Path(__file__).with_name("FINITE_BAND_ERROR_TRANSPORT.md").read_bytes()).hexdigest(),
        "alpha_star_rational": str(alpha),
        "beta_star_rational": str(beta),
        "truncation_only_bounds": rows,
        "all_time_interval": ["-9/2", "-7/2"],
        "measure_constant": "A single fixed positive Pi; rational upper bounds use Pi>=3",
        "scope": "Exact symbolic flow/contact/kernel identities and rational conditional analytic Taylor-truncation bounds only. No source, state or numerical full-ledger evaluation.",
        "premises_not_formalized_by_this_checker": ["Cauchy analytic remainder proof", "all original full W0/W2/W4 subtraction inventory identities", "integral inequality and envelope symmetry argument", "incoming-state enclosures", "numerical residual/contact/quadrature enclosures"],
        "full_twelve_case_pressure_contact_certificate": "UNRESOLVED",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({key: receipt[key] for key in ("status", "identity_and_rational_check_count", "mutation_count", "source_evaluations", "saved_array_values_loaded")}))


if __name__ == "__main__":
    main()
