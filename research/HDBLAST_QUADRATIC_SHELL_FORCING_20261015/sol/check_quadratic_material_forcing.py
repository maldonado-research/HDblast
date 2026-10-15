#!/usr/bin/env python3
"""Independent exact checks and rational consequences of retained H cells.

No archived scientific program is imported or executed. SymPy supplies exact
polynomial algebra; fractions supplies all enclosure and decision arithmetic.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from fractions import Fraction as F

import sympy as sp


HERE = Path(__file__).resolve().parent
MANIFEST_HASH = "9b42295d819f0ed5caf94d0a46af38274a27fb11648dd6bff91f5d992385f140"
PUBLIC_MAIN = "6d4070dd7fc85dfc61c3e89484b083077bb3b472"
PACKAGED_PATHS = {
    "waveform_note": "prior/INERTIAL_SHELL_WAVEFORM.md",
    "waveform_report": "prior/INERTIAL_WAVEFORM_ENCLOSURE.json",
    "waveform_cells": "prior/INERTIAL_WAVEFORM_CELLS.csv",
    "shell_box": "prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json",
    "boundary_note": "reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md",
}
ERRORS = (
    "sum_vs_copy", "factorial_half", "wrong_sign", "sphere_dimension",
    "unit_norm", "wrong_fourth_moment", "omit_degree_two",
    "wrong_harmonic_constant", "wrong_norm_jacobian", "cancel_cosh_in_quartic",
    "tension_remainder_factor", "response_sign_from_source",
    "computed_cell", "computed_J4", "computed_norm", "computed_tail",
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def same(actual, expected, message):
    require(sp.simplify(actual - expected) == 0, message)


def hash_bytes(data):
    return hashlib.sha256(data).hexdigest()


def interval(lo, hi):
    lo, hi = F(lo), F(hi)
    require(lo <= hi, "reversed interval")
    return (lo, hi)


def add(a, b):
    return interval(a[0] + b[0], a[1] + b[1])


def scale(a, c):
    c = F(c)
    return interval(min(a[0] * c, a[1] * c), max(a[0] * c, a[1] * c))


def multiply(a, b):
    p = [x * y for x in a for y in b]
    return interval(min(p), max(p))


def square(a):
    return interval(0 if a[0] <= 0 <= a[1] else min(x*x for x in a),
                    max(x*x for x in a))


def cosh_minus_three(a, b):
    """Exact enclosure over 0<=a<=t<=b, using a positive-series tail.

    Terms through t^6/6! are explicit. The remaining series starts with
    t^8/8! and its later ratios are <= b^2/(10*9).
    """
    a, b = F(a), F(b)
    require(0 <= a <= b and b*b < 90, "invalid cosh interval")
    p = lambda t: sum((t**(2*k) / math.factorial(2*k) for k in range(4)), F(0))
    lower_cosh = p(a)
    upper_cosh = p(b) + b**8 / math.factorial(8) / (1 - b*b / 90)
    return interval(1 / upper_cosh**3, 1 / lower_cosh**3)


def angular_checks(error):
    x = sp.symbols("x1:5", real=True)
    s = sum(v*v for v in x)
    A, B = x[0]**2 + x[1]**2, x[0]**2 - x[1]**2
    Z = A - s/2
    H = B**2 - s*A/2 + s*s*(sp.Rational(1, 10) if error == "wrong_harmonic_constant" else sp.Rational(1, 12))
    lap = lambda p: sp.expand(sum(sp.diff(p, v, 2) for v in x))
    for p in (B, Z, H):
        same(lap(p), 0, "ambient harmonicity")
    z_coefficient = 0 if error == "omit_degree_two" else sp.Rational(1, 2)
    same(B**2, H + s*z_coefficient*Z + s*s/6, "complete harmonic decomposition")
    dimension = 3 if error == "sphere_dimension" else 4
    for p, degree in ((B, 2), (Z, 2), (H, 4)):
        hom = sum(v*sp.diff(p, v) for v in x)
        same(hom, degree*p, "homogeneity")
        # Sphere Laplacian of a homogeneous polynomial restricted to r=1.
        spherical_lap = s*lap(p) - degree*(degree+dimension-2)*p
        same(spherical_lap, (-8 if degree == 2 else -24)*p,
             "unit-S3 eigenvalue")

    def mean(poly):
        """Normalized S3 moments from Dirichlet(1/2,1/2,1/2,1/2)."""
        result = sp.S.Zero
        for powers, coefficient in sp.Poly(sp.expand(poly), *x).terms():
            if any(k % 2 for k in powers):
                continue
            half_powers = [k//2 for k in powers]
            degree = sum(half_powers)
            moment = sp.prod(sp.rf(sp.Rational(1, 2), k) for k in half_powers) / sp.rf(2, degree)
            result += coefficient*moment
        return sp.simplify(result)

    # Homogeneous extensions restrict to the intended unit-sphere functions.
    m2, m4, z2, h2 = mean(B**2), mean(B**4), mean(Z**2), mean(H**2)
    same(mean(s), 1, "sphere radius")
    same(m2, 1 if error == "unit_norm" else sp.Rational(1, 6), "Y is not unit norm")
    same(m4, sp.Rational(1, 16) if error == "wrong_fourth_moment" else sp.Rational(3, 40),
         "fourth moment")
    same(z2, sp.Rational(1, 12), "degree-two norm")
    same(h2, sp.Rational(19, 720), "degree-four norm")
    for p, q in ((1, Z), (1, H), (Z, H)):
        same(mean(p*q), 0, "orthogonal harmonic projections")
    same(mean(B**2*Z)/z2, sp.Rational(1, 2), "degree-two projection")
    same(mean(B**2*H)/h2, 1, "degree-four projection")
    powers = [sp.Rational(1, 36), z2/4, h2]
    shares = [sp.simplify(v/m4) for v in powers]
    require(shares == [sp.Rational(10, 27), sp.Rational(5, 18), sp.Rational(19, 54)], "angular power shares")
    same(sum(shares), 1, "complete angular power")
    return {
        "sphere": "unit S3 in R4; area 2*pi^2",
        "mean_Y2": str(m2), "mean_Y4": str(m4),
        "mean_Z2_squared": str(z2), "mean_H4_squared": str(h2),
        "squared_L2_norms_divided_by_pi_squared": {"Y": "1/3", "Z2": "1/6", "H4": "19/360"},
        "orthogonal_squared_power_shares_l0_l2_l4": list(map(str, shares)),
        "Y2_projection_coefficients_l0_l2_l4": ["1/6", "1/2", "1"],
    }


def constitutive_checks(error):
    eps, psi0, u1, u2 = sp.symbols("eps psi0 psi1 psi2", real=True)
    eta, delta, c = sp.symbols("eta delta c", real=True)
    sigma = sp.Rational(2, 3) + 2*eta**2 + sp.Rational(2, 3)*eta**3 + delta*(1+c+c*eta)
    same(sp.diff(sigma, eta, 3), 4, "registered third tension derivative")
    sigma1, sigma2 = sp.diff(sigma, eta), sp.diff(sigma, eta, 2)
    require(sp.diff(sigma, eta, 4) == 0, "cubic constitutive law")
    path = psi0 + eps*u1 + eps**2*u2/2
    scalar_second = sp.diff(sigma1.subs(eta, path), eps, 2).subs(eps, 0)
    rhs = (-2 if error == "sum_vs_copy" else -4)*u1**2
    if error == "factorial_half":
        rhs /= 2
    if error == "wrong_sign":
        rhs *= -1
    same(scalar_second - sigma2.subs(eta, psi0)*u2, -rhs,
         "outward SUM source and derivative convention")
    same(rhs/2, -2*u1**2, "identical-copy specialization")
    tension_second = sp.diff(sigma.subs(eta, path), eps, 2).subs(eps, 0)
    same(tension_second, sigma1.subs(eta, psi0)*u2+sigma2.subs(eta, psi0)*u1**2,
         "pure constitutive second tension coefficient")
    h0, h1, h2 = sp.symbols("h0 h1 h2", real=True)
    metric_path = h0+eps*h1+eps**2*h2/2
    stress_second = sp.diff(-sigma.subs(eta,path)*metric_path, eps, 2).subs(eps,0)
    same(stress_second, -sigma.subs(eta,psi0)*h2-sigma1.subs(eta,psi0)*u2*h0
         -2*sigma1.subs(eta,psi0)*u1*h1-sigma2.subs(eta,psi0)*u1**2*h0,
         "complete material stress product rule")
    q = sp.symbols("q", real=True)
    cubic_factor = sp.Rational(1, 3) if error == "tension_remainder_factor" else sp.Rational(2, 3)
    remainder = sigma.subs(eta, psi0+q) - sigma.subs(eta, psi0) - sigma1.subs(eta, psi0)*q - sigma2.subs(eta, psi0)*q*q/2
    same(remainder, cubic_factor*q**3, "exact constitutive cubic remainder")
    same(sigma1.subs(eta, psi0+q), sigma1.subs(eta, psi0)+sigma2.subs(eta, psi0)*q+2*q*q,
         "exact derivative-of-tension polynomial")
    # An explicit simple counterexample to inferring response sign from source
    # alone; this is not proposed as the Einstein-shell response operator.
    t = sp.symbols("t", real=True)
    positive_response, negative_response = 1-2*t*t, -1-2*t*t
    same(sp.diff(positive_response,t,2),-4,"positive response same negative source")
    same(sp.diff(negative_response,t,2),-4,"negative response same negative source")
    require(positive_response.subs(t,0)>0 and negative_response.subs(t,0)<0,
            "source-only response-sign counterexample")
    if error == "response_sign_from_source":
        require(positive_response.subs(t,0)<=0,"negative source need not give a negative response")
    return {"sigma_third": "4", "sigma_second": "4+4*eta_b",
            "sum_scalar_RHS": "-4*psi1^2", "per_copy_scalar_RHS": "-2*psi1^2",
            "scalar_RHS_projection_coefficients_times_S_squared": ["-2/3", "-2", "-4"],
            "pure_tension_second_coefficient": "sigma_second*psi1^2",
            "exact_cubic_constitutive_remainder": "(2/3)*q^3",
            "nonlinear_solution_remainder_claimed": False,
            "response_sign_deduced": False}


def norm_checks(error):
    rho, pb, ell, G, K, C, J2, J4, sig2 = sp.symbols("rho pb ell G K C J2 J4 sig2", nonzero=True)
    G_expr, K_expr = 3/(rho**3*pb), 3*sp.pi**2/(rho**2*pb**2)
    same(sp.pi**2*rho**4*G_expr**2/3, K_expr, "proper-time quadratic prefactor")
    jacobian = ell if error == "wrong_norm_jacobian" else 1/ell
    same(K*jacobian*J2, K*J2/ell, "quadratic observable Jacobian")
    same(rho**4*G_expr**4*(3*sp.pi**2/20), sp.Rational(9, 20)*K_expr*G_expr**2,
         "quartic norm angular and background prefactors")
    cosh_power = 0 if error == "cancel_cosh_in_quartic" else -3
    require(cosh_power == 3 - 4*sp.Rational(3, 2), "quartic cosh weight survives")
    return {"background_measure": "rho_b^4*cosh(t)^3*dt*dOmega",
            "C": "K/ell*integral(H^2 dx)",
            "J4": "integral(cosh(ell*x)^(-3)*H(x)^4 dx)",
            "Qphi_L1": "4*C", "Qphi_L2_squared": "(36/5)*K*G^2*ell^(-3)*J4",
            "pure_tension_L1": "sigma_second*C",
            "pure_tension_L2_squared": "(9/20)*sigma_second^2*K*G^2*ell^(-3)*J4",
            "scalar_Israel_constitutive_amplitude_L2_squared": "pure_tension_L2_squared/9",
            "positive_tensor_norm_prescribed": False}


def authenticated_inputs(args):
    pins_path = HERE / "SOURCE_PINS.json"
    pins_data = pins_path.read_bytes()
    pins = json.loads(pins_data)
    require(pins["external_manifest_sha256"] == MANIFEST_HASH, "external manifest pin")
    require(pins["prior_public_main"] == PUBLIC_MAIN, "prior public main pin")
    input_pins_data = (args.input_root / "INPUT_PINS.json").read_bytes()
    input_pins = json.loads(input_pins_data)
    require(input_pins["external_inertial_manifest_sha256"] == MANIFEST_HASH,
            "packaged external-manifest identity")
    input_rows = {r["path"]: r for r in input_pins["files"]}
    locations = {k: args.input_root/v for k,v in PACKAGED_PATHS.items()}
    contents = {}
    authenticated = {}
    for key, path in locations.items():
        data = path.read_bytes()
        expected = pins["sources"][key]
        require(len(data) == expected["bytes"] and hash_bytes(data) == expected["sha256"],
                "source-byte authentication failed: " + key)
        row = input_rows["inputs/"+PACKAGED_PATHS[key]]
        require(row["bytes"] == len(data) and row["sha256"] == hash_bytes(data),
                "packaged source does not match input pins: " + key)
        if key != "boundary_note":
            require(row["original_public_source"]["authenticated_manifest_sha256"] == MANIFEST_HASH,
                    "source manifest identity: " + key)
        contents[key] = data
        authenticated[key] = {"bytes": len(data), "sha256": hash_bytes(data)}
    authenticated["input_pins"] = {"bytes": len(input_pins_data), "sha256": hash_bytes(input_pins_data)}
    return contents, authenticated, hash_bytes(pins_data)


def retained_data_checks(contents):
    report = json.loads(contents["waveform_report"])
    shell = json.loads(contents["shell_box"])
    require(report["status"] == "PASS", "retained report status")
    require(report["state"]["ell"] == "1/1000000", "wider tangent width")
    ell = F(report["state"]["ell"])
    gain = interval(*report["physical_gain_interval"])
    prefactor = interval(*report["physical_norm_prefactor_interval"])
    sig2 = interval(*shell["boxes"]["sigma_second_b"]["exact"])
    eta = interval(*shell["boxes"]["eta_b"]["exact"])
    require(sig2 == add(interval(4, 4), scale(eta, 4)) and sig2[0] > 0,
            "retained sigma_second correlation and sign")
    require(gain[1] < 0 and prefactor[0] > 0, "gain and norm signs")
    htail = interval(*report["tail_theorem"]["actual_H_interval"])
    require(htail[0] > 0 and report["tail_theorem"]["domain"] == [2, 3], "continuous retained tail")
    cosh_tail = cosh_minus_three(2*ell, 3*ell)
    S2tail = scale(multiply(multiply(square(gain), square(htail)), cosh_tail), 1/ell**2)
    phi_mean_tail = scale(S2tail, -F(2, 3))
    phi_coefficients_tail = {"l0": phi_mean_tail, "l2": scale(S2tail, -2), "l4": scale(S2tail, -4)}
    pure_tension_mean_tail = scale(multiply(sig2, S2tail), F(1, 6))

    rows = list(csv.DictReader(io.StringIO(contents["waveform_cells"].decode("ascii"))))
    D = int(report["cell_data"]["value_denominator"])
    grid = int(report["cell_data"]["time_grid_denominator"])
    require(D == 2**112 and grid == 4096 and len(rows) == 4353, "retained cell format")
    require(hash_bytes(contents["waveform_cells"]) == report["cell_data"]["sha256"], "report cell hash")
    sums2, sums4, sums4weighted = [0, 0], [0, 0], [0, 0]
    weight_denominator = 2**160
    end = 0
    Habs = F(0)
    tail_cells = []
    for row in rows:
        left, right = int(row["left_grid"]), int(row["right_grid"])
        require(left == end and right > left, "gap, overlap, or reversed time cell")
        end = right
        lo, hi = int(row["H_lo"]), int(row["H_hi"])
        require(lo <= hi, "reversed H cell")
        lower2 = 0 if lo <= 0 <= hi else min(lo*lo, hi*hi)
        upper2 = max(lo*lo, hi*hi)
        width = right - left
        sums2[0] += width*lower2
        sums2[1] += width*upper2
        sums4[0] += width*lower2*lower2
        sums4[1] += width*upper2*upper2
        clock = cosh_minus_three(ell*F(left, grid), ell*F(right, grid))
        clock_int = [clock[0].numerator*weight_denominator//clock[0].denominator,
                     -((-clock[1].numerator*weight_denominator)//clock[1].denominator)]
        require(F(clock_int[0], weight_denominator) >= (1-(ell*F(right,grid))**2)**2
                and F(clock_int[1], weight_denominator) <= 1, "sharper clock interval")
        sums4weighted[0] += width*lower2*lower2*clock_int[0]
        sums4weighted[1] += width*upper2*upper2*clock_int[1]
        Habs = max(Habs, F(max(abs(lo), abs(hi)), D))
        if left >= 2*grid:
            require(lo > 0, "tail cell is not positive")
            tail_cells.append(interval(F(lo,D), F(hi,D)))
    require(end == 3*grid and len(tail_cells) == 256, "domain and tail coverage")
    I2 = interval(F(sums2[0], grid*D**2), F(sums2[1], grid*D**2))
    I4 = interval(F(sums4[0], grid*D**4), F(sums4[1], grid*D**4))
    weighted_I4 = interval(*(F(v, grid*D**4*weight_denominator) for v in sums4weighted))
    C_fresh = scale(multiply(prefactor, I2), 1/ell)
    C_retained = interval(*report["quadrature"]["C_interval"])
    require(C_retained[0] <= C_fresh[0] <= C_fresh[1] <= C_retained[1], "independent exact H2 consequence")
    qphi2 = scale(multiply(multiply(prefactor, square(gain)), weighted_I4), F(36,5)/ell**3)
    qsigma2 = scale(multiply(square(sig2), qphi2), F(1,16))
    Hgainmax = max(abs(g) for g in gain)*Habs/ell
    return {
        "ell": str(ell), "G": gain, "K": prefactor, "sigma_second": sig2,
        "retained_continuous_H_tail": htail,
        "cosh_minus_three_tail": cosh_tail,
        "continuous_tail_S_squared": S2tail,
        "continuous_tail_Qphi_projection_coefficients": phi_coefficients_tail,
        "continuous_tail_Qphi_mean": phi_mean_tail,
        "continuous_tail_pure_tension_mean": pure_tension_mean_tail,
        "continuous_tail_pointwise_Qphi": scale(S2tail, -4)[0:1] + (F(0),),
        "pointwise_tail_sign": "Qphi<=0, strict where Y!=0; mean strictly negative",
        "angular_zero_set": "Y=0 iff n1^2=n2^2",
        "H2_integral_fresh": I2, "H4_integral_fresh": I4,
        "weighted_H4_integral_fresh": weighted_I4,
        "C_fresh": C_fresh, "C_retained": C_retained,
        "Qphi_L1": scale(C_fresh, 4), "Qphi_L2_squared": qphi2,
        "pure_tension_L1": multiply(sig2, C_fresh), "pure_tension_L2_squared": qsigma2,
        "H_abs_sup_upper": Habs, "psi1_abs_sup_upper": Hgainmax,
        "cubic_remainder_relative_to_half_epsilon2_pure_tension_upper_per_abs_epsilon":
            F(4,3)*Hgainmax/sig2[0],
        "tail_cell_count": len(tail_cells), "full_cell_count": len(rows),
        "sharpness_scope": "Outward interval images of retained H,G and sigma_second boxes; no new background correlation or extremum claim.",
    }


def check_computed_files(paths, contents, independently_derived, error):
    """Recompute new output coefficients and rational integrals from old data.

    The producer is not imported. Byte equality checks every new CSV entry;
    numeric report fields are compared to separately accumulated quantities.
    """
    report_path, cells_path = paths
    report_raw, cell_raw = report_path.read_bytes(), cells_path.read_bytes()
    report = json.loads(report_raw)
    prior = json.loads(contents["waveform_report"])
    sigma2 = independently_derived["sigma_second"]
    G, K = independently_derived["G"], independently_derived["K"]
    ell, D, grid = F(1, 10**6), 2**112, 4096
    source_rows = list(csv.DictReader(io.StringIO(contents["waveform_cells"].decode("ascii"))))
    header = ["left_grid", "right_grid"]
    for kind in ("Q", "T"):
        for degree in (0, 2, 4):
            header.extend([f"ell2_{kind}{degree}_lo", f"ell2_{kind}{degree}_hi"])
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(header)
    j2_acc, j4_acc = [F(0), F(0)], [F(0), F(0)]
    for index, row in enumerate(source_rows):
        ln, rn = int(row["left_grid"]), int(row["right_grid"])
        h = interval(F(int(row["H_lo"]),D), F(int(row["H_hi"]),D))
        e = interval(F(int(row["ell_S_lo"]),D), F(int(row["ell_S_hi"]),D))
        h2, e2 = square(h), square(e)
        h4 = square(h2)
        width = F(rn-ln,grid)
        clock_lo = (1-(ell*F(rn,grid))**2)**2
        j4_acc[0] += width*clock_lo*h4[0]
        j4_acc[1] += width*h4[1]
        # The inherited H2 certificate rounds the square outward to 112 bits
        # before summation; the fresh unrounded certificate is a subset.
        low_fixed = (h2[0]*D).numerator//(h2[0]*D).denominator
        high_fixed = -((-(h2[1]*D).numerator)//(h2[1]*D).denominator)
        j2_acc[0] += width*F(low_fixed,D)
        j2_acc[1] += width*F(high_fixed,D)
        tension_squared = multiply(sigma2,e2)
        factors = [scale(e2, -F(2,3)), scale(e2,-2), scale(e2,-4),
                   scale(tension_squared,F(1,6)), scale(tension_squared,F(1,2)), tension_squared]
        entries = [ln, rn]
        for lo,hi in factors:
            entries.extend([(lo*D).numerator//(lo*D).denominator,
                            -((-(hi*D).numerator)//(hi*D).denominator)])
        if error == "computed_cell" and index == 100:
            entries[2] += 1
        writer.writerow(entries)
    expected_csv = buffer.getvalue().encode("ascii")
    require(expected_csv == cell_raw, "new coefficient CSV differs from independent cell reconstruction")
    j4 = interval(*j4_acc)
    if error == "computed_J4":
        j4 = scale(j4, 2)
    j2 = interval(*j2_acc)
    C = scale(multiply(K,j2),1/ell)
    qphi2 = scale(multiply(multiply(K,square(G)),j4),F(36,5)/ell**3)
    if error == "computed_norm":
        qphi2 = scale(qphi2,F(1,2))
    tension_norm = scale(multiply(square(sigma2),qphi2),F(1,16))
    Htail = interval(*prior["tail_theorem"]["actual_H_interval"])
    E_tail = interval(G[0]*Htail[1], G[1]*(1-(3*ell)**2)*Htail[0])
    tail = scale(square(E_tail),-F(2,3))
    if error == "computed_tail":
        tail = scale(tail,-1)
    expected_numeric = {
        "quadrature.H_squared_integral_interval": j2,
        "quadrature.cosh_minus3_H_fourth_integral_interval": j4,
        "quadrature.inherited_C_interval": C,
        "quadrature.scalar_signed_integral_interval": scale(C,-4),
        "quadrature.scalar_L1_interval": scale(C,4),
        "quadrature.scalar_L2_squared_interval": qphi2,
        "quadrature.constitutive_tension_integral_interval": multiply(sigma2,C),
        "quadrature.constitutive_tension_L2_squared_interval": tension_norm,
        "tail_theorem.ell2_scalar_monopole_interval": tail,
        "parameters.sigma_second_interval": sigma2,
        "parameters.physical_gain_interval": G,
        "parameters.physical_norm_prefactor_interval": K,
    }
    for path, expected in expected_numeric.items():
        actual = report
        for key in path.split("."):
            actual = actual[key]
        require(interval(*actual) == expected, "new report numeric discrepancy: "+path)
    require(j4[0] <= independently_derived["weighted_H4_integral_fresh"][0]
            <= independently_derived["weighted_H4_integral_fresh"][1] <= j4[1],
            "independent sharper clock integral is not inside new coarse-clock bracket")
    require(qphi2[0] <= independently_derived["Qphi_L2_squared"][0]
            <= independently_derived["Qphi_L2_squared"][1] <= qphi2[1],
            "independent sharper L2 interval is not inside new bracket")
    require(report["schema"] == "hdblast-quadratic-material-tension-forcing-v1" and report["status"] == "PASS", "new schema/status")
    require(report["state"]["ell"] == str(ell) and report["state"]["domain_x"] == [0,3]
            and report["state"]["derivative_convention"] == "Ueps=U0+eps U1+eps^2 U2/2"
            and report["state"]["normal_convention"] == "outward normal SUM over both copies", "new state/convention")
    expected_angular = {"sphere_volume_over_pi_squared":"2", "mean_Y2":"1/6", "mean_Y4":"3/40",
                        "mean_Z2_squared":"1/12", "mean_H4_squared":"19/720",
                        "orthogonal_squared_power_shares_l0_l2_l4":["10/27","5/18","19/54"]}
    require(report["angular"] == expected_angular, "new exact angular report")
    for key, expected in {
        "Y_squared_decomposition":["1/6","1/2","1"],
        "scalar_coefficient_multipliers_of_S_squared":["-2/3","-2","-4"],
        "tension_coefficient_multipliers_of_sigma_second_S_squared":["1/6","1/2","1"],
    }.items():
        require(report["definitions"][key] == expected, "new definition: "+key)
    require(report["arithmetic"]["value_denominator"] == D
            and report["arithmetic"]["time_grid_denominator"] == grid, "new interval denominators")
    require(report["cell_data"]["bytes"] == len(cell_raw)
            and report["cell_data"]["sha256"] == hash_bytes(cell_raw)
            and report["cell_data"]["rows"] == len(source_rows)
            and report["cell_data"]["header"] == header, "new cell metadata")
    require(F("-2.3e-19") < tail[0] < tail[1] < F("-6.7e-20")
            and report["tail_theorem"]["published_strict_interval"] == ["-2.3e-19","-6.7e-20"], "new strict tail display")
    require(report["scope"] == {"explicit_material_tension_contribution":True,
                "complete_bulk_coordinate_quadratic_source":False, "full_U2_response_computed":False,
                "total_tension_or_stress_second_variation":False, "nonlinear_remainder_bound":False,
                "novelty_claimed":False,"cosmological_evidence":False}, "new output scope")
    return {"status":"PASS", "all_new_cells_recomputed":len(source_rows),
            "all_new_csv_bytes_identical":True, "numeric_report_paths_checked":list(expected_numeric),
            "report_sha256":hash_bytes(report_raw), "csv_sha256":hash_bytes(cell_raw),
            "sharper_independent_clock_norm_intervals_inside_new_report":True,
            "producer_imported_or_executed":False,
            "check_groups":[
                "new_csv_all_4353_cells_reconstructed_with_exact_byte_equality",
                "new_report_twelve_numeric_moment_parameter_integral_norm_and_tail_paths",
                "new_report_state_angular_conventions_scope_metadata_and_strict_tail",
                "independent_sharper_clock_integral_and_norm_inclusion"]}


def serialize(obj):
    if isinstance(obj, F):
        return str(obj)
    if isinstance(obj, dict):
        return {k: serialize(v) for k,v in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [serialize(v) for v in obj]
    return obj


def display_intervals(obj):
    if isinstance(obj, F):
        return float(obj)
    if isinstance(obj, dict):
        return {k: display_intervals(v) for k,v in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [display_intervals(v) for v in obj]
    return obj


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=HERE.parent/"inputs")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-computed", type=Path, nargs=2, metavar=("REPORT", "CELLS"))
    parser.add_argument("--inject-error", choices=ERRORS)
    args = parser.parse_args()
    contents, authenticated, pins_hash = authenticated_inputs(args)
    groups = ["portable_input_byte_authentication_and_external_manifest_identity"]
    angular = angular_checks(args.inject_error)
    groups.append("ambient_harmonics_unit_S3_spectrum_exact_moments_norms_and_power_shares")
    constitutive = constitutive_checks(args.inject_error)
    groups.append("material_sum_copy_factorial_stress_product_cubic_algebra_and_source_sign_counterexample")
    norms = norm_checks(args.inject_error)
    groups.append("proper_time_background_measure_quadratic_quartic_norm_identities_and_Jacobians")
    computed = retained_data_checks(contents)
    groups.append("retained_4353_continuous_cells_H2_H4_clock_norm_and_tail_projection_consequences")
    if args.inject_error and args.inject_error.startswith("computed_"):
        require(args.check_computed is not None, "computed fault control requires --check-computed")
    comparison = (check_computed_files(args.check_computed, contents, computed, args.inject_error)
                  if args.check_computed else {"status":"NOT_REQUESTED"})
    groups.extend(comparison.get("check_groups",[]))
    result = {
        "schema": "hdblast-quadratic-sol-independent-v1", "status": "PASS",
        "checks": [{"name":name,"status":"PASS"} for name in groups],
        "check_names": groups,
        "check_group_count": len(groups),
        "runtime": {"python": sys.version.split()[0], "sympy": sp.__version__, "optimization": sys.flags.optimize},
        "external_manifest_sha256": MANIFEST_HASH, "prior_public_main": PUBLIC_MAIN,
        "authenticated_sources": authenticated, "source_pins_sha256": pins_hash,
        "checker_sha256": hash_bytes(Path(__file__).read_bytes()),
        "derivation_note_sha256": hash_bytes((HERE/"QUADRATIC_MATERIAL_FORCING.md").read_bytes()),
        "angular": angular, "constitutive": constitutive, "proper_time_norm_identities": norms,
        "exact_rational_consequences": serialize(computed),
        "new_production_comparison": comparison,
        "display_only": display_intervals(computed),
        "scope": "Known local quadratic material constitutive terms of the retained first variation only. Full U2, metric, normal, clock and bulk reconstruction are not computed; source sign gives no response sign. Exact cubic remainder is constitutive algebra, not an Einstein-family remainder.",
        "archived_scientific_code_executed": False,
    }
    data = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(data)
    else:
        print(data, end="")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(json.dumps({"status": "FAIL", "error": type(exc).__name__+": "+str(exc)}, sort_keys=True), file=sys.stderr)
        raise SystemExit(1)
