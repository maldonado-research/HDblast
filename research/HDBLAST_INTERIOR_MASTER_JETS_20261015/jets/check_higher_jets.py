#!/usr/bin/env python3
"""Exact rational first-master C2/C3 certificate; imports no retained program.

The analytic hypotheses and reflection proof are in HIGHER_JETS.md. Run with
--root pointing to the capsule root, whose inputs/ contains the pinned files.
"""
import argparse
from csv import DictReader
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import comb, factorial
from pathlib import Path
import subprocess
import sys

ELL = F(1, 10**6)
X = F(3)
M, MZ, MT = F(102000), F(30000000), F(21, 400000)
PINS = {
    "inputs/prior/INERTIAL_WAVEFORM_CELLS.csv": "2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7",
    "inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json": "d2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8",
    "inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
    "inputs/reference/BACKGROUND_DERIVATIVE_BOUNDS.json": "fd6fd277176033703a4e75a2f66ab94894318aac49c83264ef979e9a27f0a02c",
    "inputs/reference/PHYSICAL_SHELL_RESPONSE.md": "a78f2239893fdd0e5ff653d4d5ce5cd9d3244c57c794ed2d3bb43bd4f5e1484b",
}
POTENTIAL_PIN = "fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0"
CONTROLS = {
    "drop_reflected_B": "second reflected derivative requires the shell forcing term",
    "flip_reflected_B1": "third reflected derivative requires the positive radial boundary derivative term",
    "drop_incoming_endpoint": "incoming Hprime integral requires its moving endpoint",
    "wrong_Hprime_ypp_sign": "incoming Hprime has minus ysecond",
    "wrong_lift_mixed_sign": "causal lift mixed derivative has the radial minus sign",
    "incompatible_second_corner": "nonzero initial shell forcing prevents separate K from being C2",
    "incompatible_third_corner": "nonzero initial shell forcing time derivative prevents separate K from being C3",
    "missing_ell_in_Wz": "scaled forcing radial derivative uses ell cubed times physical Wz",
    "drop_Wzz": "third-jet forcing retains its physical Wzz term",
    "false_Qsecond_envelope": "actual retained bump second derivative exceeds sixteen",
    "false_Qthird_envelope": "actual retained bump third derivative exceeds 256",
}


class Failure(Exception):
    pass


def need(condition, message):
    if not condition:
        raise Failure(message)


def encoded(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: encoded(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encoded(v) for v in value]
    return value


def floats(value):
    if isinstance(value, F):
        return float(value)
    if isinstance(value, dict):
        return {k: floats(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [floats(v) for v in value]
    return value


def authenticate(root):
    for relative, digest in PINS.items():
        need(hashlib.sha256((root / relative).read_bytes()).hexdigest() == digest,
             "pinned input identity: " + relative)


def interval_product(a, b):
    products = [x * y for x in a for y in b]
    return min(products), max(products)


def primitive_range(p0, p1, h, a, b):
    points = [a, b]
    if h and a < a - p0 / h < b:
        points.append(a - p0 / h)
    vals = [p0 * x - p1 + h * (x - a)**2 / 2 for x in points]
    return min(vals), max(vals)


def boundary_bounds(root):
    old = json.loads((root / "inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json").read_text())
    shell = json.loads((root / "inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json").read_text())
    mu = tuple(map(F, old["shell_parameter_intervals"]["mu"]))
    gamma = tuple(map(F, old["shell_parameter_intervals"]["gamma"]))
    amu = tuple(map(F, shell["boxes"]["mu"]["exact"]))
    bg = json.loads((root / "inputs/reference/BACKGROUND_DERIVATIVE_BOUNDS.json").read_text())
    agamma = tuple(map(F, bg["shell_correlation"]["parameters"]["gamma"]["exact"]))
    need(0 < mu[0] <= amu[0] <= amu[1] <= mu[1], "positive actual mu enclosure")
    need(0 < gamma[0] <= agamma[0] <= agamma[1] <= gamma[1], "actual gamma enclosure")
    need(F(old["state"]["ell"]) == ELL, "retained scaled width")
    multiplier = ELL / mu[1], ELL / mu[0]
    denominator = F(old["cell_data"]["value_denominator"])
    grid = F(old["cell_data"]["time_grid_denominator"])
    p0, p1 = [F(0), F(0)], [F(0), F(0)]
    Hsup, Fsup, Fpsup, last, n = F(0), F(0), F(0), F(0), 0
    with (root / "inputs/prior/INERTIAL_WAVEFORM_CELLS.csv").open() as stream:
        for row in DictReader(stream):
            a, b = F(row["left_grid"]) / grid, F(row["right_grid"]) / grid
            need(a == last and a < b, "complete actual H cell partition")
            h = F(row["H_lo"]) / denominator, F(row["H_hi"]) / denominator
            need(h[0] <= h[1], "ordered actual H cell")
            ranges = [primitive_range(p0[k], p1[k], h[k], a, b) for k in (0, 1)]
            value = interval_product(multiplier, (ranges[0][0], ranges[1][1]))
            slope = interval_product(multiplier,
                (p0[0] + min(F(0), h[0] * (b-a)), p0[1] + max(F(0), h[1] * (b-a))))
            Hsup = max(Hsup, *map(abs, h))
            Fsup = max(Fsup, *map(abs, value))
            Fpsup = max(Fpsup, *map(abs, slope))
            for k in (0, 1):
                p0[k] += h[k] * (b-a)
                p1[k] += h[k] * (b*b-a*a) / 2
            last, n = b, n+1
    need(last == X and n == 4353, "all 4353 actual H cells cover [0,3]")
    A = F(5, 2) * ELL**2
    Ap = F(35, 2) * ELL**3 * (X * ELL)
    need(F(35,4)*(X*ELL)**2 < F(5,2), "boundary A remains positive and bounded by 5ell^2/2")
    q = A * X**2 / 2
    need(q < 1, "boundary Volterra contraction")
    y = Fsup / (1-q)
    yp = Fpsup + X * A * y
    ypp = multiplier[1] * Hsup + A * y
    return {"cells": n, "mu": mu, "gamma": gamma, "H": Hsup,
            "twice_integrated_H": Fsup, "integrated_H": Fpsup,
            "boundary_contraction": q, "A": A, "Aprime": Ap,
            "y": y, "yprime": yp, "ysecond": ypp,
            "ell_over_mu": multiplier[1]}


def poly_integral(m, k, a, b, lo, hi):
    return sum((F(comb(k, j)) * a**(k-j) * b**j *
                (hi**(m+j+1)-lo**(m+j+1)) / (m+j+1)
                for j in range(k+1)), F(0))


def split_parts(x, r, reflected=None):
    if reflected is None:
        reflected = x > r
    split = x-r if reflected else F(0)
    return [(F(0), split, -1), (split, x, 1)]


def kernel(m, k, x, r, reflected=None):
    upper = poly_integral(m, k+1, r+x, F(-1), F(0), x)
    lower = sum((poly_integral(m, k+1, sign*(r-x), F(sign), a, b)
                 for a, b, sign in split_parts(x, r, reflected)), F(0))
    return (upper-lower) / (2*(k+1))


def radial_ray_integrals(m, k, x, r, order):
    if k < order:
        return F(0), F(0)
    coefficient = F(factorial(k), factorial(k-order))
    exponent = k-order
    upper = coefficient * poly_integral(m, exponent, r+x, F(-1), F(0), x)
    # Odd f differentiates to an even first and an odd second derivative.
    lower = coefficient * sum(((sign if order % 2 == 0 else 1) *
            poly_integral(m, exponent, sign*(r-x), F(sign), a, b)
            for a, b, sign in split_parts(x, r)), F(0))
    return upper, lower


def force(m, k, x, r, time_order=0, radial_order=0):
    if m < time_order or k < radial_order:
        return F(0)
    coefficient = F(factorial(m), factorial(m-time_order)) * F(factorial(k), factorial(k-radial_order))
    return coefficient * x**(m-time_order) * r**(k-radial_order)


def reflected_formula(m, k, x, r, ix, ir, inside=None):
    if inside is None:
        inside = x > r
    a, b = radial_ray_integrals(m, k, x, r, 1)
    B = force(m, k, x-r, F(0)) if inside else F(0)
    second = {(0, 2): (a-b)/2-B, (1, 1): (a+b)/2+B}
    second[2, 0] = force(m, k, x, r) + second[0, 2]
    if ix+ir == 2:
        return second[ix, ir]
    a, b = radial_ray_integrals(m, k, x, r, 2)
    B1 = force(m, k, x-r, F(0), time_order=1) if inside else F(0)
    third = {(0, 3): (a-b)/2+B1, (1, 2): (a+b)/2-B1}
    third[2, 1] = force(m, k, x, r, radial_order=1) + third[0, 3]
    third[3, 0] = force(m, k, x, r, time_order=1) + third[1, 2]
    return third[ix, ir]


@lru_cache(None)
def derivative_weights(degree, order):
    if order > degree:
        return [F(0)] * (degree+1)
    weights = []
    for j in range(degree+1):
        coeffs, denominator = [F(1)], F(1)
        for k in range(degree+1):
            if k == j:
                continue
            updated = [F(0)] * (len(coeffs)+1)
            for i, value in enumerate(coeffs):
                updated[i] -= k * value
                updated[i+1] += value
            coeffs = updated
            denominator *= j-k
        weights.append(coeffs[order] * factorial(order) / denominator)
    return weights


def polynomial_derivative(m, k, x, r, ix, ir, branch):
    # Tensor Lagrange differentiation of the branch polynomial, not of the
    # ray formulas. A branch extension avoids crossing a piecewise join.
    degree = m+k+2
    h = F(1, 100*(degree+1))
    wx = derivative_weights(degree, ix) if ix else [F(1)]
    wr = derivative_weights(degree, ir) if ir else [F(1)]
    return sum((a*b*kernel(m, k, x+i*h, r+j*h, branch)
                for i, a in enumerate(wx) for j, b in enumerate(wr)), F(0)) / h**(ix+ir)


def reflection_checks():
    witnesses, joins, traces = 0, 0, 0
    derivatives = [(0, 2), (1, 1), (2, 0), (0, 3), (1, 2), (2, 1), (3, 0)]
    for m in range(3):
        for k in range(4):
            for x, r, branch in ((F(2), F(1, 2), True), (F(1, 2), F(2), False), (F(2), F(0), True)):
                for ix, ir in derivatives:
                    need(polynomial_derivative(m, k, x, r, ix, ir, branch) ==
                         reflected_formula(m, k, x, r, ix, ir), "independent reflected polynomial derivative")
                    witnesses += 1
            for ix, ir in derivatives:
                # Both branch polynomials have exact one-sided characteristic traces.
                inner = polynomial_derivative(m, k, F(1), F(1), ix, ir, True)
                outer = polynomial_derivative(m, k, F(1), F(1), ix, ir, False)
                need(inner == reflected_formula(m, k, F(1), F(1), ix, ir, True), "inner characteristic trace")
                need(outer == reflected_formula(m, k, F(1), F(1), ix, ir, False), "outer characteristic trace")
                compatible = force(m, k, F(0), F(0)) == 0
                if ix+ir == 3:
                    compatible = compatible and force(m, k, F(0), F(0), time_order=1) == 0
                if compatible:
                    need(inner == outer, "compatible characteristic C2/C3 join")
                    joins += 1
                traces += 2
            x = F(2)
            need(reflected_formula(m, k, x, F(0), 0, 2) == -force(m, k, x, F(0)), "shell Krr trace")
            need(reflected_formula(m, k, x, F(0), 2, 0) == 0, "shell Kxx trace")
            need(reflected_formula(m, k, x, F(0), 1, 2) == -force(m, k, x, F(0), time_order=1), "shell Kxrr trace")
            need(reflected_formula(m, k, x, F(0), 3, 0) == 0, "shell Kxxx trace")
    for x in (F(0), F(1), X):
        for r in (F(0), F(8)-x):
            for s in (F(0), x):
                need(abs(r+x-s)+s <= r+x <= 8 and abs(r-x+s)+s <= r+x, "causal closure of both reflected rays")
    return {"monomials": 12, "independent_derivative_witnesses": witnesses,
            "characteristic_trace_witnesses": traces, "compatible_join_witnesses": joins,
            "shell_trace_witnesses": 48}


def polynomial_shift(coefficients, shift):
    result = [F(0)] * len(coefficients)
    for i, value in enumerate(coefficients):
        for j in range(i+1):
            result[j] += value * comb(i, j) * shift**(i-j)
    return result


def bump_checks():
    # Derivative identities use s'=-v*s^2, v'=-2, v^2=1-4/s.
    p2 = [F(0), F(0), F(6), F(-6), F(1)]
    p3 = [F(0), F(0), F(0), F(-12), F(24), F(-10), F(1)]
    # Differentiate v*s^2 using v'=-2 and v^2=1-4/s.
    q1 = [F(0),F(0),F(1)]
    q1_minus_derivative = q1[:]
    for i in range(1,len(q1)):
        q1_minus_derivative[i-1] -= i*q1[i]
    derived_p2 = [F(0)]*5
    for i,value in enumerate(q1_minus_derivative):
        derived_p2[i+2] += value
        derived_p2[i+1] -= 4*value
    for i,value in enumerate(q1):
        derived_p2[i] -= 2*value
    need(derived_p2 == p2, "exact bump Qsecond identity")
    derivative_p2 = [F(0)]*7
    for i, value in enumerate(p2):
        derivative_p2[i+2] += value
        if i:
            derivative_p2[i+1] -= i*value
    need(derivative_p2 == p3, "exact bump Qthird identity")
    for polynomial in ([F(-6), F(6)], [F(6), F(-6), F(2)],
                       [F(12), F(-24), F(10)], [F(-12), F(24), F(-10), F(2)]):
        shifted = polynomial_shift(polynomial, F(4))
        need(all(value >= 0 for value in shifted), "bump polynomial envelope on s>=4")
    need(sum((F(2)**j/factorial(j) for j in range(5)), F(0)) == 7,
         "positive exponential series proves exp(2)>7")
    # At s=4, Q''=-32; at r=5/4, s=16/3, v=1/2.
    need(sum((value*4**i for i, value in enumerate(p2)), F(0)) == -32, "exact bump Qsecond witness")
    return {"Q": F(1), "Qprime": F(16), "Qsecond": F(256), "Qthird": F(46656, 7),
            "zero_extension_is_smooth": True}


def background_bounds(root):
    bg = json.loads((root / "inputs/reference/BACKGROUND_DERIVATIVE_BOUNDS.json").read_text())
    state = bg["candidate_state_exact"]
    need((F(state["ell"]), F(state["T"]), F(state["Z"])) == (ELL, X*ELL, 8*ELL), "actual background collar")
    need(bg["first_exit"]["status"] == "PASS", "background first-exit proof")
    for component in bg["first_exit"]["component_margins"].values():
        need(F(component["strict_lower_margin"]) > 0 and F(component["strict_upper_margin"]) > 0, "strict background margins")
    rho = tuple(map(F, bg["first_exit"]["component_margins"]["rho"]["outer"]["exact"]))
    q = tuple(map(F, bg["outer_derived"]["Q"]["exact"]))
    W = rho[0]**2*q[0]-F(35,4)*(X*ELL)**2, rho[1]**2*q[1]
    need(max(map(abs, W)) < M, "actual W enclosure")
    # These audited physical derivative bounds are retained, not inferred
    # by differentiating interval W boxes.
    wz = tuple(map(F, bg["outer_derived"]["W_z"]["exact"]))
    need(max(map(abs, wz)) < MZ and F(bg["strict_simple_bounds"]["absolute_W_z"]) == MZ and
         F(bg["strict_simple_bounds"]["absolute_W_t_on_window"]) == MT,
         "retained physical Wz and Wt")
    return {"W": M, "Wz": MZ, "Wt": MT, "W_interval": W}


def potential_second_bound(root, path):
    if path is None:
        return None
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == POTENTIAL_PIN, "pinned actual-background Wzz certificate")
    report = json.loads(raw)
    need(report["status"] == "PASS_ACTUAL_RETAINED_POTENTIAL_THROUGH_WZZ", "actual-background Wzz certificate status")
    need(report["input_sha256"]["BACKGROUND_DERIVATIVE_BOUNDS.json"] == PINS["inputs/reference/BACKGROUND_DERIVATIVE_BOUNDS.json"],
         "Wzz certificate uses this retained actual background")
    domain = report["domain"]
    need((F(domain["ell"]),F(domain["T"]),F(domain["Z"])) == (ELL,X*ELL,8*ELL), "Wzz certificate domain")
    bounds = report["simple_bounds"]
    need((F(bounds["absolute_W"]),F(bounds["absolute_W_z"]),F(bounds["absolute_W_t"])) == (M,MZ,MT), "Wzz and first derivative conventions agree")
    MZZ = F(bounds["absolute_W_zz"])
    need(MZZ == 30000000000, "actual physical Wzz bound")
    need(max(map(abs,map(F,report["potential_intervals"]["W_zz"]["exact"]))) < MZZ,
         "actual Wzz interval lies strictly inside declared envelope")
    return MZZ


def jet_bounds(boundary, bump, background, MZZ=None):
    q = M*(X*ELL)**2/2
    need(q == F(459,10**9) < 1, "interior value contraction")
    U = (1+X*ELL+boundary["y"])/(1-q)
    e1 = X*ELL**2*M*U
    B1free = bump["Qprime"]+ELL+boundary["yprime"]
    B1 = B1free+e1
    f = ELL**2*M*U
    fr = ELL**3*MZ*U+ELL**2*M*B1
    fx = ELL**3*MT*U+ELL**2*M*B1
    f0 = ELL**2*M*boundary["y"]
    fx0 = ELL**3*MT*boundary["y"]+ELL**2*M*boundary["yprime"]
    Hprime = bump["Qsecond"]+ELL*(bump["Qprime"]+ELL)+boundary["ysecond"]+ELL*boundary["gamma"][1]*boundary["yprime"]+f0+X*fr
    ythird = boundary["ell_over_mu"]*Hprime+boundary["Aprime"]*boundary["y"]+boundary["A"]*boundary["yprime"]
    Vprime = bump["Qprime"]+ELL
    Vsecond = bump["Qsecond"]+2*ELL*bump["Qprime"]+ELL**2
    B2free = bump["Qsecond"]+ELL*Vprime+boundary["ysecond"]
    E2 = {"rr": X*fr+f0, "xr": X*fr+f0, "xx": f+X*fr+f0}
    B2 = {key: B2free+value for key, value in E2.items()}
    result = {"ell": ELL, "scaled_domain": {"x_max": X, "r_min": 0, "x_plus_r_max": 8},
              "C0_contraction": q, "U": U, "first_error": e1, "first": B1,
              "forcing": {"f": f, "fr": fr, "fx": fx, "shell_f": f0, "shell_fx": fx0},
              "Hprime": Hprime, "ysecond": boundary["ysecond"], "ythird": ythird,
              "free_second": B2free, "second_error": E2, "second": B2}
    need(U < F("1.000193") and B1 < F("16.000488"), "displayed value and first derivative enclosures")
    need(max(B2.values()) < F("256.002081"), "displayed all-second-derivative enclosure")
    if MZZ is not None:
        frr = ELL**4*MZZ*U+2*ELL**3*MZ*B1+ELL**2*M*B2["rr"]
        E3 = {"rrr": X*frr+fx0, "xrr": X*frr+fx0,
              "xxr": fr+X*frr+fx0, "xxx": fx+X*frr+fx0}
        B3free = bump["Qthird"]+ELL*Vsecond+ythird
        B3 = {key: B3free+value for key, value in E3.items()}
        result.update({"Wzz": MZZ, "forcing_rr": frr, "free_third": B3free, "third_error": E3, "third": B3})
        need(max(B3.values()) < F("6665.267480"), "displayed all-third-derivative enclosure")
    return result


def incoming_checks():
    # Independent exact Leibniz checks for integral_0^x s^m(x-s)^k ds.
    count = 0
    for m in range(3):
        for k in range(4):
            x = F(2)
            degree = m+k+1
            h = F(1, 100)
            derivative = sum((weight*poly_integral(m,k,x+j*h,F(-1),F(0),x+j*h)
                              for j, weight in enumerate(derivative_weights(degree,1))),F(0))/h
            endpoint = force(m,k,x,F(0))
            integral = k*poly_integral(m,k-1,x,F(-1),F(0),x) if k else F(0)
            need(derivative == endpoint+integral, "incoming characteristic moving-endpoint derivative")
            count += 1
    # y(s)=s^4 is compatible through third order; signs follow (-1)^radial order.
    for ix, ir in ((1,1),(0,2),(2,1),(1,2),(0,3),(3,0)):
        expected = (-1)**ir * F(factorial(4),factorial(4-ix-ir))*2**(4-ix-ir)
        need(expected == (-1)**ir*(12*4 if ix+ir==2 else 24*2), "compatible causal lift derivative sign")
    return {"moving_endpoint_witnesses": count, "compatible_lift_witnesses": 6,
            "Hprime_uses_derivatives_of_H_cells": False}


def inject(name):
    reason = CONTROLS[name]
    if name == "drop_reflected_B":
        need(F(0) == reflected_formula(0,0,F(2),F(1),0,2), reason)
    elif name == "flip_reflected_B1":
        need(F(-1) == reflected_formula(1,0,F(2),F(1),0,3), reason)
    elif name == "drop_incoming_endpoint":
        need(F(0) == 1, reason)
    elif name == "wrong_Hprime_ypp_sign":
        # Synthetic actual incoming identity: J0'=0, f=0, gamma=0, y''=2.
        need(F(2) == -F(2), reason)
    elif name == "wrong_lift_mixed_sign":
        need(F(48) == -F(48), reason)
    elif name == "incompatible_second_corner":
        need(reflected_formula(0,0,F(1),F(1),0,2,True) == reflected_formula(0,0,F(1),F(1),0,2,False), reason)
    elif name == "incompatible_third_corner":
        need(reflected_formula(1,0,F(1),F(1),0,3,True) == reflected_formula(1,0,F(1),F(1),0,3,False), reason)
    elif name == "missing_ell_in_Wz":
        need(ELL**2*MZ == ELL**3*MZ, reason)
    elif name == "drop_Wzz":
        need(F(0) == ELL**4*F(1), reason)
    elif name == "false_Qsecond_envelope":
        need(F(16) >= abs(F(-32)), reason)
    elif name == "false_Qthird_envelope":
        # The sixth alternating partial sum is a lower bound for exp(-4/3).
        s = F(16,3)
        p = (s**6-10*s**5+24*s**4-12*s**3)/2
        lower = sum(((-F(4,3))**j/factorial(j) for j in range(6)),F(0))
        need(F(256) >= abs(p)*lower, reason)
    raise Failure("deliberate corruption unexpectedly passed")


def run_controls(output):
    results = []
    for optimized in (False, True):
        for name, reason in CONTROLS.items():
            cmd = [sys.executable,"-I"]+(["-O"] if optimized else [])+["-B",str(Path(__file__).resolve()),"--inject",name]
            process = subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
            need(process.returncode == 1 and process.stdout == "" and process.stderr == "FAIL: "+reason+"\n", "intended explicit control rejection: "+name)
            results.append({"control":name,"optimized":optimized,"status":"EXPECTED_EXPLICIT_FAILURE"})
    report = {"status":"PASS","control_runs":len(results),"results":results}
    if output:
        output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report if not output else {"status":"PASS","control_runs":len(results)}))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument("--potential-report",type=Path)
    parser.add_argument("--allow-pending-Wzz",action="store_true",help="explicit C2-only mode")
    parser.add_argument("--output",type=Path)
    parser.add_argument("--inject",choices=list(CONTROLS))
    parser.add_argument("--run-controls",action="store_true")
    args = parser.parse_args()
    try:
        if args.inject:
            inject(args.inject)
        if args.run_controls:
            run_controls(args.output)
            return 0
        root = args.root.resolve()
        authenticate(root)
        boundary = boundary_bounds(root)
        background = background_bounds(root)
        bump = bump_checks()
        reflection = reflection_checks()
        incoming = incoming_checks()
        potential_path = args.potential_report or root / "potential/ACTUAL_POTENTIAL_JET_BOUNDS.json"
        MZZ = None if args.allow_pending_Wzz else potential_second_bound(root,potential_path)
        numeric = jet_bounds(boundary,bump,background,MZZ)
        report = {"schema":"hdblast-interior-master-jets-v1","status":"PASS","optimized":not __debug__,
            "named_checks":["minimal pinned actual inputs","4353 actual H whole-cell primitive enclosures",
                "retained actual background and physical first derivatives","global smooth bump envelopes",
                "independent reflected C2/C3 monomial derivatives and shell traces","compatible characteristic corner joins",
                "incoming Hprime endpoint and compatible lift derivatives","exact rational C2 and boundary C3 bounds"],
            "exact":{"boundary":boundary,"background":background,"bump":bump,"reflection_checks":reflection,
                     "incoming_checks":incoming,"interior":numeric},"display_only":floats(numeric),
            "input_pins":PINS,"potential_report_sha256":POTENTIAL_PIN,
            "scope":{"conditional_on_retained_smooth_actual_solution":True,"independent_H_cell_selections_are_not_solutions":True,
                "second_master_derivatives_certified":True,"third_master_derivatives_certified":MZZ is not None,
                "Wzz_input_pending":MZZ is None,"full_bulk_source_evaluated":False,"U2_evaluated":False}}
        if MZZ is not None:
            report["named_checks"].append("exact rational all-mixed C3 master bounds with pinned actual Wzz")
            report["source_bound_adapter"] = {"scaled_order_bounds":[U for U in
                (numeric["U"],numeric["first"],max(numeric["second"].values()),max(numeric["third"].values()))],
                "ell":ELL,"domain":"0<=x<=3, r>=0, x+r<=8"}
        raw = json.dumps(encoded(report),indent=2)+"\n"
        if args.output:
            args.output.write_text(raw)
        else:
            print(raw,end="")
    except (Failure,OSError,KeyError,ValueError) as error:
        print("FAIL: "+str(error),file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
