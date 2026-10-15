#!/usr/bin/env python3
"""Continuous interval enclosures for the exact constant-boundary free wave.

Only the standard library is used. All enclosure decisions use integer or
Fraction arithmetic, and all checks remain active with python -O. This is a
linear comparison-wave calculation, conditional on the stated PDE setting.
It does not execute an archived producer or import prior scientific code.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import sys

BITS = 112
S = 1 << BITS
ELL = F(1, 10**6)
T = 3 * ELL
KAPPA = F(1000)
DEGREE = 8
SHELL_PATH = "inputs/previous/archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json"
SHELL_PIN = "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8"


def ceil_div(n, d):
    if d <= 0:
        raise ValueError("positive divisor required")
    return -((-n) // d)


def iv(lo, hi=None):
    lo, hi = F(lo), F(lo if hi is None else hi)
    if lo > hi:
        raise ValueError("reversed interval")
    return lo.numerator * S // lo.denominator, ceil_div(hi.numerator * S, hi.denominator)


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def neg(a):
    return -a[1], -a[0]


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    p = [x * y for x in a for y in b]
    return min(p) // S, ceil_div(max(p), S)


def div_int(a, n):
    if not isinstance(n, int) or n <= 0:
        raise ValueError("positive integer divisor required")
    return a[0] // n, ceil_div(a[1], n)


def recip(a):
    if a[0] <= 0 <= a[1]:
        raise ValueError("division by an interval containing zero")
    if a[1] < 0:
        return neg(recip(neg(a)))
    return S * S // a[1], ceil_div(S * S, a[0])


def square(a):
    lo = 0 if a[0] <= 0 <= a[1] else min(a[0]**2, a[1]**2)
    return lo // S, ceil_div(max(a[0]**2, a[1]**2), S)


def widen(a, radius):
    radius = iv(radius)[1]
    return a[0] - radius, a[1] + radius


def record(a):
    return [str(F(a[0], S)), str(F(a[1], S))]


def expneg(x):
    """Integer enclosure of exp(-x), x rational and nonnegative."""
    x = F(x)
    if x < 0:
        raise ValueError("negative expneg argument")
    if not x:
        return S, S
    if x >= 128:
        # e>2 gives exp(-x)<=2^-floor(x); retain a nonzero upper tail.
        exponent = min(x.numerator // x.denominator, BITS)
        return 0, 1 << (BITS - exponent)
    reduced, halvings = x, 0
    while reduced > F(1, 8):
        reduced /= 2
        halvings += 1
    z = neg(iv(reduced))
    term = total = (S, S)
    for n in range(1, 29):
        term = div_int(mul(term, z), n)
        total = add(total, term)
    total = widen(total, F(2, 8**29 * factorial(29)))
    for _ in range(halvings):
        total = square(total)
    return max(0, total[0]), min(S, total[1])


def expneg_cell(lo, hi):
    if not 0 <= lo <= hi:
        raise ValueError("nonnegative ordered exponential cell required")
    return expneg(hi)[0], expneg(lo)[1]


def qvalue(v):
    if v is None:
        return 0, 0
    if v < 4:
        raise ValueError("bump reciprocal below four")
    return expneg(v - 4)


def q_cell(left, right):
    """Outward ranges of Q and Q' for every point of a closed cell."""
    left, right = F(left), F(right)
    if not 0 <= left <= right <= 3:
        raise ValueError("profile cell outside [0,3]")
    for split in (F(1), F(2)):
        if left < split < right:
            a, b = q_cell(left, split), q_cell(split, right)
            return tuple((min(u[0], v[0]), max(u[1], v[1])) for u, v in zip(a, b))
    if right <= 1 or left >= 2:
        return (0, 0), (0, 0)
    wl, wr = (left - 1) * (2 - left), (right - 1) * (2 - right)
    wmin = min(wl, wr)
    wmax = F(1, 4) if left <= F(3, 2) <= right else max(wl, wr)
    vmin, vmax = 1 / wmax, None if wmin == 0 else 1 / wmin
    q = qvalue(vmax)[0], qvalue(vmin)[1]
    envlo = 0 if vmax is None else mul(iv(vmax**2), qvalue(vmax))[0]
    envhi = mul(iv(vmin**2), qvalue(vmin))[1]
    # v^2 exp(4-v) decreases for v>=4.
    qp = mul(iv(3 - 2 * right, 3 - 2 * left), (envlo, envhi))
    return q, qp


class BumpMass:
    """Continuous primitive enclosure from outward interval rectangles."""
    def __init__(self, count=4096):
        if count < 16:
            raise ValueError("at least sixteen bump cells required")
        self.count = count
        self.prefix = [(0, 0)]
        for j in range(count):
            cell = q_cell(F(1) + F(j, count), F(1) + F(j + 1, count))[0]
            self.prefix.append(add(self.prefix[-1], div_int(cell, count)))

    def at(self, x):
        x = F(x)
        if x <= 1:
            return 0, 0
        if x >= 2:
            return self.prefix[-1]
        index = (x - 1) * self.count
        j = index.numerator // index.denominator
        left = F(1) + F(j, self.count)
        return add(self.prefix[j], mul(iv(x - left), q_cell(left, x)[0]))

    def cell(self, left, right):
        # Q>=0 proves monotonicity of its primitive.
        return self.at(left)[0], self.at(right)[1]


def kernel(r, mu, gamma, fault=None):
    """Enclose k,k' and P=-(k''+gamma*k'+exp(-s)*(k'+gamma*k)).

    The returned first two components are the Green kernel of
    mu*k''+k'+(gamma+5mu/2)*k=0, k(0)=0, k'(0)=1/mu.
    P additionally requires the incoming-coordinate exponential below.
    """
    if not 0 <= r[0] <= r[1] <= iv(T)[1]:
        raise ValueError("kernel interval outside certified short time")
    a = recip(mu)
    c = add(mul(gamma, a), iv(F(5, 2)))
    # In y=(KAPPA*k,k'), ||A||_infinity<=KAPPA and ||y(0)||<=a_hi.
    if F(a[1], S) + F(c[1], S) / KAPPA > KAPPA:
        raise ValueError("weighted matrix-exponential norm gate failed")
    tk = total_k = (0, 0)
    tp = total_p = neg(a) if fault == "reverse_kernel_initial" else a
    for n in range(1, DEGREE + 1):
        tk, tp = div_int(mul(tp, r), n), div_int(mul(neg(add(mul(c, tk), mul(a, tp))), r), n)
        total_k, total_p = add(total_k, tk), add(total_p, tp)
    # exp(KAPPA*T)<2 follows from KAPPA*T=3/1000<1/2 and
    # sum x^n/n! <=sum x^n=1/(1-x)<2.
    radius_p = F(a[1], S) * 2 * (KAPPA * F(r[1], S))**(DEGREE + 1) / factorial(DEGREE + 1)
    total_k = widen(total_k, radius_p / KAPPA)
    total_p = widen(total_p, radius_p)
    # A global short-time bound, independent of the supplied subinterval:
    # |k'(r)-k'(0)| <= a_hi*(exp(KAPPA*T)-1)
    # <= a_hi*KAPPA*T/(1-KAPPA*T). Thus k'>0 on [0,T], and k>=0.
    global_positive = F(a[0], S) > F(a[1], S) * (KAPPA * T) / (1 - KAPPA * T)
    if fault != "reverse_kernel_initial" and global_positive:
        total_k = max(0, total_k[0]), total_k[1]
    return total_k, total_p, a, c, radius_p


def positive_kernel(rlo, rhi, slo, shi, mu, gamma, fault=None):
    k, kp, a, c, _ = kernel(iv(rlo, rhi), mu, gamma, fault)
    e = expneg_cell(slo, shi)
    # Exact ODE reduction; it avoids subtracting two large accelerations.
    p = add(mul(sub(sub(a, gamma), e), kp), mul(sub(c, mul(gamma, e)), k))
    return neg(p) if fault == "reverse_kernel_sign" else p


def comparator_cell(left, right, mu, gamma, mass, fault=None):
    """Enclose H0 at every x in a closed [left,right] subset of [0,3]."""
    left, right = F(left), F(right)
    q, qp = q_cell(left, right)
    if right <= 1:
        return 0, 0
    a = recip(mu)
    local = expneg_cell(ELL * left, ELL * right)
    if fault != "omit_local_inertia":
        local = sub(local, a)
    if fault == "reverse_q_derivative":
        qp = neg(qp)
    base = add(qp, mul(iv(ELL), mul(local, q)))
    p = positive_kernel(F(0), ELL * (right - 1), ELL, ELL * min(right, F(2)), mu, gamma, fault)
    correction = mul(iv(ELL**2), mul(p, mass.cell(left, right)))
    return add(base, correction)


def load_shell(source_root):
    data = (Path(source_root) / SHELL_PATH).read_bytes()
    digest = sha256(data).hexdigest()
    if digest != SHELL_PIN:
        raise RuntimeError("shell input authentication failed")
    source = json.loads(data)
    if source["parameters_exact"] != {"delta": "1/1000", "c": "96714321425692/161839258930731"}:
        raise RuntimeError("retained exact parameters changed")
    mu = iv(*source["boxes"]["mu"]["exact"])
    gamma = mul(iv(*source["boxes"]["C_b"]["exact"]), iv(F(1, 2)))
    return mu, gamma, {"path": SHELL_PATH, "sha256": digest, "bytes": len(data)}


def run(source_root, count=4096, fault=None):
    checks = []

    def require(name, predicate):
        if not predicate:
            raise RuntimeError(name)
        checks.append(name)

    mu, gamma, pin = load_shell(source_root)
    require("mu lower positive", mu[0] > 0)
    require("gamma upper below 250", gamma[1] < 250 * S)
    require("weighted norm exponential argument below one half", KAPPA * T < F(1, 2))
    mass = BumpMass(count)
    require("bump integral positive and below one", 0 < mass.prefix[-1][0] <= mass.prefix[-1][1] < S)
    k, kp, _, _, radius = kernel(iv(0, T), mu, gamma, fault)
    require("kernel derivative strictly positive on complete interval", kp[0] > 0)
    require("kernel nonnegative on complete interval", k[0] >= 0)
    p = positive_kernel(0, 2 * ELL, ELL, 2 * ELL, mu, gamma, fault)
    require("positive postpulse integrand on full continuous rectangle", p[0] > 0)
    tail = comparator_cell(2, 3, mu, gamma, mass, fault)
    center = comparator_cell(F(3, 2), F(3, 2), mu, gamma, mass, fault)
    left = comparator_cell(F(5, 4), F(5, 4), mu, gamma, mass, fault)
    right = comparator_cell(F(7, 4), F(7, 4), mu, gamma, mass, fault)
    require("postpulse comparator strictly positive continuously", tail[0] > 0)
    require("center retains negative local inertia correction", center[1] < -iv(F(4, 10**4))[1])
    require("positive incoming lobe", left[0] > 3 * S)
    require("negative outgoing lobe", right[1] < -3 * S)
    require("prefront comparison exactly zero", comparator_cell(0, 1, mu, gamma, mass, fault) == (0, 0))
    require("short-time remainder below 1e-20", radius < F(1, 10**20))
    return {
        "status": "PASS_EXACT_FREE_INERTIAL_COMPARATOR",
        "python_optimization": sys.flags.optimize,
        "arithmetic": "112-bit outward integer fixed point plus exact Fraction decisions",
        "source": pin, "ell": str(ELL), "T": str(T), "mu": record(mu), "gamma": record(gamma),
        "kernel_weight": str(KAPPA), "kernel_Taylor_degree": DEGREE,
        "kernel_k_global": record(k), "kernel_kprime_global": record(kp),
        "kernel_kprime_remainder_upper": str(radius),
        "positive_tail_kernel_global": record(p), "bump_mass": record(mass.prefix[-1]),
        "H0_tail_continuous_2_3": record(tail), "H0_center": record(center),
        "H0_at_5_4": record(left), "H0_at_7_4": record(right),
        "mass_mesh_cells": count, "checks": checks, "check_count": len(checks),
        "scope": "Exact free-wave comparator only. Actual bulk and time-dependent-boundary residuals require separate analytic bounds. No nonlinear or physical-unit assertion.",
        "own_source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--mass-cells", type=int, default=4096)
    parser.add_argument("--fault", choices=("reverse_kernel_initial", "reverse_kernel_sign", "omit_local_inertia", "reverse_q_derivative"))
    args = parser.parse_args()
    try:
        print(json.dumps(run(args.source_root, args.mass_cells, args.fault), indent=2, sort_keys=True))
    except (RuntimeError, ValueError, OSError, KeyError) as exc:
        print("FAIL: " + str(exc), file=sys.stderr)
        raise SystemExit(1)
