#!/usr/bin/env python3
"""Exact fabricated fixtures for direct integration; no physical source calls."""
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from math import comb
import argparse
import json
from pathlib import Path


def need(value, label):
    if not value:
        raise RuntimeError(label)


@dataclass(frozen=True)
class CQ:
    re: Q = Q(0)
    im: Q = Q(0)

    def __add__(self, other):
        if not isinstance(other, CQ):
            other = CQ(Q(other))
        return CQ(self.re + other.re, self.im + other.im)

    __radd__ = __add__

    def __neg__(self):
        return CQ(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-other if isinstance(other, CQ) else -Q(other))

    def __mul__(self, other):
        if not isinstance(other, CQ):
            other = CQ(Q(other))
        return CQ(self.re * other.re - self.im * other.im,
                  self.re * other.im + self.im * other.re)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return CQ(self.re / other, self.im / other)

    def abs_upper(self):
        return abs(self.re) + abs(self.im)


ZERO = CQ()


def poly_add(a, b):
    return [(a[i] if i < len(a) else ZERO) +
            (b[i] if i < len(b) else ZERO)
            for i in range(max(len(a), len(b)))]


def poly_mul(a, b):
    out = [ZERO] * (len(a) + len(b) - 1)
    for i, av in enumerate(a):
        for j, bv in enumerate(b):
            out[i + j] = out[i + j] + av * bv
    return out


def poly_integral(a, h):
    return sum((v * h ** (i + 1) / (i + 1)
                for i, v in enumerate(a)), ZERO)


def poly_value(a, x):
    return sum((v * x ** i for i, v in enumerate(a)), ZERO)


def interval_disk_inside(diff, reference_radius, production_radius, label):
    margin = production_radius - reference_radius
    need(margin >= 0 and diff.re ** 2 + diff.im ** 2 <= margin ** 2, label)


def shifted_polynomial(coeff, shift):
    out = [Q(0)] * len(coeff)
    for n, av in enumerate(coeff):
        for j in range(n + 1):
            out[j] += av * comb(n, j) * shift ** (n - j)
    return out


def inv_poly(p, degree, c, h):
    coeff = [Q((-1) ** n * comb(p + n - 1, n)) / c ** (p + n)
             for n in range(degree + 1)]
    return shifted_polynomial(coeff, -h / 2)


def geom_tail(p, degree, c, h):
    r = h / (2 * abs(c))
    need(Q(0) <= r < 1, "geometry expansion excludes its pole")
    m = degree + 1
    if p == 1:
        tail = r ** m / (1 - r)
    elif p == 2:
        tail = r ** m * (Q(m + 1) / (1 - r) + r / (1 - r) ** 2)
    elif p == 3:
        tail = r ** m * (Q(comb(m + 2, 2)) / (1 - r) +
                        Q(2 * m + 3, 2) * r / (1 - r) ** 2 +
                        r * (1 + r) / (2 * (1 - r) ** 3))
    else:
        raise ValueError(p)
    return tail / abs(c) ** p


def reconstruct(p, omega, eps, u0, w0, degree):
    need(degree >= len(p) - 1, "state degree covers forcing polynomial")
    iw = CQ(Q(0), omega)
    w = [w0]
    for j in range(degree):
        aj = p[j] if j < len(p) else Q(0)
        w.append((iw * w[-1] - eps * aj) / (j + 1))
    u = [u0] + [v / (j + 1) for j, v in enumerate(w)]
    residual = [ZERO] * max(len(w), len(p))
    for j, v in enumerate(w):
        residual[j] = residual[j] - iw * v
        if j:
            residual[j - 1] = residual[j - 1] + j * v
    for j, v in enumerate(p):
        residual[j] = residual[j] + eps * v
    need(all(v == ZERO for v in residual[:degree]),
         "all lower ODE-defect coefficients vanish exactly")
    need([j * u[j] for j in range(1, len(u))] == w,
         "U prime equals W exactly")
    return u, w, [v.abs_upper() for v in residual]


def action(p, omega, eps, k, c, h, u0, w0, state_degree,
           geometry_degree, du=Q(0), dw=Q(0), forcing_error=()):
    u, w, defect = reconstruct(p, omega, eps, u0, w0, state_degree)
    q1 = inv_poly(1, geometry_degree, c, h)
    q2 = inv_poly(2, geometry_degree, c, h)
    q3 = inv_poly(3, geometry_degree, c, h)
    aa = [CQ(-3 * v / (k * eps)) for v in q3]
    dd = [CQ(v2 / (k * eps), v1 / eps) for v1, v2 in zip(q1, q2)]
    jp = poly_integral(poly_add(poly_mul(aa, u), poly_mul(dd, w)), h).re
    rho_a = 3 * geom_tail(3, geometry_degree, c, h) / abs(k * eps)
    rho_d = (geom_tail(2, geometry_degree, c, h) / abs(k * eps) +
             geom_tail(1, geometry_degree, c, h) / abs(eps))
    t0 = c - h / 2
    t1 = c + h / 2
    need(t0 * t1 > 0, "cell does not contain geometry pole")
    distance = min(abs(t0), abs(t1))
    apbar = 3 / (abs(k * eps) * distance ** 3) + rho_a
    dpbar = (1 / (abs(k * eps) * distance ** 2) +
             1 / (abs(eps) * distance) + rho_d)
    count = max(len(defect), len(forcing_error))
    cc = [(defect[j] if j < len(defect) else Q(0)) +
          abs(eps) * (forcing_error[j] if j < len(forcing_error) else Q(0))
          for j in range(count)]
    need(all(v >= 0 for v in forcing_error), "nonnegative forcing envelope")
    iwerr = h * dw + sum(v * h ** (j + 2) / ((j + 1) * (j + 2))
                         for j, v in enumerate(cc))
    iuerr = h * du + h ** 2 * dw / 2 + sum(
        v * h ** (j + 3) / ((j + 1) * (j + 2) * (j + 3))
        for j, v in enumerate(cc))
    state_bound = apbar * iuerr + dpbar * iwerr
    count = max(len(p), len(forcing_error))
    gg = [(abs(p[j]) if j < len(p) else Q(0)) +
          (forcing_error[j] if j < len(forcing_error) else Q(0))
          for j in range(count)]
    wtrue = w0.abs_upper() + dw + abs(eps) * sum(
        v * h ** (j + 1) / (j + 1) for j, v in enumerate(gg))
    utrue = u0.abs_upper() + du + h * (w0.abs_upper() + dw) + abs(eps) * sum(
        v * h ** (j + 2) / ((j + 1) * (j + 2)) for j, v in enumerate(gg))
    geometry_bound = h * (rho_a * utrue + rho_d * wtrue)
    ew = dw + sum(v * h ** (j + 1) / (j + 1) for j, v in enumerate(cc))
    eu = du + h * dw + sum(v * h ** (j + 2) / ((j + 1) * (j + 2))
                           for j, v in enumerate(cc))
    return jp, state_bound + geometry_bound, ew, eu


def approximate(q):
    with localcontext() as ctx:
        ctx.prec = 18
        return str(Decimal(q.numerator) / Decimal(q.denominator))


def check_geometry_tails():
    count = 0
    for c in (Q(-2), Q(3)):
        h = Q(1, 2)
        for p in (1, 2, 3):
            for degree in (0, 1, 4, 9):
                polynomial = inv_poly(p, degree, c, h)
                radius = geom_tail(p, degree, c, h)
                for j in range(41):
                    x = h * j / 40
                    exact = 1 / (c + x - h / 2) ** p
                    candidate = sum(v * x ** n for n, v in enumerate(polynomial))
                    need(abs(exact - candidate) <= radius, "geometry tail fixture")
                    count += 1
    return count


def check_source_disk():
    r = Q(5, 8)
    invd = Q(64, 39)
    bb = Q(4, 3)
    ll = Q(8, 27)
    exponent = Q(25, 89)
    exp_upper = 1 + exponent + exponent ** 2 / (2 * (1 - exponent / 3))
    need(exp_upper == Q(57051, 43076) < bb, "rational exponential upper bound")
    b1 = 2 * r * bb * invd ** 2
    b2 = bb * (4 * r ** 2 * invd ** 4 + 2 * invd ** 2 +
               8 * r ** 2 * invd ** 3)
    even = 4 * ll ** 2 * bb + 2 * ll * b1 + b2
    signed = 4 * ll ** 2 * r * bb + 2 * ll * (bb + r * b1) + 2 * b1 + r * b2
    need(even == Q(2737816576, 62462907) < 44 < 64, "even analytic source bound")
    need(signed == Q(2321190208, 62462907) < 38 < 64, "signed analytic source bound")
    s = Q(25, 64)
    count = 0
    for xj in range(-16, 17):
        for yj in range(-16, 17):
            vx = s * xj / 16
            vy = s * yj / 16
            if vx * vx + vy * vy <= s * s:
                real_inverse = (1 - vx) / ((1 - vx) ** 2 + vy ** 2)
                need(real_inverse >= Q(64, 89), "complex disk rational fixture")
                count += 1
    ratio = Q(1, 16)
    cauchy = 64 * ratio ** 25 / (1 - ratio)
    need(cauchy == Q(1, 15 * 2 ** 90), "degree24 Cauchy remainder")
    return {"fabricated_disk_points": count,
            "even_bound_exact": str(even), "signed_bound_exact": str(signed),
            "degree24_cauchy_remainder_exact": str(cauchy)}


def check_integrals():
    p = [Q(1, 7), Q(-2, 5), Q(3, 11)]
    eps = Q(2, 3)
    k = Q(5, 7)
    u0 = CQ(Q(2, 9), Q(-1, 13))
    w0 = CQ(Q(-1, 5), Q(3, 17))
    fixtures = [("zero_phase", Q(0), Q(-2), Q(1, 4), 24),
                ("low_phase", Q(3), Q(-2), Q(1, 4), 24),
                ("high_phase", Q(256), Q(-2), Q(1, 16), 64)]
    results = []
    for name, omega, center, h, degree in fixtures:
        prod, bound, ew, eu = action(p, omega, eps, k, center, h, u0, w0,
                                    degree, 12)
        ref, ref_bound, ewr, eur = action(p, omega, eps, k, center, h, u0, w0,
                                         144, 48)
        need(abs(prod - ref) + ref_bound <= bound,
             name + " high-order reference interval is inside production enclosure")
        up, wp, _ = reconstruct(p, omega, eps, u0, w0, degree)
        ur, wr, _ = reconstruct(p, omega, eps, u0, w0, 144)
        interval_disk_inside(poly_value(up, h) - poly_value(ur, h), eur, eu,
                             name + " inherited-u carry contains endpoint reference disk")
        interval_disk_inside(poly_value(wp, h) - poly_value(wr, h), ewr, ew,
                             name + " inherited-w carry contains endpoint reference disk")
        results.append({"name": name, "phase_exact": str(abs(omega) * h),
                        "local_bound_approx": approximate(bound),
                        "reference_bound_approx": approximate(ref_bound),
                        "reference_interval_inside": True,
                        "endpoint_reference_disks_inside": True})
    # A fabricated linear forcing error models an inherited nested moment.
    omega = Q(7)
    center = Q(3)
    h = Q(1, 8)
    du = Q(1, 1000000)
    dw = Q(1, 700000)
    error = [Q(1, 9000000), Q(1, 123000)]
    prod, bound, ew, eu = action(p, omega, eps, k, center, h, u0, w0,
                                32, 16, du, dw, error)
    pp = [p[j] + (error[j] if j < len(error) else 0) for j in range(len(p))]
    for su in (-1, 1):
        for swr in (-1, 1):
            for swi in (-1, 1):
                actual_u0 = u0 + CQ(su * du)
                actual_w0 = w0 + CQ(swr * dw / 2, swi * dw / 2)
                ref, rb, ewr, eur = action(pp, omega, eps, k, center, h,
                                           actual_u0, actual_w0, 144, 48)
                need(abs(prod - ref) + rb <= bound,
                     "inherited mode/nested forcing reference interval contained")
                up, wp, _ = reconstruct(p, omega, eps, u0, w0, 32)
                ur, wr, _ = reconstruct(pp, omega, eps, actual_u0, actual_w0, 144)
                interval_disk_inside(poly_value(up, h) - poly_value(ur, h), eur, eu,
                                     "uncertain-mode u carry contains endpoint reference disk")
                interval_disk_inside(poly_value(wp, h) - poly_value(wr, h), ewr, ew,
                                     "uncertain-mode w carry contains endpoint reference disk")
    need(ew >= dw and eu >= du + h * dw,
         "inherited cell uncertainty is retained in endpoint carry")
    results.append({"name": "inherited_modes_and_nested_moment", "corner_count": 8,
                    "local_bound_approx": approximate(bound),
                    "endpoint_w_radius_approx": approximate(ew),
                    "endpoint_u_radius_approx": approximate(eu),
                    "reference_intervals_inside": True,
                    "endpoint_reference_disks_inside": True})
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"status": "PASS_FABRICATED_DIRECT_INTEGRATION_BOUND_CHECKS",
              "arithmetic": "exact Python fractions; decimal strings are display only",
              "physical_source_calls": 0, "scientific_input_arrays": 0,
              "geometry_tail_point_checks": check_geometry_tails(),
              "analytic_source_inequalities": check_source_disk(),
              "integral_fixtures": check_integrals()}
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
