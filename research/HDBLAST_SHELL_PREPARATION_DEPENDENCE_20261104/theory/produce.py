#!/usr/bin/env python3
"""Independent exact-rational sufficient certificate; no scientific imports/oracles.

This reuses two pinned published background/response artifacts. It verifies broad
background bounds, then uses analytic bump derivative bounds on the causal
initial segment only. It is an independent method, not an independent PDE solve.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
from math import factorial
import sys


class CertificateRejected(Exception):
    pass


class QuietParser(argparse.ArgumentParser):
    def error(self, message):
        raise CertificateRejected("invalid arguments")


def need(ok, message):
    if not ok:
        raise CertificateRejected(message)


def box(value):
    return tuple(map(F, value["exact"]))


def subbox(value, lower, upper, label):
    lo, hi = box(value)
    need(lower <= lo <= hi <= upper, "background enclosure: " + label)


def record(value):
    if isinstance(value, tuple):
        return {"exact": list(map(str, value)), "display_only": list(map(float, value))}
    return {"exact": str(value), "display_only": float(value)}


def pinned_read(root, relative, expected):
    try:
        payload = (root / relative).read_bytes()
    except OSError:
        raise CertificateRejected("source unreadable: " + relative) from None
    need(hashlib.sha256(payload).hexdigest() == expected, "source hash: " + relative)
    return json.loads(payload)


def polynomial_checks():
    """Exact polynomial identities used in the continuous analytic bump proof."""
    def trim(poly):
        while len(poly) > 1 and poly[-1] == 0:
            poly.pop()
        return poly
    def add(left, right):
        return trim([(left[j] if j < len(left) else F())
                     + (right[j] if j < len(right) else F())
                     for j in range(max(len(left), len(right)))])
    def scale(poly, factor):
        return trim([factor * value for value in poly])
    def mul(left, right):
        out = [F()] * (len(left) + len(right) - 1)
        for j, u in enumerate(left):
            for k, v in enumerate(right):
                out[j + k] += u * v
        return trim(out)
    def derivative(poly):
        return [j * poly[j] for j in range(1, len(poly))] or [F()]
    def evaluate(poly, value):
        return sum((coefficient * value**j for j, coefficient in enumerate(poly)), F())
    # D'=1-2t, (D')^2=1-4D, D''=-2. These identities reduce Q''/Q.
    P2 = add(mul([F(1), F(-4)], [F(1), F(-2)]), [F(), F(), F(-2)])
    need(P2 == [1, -6, 6], "Q second derivative polynomial")
    # Differentiate Q P2(D)/D^4 and factor Q D'/D^6.
    P3 = add(add(P2, mul([F(), F(), F(1)], derivative(P2))),
             scale(mul([F(), F(1)], P2), -4))
    need(P3 == [1, -10, 24, -12], "Q third derivative polynomial")
    endpoint = F(3, 16)
    need(evaluate(P2, endpoint) == F(11, 128), "positive P2 endpoint")
    need(evaluate(P3, endpoint) == -F(113, 1024), "P3 endpoint")
    # Upper bounds for their derivatives on [0,3/16], dropping negative terms.
    need(-6 + 12 * endpoint < 0, "P2 is decreasing")
    need(-10 + 48 * endpoint < 0, "P3 is decreasing")
    return {"P2_coefficients": list(map(str, P2)), "P3_coefficients": list(map(str, P3)),
            "P2_endpoint": "11/128", "P3_endpoint": "-113/1024",
            "continuous_argument": "D in [0,3/16]; P2 decreases from 1 to 11/128; P3 decreases from 1 to -113/1024."}


def main():
    parser = QuietParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    need(args.output.is_absolute(), "output path must be absolute")
    need(str(args.output) == str(args.output.resolve()), "output path must be normalized")
    need(args.output.parent.is_dir(), "output parent missing")
    need(not args.output.exists() and not args.output.is_symlink(), "output already exists")
    static_path = "research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json"
    response_path = "research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json"
    static_sha = "3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22"
    response_sha = "1a9ebffb6bee9332bc0fd0df0f5d1d9f7f2a796825c2a563e4c644974ba4d8a3"
    st = pinned_read(args.repo_root, static_path, static_sha)
    old = pinned_read(args.repo_root, response_path, response_sha)
    bump_identity_checks = polynomial_checks()
    rho_lo, rho_hi = F(129), F(131)
    p_lo, p_hi = F(3, 25000), F(1, 7000)
    a_hi, ay_hi, py_hi = F(14, 125), F(1, 10000), F(1, 4000)
    subbox(st["y_jets"]["rho"][0], rho_lo, rho_hi, "rho")
    subbox(st["y_jets"]["p"][0], -p_hi, -p_lo, "p")
    subbox(st["y_jets"]["p"][1], -py_hi, py_hi, "p_y")
    subbox(st["y_jets"]["a"][0], F(111, 1000), a_hi, "a")
    subbox(st["y_jets"]["a"][1], -ay_hi, ay_hi, "a_y")
    subbox(old["first"]["rho"], rho_lo, rho_hi, "rho_b")
    subbox(old["first"]["Wprime"], -F(1, 1000), -F(1, 1200), "F_prime")

    ell, target, delta = F(1, 10**6), F(5, 4), F(1, 4)
    need(F(st["domain"]["t"]) == 0 and F(st["domain"]["ell"]) == ell,
         "static initial slice and scale")
    need(F(st["domain"]["z_min"]) <= ell and F(st["domain"]["z_max"]) >= ell * target,
         "static continuous interval covers the causal initial segment")
    q, eta = F(1201, 3200000), F(1, 10**11)
    need(F(old["q"]["exact"]) == q, "retained contraction q")
    # Analytic Taylor bounds for exp(4/3): 15/4 < exp(4/3) < 4.
    x = F(4, 3)
    lower_exp = sum((x**j / factorial(j) for j in range(6)), F())
    upper_exp = lower_exp + x**6 / factorial(6) / (1 - x / 7)
    need(lower_exp > F(15, 4) and upper_exp < 4, "bump exponential enclosure")
    # exp(2) > 6^6/6400, proving exp(-2)*6^6 < 6400.
    lower_exp2 = sum((F(2)**j / factorial(j) for j in range(7)), F())
    need(lower_exp2 > F(6**6, 6400), "third bump derivative exponential bound")
    need(F(4, 15) * F(16, 3)**2 < 8, "first bump derivative bound")
    need(F(4, 15) * F(16, 3)**4 < 216, "second bump derivative bound")
    # The analytic proof in THEOREM.md establishes these whole-segment bounds.
    Q0, Q1, Q2, Q3 = F(1, 3), F(8), F(216), F(6400)
    R_hi = (rho_hi / rho_lo)**3
    b_hi = ell * rho_hi * a_hi
    d_hi = ell * rho_hi * py_hi / p_lo
    c_lo = 3 * rho_lo**3 / (2 * ell * rho_hi**4 * p_hi)
    c_hi = 3 * rho_hi**3 / (2 * ell * rho_lo**4 * p_lo)
    h_hi = R_hi * Q1
    hr_hi = R_hi * (Q2 + 3 * b_hi * Q1)
    velocity_hi = R_hi * (Q2 + 4 * b_hi * Q1)
    bracket_hi = Q2 + 4 * b_hi * Q1 + 3 * ell**2 * Q0
    s_hi = c_hi * bracket_hi
    sr_hi = c_hi * (
        Q3 + 4 * b_hi * Q2 + 3 * ell**2 * Q1
        + (4 * b_hi + d_hi) * bracket_hi
        + 4 * ell**2 * rho_hi**2 * (a_hi**2 + ay_hi) * Q1
    )
    # Triangle inequality for H=(1,2,1,3,1), deliberately avoiding square roots.
    D0 = 2 * h_hi + s_hi
    Dr = 2 * hr_hi + sr_hi
    Dv = 2 * velocity_hi  # sqrt(2) <= 2
    N = max((D0 + target * Dv) / 3, Dr + Dv)
    E = q * N / (1 - q)
    scalar_lo = c_lo * F(1408, 81) - delta * E
    scalar_hi = s_hi + delta * E
    # s>0; exact kernel ratios >=1-eta; endpoint primitive moments of Q''.
    volume_history_lo = (1 - eta) * c_lo / 4 - delta**3 * E / 6
    expansion_history_lo = (1 - eta) * c_lo * F(32, 9) - delta**2 * E / 2
    volume_history_hi = delta**2 * s_hi / 2 + delta**3 * E / 6
    expansion_history_hi = delta * s_hi + delta**2 * E / 2
    need(min(scalar_lo, volume_history_lo, expansion_history_lo) > 0,
         "all three homogeneous difference sign gates")
    abs_v_lo = (rho_lo * ell)**2 * F(1, 2400) * volume_history_lo
    abs_v_hi = (rho_hi * ell)**2 * F(1, 2000) * volume_history_hi
    abs_t_lo = rho_lo * ell * F(1, 2400) * expansion_history_lo
    abs_t_hi = rho_hi * ell * F(1, 2000) * expansion_history_hi
    slope_M = (scalar_lo, scalar_hi)
    slope_v = (-abs_v_hi, -abs_v_lo)
    slope_t = (-abs_t_hi, -abs_t_lo)
    base_v = (F("-5.4520e16"), F("-2.6521e16"))
    base_t = (F("-9.0208e21"), F("-4.4458e21"))
    subbox(old["proper_v2"], *base_v, "base v2 headline")
    subbox(old["proper_Theta2"], *base_t, "base Theta2 headline")
    threshold = max(-base_v[0] / abs_v_lo, -base_t[0] / abs_t_lo)
    exponent = 0
    while F(10)**exponent <= threshold:
        exponent += 1
    scale = F(10)**exponent
    new_v = (base_v[0] + scale * abs_v_lo, base_v[1] + scale * abs_v_hi)
    new_t = (base_t[0] + scale * abs_t_lo, base_t[1] + scale * abs_t_hi)
    need(new_v[0] > 0 and new_t[0] > 0, "simultaneous reversal")
    result = {
        "status": "PASS_COARSE_ANALYTIC_PREPARATION_SIGN_CERTIFICATE",
        "scope": "Conditional theorem premises retained; exact rational sufficient bounds; no new oracle or PDE observation.",
        "profile": "chi=ell*rho_b^3*Q(r), h=-(rho_b/rho)^3*Q_r(r)",
        "initial_segment": ["1", "5/4"],
        "retained_source_pins": {static_path: static_sha, response_path: response_sha},
        "constants": {"q": record(q), "eta": record(eta), "c_lo": record(c_lo), "c_hi": record(c_hi)},
        "analytic_bump_caps": {"Q": "1/3", "Q_r": "8", "Q_rr": "216", "Q_rrr": "6400"},
        "bump_polynomial_identity_checks": bump_identity_checks,
        "initial_norm_caps": {"D0": record(D0), "Dr": record(Dr), "Dv": record(Dv), "N": record(N)},
        "full_resolvent_error": record(E),
        "positive_history_lower_bounds": {"v": record(volume_history_lo), "Theta": record(expansion_history_lo)},
        "slopes": {"M2": record(slope_M), "v2": record(slope_v), "Theta2": record(slope_t)},
        "deterministic_choice": {"threshold": record(threshold), "n": exponent, "lambda": str(-scale)},
        "new_coefficients": {"v2": record(new_v), "Theta2": record(new_t)},
        "limits": {"effective_epsilon_radius": False, "physical_preparation_selection": False,
                   "full_angular_waveform": False, "independent_PDE_observation": False},
    }
    payload = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode()
    with args.output.open("xb") as stream:
        stream.write(payload)
    print(json.dumps({"status": result["status"], "output_sha256": hashlib.sha256(payload).hexdigest()},
                     sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    try:
        main()
    except CertificateRejected as error:
        print(json.dumps({"status": "REJECTED", "reason": str(error)}, sort_keys=True, separators=(",", ":")))
        sys.exit(2)
    except (OSError, RuntimeError):
        print('{"reason":"filesystem error","status":"REJECTED"}')
        sys.exit(2)
