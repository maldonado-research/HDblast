#!/usr/bin/env python3
"""Planning-only rational shape feasibility gate; no waveform producer is run.

Inputs ell and eta come from the parent majorant.  This independent calculation
checks their consequences for the analytic bump, not the PDE majorant itself.
All checks remain active under python -O.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import isqrt
import json


class GateFailure(ValueError):
    pass


def require(condition, reason):
    if not condition:
        raise GateFailure(reason)


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __init__(self, lo, hi=None):
        object.__setattr__(self, "lo", F(lo))
        object.__setattr__(self, "hi", F(lo if hi is None else hi))
        require(self.lo <= self.hi, "reversed interval")

    def __add__(self, other):
        other = box(other)
        return Interval(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -box(other)

    def __rsub__(self, other):
        return box(other) + -self

    def __mul__(self, other):
        other = box(other)
        values = [a*b for a in (self.lo, self.hi)
                  for b in (other.lo, other.hi)]
        return Interval(min(values), max(values))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = box(other)
        require(other.lo > 0 or other.hi < 0, "zero in interval divisor")
        return self*Interval(1/other.hi, 1/other.lo)

    def __rtruediv__(self, other):
        return box(other)/self


def box(x):
    return x if isinstance(x, Interval) else Interval(x)


def sqrt_box(x):
    x = box(x)
    require(x.lo >= 0, "negative square root")
    scale = 10**30
    lo = F(isqrt(x.lo.numerator*scale**2//x.lo.denominator), scale)
    hi = F(isqrt(x.hi.numerator*scale**2//x.hi.denominator)+1, scale)
    require(lo*lo <= x.lo and x.hi < hi*hi, "square root certificate")
    return Interval(lo, hi)


def exp_minus_point(t):
    require(0 <= t <= 2, "exponential argument outside certified range")
    term = total = F(1)
    for k in range(1, 49):
        term *= -t/k
        total += term
    # At truncation 48 the tail terms decrease; S49 <= exp(-t) <= S48.
    return Interval(total+term*(-t/49), total)


def exp_minus(t):
    t = box(t)
    return Interval(exp_minus_point(t.hi).lo, exp_minus_point(t.lo).hi)


def qprime(x):
    x = F(x)
    z = 2*x-3
    d = 1-z*z
    return (-16*z/(d*d))*exp_minus(4*z*z/d)


def require_enclosure(value, lo, hi):
    require(F(lo) < value.lo and value.hi < F(hi), "target enclosure failed")


def reject(name, operation):
    try:
        operation()
    except GateFailure:
        return {"name": name, "rejected": True}
    raise GateFailure("negative control accepted: "+name)


def decimal_integer(n, digits):
    sign = "-" if n < 0 else ""
    text = str(abs(n)).zfill(digits+1)
    return sign+text[:-digits]+"."+text[-digits:]


def outward(value, digits=18):
    value = box(value)
    scale = 10**digits
    lo = value.lo.numerator*scale//value.lo.denominator
    hi = (value.hi.numerator*scale+value.hi.denominator-1)//value.hi.denominator
    return [decimal_integer(lo, digits), decimal_integer(hi, digits)]


def main():
    ell = F(1, 10**6)
    eta_exact = F(2007523995553919, 1999981973040500000)
    eta = F(251, 250000)
    require(0 < eta_exact < eta, "parent majorant rounding")

    s = sqrt_box(3)
    a = sqrt_box(2/s-1)
    m = 6*(2+s)*a*exp_minus(s-1)
    require_enclosure(m, "4.23564042213488", "4.23564042213490")
    x_left = (3-a)/2
    x_right = (3+a)/2
    require_enclosure(x_left, "1.29", "1.32")
    require_enclosure(x_right, "1.68", "1.71")

    signed_sample = qprime(F(5, 4))
    require_enclosure(signed_sample, "3.74893707542366", "3.74893707542368")
    require(qprime(F(7, 4)).hi < F("-3.7"), "negative Qprime sample")
    positive_sample = signed_sample+Interval(-eta, ell+eta)
    negative_sample = -signed_sample+Interval(-eta, ell+eta)
    require(positive_sample.lo > 0 and negative_sample.hi < 0,
            "actual waveform sign gate")

    gap_left = m-qprime(F("1.29"))
    gap_right = m-qprime(F("1.32"))
    gap_lower = min(gap_left.lo, gap_right.lo)
    require(gap_lower > F(1, 50), "simple rational endpoint gap")
    require(gap_lower > 2*eta+ell, "actual peak-window exclusion")
    require(F(1, 50) > 2*eta+ell, "coarse rational planning gate")

    max_height = m+Interval(-eta, ell+eta)
    min_height = -m+Interval(-eta, ell+eta)
    require_enclosure(max_height, "4.2346", "4.2367")
    require_enclosure(min_height, "-4.2367", "-4.2346")

    result = {
        "status": "PASS",
        "scope": "planning feasibility only; imported majorant not rederived",
        "arithmetic": "exact Fraction intervals and alternating Taylor bounds",
        "ell": str(ell), "eta_exact": str(eta_exact), "eta_rounded": str(eta),
        "enclosures": {
            "M": outward(m),
            "Qprime_positive_peak_location": outward(x_left),
            "Qprime_negative_peak_location": outward(x_right),
            "Qprime(5/4)": outward(signed_sample),
            "actual_H(5/4)": outward(positive_sample),
            "actual_H(7/4)": outward(negative_sample),
            "M-Qprime(1.29)": outward(gap_left),
            "M-Qprime(1.32)": outward(gap_right),
            "actual_global_maximum_height": outward(max_height),
            "actual_global_minimum_height": outward(min_height),
            "actual_supremum_absolute_height": outward(max_height),
        },
        "coarse_gate": {
            "certified_endpoint_gap_greater_than": "1/50",
            "required_endpoint_gap_greater_than": str(2*eta+ell),
            "simple_eta_requirement": "eta < 19999/2000000 = 0.0099995",
            "positive_global_maximum_window": ["1.29", "1.32"],
            "negative_global_minimum_window": ["1.68", "1.71"],
            "window_endpoints_excluded": True,
            "absolute_global_maximum_window": "union of the two windows",
        },
        "limitations": {
            "actual_extrema_uniqueness_certified": False,
            "actual_zero_uniqueness_certified": False,
            "absolute_maximum_lobe_identified": False,
        },
    }
    result["negative_controls"] = [
        reject("missing derivative factor 16", lambda: require_enclosure(
            m/16, "4.23564042213488", "4.23564042213490")),
        reject("reversed derivative sign", lambda: require(
            (-signed_sample).lo > F("3.7"), "wrong signed sample")),
        reject("excessive eta accepted for same windows", lambda: require(
            gap_lower > F(1, 5)+ell, "wrong exclusion gate")),
    ]
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
