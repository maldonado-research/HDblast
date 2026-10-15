#!/usr/bin/env python3
"""Certify bump-waveform peak constants using Python exact rational arithmetic.

All arithmetic affecting verification is fractions.Fraction arithmetic.  Decimal
strings are exact rationals.  Exponentials are enclosed by alternating Taylor
polynomials S_81 <= exp(-u) <= S_80 for 0 <= u <= 2.  Square roots use integer
square roots after exact rational scaling.  No floating point is used.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import isqrt
import json


class CertificateFailure(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise CertificateFailure(message)


@dataclass(frozen=True)
class I:
    lo: F
    hi: F

    def __init__(self, lo, hi=None):
        object.__setattr__(self, "lo", F(lo))
        object.__setattr__(self, "hi", F(lo if hi is None else hi))
        require(self.lo <= self.hi, "reversed interval")

    def __add__(self, other):
        other = interval(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -interval(other)

    def __rsub__(self, other):
        return interval(other) + -self

    def __mul__(self, other):
        other = interval(other)
        products = [a * b for a in (self.lo, self.hi)
                    for b in (other.lo, other.hi)]
        return I(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = interval(other)
        require(other.lo > 0 or other.hi < 0, "interval divisor contains zero")
        return self * I(1 / other.hi, 1 / other.lo)

    def __rtruediv__(self, other):
        return interval(other) / self

    def __pow__(self, n):
        require(isinstance(n, int) and n >= 0, "invalid interval exponent")
        out = I(1)
        for _ in range(n):
            out *= self
        return out


def interval(x):
    return x if isinstance(x, I) else I(x)


def sqrt_interval(x, digits=50):
    x = interval(x)
    require(x.lo >= 0, "negative square-root argument")
    scale = 10**digits
    a = isqrt(x.lo.numerator * scale**2 // x.lo.denominator)
    b = isqrt(x.hi.numerator * scale**2 // x.hi.denominator) + 1
    require(F(a, scale)**2 <= x.lo and F(b, scale)**2 > x.hi,
            "square-root endpoint verification")
    return I(F(a, scale), F(b, scale))


def exp_neg_point(u):
    require(0 <= u <= 2, "exponential argument outside proved range")
    term = total = F(1)
    for k in range(1, 81):
        term *= -u / k
        total += term
    upper = total
    lower = upper + term * (-u / 81)
    return I(lower, upper)


def exp_neg(x):
    x = interval(x)
    return I(exp_neg_point(x.hi).lo, exp_neg_point(x.lo).hi)


ELL = F(1, 10**10)
ACTUAL_ERROR = F(8000003000800720779892143999973,
                 79999999999992799998920000162000000000)
ROUNDED_ERROR = F(1000001, 10**13)


def qprime(x):
    x = interval(x)
    z = 2*x - 3
    d = 1 - z**2
    return -16*z*exp_neg(4*z**2/d)/d**2


def waveform(x):
    x = interval(x)
    z = 2*x - 3
    d = 1 - z**2
    q = exp_neg(4*z**2/d)
    return q * (-16*z/d**2 + ELL*exp_neg(ELL*x))


def critical_polynomial(x, correction_sign=1):
    """Hhat'/Q = -32 P(z)/(1-z^2)^4, z=2x-3."""
    x = interval(x)
    z = 2*x - 3
    d = 1 - z**2
    a = ELL*exp_neg(ELL*x)
    return (1 - 6*z**2 - 3*z**4
            + correction_sign*(a*z*d**2/2 + a*ELL*d**4/32))


def require_enclosure(name, value, lo, hi):
    require(F(lo) < value.lo and value.hi < F(hi), name)


def require_rejection(name, operation):
    try:
        operation()
    except CertificateFailure:
        return {"name": name, "rejected": True}
    raise CertificateFailure("negative control was accepted: " + name)


def main():
    result = {"status": "PASS", "arithmetic": "exact rational intervals",
              "enclosures": {}, "negative_controls": []}

    def certify(name, value, lo, hi):
        require_enclosure(name, value, lo, hi)
        result["enclosures"][name] = [lo, hi]

    s = sqrt_interval(3)
    r = 2/s - 1
    a = sqrt_interval(r)
    x_left = (3-a)/2
    x_right = (3+a)/2
    m = 6*(2+s)*a*exp_neg(s-1)
    certify("Qprime positive peak location", x_left,
                     "1.303340053404835680417725165684",
                     "1.303340053404835680417725165685")
    certify("Qprime negative peak location", x_right,
                     "1.696659946595164319582274834315",
                     "1.696659946595164319582274834316")
    certify("M = max |Qprime|", m,
                     "4.235640422134888218565987509535",
                     "4.235640422134888218565987509536")

    # These are the two unique critical points by the analytic argument in
    # peak_analysis.txt.  Sign changes below certify the rational root boxes.
    left = I("1.303340053406124851570706028048",
             "1.303340053406124851570706028049")
    right = I("1.696659946596453490735228079422",
              "1.696659946596453490735228079423")
    require(critical_polynomial(left.lo).hi < 0, "left peak lower sign")
    require(critical_polynomial(left.hi).lo > 0, "left peak upper sign")
    require(critical_polynomial(right.lo).lo > 0, "right peak lower sign")
    require(critical_polynomial(right.hi).hi < 0, "right peak upper sign")
    result["enclosures"]["Hhat positive peak location"] = [
        "1.303340053406124851570706028048",
        "1.303340053406124851570706028049"]
    result["enclosures"]["Hhat negative peak location"] = [
        "1.696659946596453490735228079422",
        "1.696659946596453490735228079423"]
    certify("Hhat positive peak height", waveform(left),
                     "4.235640422182980388580255",
                     "4.235640422182980388580257")
    certify("Hhat negative peak height", waveform(right),
                     "-4.235640422086796048554158",
                     "-4.235640422086796048554156")
    certify("Hhat(5/4)", waveform(F(5, 4)),
                     "3.748937075450029332710093",
                     "3.748937075450029332710095")
    certify("Hhat(7/4)", waveform(F(7, 4)),
                     "-3.748937075397309905094858",
                     "-3.748937075397309905094855")
    require(waveform(F(3, 2)).lo > 0, "zero bracket lower sign")
    require(waveform(F(3, 2) + ELL/32).hi < 0, "zero bracket upper sign")
    result["enclosures"]["Unique interior zero"] = [
        "1.5", "1.500000000003125"]

    require(ACTUAL_ERROR < ROUNDED_ERROR, "actual error rounding")
    gap_left = m-qprime(F("1.3033"))
    gap_right = m-qprime(F("1.3034"))
    certify("M-Qprime(1.3033)", gap_left,
            "0.000000263559678293613961962519",
            "0.000000263559678293613961962520")
    certify("M-Qprime(1.3034)", gap_right,
            "0.000000590303334992119176697884",
            "0.000000590303334992119176697885")
    require(min(gap_left.lo, gap_right.lo) > 2*ROUNDED_ERROR+ELL,
            "actual-waveform peak exclusion gap")
    result["actual_waveform_corollary"] = {
        "assumption": "continuous H with global |H-Hhat| <= ACTUAL_ERROR",
        "positive_maximum_location": ["1.3033", "1.3034"],
        "negative_minimum_location": ["1.6966", "1.6967"],
        "absolute_maximum_location": "union of the preceding intervals",
        "absolute_maximum_lobe_identified": False,
        "critical_point_uniqueness_certified": False,
        "actual_zero_uniqueness_certified": False,
    }
    certify("Actual waveform maximum height, conditional",
            waveform(left)+I(-ROUNDED_ERROR, ROUNDED_ERROR),
            "4.235640322182880388580255",
            "4.235640522183080388580257")
    certify("Actual waveform minimum height, conditional",
            waveform(right)+I(-ROUNDED_ERROR, ROUNDED_ERROR),
            "-4.235640522086896048554158",
            "-4.235640322086696048554156")

    result["negative_controls"].append(require_rejection(
        "missing factor 16 in Qprime", lambda: require_enclosure(
            "wrong M", m/16, "4.235640422134888218565987509535",
            "4.235640422134888218565987509536")))
    result["negative_controls"].append(require_rejection(
        "reversed Qprime sign", lambda: require(
            (-qprime(F(5, 4))).lo > 0, "wrong signed sample")))
    result["negative_controls"].append(require_rejection(
        "reversed Hhat correction sign", lambda: require(
            critical_polynomial(left.lo, -1).hi < 0,
            "wrong correction root bracket")))
    result["negative_controls"].append(require_rejection(
        "inflated actual waveform error", lambda: require(
            min(gap_left.lo, gap_right.lo) > 20*ROUNDED_ERROR+ELL,
            "wrong actual-waveform exclusion")))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
