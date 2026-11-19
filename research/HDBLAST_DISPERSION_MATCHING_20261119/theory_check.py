"""Checks for the finite-band calibrated inverse response and held-out bound."""
import json
import math
from fractions import Fraction

from scipy.integrate import quad
from scipy.optimize import brentq

NORMALIZATION = 8 * math.pi**2


def divided_response(q2, lo, hi, a, b, x):
    """Positive three-denominator spectral integral; unit squared coupling."""
    def density(t):
        return (t - q2) / ((t + a) * (t + b) * (t + x)) / NORMALIZATION
    return quad(density, lo, hi, epsabs=3e-16, epsrel=1e-12)[0]


def raw_spectral(q2, upper, x):
    # Closed-form response from the unsubtracted Dirac spectral integral.
    return (upper - q2 - (q2 + x) * math.log((upper + x) / (q2 + x))) / NORMALIZATION


def secant_tests():
    q2, upper = 16., 400.
    baseline, bare_slope = 1.7, 0.8
    def inverse(x):
        return baseline + bare_slope*x - raw_spectral(q2, upper, x)
    a, b = 0., 4.
    slope = (inverse(b) - inverse(a)) / (b - a)
    checks = {}
    for x in (2., 9., 25.):
        line = inverse(a) + (x - a)*slope
        integral = divided_response(q2, q2, upper, a, b, x)
        if not math.isclose(line - inverse(x), (x-a)*(x-b)*integral,
                            rel_tol=2e-10, abs_tol=2e-14):
            raise ArithmeticError("Secant identity")
        checks[str(x)] = {"inverse_minus_secant": inverse(x)-line,
                          "positive_curvature_residual": integral}
    # Static inverse and kinetic slope, independently imposing a repeated anchor.
    static_shift = (math.log(upper/q2)-1.+q2/upper)/NORMALIZATION
    x = 25.
    measured_slope = bare_slope + static_shift
    analytic = inverse(0.) + measured_slope*x - inverse(x)
    curvature = x*x*divided_response(q2, q2, upper, 0., 0., x)
    if not math.isclose(analytic, curvature, rel_tol=1e-12):
        raise ArithmeticError("Static-plus-slope identity")
    return {"distinct": checks, "static_plus_slope_error": analytic-curvature}


def tail_examples():
    q2, a, b, c, x = 16., 0., 4., 9., 25.
    rows = {}
    for cutoff in (12., 20., 40., 80., 160.):
        lo = cutoff**2
        measured = divided_response(q2, lo, math.inf, a, b, c)
        held_out = divided_response(q2, lo, math.inf, a, b, x)
        low, high = sorted([1., (lo+c)/(lo+x)])
        ratio = held_out / measured
        bound = (1.-q2/lo)/(1.+a/lo)/(1.+b/lo)/(1.+x/lo)
        normalized = NORMALIZATION*lo*held_out
        if not (low < ratio < high and bound < normalized < 1):
            raise ArithmeticError("Positive unknown-tail and Dirac-tail bounds")
        rows[str(int(cutoff))] = {"ratio_lower":low, "ratio":ratio, "ratio_upper":high,
                                  "dirac_scaled_tail":normalized, "dirac_scaled_tail_lower":bound,
                                  "known_missing_band_inverse_shift":(x-a)*(x-b)*held_out}
    return rows


def atomic_test():
    # Use exact rational arithmetic and allow a completely different spectrum.
    t1, t2, a, b, c, x = map(Fraction, (36, 400, 0, 4, 9, 25))
    def residual(t, w, z):
        return w/((t+a)*(t+b)*(t+z))
    vals = []
    for weight in (Fraction(0), Fraction(2,3), Fraction(1)):
        first = residual(t1, weight, c) + residual(t2, 1-weight, c)
        second = residual(t1, weight, x) + residual(t2, 1-weight, x)
        vals.append(second/first)
    if not (vals[0] == (t2+c)/(t2+x) and vals[2] == (t1+c)/(t1+x)
            and vals[2] < vals[1] < vals[0]):
        raise ArithmeticError("Sharp two-atom bound")
    return [str(v) for v in vals]


def kinetic_test():
    # A finite surrogate example, not a renormalization condition imposed on the universe.
    q2, slope, strength = 16., 0.03, 1.
    def added_slope(L):
        return strength*(math.log(L/q2)-1.+q2/L)/NORMALIZATION
    critical = brentq(lambda L: added_slope(L)-slope, q2, 1e4)
    test = 2*critical
    bare_slope = slope - added_slope(test)
    def inverse(x):
        return 1 + bare_slope*x - raw_spectral(q2,test,x) + raw_spectral(q2,test,0)
    pole = brentq(inverse, 0, 1e6)
    if not (bare_slope < 0 and abs(inverse(pole)) < 1e-10 and slope > 0):
        raise ArithmeticError("Positive-slope realizability ceiling")
    return {"static_inverse":1, "measured_kinetic_slope":slope,"squared_coupling":strength,
            "q":math.sqrt(q2), "cutoff_ceiling":math.sqrt(critical),
            "example_cutoff":math.sqrt(test),"negative_bare_slope":bare_slope,
            "unstable_euclidean_zero":pole}


if __name__ == '__main__':
    report = {"secant":secant_tests(),"tail":tail_examples(),
              "endpoint_saturation_exact":atomic_test(),"canonical_kinetic":kinetic_test()}
    print(json.dumps(report,indent=2,allow_nan=False))
