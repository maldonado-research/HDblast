"""Rational enclosures for the two-window 4D forced Klein–Gordon comparison."""
import json
import math
from fractions import Fraction as F


ELL = F(1, 10**6)
HORIZON = F(25, 16)
SOURCE_CENTER = F(129, 128)
SOURCE_RADIUS = F(1, 128)
WINDOWS = [(F(5, 4), F(11, 8)), (F(49, 32), F(25, 16))]


def require(condition, problem):
    if not condition:
        raise ValueError(problem)


def sine(x):
    """Taylor polynomial of degree 23, with Lagrange remainder."""
    center = sum((-1)**k * x**(2*k+1) / math.factorial(2*k+1) for k in range(12))
    error = abs(x)**24 / math.factorial(24)
    return center-error, center+error


def product(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def outward_thousandths(pair):
    return [math.floor(pair[0]*1000), math.ceil(pair[1]*1000)]


def main():
    length = HORIZON - 1
    # For mu in [6,7], the free sine kernel has absolute value <=1/6.
    pert = 7*ELL**2 * length/6
    u_bound = (F(1, 6))/(1-pert)
    change = (pert + F(3, 4)*(ELL*HORIZON)**2) * u_bound
    cos_lower = 1-(7*SOURCE_RADIUS)**2/2
    error = 7*change/cos_lower  # divide by the common pulse factor C(mu)/mu
    require(change < F(1, 10**12) and error < F(1, 10**11), 'curvature budget')
    require(SOURCE_CENTER-SOURCE_RADIUS == 1 and SOURCE_CENTER+SOURCE_RADIUS == F(65, 64), 'fixed pulse')
    # At every intermediate mass the denominator stays positive. Sine is
    # decreasing from 39*6/128 to 39*7/128, within (pi/2,pi).
    lag = sum(WINDOWS[0])/2-SOURCE_CENTER
    halfwidth = (WINDOWS[0][1]-WINDOWS[0][0])/2
    require(F(11, 7) < 6*lag <= 7*lag < 3, 'sine monotonicity interval')
    denominator_lower = sine(7*lag)[0] * (1-(7*halfwidth)**2/6) - error
    require(denominator_lower > F(4, 5), 'positivity throughout mass sweep')
    results = {}
    for mu in (6, 7):
        means = []
        for left, right in WINDOWS:
            phase = mu*((left+right)/2-SOURCE_CENTER)
            spread = mu*(right-left)/2
            sinc = [x/spread for x in sine(spread)]
            low, high = product(sine(phase), sinc)
            means.append((low-error, high+error))
        denominator, numerator = means
        require(denominator[0] > 0, 'denominator at mass witness')
        ratio = product(numerator, (1/denominator[1], 1/denominator[0]))
        results[str(mu)] = outward_thousandths(ratio)
    require(results['6'] == [-98, -97] and results['7'] == [-720, -719], 'witness enclosures')
    require(results['7'][1] < -406 < -330 < results['6'][0], 'coverage of target band')
    print(json.dumps({'scaled_mass_witnesses_and_ratio_bounds_in_thousandths': results,
                      'normalized_mean_error_below': '1/100000000000',
                      'normalized_first_window_mean_above': '4/5'}, indent=2))


if __name__ == '__main__':
    main()
