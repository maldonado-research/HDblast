"""Positive Gaussian memory on a compact frequency band, independently of shell evolution."""
import argparse
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import roots_legendre, spherical_jn


LOW, HIGH = 60, 120
ORDERS = (1, 2, 4, 8, 16, 32)
HORIZONS = (.1, .25, .7, 2.)
# Smallest and largest modes from the archived 24-point native momentum rule.
MOMENTA = (.019251120011914757, 7.980748879988085)


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def recurrence(q, order):
    """Exact-weight monic recurrence in y=(LOW/omega)^2 using Decimal moments."""
    lo, hi = Decimal(LOW), Decimal(HIGH)
    q2 = Decimal(str(q))**2
    a = (lo/hi)**2

    def power_integral(k):
        if k == -1:
            return -a.ln()
        return (1-a**(k+1))/(k+1)

    moments = [lo**2*power_integral(k-2)-q2*power_integral(k-1)
               for k in range(2*order+1)]

    def inner(poly, shift):
        return sum(v*w*moments[i+j+shift]
                   for i, v in enumerate(poly) for j, w in enumerate(poly))

    diagonal, subdiagonal_squared = [], []
    old, polynomial = [], [Decimal(1)]
    last_norm = None
    for n in range(order):
        norm = inner(polynomial, 0)
        require(norm > 0, 'The exact moment Gram matrix lost positivity')
        alpha = inner(polynomial, 1)/norm
        beta = norm/last_norm if n else Decimal(0)
        diagonal.append(alpha)
        if n:
            subdiagonal_squared.append(beta)
        following = [Decimal(0), *polynomial]
        for i, coefficient in enumerate(polynomial):
            following[i] -= alpha*coefficient
        for i, coefficient in enumerate(old):
            following[i] -= beta*coefficient
        old, polynomial, last_norm = polynomial, following, norm
    return diagonal, subdiagonal_squared, moments


def decimal_resolvent(zeta, q, diagonal, beta, mass):
    """Continuum normalized W and N-mode W via a high-precision continued fraction."""
    z, q2, low2, high2 = zeta*LOW**2, Decimal(str(q))**2, Decimal(LOW**2), Decimal(HIGH**2)
    reference = (high2-low2-(q2+z)*((high2+z)/(low2+z)).ln())/mass
    denominator = 1+zeta*diagonal[-1]
    for a, b in zip(reversed(diagonal[:-1]), reversed(beta), strict=True):
        denominator = 1+zeta*a-zeta**2*b/denominator
    return reference, 1/denominator


def vacuum_kernel(q, lag):
    """Exact E=integral(g²/omega exp(-i*omega*lag)) using spherical Bessel integrals."""
    center, half = (HIGH+LOW)/2, (HIGH-LOW)/2
    u = half*lag
    polynomial = ((center**2+half**2/3-q*q)*spherical_jn(0, u)
                  - 2j*center*half*spherical_jn(1, u)
                  - 2*half**2/3*spherical_jn(2, u))
    return 2*half*np.exp(-1j*center*lag)*polynomial/(4*np.pi**2)


def ellipse_bounds():
    """Conservative upper bounds on the ellipse amplitude and phase growth."""
    a, b = (LOW/HIGH)**2, 1.
    m, h = (a+b)/2, (b-a)/2
    # Re(z) and |Im(z)| are even about the real axis. Trapezium-endpoint
    # bounds on each half-ellipse arc avoid optimizing an oscillatory kernel.
    angles = np.linspace(0, np.pi, 257)
    left, right = angles[:-1], angles[1:]
    sine_max = np.maximum(np.sin(left), np.sin(right))
    sine_max[(left <= np.pi/2) & (np.pi/2 <= right)] = 1
    # Each candidate proves a bound; minimizing over these candidates is valid.
    rho = np.linspace(1.02, 2.98, 197)
    c = (rho+1/rho)[:, None]/2
    s = (rho-1/rho)[:, None]/2
    min_real = m+h*c*np.cos(right)
    max_imag = h*s*sine_max
    require(np.min(min_real) > 0, 'Ellipse crossed the branch cut')
    amplitude = LOW/(HIGH*np.sqrt(min_real))
    phase_growth = LOW*np.minimum(1/np.sqrt(8*min_real), max_imag/(2*min_real**1.5))
    return rho, amplitude, phase_growth


def time_envelope(order, horizon, ellipse):
    rho, amplitude, phase_growth = ellipse
    log_M = np.max(np.log(amplitude)+horizon*phase_growth, axis=1)
    log_bound = math.log(4)+log_M-(2*order-1)*np.log(rho)-np.log(rho-1)
    idx = np.argmin(log_bound)
    # Both continuum and discrete kernels separately obey |E| <= mass*HIGH.
    return min(2., math.exp(log_bound[idx])), float(rho[idx])


def measure_rule(diagonal, beta, mass):
    y, vec = eigh_tridiagonal(np.array(diagonal, float), np.sqrt(np.array(beta, float)),
                              lapack_driver='stev')
    return y, mass*vec[0]**2


def source_comparison(rules):
    """Cross-check the complete published source rules through order 16."""
    worst_frequency = worst_coupling = 0.
    with np.load(rules, allow_pickle=False) as archive, localcontext() as context:
        context.prec = 180
        momenta = archive['q']
        require(len(momenta) == 24 and (min(momenta),max(momenta)) == MOMENTA,
                'Unexpected native calibration grid')
        for row, q in enumerate(momenta):
            diagonal, beta, moments = recurrence(q, 16)
            mass = float(moments[0])/(8*np.pi**2)
            for n in ORDERS[:5]:
                y, mu = measure_rule(diagonal[:n], beta[:n-1], mass)
                omega = (LOW/np.sqrt(y))[::-1]
                g2 = mu[::-1]*omega**2
                om_rel = float(np.max(abs(archive[f'omega_N{n}'][row]/omega-1)))
                g2_rel = float(np.max(abs(archive[f'g2_N{n}'][row]/g2-1)))
                worst_frequency = max(worst_frequency, om_rel)
                worst_coupling = max(worst_coupling, g2_rel)
    require(worst_frequency < 1e-11 and worst_coupling < 1e-9, 'Published memory rules disagree')
    return dict(momenta=len(momenta), max_order=16, max_relative_frequency_error=worst_frequency,
                max_relative_squared_coupling_error=worst_coupling)


def experiment(rules):
    output = {'band': [LOW, HIGH], 'normalization': 'absolute complex-kernel error / (C * HIGH)',
              'lags_per_horizon': 1401, 'source_rule_comparison': source_comparison(rules),
              'momentum_extremes': {}}
    ellipse = ellipse_bounds()
    with localcontext() as context:
        context.prec = 320
        for q in MOMENTA:
            diagonal, beta, moments = recurrence(q, max(ORDERS))
            mass = float(moments[0])/(8*np.pi**2)
            node, weight = roots_legendre(160)
            freq = (HIGH+LOW+(HIGH-LOW)*node)/2
            lags = np.array([0, .113, .337, .681])
            density = weight*(HIGH-LOW)*(freq**2-q**2)/(8*np.pi**2)
            direct = density @ np.exp(-1j*freq[:, None]*lags)
            np.testing.assert_allclose(direct, vacuum_kernel(q, lags), rtol=2e-11, atol=1e-9)
            by_order = {}
            for n in ORDERS:
                a, b = diagonal[:n], beta[:n-1]
                axis_results = []
                for z in ('.001', '.1', '1', '100', '10000'):
                    zeta = Decimal(z)
                    reference, approx = decimal_resolvent(zeta, q, a, b, moments[0])
                    error = (reference-approx)/reference
                    center, width = .625, .375
                    bound = 1/math.cosh(n*math.acosh((1/float(zeta)+center)/width))**2
                    require(0 < error <= Decimal(str(bound))*(1+Decimal('1e-12')),
                            'Euclidean sign or minimax bound failed')
                    axis_results.append(dict(z_over_low_squared=float(zeta),
                                             rel_shortfall=float(error), upper_bound=bound))
                values, mu = measure_rule(a,b,mass)
                omega = LOW/np.sqrt(values)
                require(np.min(mu)>0 and np.all((values > .25) & (values < 1.)), 'Invalid positive rule')
                exact = np.array(moments[:2*n], float)/float(moments[0])
                sampled = (mu/mass) @ (values[:, None]**np.arange(2*n))
                require(np.max(abs(sampled-exact)) < 3e-14, 'Polynomial moment mismatch')
                time = []
                for horizon in HORIZONS:
                    lags = np.linspace(0, horizon, output['lags_per_horizon'])
                    actual = (mu*omega) @ np.exp(-1j*omega[:, None]*lags)
                    delta = (actual-vacuum_kernel(q, lags))/(mass*HIGH)
                    error = float(np.max(abs(delta)))
                    envelope, rho = time_envelope(n, horizon, ellipse)
                    require(error <= envelope+8e-14, 'Compact-time kernel envelope failed')
                    step = horizon/(len(lags)-1)
                    # Linear interpolation of a vector-valued function misses
                    # at most max|f''| * step²/8. Both positive spectral
                    # measures have mass C and frequencies at most HIGH.
                    allowance = (HIGH*step)**2/4
                    time.append(dict(conformal_lag_max=horizon, sampled_complex_error=error,
                                     sampled_retarded_error=float(np.max(abs(delta.imag))),
                                     sampled_twice_symmetric_noise_error=float(np.max(abs(delta.real))),
                                     analytic_upper_bound=envelope, ellipse_parameter=rho,
                                     grid_interpolation_allowance=allowance,
                                     grid_conditioned_upper_bound=error+allowance))
                by_order[str(n)] = {'imaginary_axis': axis_results, 'finite_lag': time}
            output['momentum_extremes'][str(q)] = by_order
    return output


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--rules', type=Path, default=Path(__file__).resolve().parents[1]/'source/RULES.npz')
    ap.add_argument('--output', type=Path, default=Path(__file__).with_name('CHECK.json'))
    args = ap.parse_args()
    args.output.write_text(json.dumps(experiment(args.rules), indent=2, allow_nan=False)+'\n')
