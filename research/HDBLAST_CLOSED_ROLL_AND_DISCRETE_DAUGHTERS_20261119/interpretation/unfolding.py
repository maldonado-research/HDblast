"""Angular quadrupole splitting of a first-contact two-sphere."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre


def axisymmetric(s):
    if s < 0:
        raise ValueError('s must be nonnegative')
    if s < .002:
        # Small-s expansion of the closed expression avoids cancellation.
        return s**1.5*sum(s**j/((2*j+1)*(2*j+3)) for j in range(8))
    if s == 1:
        return .5
    root = np.sqrt(s)
    if s < 1:
        return .5*(root-(1-s)*np.arctanh(root))
    return .5*(root+(s-1)*np.arctanh(1/root))


def angular_mean(fn, eta, s):
    # The only interior transition is where h(phi)=s. The integrand is
    # continuous with an integrable logarithmic derivative there.
    splits = []
    if eta < s < 1:
        splits.append(np.arcsin(np.sqrt((s-eta)/(1-eta))))
    val, _ = quad(lambda phi: fn(eta+(1-eta)*np.sin(phi)**2),
                  0., np.pi/2, points=splits, epsabs=1e-22, epsrel=5e-13, limit=120)
    return 2*val/np.pi


def reduced(eta, s):
    return angular_mean(lambda h: np.sqrt(h)*axisymmetric(s/h), eta, s)


def direct(eta, s, order):
    t, wt = roots_legendre(order)
    t, wt = (1+t)/2, wt/2
    def meridian(h):
        # Integrate over u=|cos(theta)|. For a polar cap map the root to
        # t=0 and factor the radicand to keep tiny cap volumes accurate.
        if s >= h:
            return wt @ np.sqrt((s-h)+h*t*t)
        lo = np.sqrt(1-s/h)
        width = (s/h)/(1+lo)
        du = width*t*t
        return wt @ (2*width*t*np.sqrt(h*du*(2*lo+du)))
    return angular_mean(meridian, eta, s)


def calculate():
    eta_values = (1., .5)
    s_values = (1e-6, 1e-5, .001, .01, .1, .5, .75, 1., 1.25, 2., 10., 100., 1000.)
    cases = {}
    for eta in eta_values:
        rows = []
        for s in s_values:
            value = reduced(eta, s)
            direct64, direct128 = (direct(eta, s, order) for order in (64, 128))
            rows.append(dict(s=s, scaling_function=value, relative_to_simultaneous=value/np.sqrt(s),
                             over_early_limit=value/(s**1.5/(3*np.sqrt(eta))),
                             over_late_expansion=value/(np.sqrt(s)*(1-(1+eta)/(6*s))) if s>=2 else None,
                             numerical_64=direct64, numerical_128=direct128,
                             numerical_128_relative_difference=direct128/value-1))
        early = np.log(rows[1]['scaling_function']/rows[0]['scaling_function'])/np.log(10.)
        late = np.log(rows[-1]['scaling_function']/rows[-2]['scaling_function'])/np.log(10.)
        maxerr = max(abs(r['numerical_128_relative_difference']) for r in rows)
        maxrefine = max(abs(r['numerical_128']/r['numerical_64']-1) for r in rows)
        if maxerr>2e-9 or maxrefine>1e-7 or abs(early-1.5)>1e-5 or abs(late-.5)>.002:
            raise ArithmeticError((eta, maxerr, maxrefine, early, late))
        cases[str(eta)] = dict(rows=rows, early_log_slope=early, late_log_slope=late,
                               maximum_relative_error_128=maxerr, maximum_meridian_refinement=maxrefine)
    return dict(model='prescribed angular onset times; local volume geometry only', cases=cases)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.output.write_text(json.dumps(calculate(), indent=2, allow_nan=False)+'\n')
