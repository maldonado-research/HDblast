"""Flat four-dimensional energy-condition bounds and scalar reconstruction."""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicHermiteSpline, CubicSpline, PPoly, PchipInterpolator

START, END = 1.35, 2.10
BOUND_START = 1.60


def inside(values, start, end):
    return values[np.isfinite(values) & (values > start) & (values < end)]


def extrema(piece, start=START, end=END):
    roots = inside(piece.derivative().roots(extrapolate=False), start, end)
    times = np.r_[start, inside(piece.x, start, end), roots, end]
    values = piece(times)
    return float(min(values)), float(max(values))


def potential_polynomial(h):
    # Coefficients in descending powers of s-s_i for U=3h²+h' on each span.
    coefficients = np.zeros((7, len(h.x)-1))
    for i, a in enumerate(h.c):
        for j, b in enumerate(h.c):
            coefficients[i+j] += 3*a*b
    coefficients[-3:] += h.derivative().c
    return PPoly(coefficients, h.x, extrapolate=False)


def forward(h, subdivisions):
    # Include interpolation margins so the 4D solver stays inside the scalar
    # potential throughout the comparison. Fix the field origin at START.
    left = max(h.x[0], START-.003)
    right = min(h.x[-1], END+.003)
    knots = np.unique(np.r_[left, inside(h.x, left, right), START, END, right])
    mids = (knots[1:]+knots[:-1])/2
    widths = np.diff(knots)/2
    x, w = leggauss(8)
    ss = mids[:, None]+widths[:, None]*x
    dfield = widths * (np.sqrt(-2*h(ss, 1)) @ w)
    phi = np.r_[0., np.cumsum(dfield)]
    phi -= phi[np.searchsorted(knots, START)]
    hp = h(knots, 1)
    v = np.sqrt(-2*hp)
    u = 3*h(knots)**2 + hp
    du = (6*h(knots)*hp + h(knots, 2))/v
    potential = CubicHermiteSpline(phi, u, du, extrapolate=False)

    # Use the strictly rolling field as independent variable. Each RK4 step
    # stays in one interval of the fitted potential, resolving its spline knots.
    # The Cauchy problem below receives only U(phi) and the four initial data.
    first = np.searchsorted(knots, START)
    last = np.searchsorted(knots, END)
    state = np.array([START, v[first], h(START), 0.])
    trajectory = [state.copy()]
    for index in range(first, last):
        dp = (phi[index+1] - phi[index]) / subdivisions
        c3, c2, c1, _ = potential.c[:, index]

        def equations(offset, values):
            _time, speed, rate, _loga = values
            slope = (3*c3*offset + 2*c2)*offset + c1
            return np.array([1/speed, -3*rate-slope/speed, -.5*speed, rate/speed])

        for k in range(subdivisions):
            offset = k*dp
            k1 = equations(offset, state)
            k2 = equations(offset + dp/2, state+dp*k1/2)
            k3 = equations(offset + dp/2, state+dp*k2/2)
            k4 = equations(offset + dp, state+dp*k3)
            state += dp*(k1+2*k2+2*k3+k4)/6
        trajectory.append(state.copy())
    t, speed, rate, loga = np.asarray(trajectory).T
    fields = phi[first:last+1]
    primitive = h.antiderivative()
    constraint = 3*rate**2 - .5*speed**2 - potential(fields)
    return dict(max_abs_h_error=float(max(abs(rate-h(t)))),
                max_abs_time_error=float(max(abs(t-knots[first:last+1]))),
                max_abs_log_scale_error=float(max(abs(loga-primitive(t)+primitive(START)))),
                max_abs_friedmann_constraint=float(max(abs(constraint))),
                field_excursion=float(fields[-1]), initial_field=0.,
                final_h=float(rate[-1]), input_h=float(h(END)))


def inspect(trajectory, run_forward):
    with np.load(trajectory, allow_pickle=False) as data:
        s, rate, loga = (data[name].copy() for name in ('s', 'h', 'A'))
    if min(s) >= START or max(s) <= END or min(np.diff(s)) <= 0:
        raise ValueError('Input must cover the complete proper-time interval')
    # Shape and sign checks for the interpolated rate. PCHIP supplies a second
    # differentiable reconstruction. The endpoint test also works on raw samples.
    h = CubicSpline(s, rate, extrapolate=False)
    shape = PchipInterpolator(s, rate, extrapolate=False)
    a = PchipInterpolator(s, loga, extrapolate=False)
    hp_limits = extrema(h.derivative())
    if extrema(h)[0] <= 0 or hp_limits[1] >= 0:
        raise ValueError('This history lacks the strictly rolling scalar reconstruction')
    u = potential_polynomial(h)
    roots = inside(u.roots(extrapolate=False), START, END)
    if len(roots) != 1:
        raise ValueError('Unexpected number of potential zeros')
    t = np.array([BOUND_START, END])
    positions = np.array([np.searchsorted(s, BOUND_START), np.searchsorted(s, END)-1])
    if not (0 <= positions[0] < positions[1] < len(s)):
        raise ValueError('Missing bracket for endpoint check')

    def endpoints(expansion, field):
        da = float(field[1]-field[0])
        if da <= 0 or min(expansion) <= 0:
            raise ValueError('Endpoint bound assumes expansion')
        logratio = float(np.log(expansion[1]/expansion[0]))
        return dict(delta_log_scale=da, rate_at_endpoints=expansion.tolist(),
                    ratio=float(np.exp(logratio+3*da)),
                    effective_w_averaged_over_log_scale=-1-2*logratio/(3*da))

    row = dict(hprime_extrema_spline=list(hp_limits), hprime_extrema_pchip=list(extrema(shape.derivative())),
               potential_extrema=list(extrema(u)), potential_zero=float(roots[0]),
               proper_time_endpoints=t.tolist(), endpoints=endpoints(shape(t), a(t)),
               raw_sample_times=s[positions].tolist(), raw_sample_endpoints=endpoints(rate[positions], loga[positions]))
    if run_forward:
        row['forward'] = {str(step):forward(h, step) for step in (1, 2)}
        for trial in row['forward'].values():
            if max(trial['max_abs_h_error'], trial['max_abs_log_scale_error'],
                   trial['max_abs_friedmann_constraint'], trial['max_abs_time_error']) > 2e-8:
                raise ValueError(('Scalar IVP mismatch',trial))
    return row


def calculate(runs):
    report = dict(interval=[START, END], histories={})
    for grid in ('coarse', 'fine'):
        for arm in ('quiet','plus','minus','far'):
            for suffix in ('zero','120_static_p0.01'):
                name = f'{grid}_{arm}_{suffix}'
                row = inspect(runs/(name+'.npz'), grid=='fine' and arm=='quiet')
                if not row['raw_sample_endpoints']['ratio'] < .60:
                    raise ValueError(('Dominant-energy endpoint check',name))
                report['histories'][name] = row
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    payload = json.dumps(calculate(args.runs), indent=2, allow_nan=False)+'\n'
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload, end='')
