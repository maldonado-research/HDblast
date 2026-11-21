"""Causal transverse-force limits from the age distribution of unscattered daughters."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import kve
from scipy.integrate import quad


TUBE = math.sqrt(1.05**2 - 1.)
TOLERANCES = (0., 1.e-4, .001, .01, .05, .1)


def pair_kernel_upper():
    mass2, width = 8300.**2, .30
    cutoff = mass2 / 2.
    low = 24. * math.pi * width**2 * (cutoff + 2.)**2 / (cutoff * (mass2 - cutoff)**2)
    d = width**2 * (1. - 8. / mass2)
    high = 24. * math.pi / mass2 * width**2 * (1. + 4. / mass2)**2 * (math.sqrt(1. + d) + 1.) / (2. * d)
    return max(low, high)


def thermal_cdf(y, theta):
    """Probability of p_perp/m <= y in the longitudinally boosted thermal stream."""
    z = 1. / theta
    mass = np.hypot(1., y)
    x = z * mass
    log_tail = 2. * np.log(mass) + np.log(kve(2, x) / kve(2, z)) + z - x
    return np.maximum(0., -np.expm1(np.minimum(log_tail, 0.)))


def inverse_cdf(u, theta):
    lo = np.zeros_like(u)
    hi = np.full_like(u, max(16., 16. * theta))
    while np.any(thermal_cdf(hi, theta) < u):
        hi *= 2.
    for _ in range(47):
        mid = (lo + hi) / 2.
        low = thermal_cdf(mid, theta) < u
        lo = np.where(low, mid, lo)
        hi = np.where(low, hi, mid)
    return lo


def backward_sum(data):
    return np.cumsum(data[::-1])[::-1]


def enclosed_suffix(data, time, density_prefactor):
    """Discard terminal-cell births, overcount exposure and available time to the right knot."""
    j = int(np.searchsorted(data['t'], time))
    if not 2 <= j < len(data['t']):
        raise ValueError("Marker lacks a bracketing native interval")
    t, a, rest = (data[k][:j + 1] for k in ('t', 'a', 'D'))
    dt = np.diff(t)
    if np.any(dt <= 0.) or np.any(np.diff(a) < 0.) or np.any(np.diff(rest) < 0.) or rest[0] < 0.:
        raise ValueError("Nonmonotone birth or expansion history")
    exposure = density_prefactor * rest[1:] * dt / a[:-1]**3
    unhit = np.diff(rest) * np.exp(-backward_sum(exposure)) / rest[-1]
    # Only the completed cells can contribute guaranteed pre-marker births.
    unhit[-1] = 0.
    suffix = backward_sum(unhit)
    # a(t_birth) <= a(marker) during expansion. Use right endpoint over left marker neighbor.
    redshifted_age = backward_sum(dt * np.minimum(1., a[1:] / a[-2]))
    theta = min(data['theta'][j - 1:j + 1])
    if theta <= 0 or np.any(data['R'][:j + 1] > 1.05 + 1.e-13):
        raise ValueError("Birth momentum or thermal-target assumptions changed")
    return t[:-1], redshifted_age, suffix, float(theta), j


def study(data, metadata, spec, time):
    factor = spec['channel_fraction'] * metadata['action'] * pair_kernel_upper() / (6. * metadata['daughter_mass_L']**3)
    births, ages, suffix, theta, j = enclosed_suffix(data, time, factor)
    q = metadata['charge_to_mass']
    upper_x = float(np.abs(data['X'][:j + 1]).max())
    upper_field = float(np.hypot(data['X'][:j + 1], data['B'][:j + 1]).max())
    peak_tests = {}
    for tol in TOLERANCES:
        goal = np.maximum(0., suffix - tol)
        required = np.maximum(0., inverse_cdf(goal, theta) - TUBE) / ages
        i = int(np.argmax(required))
        force = float(required[i])
        cutoff = TUBE + force * ages[i]
        loss = float(suffix[i] - thermal_cdf(cutoff, theta))
        if force > 0. and abs(loss - tol) > 2.e-10:
            raise ValueError("CDF and inverted-quantile witnesses disagree")
        var_check = float(np.max(suffix - thermal_cdf(TUBE + force * ages, theta)))
        if var_check > tol + 2.e-10:
            raise ValueError("Witness fails its forward CDF check")
        # Independent positive-density integral at the reported witness.
        # x is p_perp/(m*theta), z=1/theta.
        z = 1. / theta
        def density(x):
            e = math.hypot(x, z)
            return x * e * math.exp(z - e) * kve(1, e) / (z*z*kve(2, z))
        independent = quad(density, 0., cutoff / theta, epsabs=2.e-13, epsrel=1.e-10)[0]
        if abs(independent - float(thermal_cdf(cutoff, theta))) > 2.e-9:
            raise ValueError("Independent transverse density integral failed")
        peak_tests[str(tol)] = dict(
            peak_force_lower_L_over_mc=force, equivalent_electric_field_lower=force / q,
            earliest_tag_time_for_witness=float(births[i]), lookback_L_at_witness=float(time - births[i]),
            redshift_weighted_time_upper=float(ages[i]), minimum_unhit_recent_share=float(suffix[i]),
            momentum_cutoff_over_mc=float(cutoff), thermal_share_below_cutoff=float(thermal_cdf(cutoff, theta)),
            em_field_energy_peak_lower_over_peak_background=.5 * (force / (q * upper_field))**2,
            force_lower_relative_to_peak_parallel_electric=force / (q * upper_x))
    # A bounded extra Lorentz force of Q times the largest original axial electric field.
    max_gap = float(np.max(suffix - thermal_cdf(TUBE + q * upper_x * ages, theta)))
    return dict(marker_time=time, enclosing_time=[float(data['t'][j - 1]), float(data['t'][j])],
                guaranteed_unhit_share=float(suffix[0]), min_endpoint_theta=theta,
                sampled_peak_parallel_electric=upper_x, sampled_peak_background_EM_norm=upper_field,
                minimum_tv_for_force_at_peak_parallel_electric=max(0., max_gap),
                peak_force_tests=peak_tests)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    meta = json.loads((args.history / 'markers.json').read_text())
    rows = {}
    for name, spec in meta['runs'].items():
        with np.load(args.history / (name + '.npz'), allow_pickle=False) as source:
            data = {k: source[k] for k in source.files}
        rows[name] = {f: study(data, meta, spec, marker['time']) for f, marker in spec['markers'].items()}
    args.output.write_text(json.dumps({'native_cell_enclosures': rows}, sort_keys=True, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
