"""Construct matched B1 initial fields on an identity-chart conformal slice."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

SUMMARY_SHA256 = 'bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f'
PULSE_CENTER = -0.12
PULSE_WIDTH = 0.02
AMPLITUDE = 0.001
# All positions and coefficients use the original conformal chart coordinates.
PROFILES = {
    'delayed': {
        'centers': (-0.17, -0.22, -0.28), 'width': 0.018, 'cut': -0.32,
        'half_positive': (-0.30274771443373377, 0.5559923254887478, -0.24580011663492168),
        'half_negative': (0.5473081970670742, -0.9950587848179364, 0.43446032156457476),
        'negative_coefficients': (1.3579490275061388, -2.441325666600515, 1.0480625197638658),
        'coefficients': (-0.3705278325568672, 0.6783951490735094, -0.2972716055389422),
    },
    'interleaved': {
        'centers': (-0.05, -0.095, -0.158), 'width': 0.022, 'cut': -0.20,
        'negative_coefficients': (0.04608223945470865, -0.06677346578623844, 0.0183221735921639),
        'coefficients': (-0.008892382515266004, 0.0027532337147582486, 0.009532969691098724),
    },
}


def bump(z, center, width):
    """A unit-height smooth, compactly supported function and its derivative."""
    x = (z-center)/width
    if abs(x) >= 1:
        return 0.0, 0.0
    denominator = 1-x*x
    value = math.exp(1-1/denominator)
    return value, -2*x*value/(width*denominator**2)


def potential(phi):
    w = 1-phi+phi**3/3
    dw = phi*phi-1
    return dw*dw/2-2*w*w/3, 2*phi*dw-4*w*dw/3


def wall_values(summary_path):
    data = Path(summary_path).read_bytes()
    if hashlib.sha256(data).hexdigest() != SUMMARY_SHA256:
        raise ValueError('Wrong B1 background summary')
    summary = json.loads(data)
    phi, rho = summary['phi_b_static'], summary['rho_b']
    delta, c, d = (summary['params'][key] for key in ('delta', 'c', 'd'))
    W = 1-phi+phi**3/3
    sigma = 2*W+delta*(1+c*phi+d*phi*phi/2)
    sigma_phi = 2*(phi*phi-1)+delta*(c+d*phi)
    return math.log(rho), rho*sigma/6, phi, -rho*sigma_phi/2


def geometry_rate(z, coefficients, profile):
    return sum(c*bump(z, center, profile['width'])[0]
               for c, center in zip(coefficients, profile['centers'], strict=True))


def integrate(wall, left, amplitude, coefficients, profile, step=0.001, *,
              method='DOP853', rtol=3e-13, atol=2e-15, dense=True):
    """Integrate quiet and sourced constraints from the common wall to the left."""
    B_wall, Q_wall, phi_wall, r_wall = wall

    def rhs(z, y):
        b, q, phi, r, B, Q, u = y
        U, U_phi = potential(phi)
        exp_b = math.exp(2*b)
        f, df = bump(z, PULSE_CENTER, PULSE_WIDTH)
        trial_U = potential(phi + amplitude*f)[0]
        radial = r+amplitude*df
        momentum = -amplitude*df
        v = geometry_rate(z, coefficients, profile)
        quiet_Qz = 1-q*q-r*r/6-exp_b*U/3
        trial_Qz = u*u+u*v-Q*Q-(radial**2+momentum**2)/6-math.exp(2*B)*trial_U/3
        uz = Q*v-momentum*radial/3
        return q, quiet_Qz, r, exp_b*U_phi-3*q*r, Q, trial_Qz, uz

    result = solve_ivp(rhs, (0, left), (B_wall,Q_wall,phi_wall,r_wall,B_wall,Q_wall,1),
                       method=method, rtol=rtol, atol=atol, dense_output=dense, max_step=step)
    if not result.success:
        raise RuntimeError(result.message)
    return result


def endpoint(solution):
    b, q, _, _, B, Q, u = solution.y[:, -1]
    return np.array([B-b, Q-q, u-1])


def initial_data(z, summary_path, amplitude=AMPLITUDE, profile_name='delayed', step=0.001):
    """Return quiet, trial and delta dictionaries for (A,B,phi,A_T,B_T,phi_T).

    z is an arbitrary 1D set of finite nonpositive coordinates in the regular
    background patch. Choose ±0.001 or zero for both profiles or ±0.0005 for the delayed profile. Gluing
    uses the verified numerical root. The interval outside the cut uses the
    quiet reference, and so does the open collar nearest the wall.
    """
    z = np.asarray(z, dtype=float)
    if z.ndim != 1 or z.size == 0 or not np.isfinite(z).all() or np.max(z)>0 or step<=0:
        raise ValueError('Expected finite nonpositive coordinates and a positive step')
    profile = PROFILES[profile_name]
    keys = {AMPLITUDE: 'coefficients', -AMPLITUDE: 'negative_coefficients',
            AMPLITUDE/2: 'half_positive', -AMPLITUDE/2: 'half_negative'}
    c = np.zeros(3)
    if amplitude != 0:
        key = keys.get(amplitude)
        if key not in profile:
            raise ValueError('No matched coefficients are available for this profile and amplitude')
        c = np.asarray(profile[key])
    wall = wall_values(summary_path)
    cut = profile['cut']
    # Integrate the trial only to the matching cut, even for large quiet extensions.
    solution = integrate(wall, cut, amplitude, c, profile, step)
    if np.max(np.abs(endpoint(solution)))>1e-9:
        raise ValueError('Endpoint does not match at this precision')
    values = solution.sol(np.maximum(z, cut))
    b, _, phi, _, B, _, u = values
    outside = z < cut
    if np.any(outside):
        quiet_ext = integrate(wall,float(z.min()),0.,np.zeros(3),profile, step)
        b[outside], phi[outside] = quiet_ext.sol(z[outside])[[0, 2]]
    zero = np.zeros_like(z)
    quiet = dict(A=b.copy(), B=b.copy(), phi=phi.copy(), A_T=zero+1, B_T=zero.copy(), phi_T=zero.copy())
    f, df = np.array([bump(x, PULSE_CENTER, PULSE_WIDTH) for x in z]).T
    v = np.array([geometry_rate(x, c, profile) for x in z])
    right = max(PULSE_CENTER+PULSE_WIDTH, max(profile['centers'])+profile['width'])
    in_patch = (z>=cut)&(z<right)
    metric_delta = np.where(in_patch, B-b, 0.)
    delta = dict(A=metric_delta.copy(), B=metric_delta, phi=amplitude*f,
                 A_T=np.where(in_patch,u-1,0.), B_T=v, phi_T=-amplitude*df)
    if amplitude == 0:
        delta = {key:zero.copy() for key in quiet}
    trial = {key:quiet[key]+delta[key] for key in quiet}
    return dict(quiet=quiet, trial=trial, delta=delta)
