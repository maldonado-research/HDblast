"""Shoot the two B1 initial constraints through a compact scalar pulse on one slice."""
import argparse
import hashlib
import json
from math import exp, expm1, log
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import root

PIN = "bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f"
CENTERS = np.array([-0.17, -0.125, -0.08, -0.035])
WIDTH = 0.020
LEFT = -0.20


def potential(phi):
    w = 1 - phi + phi**3 / 3
    wp = phi**2 - 1
    return wp**2 / 2 - 2 * w**2 / 3, 2 * phi * wp - 4 * w * wp / 3


def delta_potential(phi, h):
    """Evaluate U(phi+h)-U(phi) without subtraction of order-one quantities."""
    w = 1 - phi + phi**3 / 3
    wp = phi**2 - 1
    dw = h * (wp + phi*h + h*h/3)
    dwp = h * (2*phi + h)
    return wp*dwp + dwp*dwp/2 - 2*dw*(2*w + dw)/3


def bumps(z):
    values = np.zeros(4)
    derivatives = np.zeros(4)
    for i, center in enumerate(CENTERS):
        x = (z - center) / WIDTH
        if abs(x) < 1:
            den = 1 - x*x
            values[i] = exp(1 - 1/den)
            derivatives[i] = -2*x*values[i] / (den**2 * WIDTH)
    return values, derivatives


def wall(summary):
    phi = summary["phi_b_static"]
    rho = summary["rho_b"]
    par = summary["params"]
    c, d, det = par["c"], par["d"], par["delta"]
    tension = 2*(1 - phi + phi**3/3) + det*(1 + c*phi + d*phi**2/2)
    slope = 2*(phi**2 - 1) + det*(c + d*phi)
    return [log(rho), rho*tension/6, phi, -rho*slope/2]


def integrate(summary, coeff, step=0.002, tol=2e-12, sens=True, dense=False):
    """Return the quiet fields B,Q,phi,r followed by differences dB,dQ,du.

    Initial data have A=B, B_T=0, u=A_T and p=phi_T=−F'.  The quiet
    Hamiltonian form and quiet scalar equation determine the reference.
    Additional rows are the three-by-four endpoint variational matrix.
    """
    initial = wall(summary) + [0.0]*3
    if sens:
        initial += [0.0]*12

    def rhs(z, state):
        B, Q, phi, r, d, q, v = state[:7]
        U, Up = potential(phi)
        e = exp(2*B)
        f, df = bumps(z)
        h, j = coeff @ f, coeff @ df
        hU = delta_potential(phi, h)
        ed = exp(2*d)
        # h=F, j=F'.  Then p=-j, phi_Z=r+j, and du/dZ=j(r+j)/3.
        dq = v*(2+v) - q*(2*Q+q) - j*(r+j)/3 - e*(expm1(2*d)*U + ed*hU)/3
        out = [Q, 1-Q*Q-r*r/6-e*U/3, r, e*Up-3*Q*r, q, dq, j*(r+j)/3]
        if sens:
            bd, qd, vd = state[7:].reshape(3, 4)
            qb = -2*e*ed*(U+hU)/3
            qj = -(r+2*j)/3
            phi_force = -e*ed*potential(phi+h)[1]*f/3
            out.extend(qd)
            out.extend(qb*bd - 2*(Q+q)*qd + 2*(1+v)*vd + qj*df + phi_force)
            out.extend(-qj*df)
        return out

    sol = solve_ivp(rhs, (0, LEFT), initial, rtol=tol, atol=tol*1e-3,
                    method="DOP853", max_step=step, dense_output=dense)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def outer(summary, coeff):
    end = integrate(summary, coeff, step=0.0015, tol=4e-13).y[:, -1]
    return end[4:7], end[7:].reshape(3, 4)


def shoot_fixed(summary, first, guess):
    """Hold the first lobe fixed and solve for the other three coefficients."""
    def equations(other):
        residual, jac = outer(summary, np.r_[first, other])
        return residual, jac[:, 1:]
    result = root(equations, guess, jac=True, tol=2e-10)
    coeff = np.r_[first, result.x]
    residual, jac = outer(summary, coeff)
    if np.max(np.abs(residual)) > 1e-13:
        raise RuntimeError("Shooting residual exceeds the target")
    return coeff, residual, jac


def sample_checks(summary, coeff):
    sol = integrate(summary, coeff, step=0.0004, tol=3e-14, sens=False, dense=True)
    _, Q0, _, r0, dB, dQ, du = sol.sol(np.linspace(LEFT, 0, 2001))
    # On both quiet collars the forcing vanishes. Record that the interior
    # solution is nontrivial and check both quiet collars independently.
    def collar(low, high):
        values = sol.sol(np.linspace(low, high, 41))[4:7]
        return np.max(abs(values), axis=1).tolist()
    def flux(z):
        r = sol.sol(z)[3]
        _, df = bumps(z)
        j = coeff @ df
        return r*j, j*j
    cuts = sorted(float(z+offset) for z in CENTERS for offset in [-WIDTH, WIDTH])
    terms = [quad(lambda z: flux(z)[k], LEFT, 0, epsabs=4e-17,
                  epsrel=2e-12, points=cuts)[0] for k in (0, 1)]
    return {
        "outer_finer": sol.y[4:7, -1].tolist(),
        "max_in_wall_collar_abs_dB_dQ_du": collar(-0.015, 0.0),
        "max_in_far_collar_abs_dB_dQ_du": collar(LEFT, -0.19),
        "max_patch_abs_dB_dQ_du": [float(max(abs(q))) for q in (dB,dQ,du)],
        "integrated_scalar_gradient_cross_and_square": terms,
        "integrated_momentum_residual": sum(terms),
    }


def calculate(summary):
    base = np.zeros(4)
    quiet_residual, J = outer(summary, base)
    tangent = np.linalg.svd(J)[2][-1]
    tangent *= np.sign(tangent[1])
    finite = []
    scale = 2e-7
    for direction in np.eye(4):
        endpoints = []
        for s in (-scale, scale):
            sol = integrate(summary, s*direction, step=0.0008, tol=1e-13, sens=False)
            endpoints.append(sol.y[4:7,-1])
        finite.append((endpoints[1]-endpoints[0])/(2*scale))
    jac_error = np.max(abs(np.array(finite).T-J))
    probe = np.zeros(4)
    continuation = []
    for first in (-2.5e-5, -5e-5, -7.5e-5, -1e-4):
        predictor = np.linalg.solve(J[:,1:], -J[:,0])
        guess = probe[1:] + (first-probe[0])*predictor
        probe, residual, J = shoot_fixed(summary, first, guess)
        continuation.append({"coefficients": probe.tolist(), "residual": residual.tolist()})
    second, _, _ = shoot_fixed(summary, -1e-4, [-0.00024, 0.00030, -0.000055])
    profiles = []
    for label, coeff in [("quiet_connected", probe), ("second_completion", second)]:
        residual, jac = outer(summary, coeff)
        item = {"solution":label,"coefficients":coeff.tolist(),
                "outer":residual.tolist(),
                "singular_values_fixed_first":np.linalg.svd(jac[:,1:],compute_uv=False).tolist()}
        item.update(sample_checks(summary,coeff))
        profiles.append(item)
    single = np.array([-1e-4,0,0,0])
    omit = probe.copy()
    omit[3] = 0
    quiet = outer(summary, base)[1]
    return {
        "source_summary_sha256":PIN, "chart":"identity", "dc":0,
        "left":LEFT,"wall":0,"centers":CENTERS.tolist(),"half_width":WIDTH,
        "quiet_endpoint":quiet_residual.tolist(),
        "quiet_jacobian":quiet.tolist(),
        "quiet_singular_values":np.linalg.svd(quiet,compute_uv=False).tolist(),
        "quiet_last_three_determinant":float(np.linalg.det(quiet[:,1:])),
        "quiet_tangent":tangent.tolist(),
        "centered_jacobian_max_difference":float(jac_error),
        "continuation":continuation,"profiles":profiles,
        "first_lobe_only_outer":outer(summary,single)[0].tolist(),
        "omit_fourth_lobe_outer":outer(summary,omit)[0].tolist(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary",type=Path)
    args = parser.parse_args()
    raw = args.summary.read_bytes()
    if hashlib.sha256(raw).hexdigest() != PIN:
        raise ValueError("Unexpected B1 background summary")
    print(json.dumps(calculate(json.loads(raw)),indent=2,allow_nan=False))
