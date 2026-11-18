"""A finite-band Gaussian field with amplitude damping and a massless energy reservoir."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

CUTS = (.8, 1.4, 1.55, 1.8, 2.10)
BRANCHES = {"quiet": ("classical", "quiet"), "near_plus": ("classical", "0.001"),
            "near_minus": ("classical", "-0.001"), "far_plus": ("far", "0.001")}


def background(capsule, mesh, branch):
    base = capsule / "classical/runs" / mesh
    folder, case = BRANCHES[branch]
    with np.load(base / "checkpoint.npz", allow_pickle=False) as z:
        old = z["hist"]
    with np.load(capsule / folder / "runs" / mesh / case / "rolling_timeseries.npz", allow_pickle=False) as z:
        wall = dict(zip(z["cols"], z["rec"].T, strict=True))
    event = json.loads((base / "quiet/rolling_summary.json").read_bytes())["matter"]["events"][0]
    split = wall["H0tau"][0]
    use = old[:, 0] < split - 1e-10
    s = np.r_[old[use, 0], wall["H0tau"]]
    phi = CubicSpline(s, np.r_[old[use, 1], wall["phi_b"]])
    h = CubicSpline(s, np.r_[old[use, 2], wall["H_over_H0"]])
    loga = h.antiderivative()
    normalization = wall["A_b"][0] - event["A_star"] - loga(split)
    def geometry(t):
        return np.exp(normalization + loga(t)), h(t), 100*(phi(t)-.5), 100*phi(t, 1)
    return geometry


def rule(order):
    x, w = np.polynomial.legendre.leggauss(order)
    k = 4*(1 + x)
    return k, 4*w*k*k/(2*math.pi**2)


def solve(background_fn, y, order, tol):
    k, weights = rule(order)
    a, _, mu, _ = background_fn(CUTS[0])
    energy = np.hypot(k/a, mu)
    factor = a**3*energy
    init = np.r_[.5/factor, np.zeros(order), .5*factor, np.zeros(6)]

    def terms(t, state):
        a, h, mu, dmu = background_fn(t)
        vol = a**3
        p2 = (k/a)**2
        e = np.sqrt(p2+mu*mu)
        factor = vol*e
        u, c, v = state[:3*order].reshape(3, order)
        du = u-.5/factor
        dv = v-.5*factor
        excitation = .5*(dv/vol + vol*e*e*du)
        pressure = .5*dv/vol - vol*(mu*mu/2 + p2/6)*du
        gamma = y*y*mu*mu/(8*math.pi*e)
        n = excitation/e
        source = vol*mu*dmu*du
        # The diagonal-only state in the instantaneous physical oscillator.
        diagonal_source = mu*dmu*n/e
        transfer = weights @ (gamma*excitation)
        return (u, c, v, du, dv, a, vol, h, e, gamma, excitation, pressure,
                source, diagonal_source, transfer)

    def rhs(t, state):
        u,c,v,du,dv,a,vol,h,e,g,exc,press,source,diag,loss = terms(t,state)
        Z = state[3*order]
        return np.r_[2*c/vol-g*du, v/vol-vol*e*e*u-g*c, -2*vol*e*e*c-g*dv,
                     a*loss, weights@source, -3*h*(weights@press), -h*Z/a, loss, weights@diag]

    sol = solve_ivp(rhs, (CUTS[0], CUTS[-1]), init, method="DOP853",
                    t_eval=CUTS, rtol=tol, atol=tol*1e-3, max_step=.0025)
    if not sol.success:
        raise ValueError(sol.message)
    out = {}
    for t, state in zip(sol.t,sol.y.T,strict=True):
        u,c,v,du,dv,a,vol,h,e,g,exc,press,src,diag,loss=terms(t,state)
        Z,W,WX,WR,J,WD = state[3*order:]
        U = weights@exc
        R = Z/a
        checks = (abs(U-W-WX+J), abs(R-J-WR), abs(U+R-W-WX-WR))
        det = u*v-c*c
        if min(det)<.25-1e-7 or min(exc/e)<-1e-8 or max(checks)>2e-7*(1+abs(W)+abs(WX)) or Z<-1e-9:
            raise ValueError((t, min(det), min(exc/e), checks))
        out[str(t)] = dict(chi_comoving_energy=U, daughter_comoving_energy=R,
                            daughter_physical_density=R/vol, a4_rho_daughter=Z,
                            band_comoving_number=weights@(exc/e), source_work=W,
                            diagonal_source_work=WD, coherence_source_work=W-WD,
                            chi_expansion_work=WX, daughter_expansion_work=WR,
                            decay_energy_transferred=J, min_covariance_determinant=min(det),
                            max_balance_residual=max(checks), scale_factor=a)
    return out


def calculate(capsule):
    output = {}
    for mesh in ("coarse","fine"):
        for y in (0.,1.):
            key = f"{mesh}/y{y:g}"
            output[key] = {b:solve(background(capsule,mesh,b),y,32,1.5e-9) for b in BRANCHES}
    for order in (24,48):
        key = f"refine{order}/y1"
        output[key] = {b:solve(background(capsule,"fine",b),1.,order,2e-10) for b in BRANCHES}
    output["fine/y0.5"] = {b:solve(background(capsule,"fine",b),.5,32,1.5e-9) for b in BRANCHES}
    output["fine/y2"] = {b:solve(background(capsule,"fine",b),2.,32,1.5e-9) for b in BRANCHES}
    return output


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("previous_capsule", type=Path)
    print(json.dumps(calculate(p.parse_args().previous_capsule), indent=2, allow_nan=False))
