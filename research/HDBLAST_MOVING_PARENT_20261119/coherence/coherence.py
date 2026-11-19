"""Phase-sensitive pair response for homogeneous Gaussian scalar preparations."""
import argparse
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.optimize import brentq

PRIOR = "research/HDBLAST_TAGGED_PARENT_20261119/primary/tagged.py"
PRIOR_HASH = "0acfcfbcd170c7b98a2cc21000747353db3538eb340360e8c3e448fb3a88db79"


def prior(repo):
    path = repo / PRIOR
    if hashlib.sha256(path.read_bytes()).hexdigest() != PRIOR_HASH:
        raise ValueError("Unexpected historical wall loader")
    spec = importlib.util.spec_from_file_location("rest_parent", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def spectrum(geom, start, momenta, intervals, samples, pad):
    a0, _, mass0 = geom(start)
    omega0 = np.hypot(momenta, a0*mass0)
    initial_v = 1/np.sqrt(2*omega0)
    count = len(momenta)
    def rhs(t, state):
        v, w, clock = state[:count], state[count:-1], state[-1].real
        a, h, mass = geom(clock)
        return np.r_[a*h*v+w, -(momenta**2+(a*mass)**2)*v-a*h*w, a]
    def stop(t, state):
        return state[-1].real - 2.10
    stop.terminal = True
    stop.direction = 1
    sol = solve_ivp(rhs, (0, 2*(2.10-start)/a0), np.r_[initial_v, -1j*omega0*initial_v, start],
                    events=stop, method="DOP853", dense_output=True, rtol=1e-10, atol=1e-12, max_step=.003)
    if not sol.success or len(sol.t_events[0])!=1:
        raise ValueError("Parent evolution")
    t1 = sol.t_events[0][0]
    def eta(s):
        return brentq(lambda t: sol.sol(t)[-1].real-s, 0, t1)
    af = geom(2.10)[0]
    out = {}
    nodes, weights = np.polynomial.legendre.leggauss(96)
    for window, (lo,hi) in intervals.items():
        left, right = eta(lo), eta(hi)
        delta = (right-left)/samples
        t = delta*np.arange(samples)
        state = sol.sol(left+t)
        v = state[:count]
        f = np.sin(np.pi*(state[-1].real-lo)/(hi-lo))**4
        g = f*v
        results = {key: [] for key in ("number", "energy", "screened_number", "screened_energy")}
        for q, wave in zip(momenta, g):
            # Measure frequencies from the exact light-cone edge. Opposite phases
            # from the common Fourier time origin cancel in F(+w) F(-w).
            M = samples*pad
            freq = q+2*np.pi*np.arange(M//2)/(M*delta)
            shifted = np.exp(1j*q*t)
            pos = M*delta*np.fft.ifft(wave*shifted,n=M)[:M//2]
            neg = delta*np.fft.fft(wave*shifted.conj(),n=M)[:M//2]
            def channels(A,B):
                return np.stack((abs(A)**2+abs(B)**2, A*B))
            align = (freq**2-q*q)/(8*np.pi**2)*channels(pos,neg)
            for name, pref in (("number", 1.),("energy",freq/af)):
                results[name].append(simpson(align*pref, x=freq, axis=-1))
            # An independent fixed-interval Gaussian integral ends exactly at 60.
            w = q+(nodes+1)*(60-q)/2
            kernel = np.exp(1j*w[:,None]*t)
            vals = (w**2-q*q)*weights*(60-q)/(16*np.pi**2)*channels(delta*(kernel@wave),
                                                                        delta*(kernel.conj()@wave))
            for name, pref in (("screened_number", 1.),("screened_energy",w/af)):
                results[name].append((vals*pref).sum(axis=1))
        out[window] = {key: np.asarray(val) for key,val in results.items()}
    return out, float(af)


def response_weights(q_weights, burst, quiet):
    d_b, z_b = burst[:,0].real, burst[:,1]
    d_q, z_q = quiet[:,0].real, quiet[:,1]
    if np.any(2*abs(z_q) > d_q+1e-10) or np.any(2*abs(z_b) > d_b+1e-10):
        raise ValueError("Pair kernel fails the Cauchy bound")
    n = 1.0
    cap = np.sqrt(n*(n+1))
    norm = q_weights@(n*d_q)
    d = q_weights@(n*(d_b-d_q))
    dz = z_b-z_q
    change_radius = 2*cap*(q_weights@abs(dz))
    common_radius = 2*cap*abs(q_weights@dz)
    def outcome(c):
        # Directly evaluate the Hermitian emitted-pair quadratic form. The
        # vacuum contribution cancels from the quoted excess on each wall.
        q = q_weights@(n*d_q+2*np.real(c*z_q))
        b = q_weights@(n*d_b+2*np.real(c*z_b))
        return dict(quiet=float(q), pulse=float(b), difference=float(b-q),
                    difference_percent_of_incoherent_quiet=float(100*(b-q)/norm),
                    difference_percent_of_own_quiet=float(100*(b-q)/q) if q!=0 else None)
    align = cap*np.exp(-1j*np.angle(dz))
    high, low = outcome(align), outcome(-align)
    if not np.allclose([low["difference"],high["difference"]],[d-change_radius,d+change_radius],
                       rtol=2e-9,atol=2e-13):
        raise ValueError("Explicit states fail to attain the analytic bound")
    phi = float(-np.angle(q_weights@dz))
    return dict(incoherent=outcome(0.), maximizing_state=high, minimizing_state=low,
                maximizing_constant_phase_radians=phi,
                constant_phases={str(theta):outcome(cap*np.exp(1j*theta))
                                 for theta in (-.5,0.,np.pi/2,np.pi,3*np.pi/2)},
                constant_phase_min_percent=float(100*(d-common_radius)/norm),
                constant_phase_max_percent=float(100*(d+common_radius)/norm),
                minimum_coherence_fraction_to_cancel=float(abs(d)/change_radius) if change_radius else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repo",type=Path)
    parser.add_argument("--nodes", type=int,default=24)
    parser.add_argument("--samples", type=int,default=1024)
    parser.add_argument("--padding", type=int,default=16)
    args=parser.parse_args()
    if args.nodes<4 or args.samples<256 or args.padding<8:raise ValueError("Resolution")
    mod=prior(args.repo)
    x,weights=np.polynomial.legendre.leggauss(args.nodes)
    q=4*(1+x)
    weights=4*weights*q*q/(8**3/3)
    out={"nodes":args.nodes,"samples":args.samples,"padding":args.padding,"results":{}}
    with zipfile.ZipFile(io.BytesIO(mod.checked(args.repo/mod.WALL,mod.WALL_SHA))) as z:
        for grid in ("fine","coarse"):
            ref=json.loads(mod.checked(args.repo/mod.REF.format(grid),mod.REF_SHA[grid]))["matter"]["events"][0]["A_star"]
            modes={}
            scale={}
            for arm in mod.ARMS:
                geom,start=mod.wall(z,grid,arm,ref)
                modes[arm],scale[arm]=spectrum(geom,start,q,mod.SWITCHES,args.samples,args.padding)
            for arm,windows in modes.items():
                for vals in windows.values():
                    vals["physical_number_density"]=vals["number"]*(8**3/(6*np.pi**2))/scale[arm]**3
                    vals["physical_energy_density"]=vals["energy"]*(8**3/(6*np.pi**2))/scale[arm]**3
            out["results"][grid]={arm:{win:{key:response_weights(weights,data,modes["quiet"][win][key])
                          for key,data in mode.items()} for win,mode in modes[arm].items()}
                          for arm in ("plus","minus","far")}
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
