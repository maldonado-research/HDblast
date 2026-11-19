"""Inclusive Dirac-pair spectra from moving scalar parents on prescribed B1 histories."""
import argparse
import hashlib
import importlib.util
import io
import json
import zipfile
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicHermiteSpline, CubicSpline

OLD = Path("research/HDBLAST_TAGGED_PARENT_20261119/primary")
DRIVER_HASH = "0acfcfbcd170c7b98a2cc21000747353db3538eb340360e8c3e448fb3a88db79"
REST_HASH = "a695a67bd623d8dd473b4ac326678d1570c4c19b05871fe3b8847efd2e5d6d25"
K = 8.
SCREEN = 60.
M5 = 6.49e-6**(-1/3)
SPOTS = np.array([0., 2., 4., 6., 8.])
COMPONENTS = ("ann_number", "stim_number", "ann_energy", "stim_energy")


def load_prior(repo):
    path = repo/OLD/"tagged.py"
    if hashlib.sha256(path.read_bytes()).hexdigest() != DRIVER_HASH:
        raise ValueError("Unexpected rest-study source")
    spec = importlib.util.spec_from_file_location("published_tagged", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rest = json.loads(module.checked(repo/OLD/"TAGGED.json", REST_HASH))
    return module, rest


def fields(geom, start, cut, q):
    a0, _, m0 = geom(start)
    e0 = np.hypot(q/a0, m0)
    x0 = 1/np.sqrt(2*a0**3*e0)
    n = len(q)
    def rhs(s, state):
        a, _, m = geom(s)
        x, pi = state[:n], state[n:2*n]
        return np.r_[pi/a**3, -(a**3*m**2+a*q**2)*x, 1/a]
    initial = np.r_[x0, -1j*(a0**3*e0)*x0, 0j]
    sol = solve_ivp(rhs, (start,cut), initial, method="DOP853", dense_output=True,
                    rtol=3e-11, atol=5e-13, max_step=.003)
    if not sol.success:
        raise ValueError(sol.message)
    clock = np.linspace(start,cut,4097)
    eta = sol.sol(clock)[-1].real
    s_of_eta = CubicHermiteSpline(eta,clock,geom(clock)[0],extrapolate=False)
    def sample(s):
        a,h,m = geom(s)
        z = sol.sol(s)
        v = a*z[:n]
        prime = z[n:2*n]/a + a*h*v
        return a, v, prime
    def at_eta(t):
        s = s_of_eta(t)
        return s, *sample(s)
    return lambda s: sol.sol(s)[-1].real, at_eta, sample


def spectral(q, waveform, step, padding, a_final):
    size = len(waveform)*padding
    f = np.fft.ifft(waveform, n=size)*size*step
    # Integrate both sides including the negligible Nyquist endpoint; the
    # true spectral threshold is q and generally lies between FFT nodes.
    idx = np.arange(size//2+1)
    omega = 2*np.pi*idx/(size*step)
    dw = omega[1]
    low = np.searchsorted(omega, SCREEN)+5
    full, screened = {}, {}
    removed = 0.
    for sign, channel in ((1, "ann"), (-1, "stim")):
        sq = abs(f[(sign*idx)%size])**2/(8*np.pi**2)
        number_density = (omega**2-q*q)*sq
        for name, integrand in (("number", number_density), ("energy", number_density*omega/a_final)):
            small = CubicSpline(omega[:low],integrand[:low])
            space = float(small.integrate(0,q))
            key = channel+"_"+name
            full[key] = float(simpson(integrand,dx=dw))-space
            screened[key] = float(small.integrate(q,SCREEN))
            if name == "number":
                removed -= space
            if not 0 <= screened[key] <= full[key]+1e-9:
                raise ValueError((q,key,screened[key],full[key]))
    return full,screened,removed


def measure(component):
    out = dict(component)
    out["number"] = out["ann_number"]+out["stim_number"]
    out["energy"] = out["ann_energy"]+out["stim_energy"]
    out["net_parent_loss"] = out["ann_number"]-out["stim_number"]
    return out


def read_window(old, geom, clocks, lo, hi, q, weights, samples, pad):
    eta, at_eta, proper = clocks
    a_final = float(geom(old.CUT)[0])
    left,right = eta(lo),eta(hi)
    dt = (right-left)/samples
    s, a,v,d = at_eta(left+dt*np.arange(samples))
    f, _ = old.switch(s,lo,hi)
    gs = v*f
    st = np.linspace(lo,hi,2*samples+1)
    a,v,d = proper(st)
    ff,fs = old.switch(st,lo,hi)
    loc = simpson((abs(ff*d+a*fs*v)**2-q[:,None]**2*abs(ff*v)**2)/(4*np.pi*a), x=st,axis=1)
    unit = float(np.max(abs(v*d.conj()-v.conj()*d-1j)))
    clock_error = float(np.max(abs(eta(s)-(left+dt*np.arange(samples)))))
    if unit > 2e-7 or clock_error > 5e-10 or np.min(geom(st)[1]) <= 0:
        raise ValueError(("time evolution",unit,clock_error))
    amin = float(geom(lo)[0])
    pmax = (SCREEN+K)/(2*amin)
    if pmax >= M5:
        raise ValueError("Per-daughter screen exceeds the reference M5")
    full, screened, spacelike = [], [], []
    spectrum = {}
    for j, momentum in enumerate(q):
        frow,srow,space = spectral(momentum,gs[j],dt,pad,a_final)
        full.append(frow)
        screened.append(srow)
        spacelike.append(space)
        if j < len(SPOTS):
            spectrum[f"q{momentum:g}"] = dict(full=measure(frow), screened=measure(srow),
                 local_count_without_threshold=float(loc[j]),positive_threshold_correction=space)
    error = max(abs(f["ann_number"]+f["stim_number"]-l-sp) for f,l,sp in zip(full,loc,spacelike))
    if error > 2e-7:
        raise ValueError(("Parseval identity",error))
    for label, weight in weights.items():
        norm = float(np.sum(weight))
        data = {}
        for band, rows in (("full",full),("screened",screened)):
            comp = {key:float(weight@[r[key] for r in rows[len(SPOTS):]]/norm) for key in COMPONENTS}
            result = measure(comp)
            result["number_density"] = norm/(2*np.pi**2*a_final**3)*result["number"]
            result["energy_density"] = norm/(2*np.pi**2*a_final**3)*result["energy"]
            result["isotropic_pressure"] = result["energy_density"]/3
            data[band] = result
        data["local_count_without_threshold"] = float(weight@loc[len(SPOTS):]/norm)
        data["positive_threshold_correction"] = float(weight@spacelike[len(SPOTS):]/norm)
        spectrum[label]=data
    return dict(spectrum=spectrum,a_final=a_final,min_scale_during_switch=amin,
                max_daughter_momentum_for_screen=pmax,max_parseval_error=float(error),
                max_wronskian_error=unit,max_clock_error=clock_error)


def simulate(repo, nq=16, samples=1024, padding=32):
    old,rest = load_prior(repo)
    k0, w = np.polynomial.legendre.leggauss(nq)
    k = K/2*(1+k0)
    w = w*K/2*k**2
    weights = {"uniform":w, "tapered":w*(1-(k/K)**2)**2}
    q = np.r_[SPOTS,k]
    output = dict(q_nodes=nq,samples=samples,padding=padding,upper_parent_momentum=K,
                  screen_pair_energy=SCREEN,cut_s=old.CUT,rest_driver_sha256=DRIVER_HASH,
                  rest_result_sha256=REST_HASH,
                  comoving_parent_density={s:float(sum(v)/(2*np.pi**2)) for s,v in weights.items()},
                  grids={})
    for label, analytical in (("uniform", K**3/3),("tapered",8*K**3/105)):
        if abs(sum(weights[label])-analytical)>2e-11:
            raise ValueError("Incorrect density measure")
    rest_error = 0.
    with zipfile.ZipFile(io.BytesIO(old.checked(repo/old.WALL,old.WALL_SHA))) as archive:
        for grid, pin in old.REF_SHA.items():
            ref=json.loads(old.checked(repo/old.REF.format(grid),pin))["matter"]["events"][0]["A_star"]
            cases = {}
            for arm in old.ARMS:
                geom,start = old.wall(archive,grid,arm,ref)
                clocks = fields(geom,start,old.CUT,q)
                windows = {}
                for win,(lo,hi) in old.SWITCHES.items():
                    row=read_window(old,geom,clocks,lo,hi,q,weights,samples,padding)
                    windows[win] = row
                    for band in ("full","screened"):
                        for key in ("number","energy","net_parent_loss"):
                            val = row["spectrum"]["q0"][band][key]
                            refval = rest["grids"][grid][arm][win][band][key]
                            rest_error=max(rest_error,abs(val/refval-1))
                cases[arm] = windows
            output["grids"][grid] = cases
    if rest_error>2e-7:
        raise ValueError(("rest-limit disagreement",rest_error))
    output["max_relative_error_vs_published_rest"] = rest_error
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repository",type=Path)
    parser.add_argument("--nodes",type=int,default=16)
    parser.add_argument("--samples",type=int,default=1024)
    parser.add_argument("--padding",type=int,default=32)
    args=parser.parse_args()
    if args.nodes<6 or args.samples<512 or args.samples%2 or args.padding<16:
        raise ValueError("Insufficient momentum or Fourier resolution")
    print(json.dumps(simulate(args.repository,args.nodes,args.samples,args.padding),indent=2,allow_nan=False))


if __name__=="__main__":
    main()
