"""Tagged scalar at rest coupled perturbatively to a massless Dirac pair."""
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicHermiteSpline, CubicSpline
from scipy.optimize import brentq, minimize_scalar

WALL = "research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip"
WALL_SHA = "7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15"
REF = "research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs/{}/quiet/rolling_summary.json"
REF_SHA = {"fine": "abec0887ee4de00ba0757e7d96987610f36b7157a85bc4f00696fee8064a6810",
           "coarse": "91601211fce4a617838d3a709849a138ec29a23998cc7c3ba2416ec0b0b8aa8c"}
ARMS = ("quiet", "plus", "minus", "far")
SWITCHES = {"wide": (1.15, 2.05), "narrow": (1.25, 1.95)}
CUT = 2.10
SCREEN = 60.


def checked(path, sha):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != sha:
        raise ValueError(f"Input digest mismatch for {path}")
    return data


def wall(archive, grid, arm, reference):
    def read(key):
        with archive.open(f"{grid}/y1_b1/{key}.csv") as stream:
            return np.genfromtxt(stream, delimiter=",", names=True)
    data = np.concatenate((read("pre"), read(arm)))
    clock = data["H0tau"]
    log_a = CubicHermiteSpline(clock, data["A_b"]-reference, data["H_over_H0"])
    mass = CubicHermiteSpline(clock, 100*(data["phi_b"]-.5), 100*data["v_over_H0"])
    def geometry(s):
        return np.exp(log_a(s)), log_a(s, 1), mass(s)
    return geometry, clock[0]


def field(geometry, initial):
    a0, _, mass0 = geometry(initial)
    e0 = abs(mass0)
    x0 = 1/np.sqrt(2*a0**3*e0)
    p0 = -1j*a0**3*e0*x0
    def proper_rhs(s, state):
        x, p, eta = state
        a, _, mass = geometry(s)
        return p/a**3, -a**3*mass**2*x, 1/a
    def conformal_rhs(eta, state):
        v, w, s = state
        a, h, mass = geometry(s.real)
        return a*h*v+w, -(a*mass)**2*v-a*h*w, a
    def stop(eta, state):
        return state[-1].real-CUT
    stop.terminal = True
    stop.direction = 1
    options = dict(method="DOP853", rtol=1e-11, atol=1e-13, max_step=.002, dense_output=True)
    physical = solve_ivp(proper_rhs, (initial, CUT), [x0, p0, 0], **options)
    conformal = solve_ivp(conformal_rhs, (0, 2*(CUT-initial)/a0), [a0*x0, p0/a0, initial],
                          events=stop, **options)
    if not physical.success or not conformal.success or len(conformal.t_events[0]) != 1:
        raise ValueError("Mode evolution failed")
    def eta(s):
        return brentq(lambda t: conformal.sol(t)[-1].real-s, 0, conformal.t[-1])
    def cf(t):
        v, w, s = conformal.sol(t)
        a, h, mass = geometry(s.real)
        return s.real, a, v, w+a*h*v
    def ph(s):
        x, p, t = physical.sol(s)
        a, h, mass = geometry(s)
        v = a*x
        return a, mass, v, p/a+a*h*v
    ca, _, cv, cd = ph(CUT)
    _, _, vcheck, dcheck = cf(conformal.t[-1])
    tgrid = np.linspace(0, conformal.t[-1], 2049)
    _, _, v, d = cf(tgrid)
    wronskian = float(np.max(np.abs(v*d.conj()-v.conj()*d-1j)))
    clock_error = abs(conformal.t[-1]-physical.y[-1,-1].real)
    relative_state_error = max(abs(cv-vcheck)/abs(cv),abs(cd-dcheck)/abs(cd))
    if wronskian > 1e-8 or clock_error > 1e-8 or relative_state_error > 2e-7:
        raise ValueError("Clock or Wronskian comparison failed")
    return eta, cf, ph, float(ca), dict(max_wronskian_error=wronskian,
                  final_eta_difference=float(clock_error), relative_final_mode_error=float(relative_state_error))


def switch(s, start, end):
    # The sampled field is identically zero outside the chosen support.
    t = np.clip((s-start)/(end-start), 0, 1)*np.pi
    sn = np.sin(t)
    return sn**4, 4*np.pi/(end-start)*sn**3*np.cos(t)


def spectral(g, step, padding, a_cut):
    # F(+omega) = integral g(eta) exp(+i omega eta) d eta.
    size = len(g)*padding
    F = np.fft.ifft(g, n=size)*(size*step)
    omega = 2*np.pi*np.arange(size//2)/(size*step)
    d_omega = omega[1]
    indices = np.arange(len(omega))
    full = {}
    screened = {}
    for key, multiplier in (("annihilation",1), ("stimulated",-1)):
        amplitude = F[multiplier*indices]
        dn = omega**2*abs(amplitude)**2/(8*np.pi**2)
        for name, density in (("number",dn), ("energy",dn*omega/a_cut)):
            # At DC the integrand vanishes. Nyquist is dropped and the
            # same physical frequency screen is integrated on every grid.
            full[f"{key}_{name}"] = float(np.sum(density)*d_omega)
            end = int(np.ceil(SCREEN/d_omega))+3
            interp = CubicSpline(omega[:end], density[:end])
            screened[f"{key}_{name}"] = float(interp.integrate(0, SCREEN))
    for d in (full,screened):
        d["number"] = d["annihilation_number"]+d["stimulated_number"]
        d["energy"] = d["annihilation_energy"]+d["stimulated_energy"]
        d["net_parent_loss"] = d["annihilation_number"]-d["stimulated_number"]
    return full,screened,dict(freq_spacing=float(d_omega),nyquist=float(np.pi/step))


def window(geometry, clocks, limits, samples, pad):
    eta, cf, ph, a_cut, _ = clocks
    lo, hi = limits
    left, right = eta(lo), eta(hi)
    step = (right-left)/samples
    s,a,v,d = cf(left+np.arange(samples)*step)
    f, fs = switch(s,lo,hi)
    full, screened, freq = spectral(f*v, step, pad, a_cut)
    s = np.linspace(lo,hi,2*samples+1)
    a,mu,v,d = ph(s)
    f, fs = switch(s,lo,hi)
    fe = a*fs
    def time_integral(integrand):
        return float(simpson(integrand/(a*4*np.pi), x=s))
    terms = dict(field=time_integral(abs(f*d)**2), switch=time_integral(abs(fe*v)**2),
                 interference=time_integral(2*f*fe*np.real(v.conj()*d)))
    exact = time_integral(abs(f*d+fe*v)**2)
    error = full["number"]-exact
    if abs(error) > 2e-8*exact or abs(sum(terms.values())-exact) > 1e-11*exact:
        raise ValueError("Time and spectral pair counts disagree")
    # Direct minimization includes boundary values directly. The archived
    # wall expands throughout these windows and h>0 on the control grid.
    if np.min(geometry(s)[1]) <= 0:
        raise ValueError("The screening bound needs a new minimum search")
    amin = min(a[0],a[-1], minimize_scalar(lambda t: geometry(t)[0],bounds=(lo,hi),
                       method="bounded",options={"xatol":1e-13}).fun)
    if SCREEN/(2*amin) >= 53.6107:
        raise ValueError("Daughter screen exceeds the reference UV scale")
    return dict(full=full, screened=screened,
                fraction_above_screen={k:1-screened[k]/full[k] for k in ("number","energy")},
                time_identity=dict(exact=exact,spectral_minus_time=error,terms=terms),
                instantaneous_rate_estimate=float(simpson(f**2*abs(mu)/(8*np.pi),x=s)),
                min_a=float(amin),max_screened_single_daughter_physical_momentum=float(SCREEN/(2*amin)),
                **freq)


def flat_control(samples, pad):
    # For a constant mass at rest Parseval yields the sum of the golden-rule
    # term and an analytic switching correction, also for a finite window.
    m, duration = 32., .9
    s = np.arange(samples)*duration/samples
    f,_ = switch(s,0,duration)
    g = f*np.exp(-1j*m*s)/np.sqrt(2*m)
    full, _, _ = spectral(g,duration/samples,pad,1.)
    # sin^8 has mean 35/128 and sin^6 cos^2 has mean 5/128.
    analytic = (m*(35/128)*duration+ (5/8)*np.pi**2/(m*duration))/(8*np.pi)
    if abs(full["number"]-analytic) > 1e-10 or full["annihilation_number"]<=full["stimulated_number"]:
        raise ValueError("Flat-space calibration failed")
    return dict(mass=m,duration=duration,analytic_number=analytic,**full)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("repository",type=Path)
    p.add_argument("--samples",type=int,default=2048)
    p.add_argument("--padding",type=int,default=32)
    args=p.parse_args()
    if args.samples<512 or args.padding<8 or args.samples%2:
        raise ValueError("Use an even, resolved switch grid and at least eightfold zero padding")
    out=dict(switches=SWITCHES,cut=CUT,omega_screen=SCREEN,samples=args.samples,padding=args.padding,
             flat=flat_control(args.samples,args.padding),grids={})
    with zipfile.ZipFile(io.BytesIO(checked(args.repository/WALL,WALL_SHA))) as archive:
        for grid in REF_SHA:
            ref=json.loads(checked(args.repository/REF.format(grid),REF_SHA[grid]))["matter"]["events"][0]["A_star"]
            data={}
            for arm in ARMS:
                geom,start=wall(archive,grid,arm,ref)
                clocks=field(geom,start)
                data[arm] = dict(initial_s=float(start),cut_a=clocks[3],clocks=clocks[4])
                for name,lim in SWITCHES.items():
                    data[arm][name]=window(geom,clocks,lim,args.samples,args.padding)
            for arm in ARMS:
                data[arm]["contrasts_percent"]={win:{band:{
                    k:100*(val/data["quiet"][win][band][k]-1) for k,val in data[arm][win][band].items()
                    } for band in ("full","screened")} for win in SWITCHES}
            out["grids"][grid]=data
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
