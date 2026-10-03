#!/usr/bin/env python3
"""Fixed comoving-band scalar modes on archived HDBLAST backgrounds.
This is a finite regulator/state diagnostic, not a renormalized physical source.
Inputs are read without pickle; no background evolution or decay is added.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import CubicSpline, CubicHermiteSpline
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

PINNED = "d2dea40857e1899019ec0d9872998bbede11c459"
G, PHISTAR, MASS0 = 100.0, 0.5, 0.0
END = 6.9
FINE = "main_dstar_Y0_dc1e-2_dzf5e-4"
COARSE = "main_dstar_Y0_dc1e-2_dzf1e-3"
LIMITS = dict(wronskian=1e-6, endpoint_ward=1e-5,
              clock_occupation_abs=1e-5, conformal_mode_abs=1e-6)

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (np.floating, float)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x

def load(root, tag):
    sp = root / (tag + "_summary.json")
    pp = root / (tag + "_timeseries.npz")
    summary = json.loads(sp.read_text())
    if summary.get("restart_info"):
        raise ValueError("This registered A1 test requires a complete non-restarted history")
    if float(summary["params"]["Y"]) != 0:
        raise ValueError("Production background must have Y=0")
    if abs(float(summary["params"]["delta"]) - 0.1) > 1e-12:
        raise ValueError("This test uses the tuned delta=0.1 variant")
    with np.load(pp, allow_pickle=False) as z:
        cols, rec = [str(c) for c in z["cols"]], z["rec"]
        if rec.ndim != 2 or rec.shape[1] != len(cols):
            raise ValueError("Invalid cols/rec schema")
        if len(set(cols)) != len(cols):
            raise ValueError("Duplicate columns")
        data = {c: rec[:, j].copy() for j, c in enumerate(cols)}
    required = ("H0tau", "ln_a", "H_over_H0", "phi_b", "v_over_H0",
                "Hmax_near", "Mmax_near")
    if any(k not in data for k in required):
        raise ValueError("Missing proper-time, geometry, field or constraint columns")
    t = data["H0tau"]
    if np.any(np.diff(t) < -1e-9):
        raise ValueError("Nonmonotone proper time")
    bad = np.flatnonzero((data["Hmax_near"] > 0.05) |
                         (data["Mmax_near"] > 0.05))
    last = int(bad[0]) if len(bad) else len(t)
    finite = np.logical_and.reduce([np.isfinite(data[k]) for k in required[:5]])
    failures = np.flatnonzero(~finite)
    if len(failures):
        last = min(last, int(failures[0]))
    keep = np.arange(len(t)) < last
    keep &= np.concatenate(([True], np.diff(t) > 1e-9))
    data = {k: v[keep] for k, v in data.items()}
    if len(data["H0tau"]) < 5:
        raise ValueError("Insufficient reliable records")
    if data["H0tau"][0] > 1e-10 or data["H0tau"][-1] < END:
        raise ValueError("Reliable archived history does not span [0,6.9]")
    if not np.any(np.isfinite(data["Hmax_near"])):
        raise ValueError("No finite archived constraint checks")
    meta = dict(tag=tag, summary_sha256=sha(sp), npz_sha256=sha(pp),
                summary_path=str(sp), npz_path=str(pp), params=summary["params"],
                rho_b=summary["rho_b"], H0=1/float(summary["rho_b"]),
                n_reliable=len(data["H0tau"]),
                reliable_end=float(data["H0tau"][-1]),
                first_bad_original_index=int(bad[0]) if len(bad) else None,
                stop_reason=summary.get("stop_reason"),
                archived_phi_range=[float(data["phi_b"].min()),
                                    float(data["phi_b"].max())],
                max_saved_spacing=float(np.diff(data["H0tau"]).max()))
    return data, meta

class Geometry:
    def __init__(self, data, interpolation):
        self.t = data["H0tau"]
        loga = data["ln_a"] - data["ln_a"][0]
        if interpolation == "hermite":
            self.la = CubicHermiteSpline(self.t, loga, data["H_over_H0"],
                                        extrapolate=False)
            self.ph = CubicHermiteSpline(self.t, data["phi_b"],
                                        data["v_over_H0"], extrapolate=False)
        elif interpolation == "cubic":
            self.la = CubicSpline(self.t, loga, extrapolate=False)
            self.ph = CubicSpline(self.t, data["phi_b"], extrapolate=False)
        else:
            raise ValueError(interpolation)
        self.interpolation = interpolation

    def __call__(self, t):
        a = np.exp(self.la(t))
        h, hd = self.la(t, 1), self.la(t, 2)
        ph, v = self.ph(t), self.ph(t, 1)
        return a, h, hd, ph, v

    def crossings(self):
        roots = self.ph.solve(PHISTAR, extrapolate=False)
        roots = roots[(roots >= self.t[0]) & (roots <= END)]
        roots = np.unique(np.round(roots, 13))
        return [dict(t=float(t), phi=float(self.ph(t)),
                     v=float(self.ph(t, 1)), a=float(np.exp(self.la(t))))
                for t in roots if abs(float(self.ph(t, 1))) > 1e-10]

def mass(ph):
    return MASS0**2 + G**2 * (ph - PHISTAR)**2

def initial(bg, t, k, state):
    a, h, hd, ph, v = bg(t)
    w = np.sqrt(k*k/a**2 + mass(ph))
    wd = (-h*k*k/a**2 + G*G*(ph-PHISTAR)*v) / w
    c = 1 / np.sqrt(2*a**3*w)
    if state == "wkb":
        d = (-1.5*h - wd/(2*w) - 1j*w)*c
    elif state == "hamiltonian":
        d = -1j*w*c
    else:
        raise ValueError(state)
    return c.astype(complex), d.astype(complex)

def stress(bg, t, k, c, d, xi=0.0):
    a, h, hd, ph, v = bg(t)
    f = np.abs(c)**2
    fd = 2*np.real(np.conj(c)*d)
    ds = np.abs(d)**2
    kk, mm = k*k/a**2, mass(ph)
    r4 = 6*(hd+2*h*h)
    fdd = 2*ds - 3*h*fd - 2*(kk+mm+xi*r4)*f
    rho = 0.5*(ds+(kk+mm)*f) + xi*(3*h*h*f+3*h*fd)
    p = 0.5*ds-(kk/6+mm/2)*f
    p += xi*(-(2*hd+3*h*h)*f-fdd-2*h*fd)
    current = G*G*(ph-PHISTAR)*f
    return rho, p, current

def occupation(bg, t, k, c, d):
    a, h, hd, ph, v = bg(t)
    w = np.sqrt(k*k/a**2+mass(ph))
    # beta refers to the physical Hamiltonian basis at this instant.
    scale = np.sqrt(a**3/(2*w))
    beta = scale*(w*c-1j*d)
    alpha = scale*(w*c+1j*d)
    n = np.abs(beta)**2
    coherence = 2*np.real(alpha*np.conj(beta))
    # chi^2=(1+2n+coherence)/(2 a^3 omega), up to Wronskian error.
    f_reconstruction = (np.abs(alpha)**2+n+coherence)/(2*a**3*w)
    return w, n, coherence, f_reconstruction

def integrate(bg, t0, t1, k, weights, state, rtol, grid):
    n = len(k)
    c0, d0 = initial(bg, t0, k, state)
    y0 = np.concatenate([c0, d0, [0j, 0j, 0j]])
    probe = np.unique(np.concatenate([bg.t[(bg.t >= min(t0,t1)) &
                                              (bg.t <= max(t0,t1))],
                                     np.linspace(min(t0,t1), max(t0,t1), 501)]))
    a, h, hd, ph, v = bg(probe)
    phase_bound = float(np.sqrt(k.max()**2/a**2+mass(ph)).max())
    max_step = 0.15/max(phase_bound, 1.0)
    def rhs(t, y):
        c, d = y[:n], y[n:2*n]
        a, h, hd, ph, v = bg(t)
        dd = -3*h*d-(k*k/a**2+mass(ph))*c
        rho, p, j = stress(bg, t, k, c, d)
        scalar_power = a**3*np.dot(weights, j*v)
        pressure_power = a**3*np.dot(weights, -3*h*p)
        return np.concatenate([d, dd, [scalar_power, pressure_power,
                                      abs(scalar_power)+abs(pressure_power)]])
    sol = solve_ivp(rhs, (t0,t1), y0, method="DOP853", rtol=rtol,
                    atol=rtol*0.01, max_step=max_step, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    y = sol.sol(grid)
    c, d = y[:n], y[n:2*n]
    # Evaluate vectorized geometry/stress with time axis first.
    rho, p, j = stress(bg, grid[:,None], k[None,:], c.T, d.T)
    rr, pp, jj = [z @ weights for z in (rho,p,j)]
    a, h, hd, ph, v = bg(grid)
    wr = a[None,:]**3*(c*np.conj(d)-np.conj(c)*d)/1j
    i0 = int(np.argmin(np.abs(grid-t0)))
    energy = a**3*rr
    work = y[2*n].real+y[2*n+1].real
    balance = energy-energy[i0]-work
    scale = (np.abs(energy)+abs(energy[i0])+np.abs(y[2*n].real)
             +np.abs(y[2*n+1].real)+1e-30)
    relative = np.abs(balance)/np.maximum(scale, 1e-300)
    data = dict(c=c, d=d, rho=rr, p=pp, current=jj,
                work=work, scalar_work=y[2*n].real, pressure_work=y[2*n+1].real,
                abs_power_integral=y[2*n+2].real,
                energy=energy, balance=balance, relative=relative)
    report = dict(state=state, initial_time=t0, final_time=t1, nfev=sol.nfev,
                  phase_max_step=max_step,
                  wronskian_max=float(np.max(np.abs(wr-1))),
                  ward_endpoint_relative=float(relative[-1] if t1>t0 else relative[0]),
                  ward_path_max_relative=float(relative.max()),
                  scalar_work=float(y[2*n,-1].real if t1>t0 else y[2*n,0].real),
                  pressure_work=float(y[2*n+1,-1].real if t1>t0 else y[2*n+1,0].real))
    return data, report

def clock_subset(bg, start, end, k, state, rtol, primary):
    # Each canonical integration is split at archived spline knots.
    # Add eta to X-system to construct an accurate monotonically mapped clock.
    c0, d0 = initial(bg, start, k, state)
    a0, h0, hd0, ph0, v0 = bg(start)
    x0 = a0**1.5*c0
    xd0 = a0**1.5*(d0+1.5*h0*c0)
    breaks = np.unique(np.concatenate([[start,end],
                       bg.t[(bg.t>start)&(bg.t<end)]]))
    xs, etas = [], []
    y = np.concatenate([x0,xd0,[0j]])
    nt = len(k)
    nfev = 0
    for lo, hi in zip(breaks[:-1],breaks[1:]):
        def fx(t,z):
            a,h,hd,ph,v = bg(t)
            freq2=k*k/a**2+mass(ph)-1.5*hd-2.25*h*h
            return np.concatenate([z[nt:2*nt],-freq2*z[:nt],[1/a]])
        mids=np.array([lo,(lo+hi)/2,hi])
        a,h,hd,ph,v=bg(mids)
        mx=float(np.sqrt(np.max(k)**2/a**2+mass(ph)+1.5*np.abs(hd)+2.25*h*h).max())
        sol=solve_ivp(fx,(lo,hi),y,method="DOP853",rtol=rtol,atol=rtol*.01,
                      max_step=.15/max(mx,1.),dense_output=True)
        if not sol.success: raise RuntimeError(sol.message)
        nfev+=sol.nfev
        sample=np.linspace(lo,hi,9)
        if len(etas): sample=sample[1:]
        etas.extend(sol.sol(sample)[-1].real.tolist())
        xs.extend(sample.tolist())
        y=sol.y[:,-1]
    # Integrate the conformal equation with t as an auxiliary dependent variable.
    # Per-knot conformal durations use adaptive quadrature of this same geometry.
    u0=a0*c0
    up0=a0*a0*(d0+h0*c0)
    zu=np.concatenate([u0,up0,[complex(start)]])
    clock_endpoint_error=0.
    for lo, hi in zip(breaks[:-1],breaks[1:]):
        duration=quad(lambda tt:float(np.exp(-bg.la(tt))),lo,hi,
                      epsabs=1e-13,epsrel=1e-13)[0]
        zu[-1]=complex(lo)
        def fu(e,z):
            tt=float(z[-1].real)
            a,h,hd,ph,v=bg(tt)
            potential=k*k+a*a*(mass(ph)-hd-2*h*h)
            return np.concatenate([z[nt:2*nt],-potential*z[:nt],[a]])
        samples=np.array([lo,(lo+hi)/2,hi])
        a,h,hd,ph,v=bg(samples)
        mx=float(np.sqrt(k.max()**2+a*a*(mass(ph)+np.abs(hd)+2*h*h)).max())
        sol=solve_ivp(fu,(0.,duration),zu,method="DOP853",rtol=rtol,
                      atol=rtol*.01,max_step=.15/max(mx,1.))
        if not sol.success: raise RuntimeError(sol.message)
        nfev+=sol.nfev
        zu=sol.y[:,-1]
        clock_endpoint_error=max(clock_endpoint_error,abs(float(zu[-1].real)-hi))
    ae,he,hde,phe,ve=bg(end)
    cx=y[:nt]/ae**1.5
    dx=(y[nt:2*nt]-1.5*he*y[:nt])/ae**1.5
    cu=zu[:nt]/ae
    du=zu[nt:2*nt]/ae**2-he*cu
    nref=occupation(bg,end,k,primary[0],primary[1])[1]
    nx=occupation(bg,end,k,cx,dx)[1]
    nu=occupation(bg,end,k,cu,du)[1]
    wx=(y[:nt]*np.conj(y[nt:2*nt])-np.conj(y[:nt])*y[nt:2*nt])/1j
    wu=(zu[:nt]*np.conj(zu[nt:2*nt])-np.conj(zu[:nt])*zu[nt:2*nt])/1j
    return dict(k=k,n_primary=nref,n_X=nx,n_u=nu,nfev=nfev,
                max_occupation_abs_X=float(np.max(np.abs(nx-nref))),
                max_occupation_abs_u=float(np.max(np.abs(nu-nref))),
                wronskian_X=float(np.max(np.abs(wx-1))),
                wronskian_u=float(np.max(np.abs(wu-1))),
                clock_endpoint_error=clock_endpoint_error,
                map_control="Conformal time integration, dt/deta=a, per-knot eta duration from adaptive quadrature")

def conformal_control(bg, start, end, k, rtol):
    # X=sqrt(a) exp(-ik eta)/sqrt(2k), xi=1/6, m=0.
    # Integrate the exact mapped state independently in proper time.
    a,h,hd,ph,v=bg(start)
    x=np.sqrt(a/(2*k)).astype(complex)
    xd=(h/2-1j*k/a)*x
    z=np.concatenate([x,xd,[0j]])
    nt=len(k)
    breaks=np.unique(np.concatenate([[start,end],bg.t[(bg.t>start)&(bg.t<end)]]))
    for lo,hi in zip(breaks[:-1],breaks[1:]):
        def rhs(t,y):
            a,h,hd,ph,v=bg(t)
            ff=k*k/a**2-.5*hd-.25*h*h
            return np.concatenate([y[nt:2*nt],-ff*y[:nt],[1/a]])
        sol=solve_ivp(rhs,(lo,hi),z,method="DOP853",rtol=rtol,
                      atol=rtol*.01,max_step=min(.02, .15/max(k.max(),1.)))
        if not sol.success: raise RuntimeError(sol.message)
        z=sol.y[:,-1]
    a,h,hd,ph,v=bg(end)
    u=z[:nt]/np.sqrt(a)
    up=np.sqrt(a)*(z[nt:2*nt]-h*z[:nt]/2)
    eta=z[-1].real
    target=np.exp(-1j*k*eta)/np.sqrt(2*k)
    n=np.abs(np.sqrt(1/(2*k))*(k*u-1j*up))**2
    wr=(u*np.conj(up)-np.conj(u)*up)/1j
    return dict(k=k,occupation=n,
                relative_mode_max=float(np.max(np.abs(u-target)/np.abs(target))),
                wronskian_max=float(np.max(np.abs(wr-1))),
                occupation_max=float(n.max()),
                note="Mode/state control only; no xi=0 stress assertion is made.")

def run_case(bg, reference, nk, upper, start, state, rtol, outdir, name, clocks=False, end=END):
    x,wg=leggauss(nk)
    kap=.05+(x+1)*(upper-.05)/2
    dk=(upper-.05)/2*reference["k_scale"]
    k=reference["k_scale"]*kap
    weights=wg*dk*k*k/(2*np.pi**2)
    grid=np.unique(np.concatenate([np.linspace(start,end,801),
                                  bg.t[(bg.t>=start)&(bg.t<=end)]]))
    ins,ir=integrate(bg,start,end,k,weights,state,rtol,grid)
    outs,orr=integrate(bg,end,start,k,weights,"hamiltonian",rtol,grid)
    a,h,hd,ph,v=bg(grid)
    dr,dp,dj=[ins[z]-outs[z] for z in ("rho","p","current")]
    de=ins["energy"]-outs["energy"]
    dw=(ins["work"]-ins["work"][0])-(outs["work"]-outs["work"][0])
    balance=de-de[0]-dw
    dscalar=(ins["scalar_work"]-ins["scalar_work"][0])-(outs["scalar_work"]-outs["scalar_work"][0])
    dpressure=(ins["pressure_work"]-ins["pressure_work"][0])-(outs["pressure_work"]-outs["pressure_work"][0])
    scale=np.abs(de)+abs(de[0])+np.abs(dscalar)+np.abs(dpressure)+1e-30
    rel=np.abs(balance)/np.maximum(scale,1e-300)
    w,n,coh,frecon=occupation(bg,end,k,ins["c"][:,-1],ins["d"][:,-1])
    number=float(np.dot(weights,n)/a[-1]**3)
    particle_rho=float(np.dot(weights,w*n)/a[-1]**3)
    particle_p=float(np.dot(weights,k*k/a[-1]**2*n/(3*w))/a[-1]**3)
    f=np.abs(ins["c"][:,-1])**2
    current_coherence=G*G*(ph[-1]-PHISTAR)*coh/(2*a[-1]**3*w)
    current_particles=G*G*(ph[-1]-PHISTAR)*n/(a[-1]**3*w)
    # Instantaneous zero-point removal is a distinct finite-band prescription.
    # Test its complete triplet and the inconsistent rho/p-only removal separately.
    zero_r=np.empty(len(grid)); zero_p=np.empty(len(grid)); zero_j=np.empty(len(grid))
    zero_defect=np.empty(len(grid))
    for i,t in enumerate(grid):
        om=np.sqrt(k*k/a[i]**2+mass(ph[i]))
        omdot=(-h[i]*k*k/a[i]**2+G*G*(ph[i]-PHISTAR)*v[i])/om
        zero_r[i]=np.dot(weights,om)/(2*a[i]**3)
        zero_p[i]=np.dot(weights,k*k/a[i]**2/(6*om*a[i]**3))
        zero_j[i]=np.dot(weights,G*G*(ph[i]-PHISTAR)/(2*a[i]**3*om))
        zero_defect[i]=(np.dot(weights,omdot)/(2*a[i]**3)
                        +3*h[i]*zero_p[i]-zero_j[i]*v[i])
    # The Hamiltonian zero-point triplet has zero algebraic Ward defect, but it is
    # not the stress of the exact backward out state or an action-based subtraction.
    instantaneous_delta=ins["rho"]-zero_r
    normalizer=abs(de[0])+abs(de[-1])+abs(dscalar[-1])+abs(dpressure[-1])+1e-30
    raw_scale=max(abs(ins["energy"][0]),abs(ins["energy"][-1]),
                  abs(ins["work"][-1]),abs(outs["energy"][0]),
                  abs(outs["energy"][-1]),abs(outs["work"][0]),1e-300)
    exact_alpha=1j*a[None,:]**3*(np.conj(outs["c"])*ins["d"]-np.conj(outs["d"])*ins["c"])
    exact_beta=-1j*a[None,:]**3*(outs["c"]*ins["d"]-outs["d"]*ins["c"])
    p_coherence=-(mass(ph[-1])+2*k*k/a[-1]**2/3)*coh/(2*a[-1]**3*w)
    allw,alln,allcoh,_=occupation(bg,grid[:,None],k[None,:],ins["c"].T,ins["d"].T)
    gamma2=3*h[:,None]+(-h[:,None]*k[None,:]**2/a[:,None]**2+G*G*(ph[:,None]-PHISTAR)*v[:,None])/allw**2
    gas_ward=(allw*gamma2*allcoh/2/a[:,None]**3)@weights
    max_physical_edge=reference["k_scale"]*upper/float(a.min())
    result=dict(name=name,interpolation=bg.interpolation,nk=nk,kappa_min=.05,
                kappa_max=upper,k_band=[reference["k_scale"]*.05,
                                       reference["k_scale"]*upper],
                k=k,kappa=kap,weights=weights,rtol=rtol,start=start,end=end,
                initial_state=ir,out_state=orr,spectrum_n=n,
                endpoint_omega=w,coherence=coh,
                number_density=number,particle_rho=particle_rho,particle_p=particle_p,
                paired_rho_final=float(dr[-1]),paired_p_final=float(dp[-1]),
                paired_current_final=float(dj[-1]),
                paired_scalar_work=float(dscalar[-1]),paired_pressure_work=float(dpressure[-1]),
                exact_alpha_drift_abs=float(np.max(np.abs(exact_alpha-exact_alpha[:,-1,None]))),
                exact_beta_drift_abs=float(np.max(np.abs(exact_beta-exact_beta[:,-1,None]))),
                exact_beta_endpoint_occupation_abs=float(np.max(np.abs(np.abs(exact_beta[:,-1])**2-n))),
                peak_mass=float(np.sqrt(mass(ph)).max()),
                peak_physical_band_edge=max_physical_edge,
                peak_physical_node=float(k.max()/a.min()),
                ansatz_n=np.exp(-np.pi*kap**2),
                ansatz_occupation_max_abs_difference=float(np.max(np.abs(n-np.exp(-np.pi*kap**2)))),
                paired_energy_equality_relative=abs(dr[-1]-particle_rho)/max(abs(particle_rho),1e-300),
                paired_ward_endpoint_relative=float(abs(balance[-1])/normalizer),
                paired_ward_path_relative=float(rel.max()),
                paired_ward_absolute=float(abs(balance[-1])),
                paired_ward_residual_over_raw_scale=float(abs(balance[-1])/raw_scale),
                current_particle_only=float(np.dot(weights,current_particles)),
                current_coherence=float(np.dot(weights,current_coherence)),
                pressure_coherence=float(np.dot(weights,p_coherence)),
                pressure_pair_minus_particle=float(dp[-1]-particle_p),
                particle_gas_ward_defect_max_abs=float(np.max(np.abs(gas_ward))),
                current_pair_minus_particle=float(dj[-1]-np.dot(weights,current_particles)),
                chi2_reconstruction_abs=float(np.max(np.abs(f-frecon))),
                final_H_over_omega_max=float(np.max(np.abs(h[-1])/w)),
                final_adiabatic_max=float(np.max(np.abs(
                    (-h[-1]*k*k/a[-1]**2+G*G*(ph[-1]-PHISTAR)*v[-1])/
                    (w**3)))),
                instantaneous_subtraction=dict(
                    rho_final=float(instantaneous_delta[-1]),
                    current_final=float(ins["current"][-1]-zero_j[-1]),
                    pressure_final=float(ins["p"][-1]-zero_p[-1]),
                    triplet_algebraic_ward_defect_max=float(np.max(np.abs(zero_defect))),
                    mismatched_rhop_only_subtraction_ward_defect_max=float(np.max(np.abs(zero_j*v))),
                    note="Finite-band Hamiltonian zero-point removal is not renormalization. Its algebraic Ward defect vanishes for xi=0; no exact reference state evolves with that instantaneous triplet."))
    result["controls"]=dict(
        wronskian=ir["wronskian_max"]<=LIMITS["wronskian"] and orr["wronskian_max"]<=LIMITS["wronskian"],
        raw_ward=ir["ward_endpoint_relative"]<=LIMITS["endpoint_ward"] and orr["ward_endpoint_relative"]<=LIMITS["endpoint_ward"],
        paired_ward=result["paired_ward_endpoint_relative"]<=LIMITS["endpoint_ward"])
    if clocks:
        ids=np.array([0,nk//3,2*nk//3,nk-1])
        result["clock"]=clock_subset(bg,start,end,k[ids],state,rtol,
                           (ins["c"][ids,-1],ins["d"][ids,-1]))
        cc=result["clock"]
        result["controls"]["clock"]=(max(cc["max_occupation_abs_X"],cc["max_occupation_abs_u"])<=LIMITS["clock_occupation_abs"]
                                    and max(cc["wronskian_X"],cc["wronskian_u"])<=LIMITS["wronskian"])
        result["conformal_control"]=conformal_control(bg,start,end,k[ids],rtol)
        cm=result["conformal_control"]
        result["controls"]["conformal_mode"]=(cm["relative_mode_max"]<=LIMITS["conformal_mode_abs"]
                                              and cm["wronskian_max"]<=LIMITS["wronskian"]
                                              and cm["occupation_max"]<=1e-6)
    outfile=outdir/(name+"_diagnostics.npz")
    np.savez_compressed(outfile,t=grid,a=a,H=h,phi=ph,v=v,k=k,weights=weights,
                        rho_in=ins["rho"],p_in=ins["p"],J_in=ins["current"],
                        rho_out=outs["rho"],p_out=outs["p"],J_out=outs["current"],
                        delta_rho=dr,delta_p=dp,delta_J=dj,
                        raw_balance_in=ins["balance"],raw_balance_out=outs["balance"],
                        paired_balance=balance,paired_work=dw,
                        paired_scalar_work=dscalar,paired_pressure_work=dpressure,
                        scalar_work_in=ins["scalar_work"],pressure_work_in=ins["pressure_work"],
                        scalar_work_out=outs["scalar_work"],pressure_work_out=outs["pressure_work"],
                        exact_alpha=exact_alpha,exact_beta=exact_beta,
                        particle_gas_ward_defect=gas_ward,
                        zero_rho=zero_r,zero_p=zero_p,zero_J=zero_j,
                        instantaneous_delta_rho=instantaneous_delta,
                        spectrum_n=n,coherence=coh)
    result["diagnostics_path"]=str(outfile)
    result["diagnostics_sha256"]=sha(outfile)
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-root",type=Path,default=Path("data/A1"))
    ap.add_argument("--out",type=Path,default=Path("outputs/actual_modes"))
    ap.add_argument("--registration",type=Path,required=True)
    ap.add_argument("--registration-sha256",required=True)
    ap.add_argument("--source-ref",default=PINNED)
    ap.add_argument("--audit-only",action="store_true")
    args=ap.parse_args()
    if sha(args.registration)!=args.registration_sha256:
        raise ValueError("Registration SHA-256 does not match")
    args.out.mkdir(parents=True,exist_ok=True)
    fine,fmeta=load(args.data_root,FINE)
    coarse,cmeta=load(args.data_root,COARSE)
    bgh=Geometry(fine,"hermite")
    crossings=bgh.crossings()
    if len(crossings)!=1:
        raise ValueError("Registered one-event test requires exactly one crossing in [0,6.9]")
    cr=crossings[0]
    q=G*abs(cr["v"])
    ref=dict(crossing=cr,q=q,k_scale=cr["a"]*np.sqrt(q),
             definition="k=a_star sqrt(G |dphi/d(H0tau)|_star) kappa; a(archived t0)=1",
             G_definition="G=gbar/H0, verified in Sep30 B1 registration section1.2",
             matter_units="rho/H0^4,p/H0^4,J/H0^4,chi^2/H0^2",
             shell_conversion="kappa5^2 rho/H0=b rho_hat, b=kappa5^2 H0^3")
    report=dict(status="finite comoving-band prescribed-background diagnostic",
                source_ref=args.source_ref,script_sha256=sha(__file__),
                registration_sha256=args.registration_sha256,
                registration_commit="99dd818755c493127ae4b4e0cad7bbf21fd141de",
                input_baseline="8f67197b730d4e1c43554b86f224c29cc72629eb",
                audit_only=args.audit_only,
                limits=LIMITS,model=dict(G=G,phi_star=PHISTAR,m0=MASS0,xi=0,
                                        delta=.1,Y=0,decay=False),
                sources=[fmeta,cmeta],reference=ref,
                interpolation_limits="Hermite is C1 with curvature jumps at knots; cubic is C2. Neither finite interpolation supplies a UV Hadamard certificate.",
                source_limits="Raw finite band and exact in-minus-out state stresses are regulated diagnostics. No counterterms, thermalization, sourced junctions or bulk response.",
                cases=[],errors={})
    if not args.audit_only:
        cases=[
            ("primary",fine,"hermite",32,3.,0.,"wkb",1e-9,True),
            ("tight",fine,"hermite",32,3.,0.,"wkb",1e-11,False),
            ("quadrature64",fine,"hermite",64,3.,0.,"wkb",1e-9,False),
            ("band2p5",fine,"hermite",32,2.5,0.,"wkb",1e-9,False),
            ("band3p5",fine,"hermite",32,3.5,0.,"wkb",1e-9,False),
            ("coarse",coarse,"hermite",32,3.,0.,"wkb",1e-9,False),
            ("cubic",fine,"cubic",32,3.,0.,"wkb",1e-9,False),
            ("hamiltonian_in",fine,"hermite",32,3.,0.,"hamiltonian",1e-9,False),
            ("start0p2",fine,"hermite",32,3.,.2,"wkb",1e-9,False),
            ("out_end5",fine,"hermite",32,3.,0.,"wkb",1e-9,False,5.)]
        for spec in cases:
            name,data,kind,nk,upper,start,state,rtol,clocks=spec[:9]
            end=spec[9] if len(spec)>9 else END
            try:
                result=run_case(Geometry(data,kind),ref,nk,upper,start,state,
                                rtol,args.out,name,clocks,end)
                report["cases"].append(result)
                print("HDBLAST_CASE="+json.dumps(clean(result),separators=(",",":"),allow_nan=False),flush=True)
            except Exception as e:
                report["errors"][name]=type(e).__name__+": "+str(e)
                print("HDBLAST_CASE_ERROR="+name+":"+str(e),flush=True)
        primary=next((c for c in report["cases"] if c["name"]=="primary"),None)
        if primary is not None:
            for result in (c for c in report["cases"] if c["name"]!="primary"):
                result["change_from_primary"]={
                    key:(result[key]-primary[key])/max(abs(primary[key]),1e-300)
                    for key in ("number_density","particle_rho","paired_current_final",
                                "paired_rho_final")}
                result["change_from_primary"]["interpretation"]="Descriptive; band/state/background/interpolation changes are distinct controls."
    report["registered_controls_pass"]=(None if args.audit_only else
        not report["errors"] and len(report["cases"])==10 and
        all(all(c["controls"].values()) for c in report["cases"]))
    output=args.out/"actual_modes.json"
    output.write_text(json.dumps(clean(report),indent=2,allow_nan=False)+"\n")
    print("HDBLAST_ACTUAL_MODES_JSON="+json.dumps(clean(report),separators=(",",":"),allow_nan=False),flush=True)
    failed=bool(report["errors"]) or any(not all(c["controls"].values()) for c in report["cases"])
    return int(failed)

if __name__=="__main__":
    raise SystemExit(main())
