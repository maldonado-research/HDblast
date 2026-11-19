"""Particle flux and post-switch kinetic stress of a traveling tagged scalar's daughters."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
from scipy.signal import ZoomFFT


STUDY = Path("research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip")
WALL_SHA = "7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15"
REF = "research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs/{}/quiet/rolling_summary.json"
REF_SHA = {"fine": "abec0887ee4de00ba0757e7d96987610f36b7157a85bc4f00696fee8064a6810",
           "coarse": "91601211fce4a617838d3a709849a138ec29a23998cc7c3ba2416ec0b0b8aa8c"}
PARENTS = np.array([0.,4.,8.])
ARMS = ("quiet","plus","minus","far")
WINDOWS = {"wide":(1.15,2.05),"narrow":(1.25,1.95)}
CUT = 2.10


def pinned(path, digest):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError(f"Changed historical input {path}")
    return raw


def longitudinal(beta):
    """Energy-weighted zz stress divided by energy in the shell frame."""
    b = np.asarray(beta, dtype=float)
    if np.any((b < 0) | (b > 1)):
        raise ValueError("Unphysical boost")
    result = np.ones_like(b)
    small = b < .2
    x = b[small]**2
    series = np.full_like(x, 1/3)
    for n in range(1,14):
        series += 8*x**n/((2*n-1)*(2*n+1)*(2*n+3))
    result[small] = series
    mid = (~small) & (b < 1)
    y = b[mid]
    result[mid] = 2 - 1/y**2 + (1-y**2)**2*np.arctanh(y)/y**3
    return result


def angular_check():
    beta = np.array([0.,.001,.1,.2,.5,.9,.99,.999,1.])
    mu, w = np.polynomial.legendre.leggauss(256)
    b = beta[:,None]
    # Sum both particles after a longitudinal boost. Energies are normalized
    # by the total pair energy, with the sign convention of the +q branch.
    e1,e2 = (1+b*mu)/2, (1-b*mu)/2
    p1,p2 = (b+mu)/2, (b-mu)/2
    longitudinal_direct = ((p1**2/e1+p2**2/e2) @ w)/2
    transverse_direct = ((1-b*b)*(1-mu*mu)*(1/e1+1/e2) @ w)/16
    errors = [np.max(abs(longitudinal(beta)-longitudinal_direct)),
              np.max(abs(longitudinal(beta)+2*transverse_direct-1)),
              np.max(abs((p1+p2)-b))]
    if max(errors)>4e-12:
        raise ValueError(("Angular check", errors))
    return dict(boosts=beta.tolist(), longitudinal=longitudinal_direct.tolist(),
                max_error=float(max(errors)))


def propagate(archive, grid, arm, reference):
    def rows(tag):
        with archive.open(f"{grid}/y1_b1/{tag}.csv") as stream:
            return np.genfromtxt(stream, delimiter=",", names=True)
    data = np.concatenate([rows("pre"),rows(arm)])
    t = data["H0tau"]
    A = CubicHermiteSpline(t, data["A_b"]-reference, data["H_over_H0"])
    m = CubicHermiteSpline(t, 100*(data["phi_b"]-.5),100*data["v_over_H0"])
    start = float(t[0])
    a0 = np.exp(A(start))
    e0 = np.hypot(PARENTS/a0, m(start))
    v0 = 1/np.sqrt(2*a0*e0)
    # w = a chi' in conformal time (equal to Pi/a), v = a chi.
    w0 = -1j*a0*e0*v0
    n = len(PARENTS)

    def rhs(eta, state):
        v, w, s = state[:n],state[n:2*n],state[-1].real
        a = np.exp(A(s))
        cal_h = a*A(s,1)
        return np.r_[w+cal_h*v, -(PARENTS**2+(a*m(s))**2)*v-cal_h*w, a]

    def stop(eta, state):
        return state[-1].real-CUT
    stop.terminal = True
    stop.direction = 1
    run = solve_ivp(rhs, (0,2*(CUT-start)/a0),np.r_[v0,w0,start],events=stop,
                    method="DOP853",dense_output=True,rtol=2e-11,atol=1e-13,max_step=.004)
    if not run.success or len(run.t_events[0])!=1:
        raise RuntimeError("Could not reach the readout")

    def sample(s0,s1,count):
        limits = [brentq(lambda x:run.sol(x)[-1].real-s,0,run.t[-1]) for s in (s0,s1)]
        eta,deta = np.linspace(*limits,count,retstep=True,endpoint=False)
        state = run.sol(eta)
        s = state[-1].real
        v = state[:n]
        w = state[n:2*n]
        error = float(np.max(abs(v*np.conj(w)-np.conj(v)*w-1j)))
        if error > 2e-8:
            raise ValueError(("Parent Wronskian",error))
        f = np.sin(np.pi*(s-s0)/(s1-s0))**4
        return v*f,deta
    return sample,float(np.exp(A(CUT)))


def spectrum(q, omega, forward, backward, af):
    L = np.full_like(omega,1/3)
    if q:
        L = longitudinal(q/omega)
    outputs={}
    for band,mask in [("full",omega<=omega[-1]),("screened",omega<=60)]:
        x=omega[mask]
        weight=(x*x-q*q)/(8*np.pi**2)
        ann=weight*abs(forward[mask])**2
        stim=weight*abs(backward[mask])**2
        energy=x*(ann+stim)/af
        rho = float(simpson(energy,x=x))
        longitudinal_energy=float(simpson(energy*L[mask],x=x))
        transverse_energy=float(simpson(energy*(1-L[mask])/2,x=x))
        a=float(simpson(ann,x=x))
        b=float(simpson(stim,x=x))
        pz=q*(a-b)/af
        if abs(longitudinal_energy+2*transverse_energy-rho)>1e-11*rho or abs(pz)>rho:
            raise ValueError("Massless stress trace or causal flux failed")
        outputs[band]=dict(annihilation_number=a, stimulated_number=b, energy=rho,
                           longitudinal=longitudinal_energy, transverse_each=transverse_energy,
                           signed_flux=pz, longitudinal_fraction=longitudinal_energy/rho,
                           transverse_fraction=transverse_energy/rho, flux_fraction=pz/rho)
    return outputs


def evaluate(archive, grid, arm, reference, count, dw, ceiling):
    mode,af=propagate(archive,grid,arm,reference)
    record={}
    for window,(s0,s1) in WINDOWS.items():
        g,step=mode(s0,s1,count)
        result={}
        for i,q in enumerate(PARENTS):
            # Resolve the light-cone boundary and the upper energy screen exactly.
            bins=round((ceiling-q)/dw)+1
            omega=np.linspace(q,ceiling,bins)
            if np.min(abs(omega-60))>1e-12:
                raise ValueError("Requested screen is off grid")
            transform=ZoomFFT(count,[q,ceiling],bins,fs=2*np.pi/step,endpoint=True)
            # Sign and time-origin phases cancel in the modulus squared.
            plus=np.conj(transform(np.conj(g[i])))*step
            minus=transform(g[i])*step
            result[str(int(q))]=spectrum(q,omega,plus,minus,af)
        record[window]=result
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repository",type=Path)
    args=parser.parse_args()
    scales = {g:json.loads(pinned(args.repository/REF.format(g),sha))["matter"]["events"][0]["A_star"]
              for g,sha in REF_SHA.items()}
    results={}
    with zipfile.ZipFile(io.BytesIO(pinned(args.repository/STUDY,WALL_SHA))) as data:
        for group,grid,arms,count,dw,end in [
                ("fine","fine",ARMS,2048,.125,512.),("coarse","coarse",ARMS,2048,.125,512.),
                ("refined","fine",("quiet","plus"),4096,.0625,1024.)]:
            results[group]={arm:evaluate(data,grid,arm,scales[grid],count,dw,end) for arm in arms}
    print(json.dumps(dict(readout_s=CUT,windows=WINDOWS,orientation="one tagged traveling +z parent",
                     normalization="per added parent and per y^2; energies in H0", angular=angular_check(),
                     results=results),indent=2,allow_nan=False))


if __name__ == "__main__":
    main()
