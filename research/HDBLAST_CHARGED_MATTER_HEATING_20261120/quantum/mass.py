"""Finite-band charged and neutral massive conformal scalar controls at A=0."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location("closed_wall_pairs", Path(__file__).with_name("pairs.py"))
pairs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pairs)


def evolve(hist, masses, charges, nmax, density):
    j, k = pairs.labels(nmax)
    mu = np.asarray(masses)[:, None, None]
    charge = np.asarray(charges)[None, :, None]
    def frequency(radius, amplitude):
        shift = charge * amplitude / math.sqrt(3)
        return np.sqrt(j*j + (mu*radius)**2 + shift**2 - 2*k*shift)
    wi = frequency(hist.radius(0.), 0.)
    wf = frequency(hist.radius(hist.cross), 0.)
    qmode = 1/np.sqrt(2*wi) + 0j
    pmode = -1j*np.sqrt(wi/2)
    mesh = np.r_[np.linspace(0., hist.early, 1+math.ceil(hist.early*density/1.5e-5))[:-1],
                 np.linspace(hist.early, hist.cross, 1+math.ceil((hist.cross-hist.early)*density/1.2e-3))]
    mid = (mesh[1:]+mesh[:-1])/2
    for radius, amplitude, step in zip(hist.radius(mid), hist.potential(mid), np.diff(mesh)):
        w = frequency(radius, amplitude)
        cs, sn = np.cos(step*w), np.sin(step*w)
        qmode, pmode = cs*qmode + (sn/w)*pmode, cs*pmode - (w*sn)*qmode
    alpha = (np.sqrt(wf)*qmode + 1j*pmode/np.sqrt(wf))/math.sqrt(2)
    beta = (np.sqrt(wf)*qmode - 1j*pmode/np.sqrt(wf))/math.sqrt(2)
    occ = abs(beta)**2
    err = float(np.max(abs(abs(alpha)**2-occ-1)))
    return j,k,wf,occ,err


def calculate(trace):
    masses = [0., 1/(math.sqrt(3)*1e8), .1, 1.]
    hist = pairs.History(trace)
    j,k,w,occ,err = evolve(hist,masses,[0.,1.],32,2)
    _,_,_,low,_ = evolve(hist,masses,[0.,1.],32,1)
    if err > 1e-10 or max(occ[0,0])>1e-20:
        raise ArithmeticError("Vacuum and symplectic controls failed")
    rows = []
    for i,m in enumerate(masses):
        neutral,charged=occ[i]
        n2=j==3
        n2data=charged[n2]
        if np.ptp(neutral[n2])>1e-16:
            raise ArithmeticError("Mass-only degeneracy broke")
        growing=2*float(np.sum(j*w[i,1]*charged))
        baseline=2*float(np.sum(j*w[i,0]*neutral))
        rows.append(dict(mass=m, equal_shell_energy=float(w[i,1,n2][0]),
                         neutral_shell2_occupation=float(neutral[n2][0]),
                         charged_shell2_occupations=n2data.tolist(),
                         charged_shell2_max_to_min=float(np.max(n2data)/np.min(n2data)),
                         charged_shell2_max_minus_min=float(np.ptp(n2data)),
                         neutral_pair_count=float(np.sum(j*neutral)),
                         charged_pair_count=float(np.sum(j*charged)),
                         neutral_conformal_excitation_energy=baseline,
                         charged_conformal_excitation_energy=growing,
                         charged_minus_neutral_excitation_energy=growing-baseline,
                         charged_shell2_timestep_relative_discrepancy=float(np.max(abs(low[i,1,n2]/n2data-1)))))
    tiny=float(np.max(abs(occ[1]-occ[0])))
    if tiny>5e-13 or max(x['charged_shell2_timestep_relative_discrepancy'] for x in rows)>1e-4:
        raise ArithmeticError("Near-zero-mass or time resolution control failed")
    return dict(nmax=32, charge=1., event_t=float(hist.time(hist.cross)), event_eta=hist.cross,
                max_wronskian_error=err, max_massive_to_tiny_massless_difference=tiny, cases=rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.trace),indent=2,sort_keys=True,allow_nan=False)+"\n")
