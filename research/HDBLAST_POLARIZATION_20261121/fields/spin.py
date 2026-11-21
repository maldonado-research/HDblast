"""Unrestricted holomorphic angular patterns within the selected curl blocks."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import kv


def cohort_costs(h, waves, marker, source):
    """Integrals of the dual E and B norms bound the kick per unit modal drive."""
    end = marker['time']
    j = h['t'].searchsorted(end)
    t = np.r_[h['t'][:j], end]

    def nodes(f):
        return np.r_[f[:j], np.interp(end, h['t'], f)]

    a = nodes(h['a'])
    dt = np.diff(t)
    # Every birth within a time cell may conservatively receive the whole cell's force.
    wt = np.diff(nodes(h['D'])) / marker['D']
    active = wt > 0.
    travel = np.cumsum((dt / a[:-1])[::-1])[::-1][active]
    pstar = np.sqrt(np.max(h['R'][:j+1])**2 - 1.)
    modes = np.unique(waves['kappa'])
    modes = modes[modes >= 40]
    gram = np.array(source['emitted_energy_gram'])
    pulse_order = source['profile_order']
    coeff = []
    for k in modes:
        ix = np.flatnonzero(waves['kappa'] == k)
        ids = [pulse_order.index(name) for name in waves['profile'][ix]]
        L = np.linalg.cholesky(gram[np.ix_(ids, ids)])
        terms = []
        for name in ('E', 'B'):
            values = np.array([nodes(row) for row in waves[name][ix]])
            dual = np.linalg.norm(solve_triangular(L, values, lower=True), axis=0)
            terms.append(np.maximum(dual[:-1], dual[1:]))
        impulse = dt * a[1:] * (terms[0] + terms[1]) * (1.e7 / marker['a'])
        envelope = np.cumsum(impulse[::-1])[::-1][active]
        coeff.append(envelope**2)
    return modes, np.array(coeff).T, wt[active], travel, pstar


def full_spin(history, waves, marker, table):
    kappas, coeff, wt, flight, birth_kick = cohort_costs(history, waves, marker, table)
    distortions = []
    for d in kappas - 2:
        # Haar L² norm bound on the local maximum of a degree-d holomorphic polynomial.
        term = np.ones_like(flight)
        gain = term.copy()
        for j in range(1, d+1):
            term = term * np.sqrt((d-j+1.)/j) * np.sin(flight)
            gain += term
        distortions.append(gain)
    basis = coeff * np.array(distortions).T**2
    theta = marker['theta']
    mean_target = .25*np.pi * np.exp(-1./theta) * (2.*theta + 6.*theta**2 + 6.*theta**3) / kv(2.,1./theta)
    trial = np.geomspace(.02, 7., 125) * theta
    z = np.hypot(1.,trial)
    thermal = 1. - z*z*kv(2.,z/theta)/kv(2.,1./theta)
    results = {}
    for low_kappa in (40,72,160,258,360,515):
        selected = kappas >= low_kappa
        rms = np.sqrt(basis[:,selected].sum(axis=1))
        gain = float(wt @ rms)
        eps_min = (mean_target-birth_kick)/gain
        row = dict(kappas=kappas[selected].tolist(), number_real_coefficients=int(sum(2*(k-1)*8 for k in kappas[selected])),
                   mean_transverse_momentum_target=mean_target,maximum_gain_in_mean_transverse_momentum=gain,
                   required_amplitude_from_mean=eps_min,
                   maximum_local_L2_sphere_factor=float(np.max(distortions)), direct_cdf_gaps={})
        for eps in (.02,.1,.5,1.):
            cut = np.maximum(np.finfo(float).eps,trial-birth_kick)
            guaranteed = wt @ np.maximum(0.,1.-(eps*rms[:,None]/cut[None,:])**2)
            pos = np.argmax(guaranteed-thermal)
            row['direct_cdf_gaps'][str(eps)] = dict(cutoff=float(trial[pos]),guaranteed=float(guaranteed[pos]),thermal=float(thermal[pos]),
                                                     excess=float(max(0.,guaranteed[pos]-thermal[pos])))
        results[str(low_kappa)] = row
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('history','scan','output-dir'):
        p.add_argument('--'+name,required=True,type=Path)
    args=p.parse_args()
    runs=json.loads((args.history/'markers.json').read_text())['runs']
    env=json.loads((args.scan/'STORAGE_SCAN.json').read_text())
    out={}
    for label,table in env.items():
        key='dedicated_native' if label=='dedicated_spatial' else label
        with np.load(args.history/(key+'.npz')) as h, np.load(args.scan/'traces'/(label+'.npz')) as w:
            out[label] = full_spin(h,w,runs[key]['markers']['0.01'],table)
        print(label,[(k, round(v['required_amplitude_from_mean'],4), [round(a['excess'],4) for a in v['direct_cdf_gaps'].values()]) for k,v in out[label].items()],flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/'SPIN.json').write_text(json.dumps(out,sort_keys=True,indent=2,allow_nan=False)+'\n')
if __name__=='__main__': main()
