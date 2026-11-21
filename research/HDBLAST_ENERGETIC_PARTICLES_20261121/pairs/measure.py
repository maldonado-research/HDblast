"""Sum instantaneous directed collision participations on every eligible charge pair."""
import argparse
import json
from pathlib import Path

import numpy as np


def measurement(path, block=32):
    with np.load(path) as data:
        momentum, charge = data["p"], data["sign"]
    radius = np.linalg.norm(momentum, axis=1)
    energy = np.hypot(1., radius)
    unit = np.divide(momentum, radius[:, None], out=np.zeros_like(momentum), where=radius[:, None] > 0)
    # Repeated tag energies are possible, so choose exactly ceil(N/10) in stable input order.
    n_hot = (len(energy)+9)//10
    hot = np.zeros(len(energy), dtype=bool)
    hot[np.argsort(-energy, kind="stable")[:n_hot]] = True
    pos, neg = np.flatnonzero(charge == 1.), np.flatnonzero(charge == -1.)
    if len(pos) != len(neg) or 2*len(pos) != len(energy):
        raise ValueError("Expected equal-weight opposite sign cohorts of equal size")
    rapidity = np.arcsinh(radius)
    totals = np.zeros((2, 2))  # indices denote kernel and total / high-energy participants
    for start in range(0, len(pos), block):
        ix = pos[start:start+block, None]
        # u1.u2 - 1 = cosh(eta1-eta2)-1 + |p1||p2| (1-n1.n2).
        radial = 2*np.sinh(.5*(rapidity[ix]-rapidity[neg]))**2
        chord2 = sum((unit[ix, k]-unit[neg, k])**2 for k in range(3))
        delta = radial + .5*radius[ix]*radius[neg]*chord2
        pair_energy = energy[ix]*energy[neg]
        s = 4+2*delta
        beta = np.sqrt(2*delta/s)
        pole = 6*np.pi*beta/(pair_energy*(((s-8300.**2)/(.30*(s+2)))**2+beta**2))
        flux = np.sqrt(delta*(delta+2))/pair_energy
        for row, kernel in enumerate((pole, flux)):
            totals[row, 0] += 2*kernel.sum()
            totals[row, 1] += kernel[hot[ix[:,0]]].sum() + kernel[:,hot[neg]].sum()
    fraction = n_hot/len(energy)
    share = totals[:,1]/totals[:,0]
    return dict(snapshot=path.stem, unordered_eligible_pairs=int(len(pos)*len(neg)),
                particles=len(energy), top_particles=n_hot, top_energy_share=energy[hot].sum()/energy.sum(),
                min_top_gamma=energy[hot].min(), pole_top_participation_share=share[0],
                constant_sigma_top_participation_share=share[1], pole_relative_rate=share[0]/fraction,
                constant_sigma_relative_rate=share[1]/fraction, mean_pair_kernel=totals[0,0]/(2*len(pos)*len(neg)),
                mean_pair_flux=totals[1,0]/(2*len(pos)*len(neg)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshots", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps([measurement(p) for p in args.snapshots], indent=2, allow_nan=False)+"\n")
