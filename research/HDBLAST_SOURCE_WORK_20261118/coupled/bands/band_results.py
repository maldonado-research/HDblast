"""Tabulate daughter radiation and the physical-energy screen for the wider bands."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

ARMS=('quiet','plus','minus','far')
LIMIT=6.49e-6**(-1/3)


def walls(path,pre):
    with np.load(path,allow_pickle=False) as f:
        values=dict(zip(f['cols'],f['rec'].T))
    return {k:np.r_[np.asarray(pre[k],dtype=float),a] for k,a in values.items() if k in pre and k!='mode'}


def summarize(root,base):
    report={}
    for edge in (12,14):
        for grid in ('coarse','fine'):
            folder=root/f'K{edge}'/grid/'k32_y1_b1_t0.700'
            with np.load(str(folder/'fork.npz')+'.history.npz',allow_pickle=False) as f:
                pre=dict(zip(f['cols'],f['rec'].T))
            Astar=json.loads((base/'classical/runs'/grid/'quiet/rolling_summary.json').read_text())['matter']['events'][0]['A_star']
            group={}
            for arm in ARMS:
                d=walls(folder/arm/'wall_timeseries.npz',pre)
                t=d['H0tau']
                cut=2.10
                at=lambda y:float(np.interp(cut,t,y))
                logs=d['A_b']-Astar
                mass=100*(d['phi_b']-.5)
                w=np.hypot(edge*np.exp(-logs),mass)
                if not np.allclose(w,d['max_physical_frequency'],rtol=1e-13):
                    raise ValueError('Physical frequency diagnostic mismatch')
                # On linearly interpolated (log a, phi), E² is convex on each segment.
                mask=t<cut
                maxw=max(float(max(w[mask])),math.hypot(edge*math.exp(-at(logs)),at(mass)))
                if maxw>=LIMIT:raise ValueError('Admitted frequency ceiling exceeded')
                fraction=d['rad']/d['H_over_H0']**2
                group[arm]=dict(rad_hat=at(d['rad_hat']),massive_hat=at(d['rho_hat']),rad_fraction=at(fraction),
                    max_rad_fraction_through_cut=max(max(fraction[mask]),at(fraction)),frequency_envelope_bound=maxw,
                    scalar_work=at(d['scalar_work']),pressure_work=at(d['pressure_work']),
                    ward_error=at(d['comoving_energy']-d['scalar_work']+d['pressure_work']))
            for arm,data in group.items():
                data['contrast_percent']=100*(data['rad_hat']/group['quiet']['rad_hat']-1)
            report[f'{grid}/K{edge}']=group
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('previous_capsule',type=Path)
    o=parser.parse_args()
    print(json.dumps(summarize(o.work,o.previous_capsule),indent=2,allow_nan=False))
