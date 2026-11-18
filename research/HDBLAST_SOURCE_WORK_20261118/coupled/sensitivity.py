"""Preparation time and quadrature sensitivity in the native, decaying, reciprocal model."""
from pathlib import Path
import argparse
import json

import numpy as np


def calculate(root):
    results={}
    groups=[('coarse',24,.65),('coarse',24,.70),('coarse',24,.75),('fine',24,.70),('fine',48,.70)]
    for grid,nodes,prep in groups:
        name=f'{grid}/k{nodes}_y1_b1_t{prep:.3f}'
        cases={}
        for arm in ('quiet','plus','minus','far'):
            data=root/name/arm/'wall_timeseries.npz'
            with np.load(data,allow_pickle=False) as f:
                c=dict(zip(f['cols'],f['rec'].T))
            if not c['H0tau'][0]<2.10<c['H0tau'][-1]:raise ValueError(data)
            interpolate=lambda x:float(np.interp(2.10,c['H0tau'],x))
            cases[arm]=dict(rad_hat=interpolate(c['rad_hat']),rho_hat=interpolate(c['rho_hat']),
                rad_fraction=interpolate(c['rad']/c['H_over_H0']**2),scalar_work=interpolate(c['scalar_work']),
                pressure_work=interpolate(c['pressure_work']),phi_b=interpolate(c['phi_b']))
        for record in cases.values():
            record['radiation_contrast_percent']=100*(record['rad_hat']/cases['quiet']['rad_hat']-1)
        results[name]=dict(preparation=json.loads((root/name/'prepare.json').read_text()),cases=cases)
    return results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('native_workdir',type=Path)
    args=parser.parse_args()
    print(json.dumps(calculate(args.native_workdir),indent=2,allow_nan=False))
