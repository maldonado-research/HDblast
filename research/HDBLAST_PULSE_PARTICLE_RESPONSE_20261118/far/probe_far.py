"""Try remote compensators on the same preproduction native B1 fork."""
import argparse
from pathlib import Path
import json
import os
import sys
import numpy as np


def trial(root, grid, centers, kind, amps, profile=False):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    sys.path.insert(0, str(root/'runs'/grid/'solver'))
    os.chdir(root/'runs'/grid/'solver')
    from rolling_data import add_pulse
    from evolve_b1 import Run
    from shell_matter import ChiGas
    fine = grid == 'fine'
    run = Run(delta=.1, c=.5975949350280132, d=-3.106933495673783, dc=.01,
              dzf=.0005 if fine else .001, dzc=.002 if fine else .004,
              L=16., kappa=10., project=None, cfl=.5, order=4, xc=1.9,
              matter=lambda rb: ChiGas(100.,.5,1.,6.49e-6,rb), log=lambda msg: None)
    cp = np.load(root/'runs'/grid/'checkpoint.npz')
    run.mode = str(cp['mode'].item())
    out = []
    for center, spread in centers:
        for amp in amps:
            try:
                _, row = add_pulse(run, float(cp['T']), cp['Yv'], cp['m'], amp,
                                   dict(control=center, spread=spread,kind=kind,probe=profile))
                row['status']='profile' if profile else 'compact'
            except ValueError as err:
                row = dict(amplitude=amp,control_center=center,control_halfwidth=spread,
                           status='failed', reason=str(err))
            out.append(row)
    return out


if __name__ == '__main__':
    par = argparse.ArgumentParser()
    par.add_argument('root', type=Path)
    par.add_argument('--grid', default='coarse')
    par.add_argument('--kind', default='incident')
    par.add_argument('--profile',action='store_true')
    par.add_argument('--amps', default='0.001,-0.001')
    par.add_argument('--centers', default='-.85,.25;-.90,.30')
    args=par.parse_args()
    centers=[tuple(map(float,pair.split(','))) for pair in args.centers.split(';')]
    print(json.dumps(trial(args.root.resolve(),args.grid,centers,args.kind,list(map(float,args.amps.split(','))),args.profile),indent=2,allow_nan=False))
