"""Recover compensated finite-amplitude initial data from a zero initial control."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import root

sys.path.insert(0, str(Path(__file__).resolve().parent))
from matched_pulse import PROFILES, endpoint, integrate, wall_values


def shoot(summary, profile_name, amplitude, step=0.001):
    wall = wall_values(summary)
    profile = PROFILES[profile_name]

    def residual(c):
        return endpoint(integrate(wall,profile['cut'],amplitude,c,profile,step,dense=False))

    zero = np.zeros(3)
    eps = 0.001
    jacobian = np.column_stack([(residual(x)-residual(-x))/(2*eps) for x in eps*np.eye(3)])
    initial = np.linalg.solve(jacobian,-residual(zero))
    solution = root(residual,initial,tol=1e-8)
    error = residual(solution.x)
    if not np.isfinite(error).all() or max(abs(error))>2e-10:
        raise ValueError(f'Root finder failed to close the endpoint residual: {error}')
    return dict(profile=profile_name,amplitude=amplitude,coefficients=list(solution.x),
                raw_endpoint=list(error),optimizer_message=solution.message,
                condition_at_zero=float(np.linalg.cond(jacobian)),max_step=step)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary',type=Path)
    parser.add_argument('profile',choices=PROFILES)
    parser.add_argument('amplitudes',type=float,nargs='+')
    parser.add_argument('--step',type=float,default=.001)
    args=parser.parse_args()
    print(json.dumps([shoot(args.summary,args.profile,a,args.step) for a in args.amplitudes],
                     indent=2,allow_nan=False))
