"""Singular values of the fixed-wall covariance-to-force map in specified L2 norms."""
import argparse
import importlib.util
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("fork_dynamics",HERE.parent/"initial/corner.py")
fork=importlib.util.module_from_spec(spec)
spec.loader.exec_module(fork)


def analyze(background, count, samples):
    loga, phi, start, cut = background
    q, w = fork.quadrature(count)
    _, _, c = fork.reconstruct((loga,phi),start,cut,q)
    # Columns evolve the derivative of the final field variance w.r.t. the
    # off-diagonal input covariance. The diagonal input entries stay fixed.
    times=np.linspace(cut,2.10,samples)
    data=solve_ivp(fork.covariance_rhs,(cut,times[-1]),np.r_[0*q,0*q,1+0*q],
                   args=(loga,phi,q),method="DOP853",t_eval=times,rtol=2e-11,atol=2e-13,max_step=.001)
    if not data.success:
        raise RuntimeError(data.message)
    weights_time=np.ones(samples)/(samples-1)
    weights_time[[0,-1]]*=.5
    # L2 input uses probability measure w/sum(w). Physical j uses w itself.
    # Output uses the mean squared force over the proper-time interval.
    kernel=(1e4*(phi(times)-.5))[:,None]*data.y[:count].T
    op=kernel*(abs(c)*np.sqrt(w*w.sum()))[None,:]*np.sqrt(weights_time)[:,None]
    sigma=np.linalg.svd(op,compute_uv=False)
    ratios=sigma/sigma[0]
    return dict(momentum_nodes=count,proper_time_samples=samples,cut=float(cut),
                leading_gain=float(sigma[0]),first_twelve_ratios=ratios[:12].tolist(),
                above_one_percent=int(np.sum(ratios>.01)),above_one_per_mille=int(np.sum(ratios>.001)),
                above_one_part_per_million=int(np.sum(ratios>1e-6)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repository",type=Path)
    args=parser.parse_args()
    report={}
    native=json.loads((HERE.parent/"native/FORK.json").read_bytes())
    checked=fork.checked(args.repository/fork.WALLS,fork.WALL_SHA)
    with ZipFile(io.BytesIO(checked)) as z:
        for grid,row in native.items():
            bg=fork.background(z,grid,row["reference_A"])
            report[grid]=[analyze(bg,n,t) for n,t in ((24,257),(48,513),(96,1025))]
    print(json.dumps(report,allow_nan=False,indent=2))


if __name__=="__main__":
    main()
