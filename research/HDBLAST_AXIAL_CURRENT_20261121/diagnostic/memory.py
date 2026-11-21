"""Attribute the axial field difference to earlier current differences with the exact discrete adjoint."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(reference, particle, background):
    with np.load(reference) as first, np.load(particle) as second:
        t = first["t"]
        if not np.array_equal(t, second["t"]) or not np.array_equal(first["incoming"], second["incoming"]):
            raise ValueError("Current histories must have the same clock and incoming wave")
        current = second["J_used"][:,2] - first["J_used"][:,2]
        field = second["E"][:,2] - first["E"][:,2]
        a = np.interp(t, background["t"], background["a"])
        if abs(field[0]) > 1.e-15:
            raise ValueError("Initial fields differ")
    current = current[1:]
    dt = np.diff(t)
    amid = (a[:-1]+a[1:])/2
    deta = dt/amid
    # Each free half-step is a damped curl rotation bracketed by exact impedance steps.
    damp = np.exp(-dt/(2*2.54))
    M = np.empty((len(dt),2,2))
    M[:,0,0] = damp**2*np.cos(deta)
    M[:,0,1] = -damp*np.sin(deta)
    M[:,1,0] = -M[:,0,1]
    M[:,1,1] = np.cos(deta)
    time_splits = (6.37, 6.45, 6.47163486, 6.57481192)
    sample_times = (6.45,6.47163486,6.57481192,t[-1])
    mid = (t[:-1]+t[1:])/2
    observations=[]
    for end in sample_times:
        n=int(np.searchsorted(t,end,side="left"))
        adj=np.array([1.,0.])
        kernels=np.zeros(n)
        for j in range(n-1,-1,-1):
            half=adj @ M[j]
            kernels[j]=half[0]
            adj=half @ M[j]
        terms=-kernels*current[:n]*deta[:n]/(2.54*a[n]**2)
        cuts=np.searchsorted(time_splits,mid[:n])
        groups=np.bincount(cuts,weights=terms,minlength=len(time_splits)+1)
        difference=abs(terms.sum()-field[n])
        if np.min(kernels)<=0. or difference>1.e-10:
            raise ValueError(("Field-response identity failed",difference))
        observations.append(dict(time=float(t[n]),reference_difference=float(field[n]),sum=float(terms.sum()),
                            positive=float(np.maximum(terms,0).sum()),negative=float(np.minimum(terms,0).sum()),
                            contributions_by_time_bin=groups.tolist(),check_residual=float(difference)))
    return observations


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs",type=Path,required=True)
    p.add_argument("--snapshots",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    with np.load(args.inputs/"dedicated_half.npz") as background:
        results={}
        for file in sorted(args.snapshots.glob("t6.67000-*-particle-pole.npz")):
            ref=file.with_name(file.name.replace("-particle-", "-fluid-"))
            results[file.stem]=compare(ref,file,background)
    if len(results)!=4: raise ValueError("Expected four paired numerical refinements")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(time_bin_edges=[6.37,6.45,6.47163486,6.57481192],runs=results),indent=2,sort_keys=True)+"\n")
