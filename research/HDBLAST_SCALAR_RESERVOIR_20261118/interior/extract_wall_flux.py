"""Reduce the read-only B1 collar outputs to two derived wall-flux CSVs.

Pass a run_instrumented.py directory, 'coarse' or 'fine', an output CSV and,
optionally, the archived original wall timeseries.npz for trajectory verification.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ARCHIVE_PINS = {
    "coarse": "3c9e732322715bc33d09eff6212066b15448da70c7ea10a177edc39e9f6ac9b2",
    "fine": "9eae8828655af2ec8cc91bb658ecf0e768ca9f54215b8927313937064f1b1646",
}
FIELDS = ["s","tension_rate","matter_rate","junction_rate"]
LABELS = [(incl,p) for incl in ("interior","with_wall") for p in (1,2,3,4)]
FIELDS += [f"{incl}_{p}" for incl,p in LABELS]


def derivative(z, degree, include_wall):
    # Fit a degree-p polynomial to p+1 physical grid coordinates. Subtracting
    # the first phi value before the dot product cancels the constant exactly.
    last = len(z) if include_wall else len(z)-1
    loc = slice(last-degree-1,last)
    dz = z[-1]-z[-2]
    nodes = (z[loc]-z[-1])/dz
    weights = np.linalg.solve(np.vander(nodes,degree+1,increasing=True).T,
                              np.eye(degree+1)[1])/dz
    for power in range(1,degree+1):
        got = np.dot(weights*dz, nodes**power-nodes[0]**power)
        if abs(got-float(power==1))>1e-10:
            raise ValueError("manufactured derivative")
    return loc,weights


def check_archive(path, grid, run_path, clock_end):
    if hashlib.sha256(path.read_bytes()).hexdigest()!=ARCHIVE_PINS[grid]:
        raise ValueError("historical trajectory hash")
    def load(p):
        with np.load(p,allow_pickle=False) as z:
            return dict(zip(z["cols"].tolist(),z["rec"].T))
    orig, run=load(path),load(run_path)
    stop=int(np.searchsorted(run["H0tau"],clock_end,"left"))+1
    dkeys=("T","H0tau","phi_b","A_b","B_b","H_over_H0","v_over_H0","R","X","J","sigma1")
    # Compare the same original coordinate-time output indices. Separate clock
    # differences are reported, so phase and amplitude errors stay visible.
    if np.max(np.abs(orig["T"][:stop]-run["T"][:stop]))>1e-10:
        raise ValueError("coordinate cadence changed")
    return {k:float(np.max(np.abs(orig[k][:stop]-run[k][:stop]))) for k in dkeys}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory",type=Path)
    parser.add_argument("grid",choices=["coarse","fine"])
    parser.add_argument("output",type=Path)
    parser.add_argument("--archive",type=Path)
    p=parser.parse_args()
    folder=p.run_directory/"run"
    summary=json.loads((folder/f"b1_{p.grid}_summary.json").read_text())
    rb=summary["rho_b"]
    files=sorted((folder/"collar").glob("sample_*.npz"))
    with np.load(files[0]) as r0:
        z=r0["z"].copy(); initial_A=r0["A"][-1]
    if z[-1] != 0 or not np.all(np.diff(z)>0):
        raise ValueError("expected interior at negative z and fixed boundary at zero")
    stencils = [derivative(z,order,label=="with_wall") for label,order in LABELS]
    rows=[]
    for f in files:
        with np.load(f,allow_pickle=False) as u:
            if not np.array_equal(z,u["z"]):
                raise ValueError("moving physical mesh")
            s=float(u["s"]); lapse_factor=math.exp(-u["B"][-1])
            scale=math.exp(3*(u["A"][-1]-initial_A))
            velocity=rb*lapse_factor*u["phi_T"][-1]
            # Rates in d/ds of a^3 kappa_5^2 rho/H0; both copies, inflow to wall.
            factor=rb*scale*velocity
            tension_rate=factor*float(u["sigma1"])
            matter_rate=factor*float(u["J"])
            row=[s,tension_rate,matter_rate,tension_rate+matter_rate]
            phi=u["phi"]
            for loc,w in stencils:
                phiz=float(np.dot(w,phi[loc]-phi[loc][0]))
                row.append(-2*factor*lapse_factor*phiz)
            rows.append(row)
    raw=np.array(rows).T
    s=raw[0]
    if not (s[0]==0 and s[-1]>4 and np.all(np.diff(s)>0)):
        raise ValueError("proper-time coverage")
    samples=np.r_[1,s[(s>1)&(s<4)],4]
    values=[samples]+[np.interp(samples,s,rate) for rate in raw[1:]]
    with p.output.open("w",newline="") as out:
        writer=csv.writer(out,lineterminator="\n")
        writer.writerow(FIELDS)
        for row in zip(*values):
            writer.writerow(format(float(x),".17g") for x in row)
    receipt={"grid":p.grid,"rho_b":rb,"start_proper":float(s[0]),"end_proper":float(s[-1]),
             "records":len(files),"samples_exported":len(samples),"z_min":float(z[0]),
             "z_spacing_at_wall":float(z[-1]-z[-2]),"sigma1_and_J_units":"unscaled model densities in junction",
             "csv_sha256":hashlib.sha256(p.output.read_bytes()).hexdigest(),
             "last_step_original_stop":summary["stop_reason"]}
    if p.archive:
        receipt["original_npz_sha256"]=ARCHIVE_PINS[p.grid]
        receipt["max_abs_wall_trajectory_discrepancy_at_original_samples_through_s4_bracket"]=check_archive(
            p.archive,p.grid,folder/f"b1_{p.grid}_timeseries.npz",4)
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
