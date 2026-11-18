#!/usr/bin/env python3
"""Replay the declared finite-collar stress balance using only derived CSVs."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
HEADER=('T','s','collar_energy','inner_flux','wall_flux','metric_work','junction_wall_flux')


def require(ok,message):
    if not ok:raise RuntimeError(message)


def trap(t,y):return math.fsum((t[i+1]-t[i])*(y[i]+y[i+1])/2 for i in range(len(t)-1))


def evaluate(path):
    with path.open(newline='') as f:
        r=csv.DictReader(f);require(tuple(r.fieldnames)==HEADER,'wrong columns')
        rows=[{k:float(v) for k,v in row.items()} for row in r]
    require(len(rows)>=3 and all(math.isfinite(v) for row in rows for v in row.values()),'bad values')
    t=[r['T'] for r in rows];require(all(a<b for a,b in zip(t,t[1:])),'bad times')
    require(rows[0]['s']==1 and rows[-1]['s']==4,'wrong selected window')
    def integ(key,indices=None):
        if indices is None:return trap(t,[r[key] for r in rows])
        return trap([t[i] for i in indices],[rows[i][key] for i in indices])
    dE=rows[-1]['collar_energy']-rows[0]['collar_energy']
    left=integ('inner_flux');wall=integ('wall_flux');work=integ('metric_work')
    prediction=left-wall+work;residual=dE-prediction
    total=abs(dE)+abs(left)+abs(wall)+abs(work)
    scale=max(abs(dE),abs(prediction))
    thin=list(range(0,len(rows),2))
    if thin[-1]!=len(rows)-1:thin.append(len(rows)-1)
    pred_thin=integ('inner_flux',thin)-integ('wall_flux',thin)+integ('metric_work',thin)
    # Deliberately report discrepancies; no outcome-dependent pass threshold.
    return dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),rows=len(rows),
        T_start=t[0],T_end=t[-1],energy_start=rows[0]['collar_energy'],energy_end=rows[-1]['collar_energy'],
        energy_change=dE,inner_influx_integral=left,outward_wall_flux_integral=wall,metric_work_integral=work,
        predicted_energy_change=prediction,residual=residual,
        residual_over_absolute_term_sum=residual/total if total else 0.,
        residual_over_energy_change_scale=residual/scale if scale else 0.,
        temporal_thinning_prediction=pred_thin,temporal_thinning_residual=dE-pred_thin,
        temporal_thinning_relative_prediction_change=(pred_thin-prediction)/scale if scale else 0.,
        wall_flux_from_junction_integral=integ('junction_wall_flux'),
        two_copy_measured_wall_flux=2*wall,
        wrong_drop_metric_work_residual=dE-(left-wall),
        wrong_reverse_wall_sign_residual=dE-(left+wall+work))


def main():
    p=argparse.ArgumentParser();p.add_argument('--grid',choices=['coarse','fine','both'],default='both');p.add_argument('--output',type=Path);a=p.parse_args()
    grids=('coarse','fine') if a.grid=='both' else (a.grid,)
    report=dict(scope='Numerical finite-collar scalar-stress diagnostic, unit-observer energy, per copy. No global positive-energy or cosmological claim.',
                formula='Delta E_coll = integral [F_left-F_wall+G_coll] dT',
                rules='Fixed Z[-.04,0], s[1,4]; derivative degree2/4; strict-interior wall stencil; spatial and temporal trapezoids.',
                results=[evaluate(HERE/'collar_data'/f'{g}_d{d}.csv') for g in grids for d in (2,4)])
    if a.output: a.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__=='__main__':main()
