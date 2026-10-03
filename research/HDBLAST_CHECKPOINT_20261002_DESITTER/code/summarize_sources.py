#!/usr/bin/env python3
"""Postprocess frozen producer and independent-review artifacts; no source solve."""
import argparse
from decimal import Decimal as D, getcontext
import hashlib
import json
from pathlib import Path

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def sci(value):return f'{D(value):.3E}'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--primary',type=Path,required=True)
    ap.add_argument('--independent',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();getcontext().prec=70
    p=json.loads(a.primary.read_text());r=json.loads(a.independent.read_text())
    if p['status']!='PASS' or r['status']!='PASS':raise ValueError('Cannot report PASS from failing inputs')
    rows=[];comparison=[]
    for v in p['points']:
        selected=[t for t in r['results'] if D(t['x'])==D(v['x_over_r']) and D(t['z'])==D(v['H2_over_r']) and t['settings']['dps']==55]
        if len(selected)!=1:raise ValueError('Missing or duplicate independently refined point')
        alt=selected[0]
        for q in ['W','rho','Q']:
            difference=abs(D(alt[q])-D(v['refined'][q]));threshold=D('1e-9')+D('1e-7')*abs(D(v['refined'][q]))
            comparison.append({'x_over_r':v['x_over_r'],'H2_over_r':v['H2_over_r'],'observable':q,'difference':str(difference),'threshold':str(threshold),'pass':difference<=threshold})
        rows.append('| '+v['x_over_r']+' | '+v['H2_over_r']+' | '+' | '.join(f'{D(v["refined"][q]):.12E}' for q in ['W','rho','Q'])+' |')
    if len(comparison)!=12 or not all(v['pass'] for v in comparison):raise ValueError('Failed cross-method comparison')
    errorrows=[]
    for q in ['W','rho','Q']:
        low=max(D(v['primary']['tail_bound'][q]) for v in p['points'])
        high=max(D(v['refined']['tail_bound'][q]) for v in p['points'])
        change=max(abs(D(v['primary'][q])-D(v['refined'][q])) for v in p['points'])
        cross=max(D(v['difference']) for v in comparison if v['observable']==q)
        errorrows.append('| '+q+' | '+' | '.join(map(sci,[low,high,change,cross]))+' |')
    summary={'status':'PASS','primary_sha256':digest(a.primary),'independent_sha256':digest(a.independent),'postprocessor_sha256':digest(__file__),'comparisons':comparison,'total_error_enclosure':False}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/'CROSS_METHOD_COMPARISON.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    md='''# Actual massive de Sitter quantum sources: four-point mathematical benchmark

The registered determinant calculation passed all 57 primary gates and detected
three deliberately wrong source prescriptions. The separately implemented
proper-time calculation agrees at all four points and passes all 12 registered
cross-method comparisons. These are one-loop Euclidean/Bunch–Davies vacuum
sources in the explicitly declared finite convention, with no physical HDBLAST
mass or gravitational coupling selected.

| x/r | H²/r | W/r² | rho/r² | Q/r |
|---:|---:|---:|---:|---:|
'''+ '\n'.join(rows)+'''

The pressure is p=−rho by de Sitter invariance. The scalar-shell source remains
J_phi=x_phi Q/2; no numerical x_phi is invented. The negative rho and Q at
(x/r,H²/r)=(2,1/2) are retained. They are finite-scheme vacuum expectation
values, not a negative occupation number or a particle-heating result.

The calculation derives both sources from the same finite action: Q=2W_x and
rho=W−(H²/2)W_(H²), with r fixed. The independent calculation instead integrates
the heat trace and its separately differentiated mass and metric kernels.
It does not import the zeta producer. The trace relation is checked only after
the metric density has been computed; it is an algebraic consistency check of
the shared action and is not independent evidence for determinant truncation.

## Numerical evidence and its limits

The table gives the largest absolute quantities across the four points, in
the same reference units as the source table. The two primary calculations
use 50 decimal digits/N=128 and 80 decimal digits/N=256, respectively.

| Observable | Primary analytic series-tail bound | Refined analytic series-tail bound | Precision/truncation change | Refined proper-time difference |
|---|---:|---:|---:|---:|
'''+ '\n'.join(errorrows)+'''

The explicit geometric bounds control the convergent zeta-series tail only.
They do not certify mpmath transcendental evaluation or floating-point roundoff.
The independent proper-time calculation has analytic harmonic and infrared
tail bounds, but its small-time asymptotic remainder and quadrature error are
tested by refinement without a proved total-error enclosure. Thus the agreement
is strong numerical evidence within the registered absolute/relative gates,
not an interval-certified claim of every displayed digit. No tolerance or grid
was changed after seeing the sources.

Independent centered action differences with Richardson extrapolation agree
with Q to 1.91e−19 and rho to 1.33e−21 at worst. The independently differenced
common-action pairing residual is at most 5.76e−21; the separate digamma Green
function differs from Q by at most 5.45e−31. Negative controls omit the metric
variation, reverse the scalar current, and drop the current paired to a local
F(x)R term; all are rejected at the original gates.

## Prospectivity and scope

The source baseline is 0205cc651bfb614c32e39dfe833d93d957264229. The local
registration and code hashes were frozen at 2026-10-02T06:17:06.731272Z,
before any primary source evaluation. The identical files were publicly
committed/pushed as 85aea9955b99bba911a0e66e5869c184dd86a660 before the
original primary run began at 06:18:58.957725Z and completed at 06:18:59.751705Z.
The independent reviewer ran under an earlier local registration; its run was
completed before this public commit. These distinct timing claims must be
preserved. The registered producer, validator, source inputs and run provenance
are recorded by SHA-256. No primary run failed or required a repair.

This is a bounded mathematical source benchmark. It does not establish a
physical shell parameter choice, solve either junction or a coupled radial
boundary-value problem, evolve the bulk or shell, prove quantum stability,
or calculate heating or particle production. Its four points are not an
extrapolation over the full mass/curvature domain.
'''
    (a.output/'NUMERICAL_RESULTS.md').write_text(md)
    print(json.dumps({'status':'PASS','cross_method_gates':len(comparison),'output':str(a.output)},indent=2))

if __name__=='__main__':main()
