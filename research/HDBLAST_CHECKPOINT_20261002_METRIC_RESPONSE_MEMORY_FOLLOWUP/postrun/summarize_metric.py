#!/usr/bin/env python3
"""Summarize saved registered metric results; perform no physical evaluation."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

QUANTITIES=('q','q_prime','q_second','rho','p','Q0','rho0','p0','current')
CORE=('q','q_prime','q_second','rho','p')

def require(ok,message):
    if not ok:raise RuntimeError(message)
def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def maximum(rows,key,quantity):
    return max((row[key][quantity] for row in rows),default=0.0)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('checkpoint','primary','modes','checks','checks-optimized','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    require(not args.output.exists(),'Fresh summary output required')
    require(not args.output.resolve().is_relative_to(args.checkpoint.resolve()),'Summary output must be external')
    primary,modes,checks,optimized=[read(path) for path in (args.primary,args.modes,args.checks,args.checks_optimized)]
    require(primary['status']=='completed' and modes['status']=='passed_internal_gates','Completed producers required')
    require(checks['status']==optimized['status']=='PASS' and checks['counts']==optimized['counts'],'Passing normal/optimized validation required')
    registration=read(args.checkpoint/'FULL_REGISTRATION.json')
    receipt=read(args.checkpoint/'FREEZE_RECEIPT.json')
    require(sha(args.checkpoint/'FULL_REGISTRATION.json')==receipt['registration_sha256'],'Registration receipt mismatch')
    indexed={(r['source'],r['eta']):r for r in modes['rows']}
    require(len(primary['rows'])==len(indexed)==12,'Response coverage must be12points')
    registered=[];points=[]
    for row in primary['rows']:
        direct=indexed[(row['source'],row['eta'])]
        view={'source':row['source'],'eta':row['eta'],'a':row['a'],'source_jet_over_epsilon':row['source_jet_over_epsilon'],'canonical_forcing_jet_over_epsilon':row['canonical_forcing_jet_over_epsilon'],'continuum_values':row['continuum']['values'],'continuum_quadrature_error_estimates':row['continuum']['quadrature_error_estimates'],'continuum_components':row['continuum']['components'],'finite_k':[]}
        other={f['K']:f for f in direct['finite_k']}
        for finite in row['finite_k']:
            mode=other[finite['K']]
            require(set(finite['values'])==set(mode['values'])==set(QUANTITIES),'Quantity coverage differs')
            item={'source':row['source'],'eta':row['eta'],'K':finite['K'],'primary_values':finite['values'],'direct_mode_values':mode['values'],'quadrature_error_estimates':finite['quadrature_error_estimates'],'mode_refinement_differences':mode['refinement_differences'],'ward':mode['ward'],'primary_to_continuum_differences':{q:abs(finite['values'][q]-view['continuum_values'][q]) for q in QUANTITIES},'direct_to_continuum_differences':{q:abs(mode['values'][q]-view['continuum_values'][q]) for q in QUANTITIES},'direct_to_primary_differences':{q:abs(mode['values'][q]-finite['values'][q]) for q in QUANTITIES}}
            require(all(math.isfinite(x) for group in ('primary_values','direct_mode_values','direct_to_primary_differences') for x in item[group].values()),'Nonfinite saved values')
            view['finite_k'].append(item);points.append(item)
        require(len(view['finite_k'])==3,'Cutoff coverage differs')
        registered.append(view)
    metrics={'counts':checks['counts'],'maximum_cross_route_difference':{q:maximum(points,'direct_to_primary_differences',q) for q in QUANTITIES},'maximum_refinement_difference':modes['maximum_refinement_differences'],'maximum_primary_quadrature_estimate':{q:max([r['continuum']['quadrature_error_estimates'][q] for r in primary['rows']]+[f['quadrature_error_estimates'][q] for r in primary['rows'] for f in r['finite_k']]) for q in QUANTITIES},'maximum_fine_ward_endpoint_difference':modes['maximum_ward_endpoint_difference'],'maximum_ward_refinement_difference':modes['maximum_ward_refinement_difference'],'primary_elapsed_seconds':primary['elapsed_seconds'],'independent_resources':modes['resources'],'maximum_K256_primary_to_continuum_difference':{q:maximum([x for x in points if x['K']==256],'primary_to_continuum_differences',q) for q in QUANTITIES}}
    tails=checks.get('continuum_tail_comparisons',[])
    require(len(tails)==324,'Registered omitted-band comparison coverage must be324')
    metrics['maximum_K256_analytic_omitted_band_envelope']={q:max(t['analytic_tail_upper_bound'] for t in tails if t['K']==256 and t['quantity']==q) for q in QUANTITIES}
    metrics['maximum_declared_amplitude_wronskian_over_epsilon']=max(run.get('wronskian_max_scaled',0) for run in modes['runs'])
    metrics['maximum_declared_physical_wronskian_over_epsilon']=max(run.get('physical_wronskian_max_scaled',0) for run in modes['runs'])
    inputs={'FULL_REGISTRATION.json':args.checkpoint/'FULL_REGISTRATION.json','FREEZE_RECEIPT.json':args.checkpoint/'FREEZE_RECEIPT.json','primary/results.json':args.primary,'independent/results.json':args.modes,'validation/CHECKS.json':args.checks,'validation_optimized/CHECKS.json':args.checks_optimized,'postrun/summarize_metric.py':Path(__file__)}
    for run in modes['runs']:
        if 'archive' in run:
            inputs['independent/'+run['archive']['path']]=args.modes.parent/run['archive']['path']
    summary={'schema_version':1,'status':'registered_metric_calibration_passed','postrun_only':True,'new_physical_evaluations':0,'public_freeze_commit':receipt['public_freeze_commit'],'registration_sha256':receipt['registration_sha256'],'frozen_source_sha256':registration['files'],'normalizations':primary['normalization'],'epsilon':primary['epsilon'],'registered_rows':registered,'metrics':metrics,'validator_counts':checks['counts'],'continuum_tail_comparisons':tails,'controls':checks['controls'],'raw_and_presentation_sha256':{n:sha(path) for n,path in inputs.items()},'limits':['This is one homogeneous linear conformal metric channel at H=1, fixed x=r=2, xi=0, with fixed incoming BD state.','Finite-K direct density and pressure comparison is the main numerical calibration.','The omitted-band envelopes are deliberately conservative, especially for high derivatives; they do not establish precise continuum stress or its sign.','Quadrature, refinement and arithmetic estimates are empirical, not interval-certified total errors.','Actual shifted-root, general lapse, bulk, boundary and state-response prerequisites remain incomplete.','No coupled shell evolution, stability, heating, particle production or Big Bang claim follows.']}
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)+'\n')
    lines=['# Registered homogeneous metric calibration','',f'Public prospective freeze: `{receipt["public_freeze_commit"]}`.','',f'Normal and optimized validation pass for 12 source/time points and 36 finite-cutoff points. Both density and pressure were evaluated directly; the Ward ledger is an independent check.','','| Quantity | Maximum primary/direct difference | Maximum primary quadrature estimate |','|---|---:|---:|']
    for q in QUANTITIES:lines.append(f'| {q} | {metrics["maximum_cross_route_difference"][q]:.6g} | {metrics["maximum_primary_quadrature_estimate"][q]:.6g} |')
    lines+=['',f'Maximum fine Ward endpoint difference: {metrics["maximum_fine_ward_endpoint_difference"]:.6g}. Maximum Ward refinement difference: {metrics["maximum_ward_refinement_difference"]:.6g}.','',f'The validator records {len(tails)} omitted-band comparisons. These envelopes bound the removed momentum band separately from finite integration and arithmetic estimates. High derivative bounds can be much larger than the observed cutoff differences.','','| Quantity | Maximum primary K=256 difference from continuum formula | Maximum omitted-band envelope |','|---|---:|---:|']
    for q in QUANTITIES:lines.append(f'| {q} | {metrics["maximum_K256_primary_to_continuum_difference"][q]:.6g} | {metrics["maximum_K256_analytic_omitted_band_envelope"][q]:.6g} |')
    lines+=['','Wronskians at saved observations are reconstructed by the validator; global run maxima remain producer declarations.','',*summary['limits'],'']
    (args.output/'NUMERICAL_RESULTS.md').write_text('\n'.join(lines))
    print(json.dumps({'status':summary['status'],'response_points':len(registered),'finite_K_points':len(points),'new_physical_evaluations':0}))

if __name__=='__main__':main()
