#!/usr/bin/env python3
"""Post-run summaries of registered saved stress results; no producer imports."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path


FREEZE='4a5dad8dda6d0a57a9cf88cc8b4a9b8feeaea45d'
REGISTRATION='4a1533fa7bf9df41cc635d76fa61a22b3aa119d8fbf215349fe09cf3ee54c893'
NAMES=('q','q_prime','q_second','rho','p','Q0','anomaly','current')


def require(condition,message):
    if not condition:raise RuntimeError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(path):return json.loads(Path(path).read_text())
def maximum_by(items,key,value):
    return {name:max(row[value] for row in items if row[key]==name) for name in sorted({r[key] for r in items})}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--primary',type=Path,required=True)
    parser.add_argument('--modes',type=Path,required=True)
    parser.add_argument('--checks',type=Path,required=True)
    parser.add_argument('--checks-optimized',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(not args.output.exists(),'Refusing to overwrite a post-run summary')
    require(sha(args.checkpoint/'FULL_REGISTRATION.json')==REGISTRATION,'Registration digest changed')
    registered=load(args.checkpoint/'FULL_REGISTRATION.json')
    for name,digest in registered['files'].items():
        require(sha(args.checkpoint/name)==digest,'Frozen source changed: '+name)
    primary,modes,checks,optimized=map(load,(args.primary,args.modes,args.checks,args.checks_optimized))
    require(primary['status']=='completed' and modes['status']=='passed_internal_gates','Physical producers did not pass')
    require(primary['provenance']['public_freeze_commit']==modes['freeze_commit']==FREEZE,'Producer freeze mismatch')
    for label,document,flag in (('normal',checks,0),('optimized',optimized,1)):
        require(document['status']=='PASS' and document['python_optimization']==flag,label+' validation status/optimization mismatch')
        require(document['registration_sha256']==REGISTRATION and document['public_freeze_commit']==FREEZE,label+' prospective pin mismatch')
        require(document['input_sha256']=={'primary':sha(args.primary),'modes':sha(args.modes)},label+' raw data hashes mismatch')
    require({k:v for k,v in checks.items() if k!='python_optimization'}=={k:v for k,v in optimized.items() if k!='python_optimization'},
            'Normal and optimized validation results differ')
    require(len(primary['rows'])==len(modes['rows'])==12 and len(modes['runs'])==4,'Registered producer coverage changed')
    mode_rows={(r['source'],r['eta']):r for r in modes['rows']}
    rows=[]
    for row in primary['rows']:
        independent=mode_rows[(row['source'],row['eta'])]
        points=[]
        for finite,direct in zip(row['finite_k'],independent['finite_k']):
            require(finite['K']==direct['K'],'Cutoff row mismatch')
            points.append({'K':finite['K'],'primary_values':finite['values'],'direct_mode_values':direct['values'],
                           'primary_to_continuum_differences':{n:abs(finite['values'][n]-row['continuum']['values'][n]) for n in NAMES},
                           'direct_to_continuum_differences':{n:abs(direct['values'][n]-row['continuum']['values'][n]) for n in NAMES},
                           'direct_to_primary_differences':{n:abs(direct['values'][n]-finite['values'][n]) for n in NAMES},
                           'quadrature_error_estimates':finite['quadrature_error_estimates'],
                           'analytic_omitted_band_bounds':{n:finite['analytic_tail_bounds']['Q0_difference' if n=='Q0' else n] for n in NAMES},
                           'mode_refinement_differences':direct['refinement_differences'],
                           'ward':direct['ward']})
        a=row['a'];eps=primary['epsilon'];v=row['continuum']['values']
        rows.append({'source':row['source'],'eta':row['eta'],'a':a,
                     'source_jet_over_epsilon':row['source_jet_over_epsilon'],
                     'continuum_values':v,'continuum_quadrature_error_estimates':row['continuum']['quadrature_error_estimates'],
                     'physical_continuum_response':{'delta_Q':eps*v['q']/a**2,'delta_rho':eps*v['rho']/a**4,
                                                    'delta_p':eps*v['p']/a**4,'delta_j':eps*v['current']/a**2},
                     'continuum_components':row['continuum']['components'],'finite_k':points})
    final=[r['finite_k'][-1] for r in rows]
    allpoints=[p for r in rows for p in r['finite_k']]
    controls=checks['controls']
    sensitivities=[c for c in controls if c.get('nonzero_algebraic_sensitivity_witness')]
    guards=[c for c in controls if c.get('rejected')]
    metrics={'counts':checks['counts'],
             'maximum_core_cross_route_difference':maximum_by(checks['core_comparisons'],'quantity','gap'),
             'maximum_extra_cross_route_difference':maximum_by(checks['baseline_anomaly_current_comparisons'],'quantity','gap'),
             'maximum_direct_stress_closed_difference':maximum_by(checks['direct_stress_closed_comparisons'],'quantity','gap'),
             'maximum_K256_direct_to_continuum_difference':{n:max(p['direct_to_continuum_differences'][n] for p in final) for n in NAMES},
             'maximum_K256_primary_to_continuum_difference':{n:max(p['primary_to_continuum_differences'][n] for p in final) for n in NAMES},
             'maximum_K256_analytic_tail_bound':{n:max(p['analytic_omitted_band_bounds'][n] for p in final) for n in NAMES},
             'maximum_primary_quadrature_estimate':{n:max([r['continuum']['quadrature_error_estimates'][n] for r in primary['rows']]+[p['quadrature_error_estimates'][n] for r in primary['rows'] for p in r['finite_k']]) for n in NAMES},
             'maximum_mode_refinement_difference':modes['maximum_refinement_differences'],
             'maximum_fine_Ward_endpoint_difference':modes['maximum_ward_endpoint_difference'],
             'maximum_Ward_refinement_difference':modes['maximum_ward_refinement_difference'],
             'maximum_finite_K_trace_difference':max(r['gap'] for r in checks['finite_K_trace_comparisons']),
             'maximum_linear_Wronskian_over_epsilon':max(r['wronskian_max_scaled'] for r in modes['runs']),
             'maximum_physical_Wronskian_over_epsilon':max(r['physical_wronskian_max_scaled'] for r in modes['runs']),
             'algebraic_or_saved_operator_sensitivity_controls':len(sensitivities),
             'sensitivity_controls_exceeding_operational_cross_gate':sum(c['operational_gate_exceeded_somewhere'] for c in sensitivities),
             'synthetic_guard_mutations_rejected':len(guards),
             'primary_elapsed_seconds':primary['elapsed_seconds'],
             'independent_elapsed_seconds':modes['resources']['elapsed_seconds'],
             'independent_peak_rss_kib':modes['resources']['peak_rss_kib']}
    hashes={'FULL_REGISTRATION.json':REGISTRATION,'outputs/primary/results.json':sha(args.primary),
            'outputs/independent/results.json':sha(args.modes),'outputs/validation/CHECKS.json':sha(args.checks),
            'outputs/validation_optimized/CHECKS.json':sha(args.checks_optimized),'postrun/summarize_stress.py':sha(__file__)}
    for run in modes['runs']:
        archive=args.modes.parent/run['archive']['path']
        require(sha(archive)==run['archive']['sha256'],'Independent archive changed')
        hashes['outputs/independent/'+archive.name]=sha(archive)
    for name in ('started.json','partial_results.json','EXECUTION.json'):
        hashes['outputs/primary/'+name]=sha(args.primary.parent/name)
    summary={'status':'registered_matched_stress_calibration_passed','public_freeze_commit':FREEZE,
             'registration_sha256':REGISTRATION,'normalizations':primary['normalization'],'epsilon':primary['epsilon'],
             'postrun_only':True,'new_physical_source_response_or_mode_evaluations':0,
             'registered_rows':rows,'metrics':metrics,'controls':controls,
             'Ward_endpoints':checks['raw_direct_history_Ward_endpoints'],'trace_comparisons':checks['finite_K_trace_comparisons'],
             'core_cross_route_comparisons':checks['core_comparisons'],
             'raw_and_presentation_sha256':hashes,'frozen_source_sha256':registered['files'],
             'uncertainty':checks['uncertainty'],'scope':checks['scope']}
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)+'\n')
    lines=['# Registered matched stress and current calibration','',
           'The prescribed-geometry calibration passed both normal and optimized validators. At the twelve registered source/observation points, independent forced canonical modes with direct minimal stress agree with the logarithmic-memory and matched closed-stress route. Independent raw-mode and complete direct-stress histories support the finite-band trace and Ward checks.','',
           'This is a linear homogeneous scalar-to-stress/current response on fixed de Sitter geometry with an unchanged incoming BD state. It does not compute metric response, coupled shell dynamics, stability, particle yield, radiation transfer, thermalization, or heating.','',
           '## Registration and scales','',
           'H=1, r=2, a=-1/eta, epsilon=1e-4, and b=1 for the current translation. The compact sources s=a²delta_x are epsilon*B(eta+4) and epsilon*(eta+4)B(eta+4), supported on -5<eta<-3. Only eta={-5.5,-4.5,-4,-3.5,-2.5,-1.5} and K={64,128,256} enter the displayed response plots.','',
           '**Every numerical table below uses the registered normalization:** q=a²deltaQ/epsilon; q_prime and q_second differentiate a²deltaQ before division by epsilon; R=a⁴delta_rho/epsilon; P=a⁴delta_p/epsilon; J=a²delta_j/epsilon. Q0K alone is a physical baseline without epsilon division. Thus delta_rho=epsilon*R/a⁴, delta_p=epsilon*P/a⁴, and delta_j=epsilon*J/a². The changing conformal factors matter when comparing physical amplitudes at different times. `SUMMARY.json` records both normalized and physical continuum values.','',
           '## Independent agreement and uncertainty','',
           '| Quantity | Max direct-mode / primary finite-K gap | Max refinement gap | Max K=256 / continuum gap | Max analytic K=256 tail |',
           '|---|---:|---:|---:|---:|']
    for n in ('q','q_prime','q_second','rho','p','current'):
        cross=(metrics['maximum_core_cross_route_difference']|metrics['maximum_extra_cross_route_difference'])[n]
        lines.append(f"| {n} | {cross:.8e} | {metrics['maximum_mode_refinement_difference'][n]:.8e} | {metrics['maximum_K256_direct_to_continuum_difference'][n]:.8e} | {metrics['maximum_K256_analytic_tail_bound'][n]:.8e} |")
    lines+=['',
            'The fourth column is an observed difference, whereas the final column is a constructive analytic enclosure of the omitted combined momentum band. The derivative/stress tail bounds include the higher source derivatives and finite-band local-contact differences. They are not the earlier variance bound reused as a stress bound. Quadrature estimates, finite-step refinement and floating-point allowances remain numerical evidence; neither tiny observed agreement nor the UV theorem certifies the complete numerical result.','',
            f"The maximum fine Ward-ledger endpoint residual is {metrics['maximum_fine_Ward_endpoint_difference']:.9e}; the maximum coarse/fine Ward-ledger difference is {metrics['maximum_Ward_refinement_difference']:.9e}. The maximum finite-band trace residual is {metrics['maximum_finite_K_trace_difference']:.9e}. These quantities use the registered normalization. The ledger independently integrates the saved direct-stress history; the density was not defined by that ledger. The finite-band checks retain Q0K and the finite-band local trace remainder.",'',
            f"The maximum linear Wronskian residual/epsilon is {metrics['maximum_linear_Wronskian_over_epsilon']:.9e}, and the independently reconstructed physical residual/epsilon is {metrics['maximum_physical_Wronskian_over_epsilon']:.9e}. Observation modes and direct-stress histories are archived; global all-time Wronskian maxima are recorded producer evidence.",'',
            '## All registered stress/current values','',
            '| Source | eta | f=s/epsilon | R continuum | P continuum | J continuum | K=256 rho UV bound | K=256 p UV bound |',
            '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        v=r['continuum_values'];tail=r['finite_k'][-1]['analytic_omitted_band_bounds']
        lines.append(f"| {r['source']} | {r['eta']:g} | {r['source_jet_over_epsilon']['f0']:.8g} | {v['rho']:.11g} | {v['p']:.11g} | {v['current']:.11g} | {tail['rho']:.5e} | {tail['p']:.5e} |")
    lines+=['',
            'Before the pulse the perturbative responses vanish. At observations after support, the source and all source jets are zero, so the local source contacts vanish while retarded coherent variance/stress terms can remain. A signed first-order stress perturbation is not a positive particle-energy yield. The Ward source at this order is Q0*delta_x_prime/2; the product deltaQ*delta_x_prime first appears at second order.','',
            '## Control accounting','',
            f"The validator records {len(controls)} controls: {len(sensitivities)} nonzero algebraic or saved-operator sensitivity witnesses and {len(guards)} explicitly rejected synthetic metadata/raw-format mutations. Of the sensitivity witnesses, {metrics['sensitivity_controls_exceeding_operational_cross_gate']} exceed their corresponding operational cross-route tolerance somewhere on the registered grid. The sensitivity witnesses are not all described as operational gate rejections. They reuse saved modes or exact contact algebra; they are not additional altered-state physical experiments.",'',
            '| Control | Classification | Recorded outcome |','|---|---|---|']
    for c in controls:
        if c.get('rejected'):
            outcome='Synthetic guard rejected';kind='Metadata/raw-format mutation'
        else:
            outcome=f"max residual {c['maximum_absolute_residual']:.7e}; operational gate exceeded: {c['operational_gate_exceeded_somewhere']}"
            kind='Saved-operator/algebraic sensitivity'
        lines.append(f"| {c['name']} | {kind} | {outcome} |")
    lines+=['','## Reproducibility and presentation','',
            f"Public prospective freeze: `{FREEZE}`. Full-registration SHA256: `{REGISTRATION}`. Primary runtime {metrics['primary_elapsed_seconds']:.6f} s; independent runtime {metrics['independent_elapsed_seconds']:.6f} s; independent peak RSS {metrics['independent_peak_rss_kib']} KiB. Timings describe this execution and are not portable performance guarantees.",'',
            '`stress_response.svg/.png/.pdf` plots normalized density, pressure and current at registered times. `stress_uv_bounds.svg/.png/.pdf` separates observed K=256 discrepancies from analytic removed-band bounds. `stress_checks.svg/.png/.pdf` summarizes registered finite-K agreement and conservation residuals. Lines only guide the eye; no denser response history was evaluated for presentation. Zero values are omitted on logarithmic axes.','',
            'This report and renderer are explicitly post-run. No source, quadrature, mode, or response evaluation was added by them, and no frozen scientific source was modified. Full-precision values, all raw archive hashes and all frozen source hashes are retained in `SUMMARY.json`; figures have their own rendering provenance.','',
            '## Exact raw and presentation hashes','','| Artifact | SHA256 |','|---|---|']
    lines += [f"| {name} | `{digest}` |" for name,digest in hashes.items()]
    (args.output/'NUMERICAL_RESULTS.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(metrics,indent=2,sort_keys=True))


if __name__=='__main__':main()
