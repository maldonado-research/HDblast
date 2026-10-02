#!/usr/bin/env python3
"""Postprocess recorded stationary evidence without evaluating sources or roots."""
import argparse
import hashlib
import json
from pathlib import Path

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--runs',type=Path,required=True)
    ap.add_argument('--checks',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    runs=json.loads(args.runs.read_text());checks=json.loads(args.checks.read_text())
    rows=checks['shift_diagnostics'];gates=checks['gates']
    failed=[g for g in gates if not g['pass']]
    controls=runs['negative_controls']
    finest=[r for r in runs['runs'] if r['resolution']=='refined']
    zref=next(r['H2'] for r in finest if r['gamma']==0)
    data={'status':checks['status'],'gate_count':len(gates),'failed_gates':failed,
          'runs_sha256':sha(args.runs),'checks_sha256':sha(args.checks),
          'postprocessor_sha256':sha(Path(__file__)),
          'public_freeze':runs['public_freeze'],
          'solver_run_count':len(runs['runs']),'control_count':len(controls),
          'solver_failures':runs['failures'],
          'rejected_newton_trials':sum(len(r.get('rejected_trials',[])) for r in runs['runs']+controls),
          'maximum_producer_C':max(abs(r['C']) for r in runs['runs']),
          'maximum_producer_B':max(abs(r['B']) for r in runs['runs']),
          'resolved_shift_count':{k:sum(s['observables'][k]['shift_resolved'] for s in rows) for k in ['eta','H2']},
          'resolved_nonlinear_count':{k:sum(s['observables'][k]['nonlinear_resolved'] for s in rows) for k in ['eta','H2']},
          'positive_gamma_diagnostic_count':len(rows),'shifts':rows,
          'baseline_sources':runs['baseline_sources'],
          'negative_controls':[{'mutation':r['mutation'],'correct_equation_residual':r['correct_equation_residual']} for r in controls],
          'validator_failure_history':'The original frozen validator exited 1 with 48 identical diagnostic-keyword TypeErrors. Its code, CHECKS.json, log and exit are retained. A separately published one-keyword metadata amendment was run on unchanged producer results; no gate, model or grid was changed.'}
    (args.output/'SUMMARY.json').write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    lines=['# Registered stationary shell results','',
      f"Numerical validator status: **{checks['status']}**, {len(gates)} gates, {len(failed)} failed. "
      f"The producer completed {len(runs['runs'])} registered integration/root cases and {len(controls)} wrong-model controls, "
      f"with {len(runs['failures'])} solver failures.",'',
      'These are stationary regular-cone roots of the declared one-field mathematical model. '
      'The sampled hierarchy screen is reported separately; it does not certify an EFT cutoff or physical viability.','',
      '| b | gamma | Delta phi | Delta H²/H0² | phi/H² shift resolved | phi/H² nonlinear resolved | hierarchy screen |',
      '|---:|---:|---:|---:|:---:|:---:|:---:|']
    yes=lambda x:'yes' if x else 'no'
    for s in rows:
        e=s['observables']['eta'];z=s['observables']['H2']
        lines.append(f"| {s['b']} | {s['gamma']:g} | {e['delta_from_own_solved_control']:.10g} | "
                     f"{z['delta_from_own_solved_control']/zref:.10g} | "
                     f"{yes(e['shift_resolved'])}/{yes(z['shift_resolved'])} | "
                     f"{yes(e['nonlinear_resolved'])}/{yes(z['nonlinear_resolved'])} | "
                     f"{'pass' if s['hierarchy']['tenfold_screen_pass'] else 'fail'} |")
    lines+=['','Shifts subtract each integration method’s own solved gamma-zero control. '
      'The registered nonlinear detection requires departure from the archived linear prediction '
      'of at least 0.1% of the measured shift and more than ten empirical refinement spreads plus '
      'the fixed floating-point floor. Detection classification is separate from root acceptance.','',
      '## Largest registered amplitude','']
    for s in rows:
        if s['gamma']!=1000000:continue
        e=s['observables']['eta'];z=s['observables']['H2']
        lines.append(f"- b={s['b']:+d}: Delta phi={e['delta_from_own_solved_control']:.12g}; "
                     f"Delta H²/H0²={z['delta_from_own_solved_control']/zref:.12g}; "
                     f"nonlinear fractions phi={e['nonlinear_fraction']:.8g}, H²={z['nonlinear_fraction']:.8g}; "
                     f"Ehat/M5={s['hierarchy']['energy_over_M5']:.8g}.")
    lines+=['','## Controls, failures and provenance','',
      f"Maximum producer residuals: |C|={data['maximum_producer_C']:.5g}, |B|={data['maximum_producer_B']:.5g}. "
      f"Rejected Newton trial steps: {data['rejected_newton_trials']}.",'']
    for r in controls:
        lines.append(f"- {r['mutation']}: correct-equation (C,B) at the wrong root = {r['correct_equation_residual']}.")
    lines+=['',data['validator_failure_history'],'',
      f"Public freeze: `{runs['public_freeze']}`. Producer results SHA-256: `{data['runs_sha256']}`. "
      f"Amended checks SHA-256: `{data['checks_sha256']}`.",'',
      'No uniqueness, dynamical or quantum stability, heating, particle production, history of exit, '
      'observational fit or physical parameter inference is established.']
    (args.output/'NUMERICAL_RESULTS.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:data[k] for k in ['status','gate_count','solver_run_count','control_count','resolved_shift_count','resolved_nonlinear_count']},indent=2))

if __name__=='__main__':
    main()
