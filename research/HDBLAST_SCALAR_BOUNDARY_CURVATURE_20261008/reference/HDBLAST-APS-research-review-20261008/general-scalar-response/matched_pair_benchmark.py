"""Additional matched variations; uses original reviewed second-order solver.

This isolates coefficient dependence but is not an independent solver implementation.
"""
import json, hashlib
from pathlib import Path
from family_benchmark import solve_case, CASES, DELTAS
OUT=Path(__file__).resolve().parent

def main():
    baseline=json.loads((OUT/'family-benchmark-results.json').read_text())
    base=baseline['cases'][0]
    original=dict(CASES[0])
    variants=[dict(original,name='matched_bulk_cubic',v=-3.),
              dict(original,name='matched_tension_curvature',f2=1.2)]
    result=dict(complete=False,scope='Matched coefficient-dependence check using the original numerical solver; not independent software or a proof.',
                baseline_script_sha256=hashlib.sha256((OUT/'family_benchmark.py').read_bytes()).hexdigest(),
                baseline_results_sha256=hashlib.sha256((OUT/'family-benchmark-results.json').read_bytes()).hexdigest(),
                cases=[])
    for variant in variants:
        case=dict(parameters=variant,rows=[])
        for delta in DELTAS:
            standard=solve_case(variant,delta,False)
            refined=solve_case(variant,delta,True)
            b=next(r['refined'] for r in base['rows'] if r['delta']==delta)
            case['rows'].append(dict(delta=delta,standard=standard,refined=refined,
                 refinement_difference=abs(standard['H2']-refined['H2']),
                 H2_difference_from_baseline=refined['H2']-b['H2'],
                 H2_difference_from_baseline_over_delta3=(refined['H2']-b['H2'])/delta**3,
                 coefficient_difference_from_baseline=refined['H2_minus_metric_over_delta2']-b['H2_minus_metric_over_delta2']))
            print(variant['name'],delta,refined['numerical_checks_pass'],flush=True)
        result['cases'].append(case)
        (OUT/'matched-pair-results.json').write_text(json.dumps(result,indent=2)+'\n')
    result['complete']=True
    result['pass_all']=all(r[t]['numerical_checks_pass'] for c in result['cases'] for r in c['rows'] for t in ['standard','refined'])
    result['integrations']=12
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (OUT/'matched-pair-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(complete=True,pass_all=result['pass_all'],integrations=12)),flush=True)

if __name__=='__main__':main()
