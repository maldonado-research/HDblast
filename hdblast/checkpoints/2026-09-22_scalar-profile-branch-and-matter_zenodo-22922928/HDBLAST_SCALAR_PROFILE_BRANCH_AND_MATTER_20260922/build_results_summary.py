"""Summarize completed calculations; does not solve shooting or evolution."""
from pathlib import Path
import hashlib,json,math

ROOT=Path(__file__).resolve().parent
names=['balanced_epsp01_h2','balanced_epsp01_h1','balanced_epsp01_wide_h2',
       'balanced_epsp01_wide_h1','balanced_epsp001_wide_h2','balanced_epsm001_wide_h2','balanced_epsp01_wide_h05']
hashes={}
def read(path):
    f=ROOT/path;hashes[path]=hashlib.sha256(f.read_bytes()).hexdigest()
    return json.loads(f.read_text())
checks=[]
def check(name,condition):
    checks.append(dict(name=name,passed=bool(condition)))
    if not condition:raise RuntimeError(name)
code_hashes={hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['frozen/evolve_balanced_v1.py','evolution/evolve_balanced.py']}
rows=[]
for name in names:
    d=read('runs/'+name+'.json');last=d['records'][-1]
    check(name+' completed requested bounded run',d['stop_reason']=='tf' and last['time'] in [.5,1.,2.5])
    check(name+' executed code version retained',d['code_sha256'] in code_hashes)
    check(name+' remains before taper can reach shell',last['time']<d['metadata']['earliest_taper_influence_on_shell'])
    rows.append(dict(name=name,epsilon=d['epsilon'],hmin=d['hmin'],stretch=d.get('stretch',.05),L=d['L'],
                     nodes=d['nodes'],dt=d['dt'],t_end=last['time'],initial=d['records'][0],final=last,
                     maximum_over_recorded_times={key:max(r[key] for r in d['records']) for key in
                         ['hamiltonian_max','momentum_max','core_H','core_M','weighted_characteristic_max']},
                     runtime_seconds=d['runtime_seconds']))
branch=read('static_branch/PLUS_BRANCH_RESULTS.json')
selected=next(r['refined'] for r in branch['rows'] if r['delta']==.001)
check('registered static branch satisfies both junctions numerically',max(abs(x) for x in selected['junction_residual'])<1e-10)
check('constant-phi benchmark violates scalar junction',selected['constant_phi_scalar_junction_residual']>1e-4)
check('new branch retains a nonconstant scalar',selected['eta_b']<0 and selected['phi_y_b']!=0)
check('full H2 identity includes scalar boundary gradient',abs(selected['H2_identity_residual'])<1e-12)
metric_H=math.sqrt(selected['metric_only_H2']);new_H=math.sqrt(selected['H2'])
correction_ppm=(new_H/metric_H-1)*1e6
check('small benchmark correction is not percent-scale',abs(correction_ppm)<10)
result=dict(status='PASS_COMPLETED_RESULT_ACCOUNTING',checks=checks,static_registered_candidate=selected,
            expansion_rate_correction_ppm=correction_ppm,evolutions=rows,input_hashes=hashes,
            physical_gate='OPEN: none of these checks certifies full-domain late-time evolution, branch attraction, quantum production or a hot universe.',
            count_scope='Completion and result-accounting checks include explicit negative physical findings; PASS is not a global nonlinear-accuracy verdict.')
(ROOT/'RESULTS_SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['# Completed calculations and their limits','',
       'Seven new finite-amplitude PDE evolutions are retained, including failed accuracy controls. The static branch has six detunings, two producer tolerances, and a separate 18-integration review. The matter extension is derived and screened, not evolved.','',
       '| Run | ε | h at shell | Stretch | L | Final t | Final φ_b | Max C_H/background scale over stored times | Max C_M/background scale over stored times |',
       '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
for r in rows:
    m=r['maximum_over_recorded_times'];v=r['final']
    lines.append(f"| {r['name']} | {r['epsilon']} | {r['hmin']} | {r['stretch']} | {r['L']} | {r['t_end']} | {v['phi_b']:.10g} | {m['hamiltonian_max']:.5g} | {m['momentum_max']:.5g} |")
lines+=['',
        'The denominator is 1+6Hc²+φs,z². The independent residual review also gives actual-term cancellation ratios. These are constraints, not error bars on φ_b. The monitor domain is z>−0.8L with the first six grid nodes excluded; a separate near-shell domain is z>−0.2. Maxima in this table cover stored times and can exceed the final-time values quoted in comparisons.',
        '',f'The corrected registered static rate differs from the old constant-scalar metric benchmark by {correction_ppm:.6f} parts per million. This small correction cannot account for the old late-time evolution gap of several percent.',
        '',f'All {len(checks)} completion/accounting checks pass. The global nonlinear-accuracy gate remains open. The JSON retains exact values, source hashes, first/final rows and maximum-over-time diagnostics.','']
(ROOT/'RESULTS_SUMMARY.md').write_text('\n'.join(lines))
print(json.dumps({'status':result['status'],'checks':len(checks),'completed_PDE_runs':len(rows),'rate_correction_ppm':correction_ppm}))
