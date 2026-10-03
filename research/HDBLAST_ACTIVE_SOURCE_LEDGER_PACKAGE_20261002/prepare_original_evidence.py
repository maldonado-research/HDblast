"""Post-freeze presentation and evidence copying; no physical evaluation."""
from pathlib import Path
from fractions import Fraction
import json, shutil

stage=Path('/workspace/hdblast-research-work/active-source-20261002')
cp=Path('/workspace/HDblast/research/HDBLAST_ACTIVE_SOURCE_LEDGER_EXECUTION_REPAIR_20261003/HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER')
run=stage/'repair-original-replay'
out=cp/'outputs/registered_repaired_run'
if out.exists(): raise RuntimeError('fresh output required')
shutil.copytree(run/'fresh',out)
ev=cp/'evidence/registered_repaired_run';ev.mkdir()
for name in ['REPLAY.json','EXECUTION.json','COMMAND_PLAN.json','INPUT_MANIFEST.json']:
 shutil.copy2(run/name,ev/name)
shutil.copytree(run/'logs',ev/'logs')
shutil.copy2(stage/'INITIAL_VS_REPAIRED_PRIMARY.json',ev/'INITIAL_VS_REPAIRED_PRIMARY.json')
for name in ['ACTUAL_PRIMARY_INTERPRETATION_002.md','INITIAL_VS_REPAIRED_PRIMARY_ADVISORY.json']:
 shutil.copy2(stage/'primary'/name,ev/name)
p=json.loads((out/'primary/diagnostic.json').read_text())
v=json.loads((out/'validation/VALIDATION.json').read_text())
rows=[]
for group in p['records']:
 control=next(c for c in group['controls'] if c['quadrature_order']==32)
 for r in control['levels']['100']:
  rows.append(dict(source=group['source'],setting=group['setting'],K=r['K'],**{k:r[k] for k in ['D_S','D_cont','E_Q','E_flow','E_momentum','E_operator','E_reconstruction']}))
maxima={k:max(abs(Fraction(r[k])) for r in rows) for k in ['D_S','D_cont','E_flow','E_momentum','E_operator','E_reconstruction']}
gap=Fraction(v['cross_route_max_scientific_gap_exact_rational'])
summary={'classification':v['classification'],'fixed_cases':v['fixed_cases'],'attribution_cases':v['attribution_cases'],'old_metric_status':'FAIL','maxima_exact_rational':{k:str(x) for k,x in maxima.items()},'cross_route_max_scientific_gap_exact_rational':str(gap),'canonical_primary_rows':rows,'scope':'Finite-momentum, fixed-background anchored active-source diagnostic; empirical quadrature controls; no blast inference.','fresh_zip_replay':'PENDING_EXTERNAL_TO_THIS_PACKAGE'}
(cp/'reports/REGISTERED_RESULT_SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
text=f'''# Active-source ledger diagnostic: registered repaired run

The two independently implemented routes agree on **LEDGER_ERROR_DEMONSTRATED** for the fixed interval [-4.5, -3.5]. Eight of twelve registered source/grid/momentum cases meet the strict attribution criterion, and all registered consistency controls pass. Both the normal and optimized validators accept the same classification. Historical metric-calibration and refinement results remain **FAIL**.

The earlier Simpson ledger produces a maximum canonical primary discrepancy |Delta R - S| of {float(maxima['D_S']):.8e}. The separately constructed matched continuity integral reduces |Delta R - I| to at most {float(maxima['D_cont']):.8e}. The largest registered cross-route scientific comparison is {float(gap):.8e}, below the fixed empirical absolute gate 2e-7. Four cases do not cross the strict discrepancy threshold; they remain in the full record.

This attributes a discrepancy in the retained finite-momentum ledger. It does not establish continuum convergence, validate the initial quantum state, solve gravitational backreaction, demonstrate reheating, or support a higher-dimensional origin of the Big Bang. MP80/100 reductions do not turn the native trajectories or promoted binary64 quadrature nodes into fully arbitrary-precision dynamics. Quadrature comparisons are empirical controls, not certified error bounds.

## Accounting and contact conventions

Define D_S = Delta R - S, D_cont = Delta R - I, and E_Q = I - S. Thus D_S = D_cont + E_Q. The primary retains its raw analytic-contact result I_A and a separately computed momentum-contact correction E_momentum, giving the matched I = I_A + E_momentum. The raw and matched discrepancies, endpoint and midpoint terms, source primitive, invariant E-J, action-derived pressure, full mode profiles, three discrete contact/baseline knots, and all source/ledger control combinations remain available in the raw outputs. No final measured density defines the continuity integral, and pressure is not defined from the Ward identity.

## Execution chronology

The original scientific methods and inputs were publicly frozen at 2e7dfddaa6b7aa2b43bcb12c7ce60a17eed39d3c, with 169 remote files verified before evaluation. The first primary completed; the first independent attempt stopped at a registration-marker guard before numerical evaluation. Its failure and primary result are preserved. After that known primary result, a disclosed registration-only repair was independently reviewed and publicly frozen at 79342726f097affa0c38d7726105d109f9b52a5c, with 205 remote registered files verified before the reruns. The numerical methods, inputs, controls and gates were unchanged. The repaired primary matches the first primary exactly across 152,524 projected scientific fields and values. This is a prospectively fixed diagnostic motivated by known earlier failures; the repaired run is not a blind first evaluation.

The complete primary route took 236.037 seconds and peaked at 193,852 KiB; the independent route took 475.934 seconds and peaked at 131,980 KiB. Each complete route stayed within its registered 900-second and 262,144-KiB limits, with non-root execution, one numerical thread and no child processes. Execution receipts and both validator outputs are included.

## Reproduction and next work

This package includes the retained source/input capsules, exact frozen executables, original failed attempt, repaired original results and public-freeze proofs. A full physical replay from a fresh extraction is required separately; its evidence belongs outside the immutable ZIP to avoid a circular package checksum. The distribution report will record that replay's actual state.

The next scientific priority is a separately registered integration/convergence study against an independent error certificate, followed by an explicit conserved stress tensor and gravitational backreaction model. A proposed higher-dimensional mechanism requires specified bulk dynamics, junction conditions, energy transfer and falsifiable observables; this diagnostic supplies none of those by itself.
'''
(cp/'reports/REGISTERED_RESULT.md').write_text(text)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(10,4.5),layout='constrained')
x=list(range(len(rows)))
ax.semilogy(x,[float(abs(Fraction(r['D_S']))) for r in rows],'o-',label='|Delta R - Simpson ledger|')
ax.semilogy(x,[max(float(abs(Fraction(r['D_cont']))),1e-20) for r in rows],'s-',label='|Delta R - matched continuity integral|')
ax.axhline(2e-6,color='gray',ls=':',label='strict discrepancy threshold')
ax.set_xticks(x,[f"{r['source']}\n{r['setting']} K={r['K']}" for r in rows],rotation=65,ha='right',fontsize=7)
ax.set_ylabel('Absolute discrepancy (registered units)')
ax.set_title('Canonical primary cases; independent route agrees within registered controls')
ax.legend(fontsize=8);ax.grid(axis='y',alpha=.25)
(cp/'figures').mkdir(exist_ok=True)
fig.savefig(cp/'figures/ACTIVE_SOURCE_LEDGER_ATTRIBUTION.png',dpi=180)
plt.close(fig)
print(json.dumps({'cases':len(rows),'classification':v['classification'],'maxima':{k:float(x) for k,x in maxima.items()},'cross_route_gap':float(gap)},indent=2))
