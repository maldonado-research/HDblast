from pathlib import Path
import importlib.util,json,time
import numpy as np
p=Path(__file__).with_name('diagnostic_primary.py');spec=importlib.util.spec_from_file_location('primary_schema',p)
P=importlib.util.module_from_spec(spec);spec.loader.exec_module(P)
n=3;k=np.asarray([P.LD(32),P.LD(96),P.LD(192)],dtype=P.LD);count=577
history=P.LD(-6)+np.arange(count,dtype=P.LD)/128
zero=lambda:np.zeros(n,dtype=P.CD)
a={'k':k,'momentum_weights':np.ones(n,dtype=P.LD),'u_1':zero(),'w_1':zero(),
   'u_2':zero(),'w_2':zero(),'u_3':zero(),'w_3':zero(),'history_eta':history,
   'history_values':np.zeros((count,3,9),dtype=P.LD),'history_ledger_integrand':np.zeros((count,3),dtype=P.LD),
   'history_baseline_contact':np.zeros((count,3),dtype=P.LD),'history_source_jet':P.synthetic_jets(history).T,
   'history_forcing_jet':P.forcing_jet(history,P.synthetic_jets(history)).T}
base=Path('/workspace/HDblast/research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP')
case=P.run_case(a,'positive_B','coarse',base/'code/metric_contact_coefficients.py',time.monotonic(),P.synthetic_jets,
                (base/'independent/metric_wkb.py',base/'independent/stable_baselines.py'))
row=case['controls'][0]['levels']['80'][0]
result={'producer_sha256':P.sha(p),'case_keys':list(case),'control_keys':list(case['controls'][0]),
        'scientific_profile_lengths':{k:len(v) for k,v in row.items() if isinstance(v,list)},
        'scientific_scalar_keys':[k for k,v in row.items() if not isinstance(v,list) and k!='K'],
        'case':case,'analysis':P.analyze_records([case]),'source_guard_authorized':P.AUTHORIZED,
        'scope':'Schema/algebra only: intentionally off-domain three-node momentum rule [32,96,192] with unit weights and zero fabricated stress history. Large numerical discrepancies are expected and have no scientific meaning.'}
path=Path(__file__).with_name('SYNTHETIC_SCHEMA_FINAL.json');path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('case','analysis')},sort_keys=True))
