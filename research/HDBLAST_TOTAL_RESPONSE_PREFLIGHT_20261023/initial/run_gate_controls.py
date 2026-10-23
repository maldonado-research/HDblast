#!/usr/bin/env python3
import argparse,copy,importlib.util,json,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--report',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
here=Path(__file__).resolve().parent;s=importlib.util.spec_from_file_location('checker',here/'check_initial_gate_ledger.py');checker=importlib.util.module_from_spec(s);s.loader.exec_module(checker)
base=json.loads(a.report.read_text());controls=[]
mutations={
 'coverage_gap':lambda r:r['partition']['rows'][0].update({'r':['129/128','130/128']}),
 'forged_small_H_error':lambda r:r['error_budget'].update({'H_error_caps':['0','0','0']}),
 'forged_success':lambda r:r['error_budget'].update({'passes':True}),
 'false_impossibility':lambda r:r.update({'mathematical_impossibility_claimed':True}),
 'changed_point_radius':lambda r:r['point_preflight'][0].update({'scalar_radial_oracle_radius':'0'}),
}
with tempfile.TemporaryDirectory() as td:
 for name,mutate in mutations.items():
  r=copy.deepcopy(base);mutate(r);path=Path(td)/'tampered.json';path.write_text(json.dumps(r))
  try:checker.run(path)
  except ValueError as e:controls.append({'name':name,'rejected':True,'reason':str(e)})
  else:raise ValueError('control accepted '+name)
a.output.write_text(json.dumps({'status':'PASS_ALL_INITIAL_GATE_LEDGER_REJECTION_CONTROLS','controls':controls},indent=2)+'\n')
print('PASS_ALL_INITIAL_GATE_LEDGER_REJECTION_CONTROLS')
