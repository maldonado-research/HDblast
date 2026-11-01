#!/usr/bin/env python3
"""Exact-witness formula controls and separate custody rejection control."""
import json,subprocess,sys,tempfile
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent))
from portable_io import arguments,authenticate,finish,HERE
def run(a):
    common=['--repo-root',str(a.repo_root),'--initial-report',str(a.initial_report),
        '--conditioned-report',str(a.conditioned_report),'--response-report',str(a.response_report)]
    prefix=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])
    tests=[('reconstruct_response.py','omit_measure','independent reconstructed mismatch: proper_M2_target'),
      ('reconstruct_response.py','halve_scalar','independent reconstructed mismatch: proper_M2_target'),
      ('reconstruct_response.py','omit_history_error','independent reconstructed mismatch: proper_v2'),
      ('reaggregate_initial.py','wrong_moment','METRIC_MOMENT'),
      ('reaggregate_initial.py','wrong_primitive_sign','CELL_XI')]
    outcomes=[]
    with tempfile.TemporaryDirectory(prefix='prepared-mean-independent-controls-') as tmp:
        temp=Path(tmp)
        for script,mutant,reason in tests:
            dest=temp/(mutant+'.json')
            cp=subprocess.run(prefix+[str(HERE/script),*common,'--mutant',mutant,'--output',str(dest)],text=True,capture_output=True)
            expected={'status':'REJECT','stage':'SCIENCE','reason':reason}
            if cp.returncode!=2 or cp.stderr or len(cp.stdout.splitlines())!=1 or dest.exists() or json.loads(cp.stdout)!=expected:
                raise ValueError('CONTROL_FAILED:'+mutant)
            outcomes.append({'case':mutant,'exit_code':2,'reason':reason,'stage':'SCIENCE','no_output':True,'final_witness':json.loads(cp.stdout),'stderr_empty':True,'stdout_one_json':True})
        bad=temp/'bad-response.json';payload=json.loads(a.response_report.read_text())
        payload['proper_v2']['exact'][1]='0';bad.write_text(json.dumps(payload))
        dest=temp/'bad-result.json';args=list(common);args[args.index('--response-report')+1]=str(bad)
        cp=subprocess.run(prefix+[str(HERE/'reconstruct_response.py'),*args,'--output',str(dest)],text=True,capture_output=True)
        expected={'status':'REJECT','stage':'AUTH_REPORT','reason':'response_report'}
        if cp.returncode!=2 or cp.stderr or len(cp.stdout.splitlines())!=1 or dest.exists() or json.loads(cp.stdout)!=expected:
            raise ValueError('CONTROL_FAILED:report_custody')
        outcomes.append({'case':'report_custody','exit_code':2,'reason':'response_report','stage':'AUTH_REPORT','no_output':True,'final_witness':json.loads(cp.stdout),'stderr_empty':True,'stdout_one_json':True})
    return {'status':'PASS_INDEPENDENT_EXACT_WITNESS_CONTROLS','positive_source_authentication':'PASS',
        'counts':{'scientific_formula_rejections':5,'separate_report_custody_rejections':1,'no_output_rejections':6},'controls':outcomes}
if __name__=='__main__':
    a=arguments();authenticate(a)
    try:finish(a,run(a))
    except (ValueError,OSError) as error:
        print(json.dumps({'status':'REJECT','stage':'CONTROLS','reason':str(error)},sort_keys=True))
        raise SystemExit(2)
