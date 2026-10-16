#!/usr/bin/env python3
"""Portable direct negative formula/input controls plus golden mode identity."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PYTHON=[sys.executable,'-I','-S']+(['-O'] if sys.flags.optimize else [])+['-B']

def need(ok,msg):
 if not ok:raise ValueError(msg)

def run():
 pairs=('ACTUAL_FIRST_JETS','JET_ADAPTER_CONTROLS','EXTERIOR_ERRORS','ACTUAL_SOURCES')
 identities=[]
 for name in pairs:
  p=ROOT/(name+'_NORMAL.json');q=ROOT/(name+'_OPTIMIZED.json');raw=p.read_bytes()
  need(raw==q.read_bytes(),'normal/optimized complete report identity: '+name)
  identities.append({'report':name,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
 controls=[]
 cases=[('check_actual_jets.py','drop_gauge_product','characteristic gauge product derivative'),
        ('check_actual_jets.py','drop_radial_conversion','variable-rho physical second derivative conversion'),
        ('check_actual_jets.py','wrong_scalar_shift','registered scalar pullback shift is +pf'),
        ('check_actual_jets.py','drop_radial_rho','registered raw gauge radial frame factor includes rho'),
        ('evaluate_actual_sources.py','comparator_only','actual solution derivative errors retained rather than comparator substitution'),
        ('evaluate_actual_sources.py','halve_hessian','actual full-cell Einstein Hessian negative witness and second-derivative normalization'),
        ('evaluate_actual_sources.py','reverse_metric_forcing','registered reduced metric forcing sign and scale'),
        ('evaluate_actual_sources.py','drop_corner_cell','complete clipped causal triangle box cover including diagonal cells'),
        ('evaluate_actual_sources.py','omit_shell_normalization','physical shell unit-normal affine Hessian and first-row elimination')]
 with tempfile.TemporaryDirectory(prefix='actual-source-controls-') as tmp:
  tmp=Path(tmp)
  for number,(script,mutation,reason) in enumerate(cases):
   output=tmp/f'wrong-{number}.json'
   result=subprocess.run(PYTHON+[str(ROOT/script),'--inject-error',mutation,'--output',str(output)],capture_output=True,text=True)
   need(result.returncode!=0 and reason in result.stderr and not output.exists(),'direct negative formula control: '+mutation)
   controls.append({'script':script,'mutation':mutation,'required_stderr':reason,'exit_nonzero':True,'no_report_written':True})
  # The inherited checker is a separately identified nested control layer;
  # direct controls above exercise this checkpoint's own formula adapters.
  nested=[]
  for mutation,reason in [('drop_dynamic_lift','dynamic boundary third-jet retained'),('half_comparator','explicit product-profile commutator'),('drop_initial_offset','stabilized third initial jet'),('drop_radial_residual','third radial recovery retains residual_r')]:
   output=tmp/('inherited-'+mutation+'.json')
   result=subprocess.run(PYTHON+[str(ROOT/'check_inherited_errors.py'),'--inject-error',mutation,'--output',str(output)],capture_output=True,text=True)
   need(result.returncode!=0 and reason in result.stderr and not output.exists(),'nested inherited control: '+mutation)
   nested.append({'mutation':mutation,'required_stderr':reason,'exit_nonzero':True,'no_report_written':True})
  copied=tmp/'corrupt-input';copied.mkdir();(copied/'inputs').mkdir()
  pins=json.loads((ROOT/'PREDECESSOR_PINS.json').read_text())
  for relative in list(pins)+['actual_first_jets.py','exterior_error_caps.py','check_actual_jets.py','PREDECESSOR_PINS.json']:
   target=copied/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/relative,target)
  bad=copied/'inputs/BUMP_DERIVATIVE_BOUNDS.json';bad.write_bytes(bad.read_bytes()+b' ')
  output=tmp/'corrupt-output.json';result=subprocess.run(PYTHON+[str(copied/'check_actual_jets.py'),'--output',str(output)],capture_output=True,text=True)
  reason='predecessor authentication: inputs/BUMP_DERIVATIVE_BOUNDS.json'
  need(result.returncode!=0 and reason in result.stderr and not output.exists(),'direct copied-input corruption control')
  controls.append({'script':'check_actual_jets.py','mutation':'copied bump input byte corruption','required_stderr':reason,'exit_nonzero':True,'no_report_written':True})
 return {'status':'PASS_ACTUAL_SOURCE_NUMERIC_CONTROLS','complete_normal_optimized_identity':identities,
         'direct_checkpoint_controls':controls,'separate_nested_inherited_controls':nested,
         'total_rejection_processes':len(controls)+len(nested),'only_standard_library':True}

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args();need(not args.output.exists(),'refuse overwrite');r=run();args.output.write_text(json.dumps(r,indent=2)+'\n');print(r['status'])
if __name__=='__main__':main()
