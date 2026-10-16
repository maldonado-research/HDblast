#!/usr/bin/env python3
"""Normal/optimized complete reports and wrong-formula/input process controls."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
CASES=[('drop_dynamic_lift','dynamic boundary third-jet retained'),('half_comparator','explicit product-profile commutator'),('drop_initial_offset','stabilized third initial jet A(0) residual(0)+residual_x(0)'),('drop_radial_residual','third radial recovery retains residual_r')]
def need(ok,msg):
 if not ok:raise ValueError(msg)
def invoke(script,out,mode,extra=(),source=ROOT):
 cmd=[sys.executable,'-B']+(['-O'] if mode else [])+[str(source/script),'--output',str(out)]+list(extra)
 return subprocess.run(cmd,capture_output=True,text=True,timeout=120)
def run():
 results=[]
 with tempfile.TemporaryDirectory(prefix='hdblast-numeric-controls-') as temp:
  td=Path(temp)
  for optimized in (False,True):
   for mutation,reason in CASES:
    out=td/(mutation+str(optimized)+'.json');p=invoke('check_approximation_errors.py',out,optimized,('--inject-error',mutation));need(p.returncode!=0 and reason in p.stderr and not out.exists(),'wrong-formula intended rejection: '+mutation)
    results.append({'mode':'optimized' if optimized else 'normal','mutation':mutation,'reason':reason,'report_absent':True})
   out=td/('empty-background'+str(optimized)+'.json');p=invoke('evaluate_background_approximation.py',out,optimized,('--cells','0'));reason='nonempty background collar partition';need(p.returncode!=0 and reason in p.stderr and not out.exists(),'empty collar intended rejection');results.append({'mode':'optimized' if optimized else 'normal','mutation':'empty_background_partition','reason':reason,'report_absent':True})
   out=td/('initial-scalar'+str(optimized)+'.json');p=invoke('evaluate_initial_profiles.py',out,optimized,('--inject-error','omit_scalar_initial','--cells','64'));reason='nonzero second scalar initial data retained';need(p.returncode!=0 and reason in p.stderr and not out.exists(),'wrong initial scalar intended rejection');results.append({'mode':'optimized' if optimized else 'normal','mutation':'omit_scalar_initial','reason':reason,'report_absent':True})
  # True input corruption: keep frozen pin and modify exactly one input byte.
  fixture=td/'fixture';fixture.mkdir();shutil.copy2(ROOT/'check_approximation_errors.py',fixture/'check_approximation_errors.py');shutil.copy2(ROOT/'INPUT_PINS.json',fixture/'INPUT_PINS.json');shutil.copytree(ROOT/'inputs',fixture/'inputs')
  pins=json.loads((fixture/'INPUT_PINS.json').read_text());wrong=fixture/pins['master']['path'];wrong.write_bytes(wrong.read_bytes()+b' ')
  for optimized in (False,True):
   out=td/('wrong-byte'+str(optimized)+'.json');p=invoke('check_approximation_errors.py',out,optimized,source=fixture);reason='input authentication: master';need(p.returncode!=0 and reason in p.stderr and not out.exists(),'wrong input intended rejection');results.append({'mode':'optimized' if optimized else 'normal','mutation':'changed_master_input_byte','reason':reason,'report_absent':True})
 # Full final reports, excluding only output paths (paths are absent from reports).
 comparisons=[]
 for name in ('APPROXIMATION_ERRORS','INITIAL_PROFILES','BACKGROUND_APPROXIMATION'):
  a=(ROOT/(name+'_NORMAL.json')).read_bytes();b=(ROOT/(name+'_OPTIMIZED.json')).read_bytes();need(a==b,'complete normal optimized report mismatch: '+name);comparisons.append({'report_family':name,'bytes':len(a),'sha256':hashlib.sha256(a).hexdigest(),'complete_byte_identity':True})
 return {'status':'PASS_NORMAL_OPTIMIZED_REPLAY_AND_14_WRONG_INPUT_FORMULA_PROCESSES','complete_report_comparisons':comparisons,'rejected_processes':results,'scope':'Process controls enforce specific independently justified identities and input pins; the mathematical proof and independent review establish continuous validity.'}
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();need(not a.output.exists(),'refuse overwrite');r=run();a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'status':r['status'],'rejected_processes':len(r['rejected_processes'])}))
if __name__=='__main__':main()
