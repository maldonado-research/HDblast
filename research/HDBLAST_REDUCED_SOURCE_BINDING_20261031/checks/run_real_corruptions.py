#!/usr/bin/env python3
"""Fresh normal/-O subprocesses against changed disk bytes, never old files."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ADAPTER_PIN = '4847c79f034a09f423cc020ab3bef86e03d729970084a9bb088e94caa3b9b18a'
GEOM = 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
SELECTED = (
 'HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/MASTER4.json',
 'HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/checks/check_master_fourth.py',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/GAUGE4.json',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/checks/check_actual_gauge_fourth.py',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/checks/check_original_bulk_fx.py',
)
MUTATIONS = (
 ('radial_leading_sign','third=-L**3','third=L**3','direct substitution witness'),
 ('omit_radial_middle','-3*av*L**2*j.derivative(0,2)','','direct substitution witness'),
 ('omit_radial_background','+(ay-av**2)*L*j.derivative(0,1)','','direct substitution witness'),
 ('wrong_radial_background','(ay-av**2)','(ay+av**2)','direct substitution witness'),
 ('omit_yyt_rho','+av*L*j.derivative(1,1)','','direct substitution witness'),
 ('wrong_ytt_sign',"'ytt':-L*j.derivative(2,1)","'ytt':L*j.derivative(2,1)",'direct substitution witness'),
 ('lost_time_ell',"'ttt':j.derivative(3,0)/ell**3","'ttt':j.derivative(3,0)/ell**2",'direct substitution witness'),
 ('omit_normalized_factorials','return math.factorial(i)*math.factorial(j)*self.c[i,j]',
  'return self.c[i,j]','independent composite derivative'),
 ('zero_master_fourth','self.global_third,self.coefficient_boxes=self._global_third()',
  'self.master_caps[4]=F(0)\n        self.global_third,self.coefficient_boxes=self._global_third()',
  'faithful authenticated master0through4 cap join'),
 ('omit_pure_radial_slot','self.global_third,self.coefficient_boxes=self._global_third()',
  "self.global_third,self.coefficient_boxes=self._global_third()\n        del self.global_third['Ayyy']",
  'complete ordinary third-jet interface'),
)

def need(ok,message):
    if not ok: raise ValueError(message)

def sha(raw): return hashlib.sha256(raw).hexdigest()

PROBE = r'''
from pathlib import Path
import hashlib,importlib.util,sys
adapter,root,pin,sentinel,action=map(str,sys.argv[1:])
p=Path(adapter)
if hashlib.sha256(p.read_bytes()).hexdigest()!=pin:raise ValueError('probe adapter pin')
spec=importlib.util.spec_from_file_location('pinned_producer_adapter',p)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
instance=module.Adapter(root)
if action=='post_import':
    q=Path(root)/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/interval_core.py'
    q.write_bytes(q.read_bytes()+b'\n# changed after authenticated import\n')
    instance.actual_jet(0,0)
raise ValueError('UNEXPECTED_ACCEPTANCE')
'''

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--adapter',type=Path,required=True)
    ap.add_argument('--checker',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    need(not args.output.exists(),'fresh receipt required')
    original=args.adapter.read_bytes()
    need(sha(original)==ADAPTER_PIN,'frozen adapter before corruptions')
    initial_checker=sha(args.checker.read_bytes())
    work=args.output.parent/'real-corruption-fixtures'
    need(not work.exists(),'fresh private corruption fixture directory required')
    work.mkdir()
    rows=[]
    for name,old,new,marker in MUTATIONS:
        source=original.decode()
        need(old in source,'mutation anchor '+name)
        modified=source.replace(old,new,1).encode()
        need(modified!=original,'real changed adapter bytes')
        path=work/(name+'.py');path.write_bytes(modified)
        for optimized in (False,True):
            output=work/(name+('-optimized' if optimized else '-normal')+'.json')
            command=[sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])+[
              str(args.checker),'--repo-root',str(args.repo_root),'--adapter',str(path),
              '--adapter-sha256',sha(modified),'--output',str(output)]
            done=subprocess.run(command,capture_output=True,text=True,timeout=120)
            need(done.returncode!=0 and not output.exists() and marker in done.stderr,
                 'real semantic corruption not rejected '+name+': '+done.stderr[-1000:])
            witness=done.stderr.split('ValueError: ',1)[-1].strip().splitlines()[-1]
            rows.append({'case':name,'optimized':optimized,'command':command,'returncode':done.returncode,
                         'changed_adapter_sha256':sha(modified),'output_absent':True,'witness':witness})
    # Authenticate before import using copied byte-corrupted original helpers.
    copied=work/'copied-repo'
    src=args.repo_root/'research'
    shutil.copytree(src/GEOM,copied/'research'/GEOM)
    for rel in SELECTED:
        target=copied/'research'/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src/rel,target)
    sentinel=work/'MUST_NOT_EXECUTE.txt'
    cases=(
      ('corrupt_actual_helper',GEOM+'/numeric/actual_first_jets.py','source member authentication','append'),
      ('corrupt_transitive_helper',GEOM+'/numeric/interval_core.py','source member authentication','append'),
      ('corrupt_gauge4_report',SELECTED[2],'selected higher-jet source','space'),
      ('corrupt_master4_report',SELECTED[0],'selected higher-jet source','space'),
      ('unexpected_source_member',GEOM+'/numeric/UNREGISTERED.py','closed geometric source tree','new'),
      ('post_import_change',GEOM+'/numeric/interval_core.py','source member authentication','post_import'),
    )
    for name,rel,marker,action in cases:
        target=copied/'research'/rel
        before=target.read_bytes() if target.exists() else None
        for optimized in (False,True):
            if before is not None:target.write_bytes(before)
            if action in ('append','new'):
                target.write_bytes((before or b'')+ ('\nfrom pathlib import Path\nPath('+repr(str(sentinel))+').write_text("executed")\n').encode())
            elif action=='space':target.write_bytes(before+b' ')
            need(not sentinel.exists(),'fresh nonexecution sentinel')
            command=[sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])+[
              '-c',PROBE,str(args.adapter),str(copied),ADAPTER_PIN,str(sentinel),action]
            done=subprocess.run(command,capture_output=True,text=True,timeout=120)
            need(done.returncode!=0 and marker in done.stderr and not sentinel.exists(),
                 'source corruption not rejected before use '+name+': '+done.stderr[-1000:])
            witness=done.stderr.split('ValueError: ',1)[-1].strip().splitlines()[-1]
            rows.append({'case':name,'optimized':optimized,'command':command,'returncode':done.returncode,
                         'side_effect_absent':True,'witness':witness})
        if before is None:target.unlink()
        else:target.write_bytes(before)
    need(sha(args.adapter.read_bytes())==ADAPTER_PIN,'frozen original adapter unchanged')
    need(sha(args.checker.read_bytes())==initial_checker,'independent checker unchanged')
    report={'status':'PASS_REAL_THIRD_JET_CORRUPTIONS_REJECTED',
            'adapter_sha256':ADAPTER_PIN,'independent_checker_sha256':initial_checker,
            'semantic_cases':len(MUTATIONS),'source_authentication_cases':len(cases),
            'fresh_subprocesses':len(rows),'executions':rows}
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(report['status'])

if __name__=='__main__':main()
