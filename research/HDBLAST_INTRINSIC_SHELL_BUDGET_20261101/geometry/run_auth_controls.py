#!/usr/bin/env python3
"""Verify failed authentication precedes third-party imports, using isolated copies."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile

def need(ok,label):
 if not ok:raise ValueError(label)
def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();here=Path(__file__).resolve().parent
 records=json.loads((here/'SOURCE_PINS.json').read_text());results=[]
 optimization=['-O']*sys.flags.optimize
 with tempfile.TemporaryDirectory(prefix='shell-geometry-auth-controls-') as tmp:
  tmp=Path(tmp);source=tmp/'source';source.mkdir()
  for row in records['sources']:
   rel=Path(row['path']);dst=source/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.repo_root/rel,dst)
  altered=source/records['sources'][0]['path'];altered.write_bytes(altered.read_bytes()+b'\n')
  # -S prevents site-package imports. An early SymPy import would fail with a
  # different error; the authentication error must occur first.
  proc=subprocess.run([sys.executable,'-I','-S','-B']+optimization+[str(here/'run_geometry.py'),'--repo-root',str(source),'--output',str(tmp/'badsource.json')],capture_output=True,text=True)
  lines=proc.stderr.rstrip().splitlines()
  need(proc.returncode!=0 and lines and lines[-1]=='ValueError: SOURCE_MEMBER:'+records['sources'][0]['path'] and proc.stdout=='' and not (tmp/'badsource.json').exists() and not (tmp/'badsource.json').is_symlink(),'SOURCE_REJECT_BEFORE_IMPORT')
  results.append({'name':'changed_pinned_source','status':'REJECTED_BEFORE_THIRD_PARTY_IMPORT','expected_final_witness':'ValueError: SOURCE_MEMBER:'+records['sources'][0]['path'],'observed_final_witness':lines[-1],'subprocess_exit':proc.returncode,'output_absent':not (tmp/'badsource.json').exists() and not (tmp/'badsource.json').is_symlink()})
  component=tmp/'component';component.mkdir()
  for name in ['run_geometry.py','COMPONENT_PINS.json']+[r['path'] for r in json.loads((here/'COMPONENT_PINS.json').read_text())['files']]:shutil.copyfile(here/name,component/name)
  altered=component/'check_geometry.py';altered.write_bytes(altered.read_bytes()+b'\n')
  proc=subprocess.run([sys.executable,'-I','-S','-B']+optimization+[str(component/'run_geometry.py'),'--repo-root',str(a.repo_root),'--output',str(tmp/'badcomponent.json')],capture_output=True,text=True)
  lines=proc.stderr.rstrip().splitlines()
  need(proc.returncode!=0 and lines and lines[-1]=='ValueError: COMPONENT_MEMBER:check_geometry.py' and proc.stdout=='' and not (tmp/'badcomponent.json').exists() and not (tmp/'badcomponent.json').is_symlink(),'COMPONENT_REJECT_BEFORE_IMPORT')
  results.append({'name':'changed_checker_source','status':'REJECTED_BEFORE_THIRD_PARTY_IMPORT','expected_final_witness':'ValueError: COMPONENT_MEMBER:check_geometry.py','observed_final_witness':lines[-1],'subprocess_exit':proc.returncode,'output_absent':not (tmp/'badcomponent.json').exists() and not (tmp/'badcomponent.json').is_symlink()})
 report={'status':'PASS','expected_rejections':results,'subprocesses':2,'oracle_calls':0,'site_packages_disabled_for_negative_controls':True}
 need(not a.output.exists(),'FRESH_OUTPUT');a.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');print(json.dumps({'status':'PASS','expected_rejections':2,'report_sha256':hashlib.sha256(a.output.read_bytes()).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
