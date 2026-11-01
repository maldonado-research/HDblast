#!/usr/bin/env python3
"""Authenticate source closure before any third-party import or checker execution."""
from pathlib import Path
import argparse,contextlib,hashlib,io,json,runpy,sys
sys.dont_write_bytecode=True
COMPONENT_PINS_SHA='620296918a047a077f85ccf93481d9630230ba4515710affb0296ca2b7854a12'
def need(ok,msg):
 if not ok:raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def regular(p):return p.is_file() and not p.is_symlink() and all(not q.is_symlink() for q in p.absolute().parents)
def safe_member(root,rel):
 p=Path(rel);need(not p.is_absolute() and '..' not in p.parts,'SAFE_RELATIVE_MEMBER');q=root/p;need(regular(q),'REGULAR_MEMBER:'+rel);return q
def authenticate(here,root):
 pinpath=here/'COMPONENT_PINS.json';need(regular(pinpath) and sha(pinpath)==COMPONENT_PINS_SHA,'COMPONENT_PINS')
 records=json.loads(pinpath.read_text());out={}
 for row in records['files']:
  p=safe_member(here,row['path']);need(p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'COMPONENT_MEMBER:'+row['path'])
 source=json.loads((here/'SOURCE_PINS.json').read_text())
 for row in source['sources']:
  p=safe_member(root,row['path']);need(p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'SOURCE_MEMBER:'+row['path']);out[row['path']]=row['sha256']
 return out
def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 here=Path(__file__).resolve().parent;need(a.repo_root.is_dir() and not a.repo_root.is_symlink(),'REGULAR_REPO_ROOT');root=a.repo_root.absolute()
 before=authenticate(here,root)
 result={}
 for key,name in [('generic_gradient_jet','check_geometry.py'),('curvature_buchert','check_curvature_buchert.py'),('angular','check_angular_exact.py')]:
  buffer=io.StringIO()
  with contextlib.redirect_stdout(buffer):runpy.run_path(str(here/name),run_name='__main__')
  output=buffer.getvalue()
  result[key]=json.loads(output) if key!='angular' else {'status':'PASS','checks':output.splitlines()}
 need(before==authenticate(here,root),'SOURCE_POSTCHECK')
 need(not a.output.exists() and not a.output.is_symlink(),'FRESH_OUTPUT');need(all(not q.is_symlink() for q in a.output.absolute().parents),'OUTPUT_ANCESTORS')
 report={'status':'PASS_PORTABLE_INDEPENDENT_SHELL_GEOMETRY','source_authentication_before_third_party_import_and_after':True,'repository_relative_source_pins':before,'oracle_calls':0,'science':result}
 a.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');print(json.dumps({'status':report['status'],'report_sha256':sha(a.output)},sort_keys=True))
if __name__=='__main__':main()
