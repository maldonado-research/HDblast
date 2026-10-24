#!/usr/bin/env python3
"""Mutate owned copies only; reject caches/symlinks before numerical imports."""
import argparse,hashlib,importlib.util,json,shutil,sys,tempfile
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--dependency-map',type=Path,default=HERE/'DEPENDENCY_MAP.json');ap.add_argument('--output',type=Path,default=HERE/'AUTHENTICATION_CONTROLS.json');ap.add_argument('--scratch-dir',type=Path,help='fresh nonexistent owned directory; removed after controls');args=ap.parse_args()
spec=importlib.util.spec_from_file_location('difference_producer',HERE/'direct_difference_pilot.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
source=json.loads(args.dependency_map.read_text())
if args.scratch_dir is None:scratch=Path(tempfile.mkdtemp(prefix='hdblast-authentication-controls-',dir='/tmp'))
else:
 scratch=args.scratch_dir.resolve()
 if scratch.exists():raise ValueError('scratch directory must be fresh and nonexistent')
 scratch.mkdir()
results=[]
try:
 for label in ('T','N'):shutil.copytree(source[label],scratch/label)
 mapped={**source,'T':str(scratch/'T'),'N':str(scratch/'N')}
 for label,kind in [('T','cache'),('T','extra'),('T','file_symlink'),('T','directory_symlink'),('N','cache'),('N','directory_symlink')]:
  root=scratch/label
  if kind=='cache':p=root/'__pycache__'/'unhashed.cpython-312.pyc';p.parent.mkdir();p.write_bytes(b'noncanonical stale bytecode')
  elif kind=='extra':p=root/'UNDECLARED.txt';p.write_text('extra file')
  elif kind=='file_symlink':p=root/'linked.py';p.symlink_to(root/'MANIFEST.json')
  else:p=root/'linked-directory';p.symlink_to(root/'initial' if label=='T' else root/'numeric',target_is_directory=True)
  calls=[];original=m.loadmod
  def logged(name,path):
   calls.append(name)
   if name!='audit':raise ValueError('numerical import attempted during authentication')
   return original(name,path)
  m.loadmod=logged
  try:
   m.authenticate(mapped)
  except ValueError as e:
   if 'symbolic' not in str(e) and 'undeclared' not in str(e):raise
   if label=='T' and calls:raise ValueError('T import preceded rejection')
   results.append({'label':label,'mutation':kind,'rejected':True,'error':str(e).replace(str(scratch),'<owned-control-copy>'),'modules_loaded_before_rejection':calls})
  else:raise ValueError('mutation accepted')
  finally:
   m.loadmod=original;p.unlink()
   if kind=='cache':p.parent.rmdir()
finally:shutil.rmtree(scratch)
r={'status':'PASS_REJECT_UNHASHED_CACHES_AND_SYMLINK_PATHS_BEFORE_NUMERICAL_IMPORT','producer_sha256':hashlib.sha256((HERE/'direct_difference_pilot.py').read_bytes()).hexdigest(),'controls':results,'original_input_trees_modified':False}
args.output.write_text(json.dumps(r,indent=2)+'\n');print(r['status'])
