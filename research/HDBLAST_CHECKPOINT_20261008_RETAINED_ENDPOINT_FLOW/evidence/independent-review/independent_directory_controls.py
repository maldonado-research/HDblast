"""Fault-injected closure scanning controls with harmless fabricated source files."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import sys
import types
BASE=Path(__file__).resolve().parent

def run(guard_path,label):
 guard_path=Path(guard_path);destination=BASE/label;destination.mkdir()
 raw=guard_path.read_bytes();guard=types.ModuleType('independent_directory_guard')
 exec(compile(raw,str(guard_path),'exec'),guard.__dict__)
 plans=['clean','unknown_empty_directory','hidden_execute_only_package','unreadable_registered_directory',
 'directory_symlink','file_symlink','hardlink','injected_scandir_error','replace_directory_before_open',
 'replace_directory_before_scandir','replace_directory_with_symlink_before_open']
 results=[]
 for case in plans:
  root=destination/case;root.mkdir();source=root/'source';source.mkdir();(source/'foo.py').write_text('FABRICATED=1\n')
  files={'source/foo.py':{}};blocked=None;injected=False;original_open=os.open;original_scan=os.scandir
  error=None;accepted=False
  if case=='unknown_empty_directory':(root/'unknown').mkdir()
  if case=='hidden_execute_only_package':
   blocked=source/'foo';blocked.mkdir();(blocked/'__init__.py').write_text('raise RuntimeError("MUST_NEVER_IMPORT")\n');blocked.chmod(0o111)
  if case=='unreadable_registered_directory':blocked=source;source.chmod(0o111)
  if case in ('directory_symlink','file_symlink','hardlink'):
   target=destination/(case+'-outside');target.mkdir();(target/'foo.py').write_text('FABRICATED=1\n')
   if case=='directory_symlink':shutil.rmtree(source);source.symlink_to(target,target_is_directory=True)
   else:
    (source/'foo.py').unlink()
    if case=='file_symlink':(source/'foo.py').symlink_to(target/'foo.py')
    else:(source/'foo.py').hardlink_to(target/'foo.py')
  source_identity=None
  if source.is_dir() and not source.is_symlink():
   meta=source.stat();source_identity=(meta.st_dev,meta.st_ino)
  def is_source(fd):
   meta=os.fstat(fd);return (meta.st_dev,meta.st_ino)==source_identity
  def scan(fd):
   nonlocal injected
   if case in ('injected_scandir_error','replace_directory_before_scandir') and isinstance(fd,int) and is_source(fd):
    injected=True
    if case=='injected_scandir_error':raise PermissionError('FABRICATED_SCANDIR_ERROR')
    source.rename(destination/(case+'-held'));source.mkdir();(source/'foo.py').write_text('FABRICATED=1\n')
   return original_scan(fd)
  def opening(name,flags,*args,**kwargs):
   nonlocal injected
   if case in ('replace_directory_before_open','replace_directory_with_symlink_before_open') and name=='source' and 'dir_fd' in kwargs:
    injected=True;held=destination/(case+'-held');source.rename(held)
    if case=='replace_directory_before_open':source.mkdir();(source/'foo.py').write_text('FABRICATED=1\n')
    else:source.symlink_to(held,target_is_directory=True)
   return original_open(name,flags,*args,**kwargs)
  try:
   os.open=opening;os.scandir=scan
   accepted=guard.source_tree(root,files)==set(files)
  except Exception as exc:error=type(exc).__name__+': '+str(exc)
  finally:
   os.open=original_open;os.scandir=original_scan
   if blocked is not None:blocked.chmod(0o755)
  needs_injection=case.startswith('injected_') or case.startswith('replace_')
  passed=(accepted and error is None) if case=='clean' else (error is not None and not accepted and (injected or not needs_injection))
  results.append({'case':case,'passed':passed,'accepted':accepted,'error':error,'fault_injected':injected})
 receipt={'status':'PASS_INDEPENDENT_DIRECTORY_SCANNER_CONTROLS' if all(x['passed'] for x in results) else 'FAIL_INDEPENDENT_DIRECTORY_SCANNER_CONTROLS',
  'guard_sha256':hashlib.sha256(raw).hexdigest(),'checks':len(results),'cases':results,
  'production_numeric_modules_imported':0,'retained_decodes':0,'physical_source_calls':0}
 (destination/'RESULT.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
 if not all(x['passed'] for x in results):raise SystemExit(1)
if __name__=='__main__':run(*sys.argv[1:])
