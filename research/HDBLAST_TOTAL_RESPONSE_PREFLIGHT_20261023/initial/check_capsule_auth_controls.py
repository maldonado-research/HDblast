#!/usr/bin/env python3
"""Reject manifest-only and changed-byte fake dependencies before importing numerics."""
import argparse,hashlib,importlib.util,json,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--dependency-map',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
here=Path(__file__).resolve().parent;s=importlib.util.spec_from_file_location('auth',here/'authenticate_inputs.py');auth=importlib.util.module_from_spec(s);s.loader.exec_module(auth)
base=json.loads(a.dependency_map.read_text());out=[]
with tempfile.TemporaryDirectory() as td:
 root=Path(td);(root/'MANIFEST.json').write_bytes((Path(base['P'])/'MANIFEST.json').read_bytes());deps={**base,'P':str(root)}
 for label,change in [('manifest_without_closure',False),('changed_declared_file',True)]:
  if change:(root/'LICENSE').write_bytes(b'changed bytes\n')
  try:auth.run(deps)
  except (ValueError,FileNotFoundError) as e:out.append({'name':label,'rejected':True,'reason':'missing declared file' if isinstance(e,FileNotFoundError) else 'declared file byte-length or SHA256 mismatch'})
  else:raise ValueError('invalid capsule accepted '+label)
a.output.write_text(json.dumps({'status':'PASS_COMPLETE_CLOSURE_AUTHENTICATION_REJECTION_CONTROLS','controls':out,'numeric_capsule_modules_imported':False},indent=2)+'\n')
print('PASS_COMPLETE_CLOSURE_AUTHENTICATION_REJECTION_CONTROLS')
