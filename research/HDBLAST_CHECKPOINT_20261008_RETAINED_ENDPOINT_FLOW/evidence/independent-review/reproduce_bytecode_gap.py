"""Harmless counterfeit timestamp cache shows that -I -B reads bytecode."""
from pathlib import Path
import hashlib
import importlib._bootstrap_external as bootstrap
import importlib.util
import json
import subprocess
import sys

BASE=Path(__file__).resolve().parent
root=BASE/'fabricated_bytecode_root'
root.mkdir()
source=root/'innocent.py'
source.write_text('VALUE = "AUTHENTICATED_SOURCE"\n')
info=source.stat()
payload=compile('VALUE = "UNREGISTERED_CACHE"\n',str(source),'exec')
cache=Path(importlib.util.cache_from_source(str(source)))
cache.parent.mkdir()
cache.write_bytes(bootstrap._code_to_timestamp_pyc(payload,int(info.st_mtime),info.st_size))
child=subprocess.run([sys.executable,'-I','-B','-c',
    'import sys;sys.path.insert(0,sys.argv[1]);import innocent;print(innocent.VALUE)',str(root)],
    capture_output=True,text=True,check=True)
result={'status':'INITIAL_BYTECODE_CLOSURE_GAP_REPRODUCED',
        'flags':['-I','-B'],'authenticated_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'unregistered_cache_sha256':hashlib.sha256(cache.read_bytes()).hexdigest(),
        'actual_imported_value':child.stdout.strip(),
        'expected_source_value':'AUTHENTICATED_SOURCE','retained_decodes':0,'physical_source_calls':0}
if result['actual_imported_value']!='UNREGISTERED_CACHE':raise RuntimeError('counterexample did not reproduce')
with (BASE/'INITIAL_BYTECODE_GAP.json').open('x') as out:json.dump(result,out,sort_keys=True,indent=2);out.write('\n')
print(json.dumps(result,sort_keys=True))
