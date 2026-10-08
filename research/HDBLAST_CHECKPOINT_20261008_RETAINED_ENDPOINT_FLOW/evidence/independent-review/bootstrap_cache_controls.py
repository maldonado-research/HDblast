"""Counterfeit harmless local-guard cache at both supported entry boundaries."""
from pathlib import Path
import hashlib
import importlib._bootstrap_external as bootstrap
import importlib.util
import json
import subprocess
import sys

BASE=Path(__file__).resolve().parent
ROOT=BASE.parent/'implementation'
def run(label):
    destination=BASE/label;destination.mkdir();cases=[]
    for entry,module in (('run_endpoint.py','registration_guard'),('bounded_launcher.py','kernel_guard')):
        root=destination/entry.removesuffix('.py');execution=root/'execution';execution.mkdir(parents=True)
        original=(ROOT/'execution'/entry).read_bytes()
        (execution/entry).write_bytes(original)
        source=execution/(module+'.py');source.write_bytes((ROOT/'execution'/(module+'.py')).read_bytes())
        marker='UNREGISTERED_BOOTSTRAP_EXECUTED'
        cache=Path(importlib.util.cache_from_source(str(source)));cache.parent.mkdir()
        metadata=source.stat()
        cache.write_bytes(bootstrap._code_to_timestamp_pyc(compile('raise RuntimeError("'+marker+'")',str(source),'exec'),
                                                         int(metadata.st_mtime),metadata.st_size))
        result=subprocess.run([sys.executable,'-I','-B',str(execution/entry),'--help'],capture_output=True,text=True)
        log=result.stdout+result.stderr
        (root/'CHILD.log').write_text(log)
        cases.append({'entry':entry,'entry_sha256':hashlib.sha256(original).hexdigest(),
                      'cached_module':module,'unregistered_cache_executed':marker in log,
                      'exit_code':result.returncode,'child_log_sha256':hashlib.sha256(log.encode()).hexdigest()})
    receipt={'status':'BOOTSTRAP_CACHE_CONTROLS','cases':cases,'retained_decodes':0,'physical_source_calls':0}
    with (destination/'RESULT.json').open('x') as out:json.dump(receipt,out,sort_keys=True,indent=2);out.write('\n')
    print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':run(sys.argv[1])
