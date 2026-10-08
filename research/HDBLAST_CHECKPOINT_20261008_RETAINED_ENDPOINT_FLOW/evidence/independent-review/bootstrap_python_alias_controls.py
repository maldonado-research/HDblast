"""Harmless pure-Python bootstrap aliases; no input or physical evaluation."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

BASE=Path(__file__).resolve().parent
def run(label,source_root):
    source_root=Path(source_root);destination=BASE/label;destination.mkdir();cases=[]
    plans=[('run_endpoint.py','registration_guard/__init__.py'),
           ('bounded_launcher.py','kernel_guard/__init__.py'),
           ('run_endpoint.py','resource.py'),('bounded_launcher.py','ctypes.py')]
    for index,(entry,alias) in enumerate(plans):
        root=destination/str(index);execution=root/'execution';execution.mkdir(parents=True)
        for name in (entry,'registration_guard.py','kernel_guard.py'):
            (execution/name).write_bytes((source_root/'execution'/name).read_bytes())
        marker='UNREGISTERED_PYTHON_ALIAS_EXECUTED'
        injected=execution/alias;injected.parent.mkdir(parents=True,exist_ok=True)
        injected.write_text('raise RuntimeError("'+marker+'")\n')
        result=subprocess.run([sys.executable,'-I','-B',str(execution/entry),'--help'],capture_output=True,text=True)
        log=result.stdout+result.stderr;(root/'CHILD.log').write_text(log)
        cases.append({'entry':entry,'alias':alias,
                      'entry_sha256':hashlib.sha256((execution/entry).read_bytes()).hexdigest(),
                      'unregistered_alias_executed':marker in log,'exit_code':result.returncode,
                      'child_log_sha256':hashlib.sha256(log.encode()).hexdigest()})
    receipt={'status':'BOOTSTRAP_PURE_PYTHON_ALIAS_CONTROLS','cases':cases,'retained_decodes':0,'physical_source_calls':0}
    with (destination/'RESULT.json').open('x') as out:json.dump(receipt,out,sort_keys=True,indent=2);out.write('\n')
    print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':run(*sys.argv[1:])
