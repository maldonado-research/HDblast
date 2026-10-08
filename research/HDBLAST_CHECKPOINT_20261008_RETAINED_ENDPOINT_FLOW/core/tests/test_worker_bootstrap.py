"""Harmless source/package/cache/stdlib bootstrap markers must never execute."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import types
sys.dont_write_bytecode=True
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'execution'))
from registration_guard import source_tree
sys.path.insert(0,str(root/'tests'))
from test_registration_guard import CONTRACT

def run():
    passed=[]
    for mutation in ('registration package','kernel package','ctypes sibling','resource sibling',
                     'guard byte changed','cached guard','nonisolated json sibling','no actual GO',
                     'unreadable source package','unreadable registered directory','empty unregistered directory'):
        with tempfile.TemporaryDirectory(prefix='endpoint-worker-bootstrap-',dir='/tmp') as temporary:
            temp=Path(temporary);candidate=temp/'candidate'
            shutil.copytree(root,candidate)
            files={}
            for relative in source_tree(candidate):
                raw=(candidate/relative).read_bytes();files[relative]={'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
            registration={'schema_version':1,'scope':'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE','files':files,'contract':CONTRACT}
            raw=(json.dumps(registration,sort_keys=True,separators=(',',':'))+'\n').encode()
            reg=temp/'registration.json';reg.write_bytes(raw);digest=hashlib.sha256(raw).hexdigest()
            marker=temp/'UNAUTHORIZED_MARKER'
            code=f'from pathlib import Path\nPath({str(marker)!r}).write_text("EXECUTED")\nraise RuntimeError("UNAUTHORIZED_MARKER")\n'
            if mutation in ('registration package','kernel package'):
                directory=candidate/'execution'/('registration_guard' if mutation.startswith('registration') else 'kernel_guard')
                directory.mkdir();(directory/'__init__.py').write_text(code)
            elif mutation.endswith('sibling'):
                name=mutation.split()[0];(candidate/'execution'/('json.py' if name=='nonisolated' else name+'.py')).write_text(code)
            elif mutation=='guard byte changed':
                (candidate/'execution/registration_guard.py').write_text(code)
            elif mutation=='cached guard':
                cache=candidate/'execution/__pycache__';cache.mkdir();(cache/'registration_guard.cpython-312.pyc').write_bytes(b'fabricated cached executable')
            elif mutation=='unreadable source package':
                directory=candidate/'source/later_source';directory.mkdir();(directory/'__init__.py').write_text(code)
                directory.chmod(0o111)
            elif mutation=='unreadable registered directory':
                directory=candidate/'engine';directory.chmod(0o111)
            elif mutation=='empty unregistered directory':
                (candidate/'unregistered_empty').mkdir()
            command=[sys.executable]+([] if mutation.startswith('nonisolated') else ['-I'])+['-B']
            if sys.flags.optimize:command.append('-O')
            command.extend([str(candidate/'execution/run_endpoint.py'),'--root',str(candidate),
                '--registration',str(reg),'--registration-sha256',digest,'--output',str(temp/'output')])
            if mutation!='no actual GO':command.append('--manufactured-only')
            else:command.extend(['--repository-root',str(temp/'not-a-repository')])
            try:result=subprocess.run(command,capture_output=True,text=True,timeout=10)
            finally:
                if mutation in ('unreadable source package','unreadable registered directory'):
                    directory.chmod(0o755)
            if result.returncode==0 or marker.exists() or 'UNAUTHORIZED_MARKER' in result.stdout:
                raise RuntimeError('bootstrap marker executed or control accepted: '+mutation)
            expected=('unregistered source directory' if mutation in ('registration package','kernel package',
                       'unreadable source package','empty unregistered directory') else
                      'source directory inspection failed' if mutation=='unreadable registered directory' else
                      'package aliases' if mutation.endswith('package') else
                      'standard-library aliases' if mutation.endswith('sibling') and not mutation.startswith('nonisolated') else
                      'source byte pin mismatch' if mutation=='guard byte changed' else
                      'cached/native executable aliases' if mutation=='cached guard' else
                      'isolated Python -I required' if mutation.startswith('nonisolated') else 'immutable public byte GO')
            if expected not in result.stderr:raise RuntimeError('wrong rejection reason: '+mutation+' '+result.stderr)
            if (temp/'output').exists():raise RuntimeError('bootstrap rejection created worker output')
            passed.append(mutation)
    # Captured bytes and module identities remain effective after the named
    # source is replaced. This harmless dependency chain uses dataclasses and
    # proves import search and the changed filesystem are never consulted.
    entry=types.ModuleType('manufactured_entry_functions')
    entry.__file__=str(root/'execution/run_endpoint.py')
    exec(compile((root/'execution/run_endpoint.py').read_bytes(),entry.__file__,'exec'),entry.__dict__)
    with tempfile.TemporaryDirectory(prefix='endpoint-captured-module-',dir='/tmp') as temporary:
        temp=Path(temporary);path=temp/'fixture.py';marker=temp/'UNAUTHORIZED_MARKER'
        trusted=b'from dataclasses import dataclass\n@dataclass\nclass ExactFixture:\n    value:int=37\n'
        path.write_bytes(trusted);pin={'bytes':len(trusted),'sha256':hashlib.sha256(trusted).hexdigest()}
        captured=entry.capture(path,1<<20,pin)
        path.write_text('from pathlib import Path\nPath('+repr(str(marker))+').write_text("EXECUTED")\n')
        before=tuple(sys.path)
        module=entry.source_module('manufactured_captured_dependency',path,captured)
        if module.ExactFixture().value!=37 or marker.exists() or tuple(sys.path)!=before:
            raise RuntimeError('captured source consulted changed filesystem or import search')
        passed.append('captured bytes survive source replacement without path exposure')
        consumer=entry.source_module('manufactured_captured_consumer',temp/'consumer.py',
            b'from manufactured_captured_dependency import ExactFixture\nresult=ExactFixture().value+5\n')
        if consumer.result!=42 or module.__file__!=str(path) or marker.exists():
            raise RuntimeError('captured dependency identity or file provenance changed')
        passed.append('captured dataclass dependency identity and file provenance')
        try:entry.capture(path,1<<20,pin)
        except ValueError:passed.append('changed named source rejected on later pin capture')
        else:raise RuntimeError('changed named source accepted with old byte pin')
        sys.modules.pop('manufactured_captured_dependency');sys.modules.pop('manufactured_captured_consumer')
    print(json.dumps({'status':'PASS_SOURCE_ONLY_WORKER_BOOTSTRAP','checks':len(passed),'cases':passed,
                      'physical_callbacks':0,'retained_array_reads':0,'untrusted_markers_executed':0},sort_keys=True))

if __name__=='__main__':run()
