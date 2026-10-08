"""Isolated source-only bootstrap authenticates before local import search."""
import sys
if not sys.flags.isolated:
    raise ValueError('isolated Python -I required before bootstrap imports')
sys.dont_write_bytecode=True

# Candidate directories remain outside import search while trusted standard
# library modules and captured guard/kernel sources are loaded.
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import types

def require(ok,text):
    if not ok:raise ValueError(text)

def unique(items):
    result={}
    for key,value in items:
        require(key not in result,'duplicate bootstrap JSON key');result[key]=value
    return result

def capture(path,limit,pin=None):
    with os.fdopen(os.open(path,os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK),'rb') as stream:
        info=os.fstat(stream.fileno())
        require(stat.S_ISREG(info.st_mode) and 0<info.st_size<=limit,'bounded regular bootstrap file required')
        raw=stream.read(limit+1)
    require(len(raw)==info.st_size and len(raw)<=limit,'bootstrap file changed or grew')
    if pin is not None:
        require(type(pin) is dict and type(pin['bytes']) is int and len(raw)==pin['bytes'] and
                hashlib.sha256(raw).hexdigest()==pin['sha256'],'captured source byte pin mismatch')
    return raw

def source_module(name,path,raw):
    module=types.ModuleType(name);module.__file__=str(path);module.__package__=''
    sys.modules[name]=module
    exec(compile(raw,str(path),'exec'),module.__dict__)
    return module

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--registration',type=Path,required=True)
    p.add_argument('--registration-sha256',required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--manufactured-only',action='store_true')
    p.add_argument('--repository-root',type=Path)
    p.add_argument('--go',type=Path)
    p.add_argument('--go-sha256')
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    require(args.root.absolute()==root,'entry must execute the declared candidate source root')
    raw=capture(args.registration,8<<20)
    require(re.fullmatch('[0-9a-f]{64}',args.registration_sha256) and
            hashlib.sha256(raw).hexdigest()==args.registration_sha256,'bootstrap registration SHA mismatch')
    registration=json.loads(raw,object_pairs_hook=unique,
                            parse_constant=lambda _:(_ for _ in ()).throw(ValueError('nonfinite bootstrap JSON')))
    guard_path=root/'execution/registration_guard.py'
    guard=source_module('registration_guard',guard_path,
                        capture(guard_path,1<<20,registration['files']['execution/registration_guard.py']))
    v=guard.authenticate(root,args.registration,args.registration_sha256,args.manufactured_only,args.go,args.go_sha256)
    if args.manufactured_only and args.repository_root is not None:
        raise ValueError('manufactured branch cannot receive retained repository root')
    if not args.manufactured_only and args.repository_root is None:
        raise ValueError('actual branch requires exact retained repository root')
    kernel_path=guard.safe_file(root,'execution/kernel_guard.py')
    kernel=source_module('kernel_guard',kernel_path,
                         capture(kernel_path,1<<20,registration['files']['execution/kernel_guard.py']))
    kernel.require_worker_limits()
    # Capture every consumed local module before executing any of them. Local
    # directories remain outside import search for the entire worker lifetime;
    # imports between these modules resolve only their captured sys.modules IDs.
    modules=(('exact_binary80','provenance_decoder/exact_binary80.py'),
             ('fabricated_fixtures','provenance_decoder/fabricated_fixtures.py'),
             ('input_spec_adapter','provenance_decoder/input_spec_adapter.py'),
             ('endpoint_engine','engine/endpoint_engine.py'),
             ('endpoint_aggregate','engine/endpoint_aggregate.py'),
             ('source_algebra_baseline','source/source_algebra_baseline.py'),
             ('later_source','source/later_source.py'),
             ('validate_outputs','execution/validate_outputs.py'),
             ('endpoint_worker','execution/endpoint_worker.py'))
    captured=[(name,guard.safe_file(root,relative),
               capture(guard.safe_file(root,relative),1<<20,registration['files'][relative]))
              for name,relative in modules]
    trusted_path=tuple(sys.path)
    require(all(not Path(path).absolute().is_relative_to(root) for path in trusted_path),
            'candidate must stay outside trusted import search')
    for name,path,data in captured:
        require(name not in sys.modules,'unregistered prior local module identity')
        source_module(name,path,data)
        require(tuple(sys.path)==trusted_path,'registered module changed trusted import search')
    import flint
    if sys.version_info[:3]!=(3,12,14) or flint.__version__!='0.9.0':
        raise ValueError('pinned Python3.12.14/python-flint0.9.0 required')
    v['kernel_denied_syscalls']=kernel.install()
    guard.reauthenticate(v)
    try:
        sys.modules['endpoint_worker'].execute(v,args.output,args.repository_root)
    except Exception as error:
        if args.output.is_dir():
            failure=args.output/'FAILURE.json'
            if not failure.exists():
                with failure.open('x') as out:
                    json.dump({'status':'FAIL_CLOSED','error_type':type(error).__name__,'message':str(error)},out,sort_keys=True)
                    out.write('\n')
        raise

if __name__=='__main__':main()
