"""Stdlib-only candidate closure/GO guard. Numeric modules import afterward."""
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import sys
import time
import resource

REQUIRED_FILES = {
    'execution/registration_guard.py','execution/run_endpoint.py','execution/endpoint_worker.py',
    'execution/kernel_guard.py','execution/bounded_launcher.py','execution/validate_outputs.py',
    'engine/endpoint_engine.py','engine/endpoint_aggregate.py','source/later_source.py',
    'source/source_algebra_baseline.py','provenance_decoder/exact_binary80.py',
    'provenance_decoder/fabricated_fixtures.py','provenance_decoder/input_spec_adapter.py',
    'provenance_decoder/INPUT_SPEC.json','provenance_decoder/UPSTREAM_HASH_RECEIPT.json',
    'provenance_decoder/LAYOUT_PROVENANCE_BASIS.json','provenance_decoder/ORIGINAL_ZIP_INVENTORIES.json',
    'PROTOCOL.md',
    'theory/SOURCE_AND_ENDPOINT_ENGINE_PROOF.md','theory/LATER_TRAJECTORY_THEOREM.md',
    'theory/MATHEMATICAL_REVIEW_RECEIPT.json',
    'theory/DYADIC_SOURCE_ERROR_UPDATE.md','theory/DYADIC_SOURCE_ERROR_UPDATE_RECEIPT.json',
}

def normal_relative(relative):
    require(type(relative) is str and relative and '\\' not in relative,'safe relative path required')
    path=Path(relative)
    require(not path.is_absolute() and '..' not in path.parts and
            relative==path.as_posix() and relative!='.','normalized relative path required')
    return path

def directory_signature(info):
    return (info.st_dev,info.st_ino,info.st_mode,info.st_size,info.st_mtime_ns,info.st_ctime_ns)

def source_tree(root,files=None):
    # Descriptor-relative traversal never suppresses a directory enumeration
    # error and never follows a link or a parent-entry replacement. Registered
    # file paths determine the exact directory roster, including no empty dirs.
    root=Path(root).absolute();found=set();directories=set()
    expected=None
    if files is not None:
        expected={str(parent) for relative in files for parent in normal_relative(relative).parents
                  if str(parent)!='.'}
    flags=os.O_RDONLY|os.O_DIRECTORY|os.O_NOFOLLOW|os.O_CLOEXEC
    def visit(fd,relative,before):
        require(directory_signature(os.fstat(fd))==directory_signature(before),'source directory changed before scan')
        with os.scandir(fd) as entries:
            for entry in entries:
                name=entry.name;path=relative/name;key=path.as_posix()
                normal_relative(key)
                info=os.stat(name,dir_fd=fd,follow_symlinks=False)
                require(not stat.S_ISLNK(info.st_mode),'symlink in candidate source tree')
                require(name!='__init__.py','candidate package aliases forbidden')
                require(name!='__pycache__' and path.suffix not in ('.pyc','.pyo','.so','.pyd','.dll','.dylib'),
                        'candidate cached/native executable aliases forbidden')
                require(not (path.suffix=='.py' and path.stem in sys.stdlib_module_names) and
                        not (stat.S_ISDIR(info.st_mode) and name in sys.stdlib_module_names),
                        'candidate standard-library aliases forbidden')
                if stat.S_ISREG(info.st_mode):
                    require(info.st_nlink==1,'source hardlink aliases forbidden')
                    found.add(key)
                else:
                    require(stat.S_ISDIR(info.st_mode),'special file in candidate source tree')
                    require(expected is None or key in expected,'unregistered source directory: '+key)
                    directories.add(key)
                    child=os.open(name,flags,dir_fd=fd)
                    try:visit(child,path,info)
                    finally:os.close(child)
                require(directory_signature(os.stat(name,dir_fd=fd,follow_symlinks=False))==
                        directory_signature(info),'source entry changed during scan: '+key)
        require(directory_signature(os.fstat(fd))==directory_signature(before),'source directory changed during scan')
    try:
        before=os.stat(root,follow_symlinks=False)
        require(stat.S_ISDIR(before.st_mode),'real source directory required')
        fd=os.open(root,flags)
        try:visit(fd,Path(),before)
        finally:os.close(fd)
        require(directory_signature(os.stat(root,follow_symlinks=False))==directory_signature(before),
                'source root changed during scan')
    except OSError as error:
        raise ValueError('source directory inspection failed: '+str(error)) from error
    require(expected is None or directories==expected,'registered source directory roster differs')
    return found

def require(ok,text):
    if not ok:raise ValueError(text)

def sha(raw):return hashlib.sha256(raw).hexdigest()

def load(raw):
    def pairs(items):
        d={}
        for k,v in items:
            require(k not in d,'duplicate JSON key');d[k]=v
        return d
    return json.loads(raw,object_pairs_hook=pairs,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('nonfinite JSON')))

def safe_file(root,relative):
    p=normal_relative(relative)
    out=root
    for component in p.parts:
        out=out/component
        require(not out.is_symlink(),'symlink forbidden')
    require(out.is_file() and stat.S_ISREG(out.stat().st_mode),'regular file required')
    require(out.resolve().is_relative_to(root.resolve()),'path escapes root')
    return out

def authenticate(root,registration_path,registration_sha,manufactured,go_path=None,go_sha=None):
    require(type(manufactured) is bool,'explicit manufactured branch required')
    root=Path(root).absolute()
    require(root.is_dir() and not root.is_symlink(),'real candidate directory required')
    raw=Path(registration_path).read_bytes()
    require(re.fullmatch('[0-9a-f]{64}',registration_sha or '') and sha(raw)==registration_sha,'registration SHA mismatch')
    registration=load(raw)
    require(registration['schema_version']==1 and registration['scope']=='LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE',
            'registration scope/schema differs')
    files=registration['files'];require(type(files) is dict and REQUIRED_FILES<=set(files),
                                     'mandatory complete semantic source closure required')
    for relative,pin in files.items():
        data=safe_file(root,relative).read_bytes()
        require(len(data)==pin['bytes'] and sha(data)==pin['sha256'],'frozen source byte mismatch: '+relative)
    require(source_tree(root,files)==set(files),'unregistered or omitted noncache source artifact')
    contract=registration['contract']
    require(contract=={'arithmetic_bits':1024,'momentum_degree':2048,'source_degree':24,
        'source_coefficient_bits':512,'source_cells':64,'export_bits':96,'export_gate':'1/1000000000000000000',
        'anchor':'-9/2','endpoints':['-4','-7/2'],'selected_decode_members':['k.npy','momentum_weights.npy',
        'observation_eta.npy','u_1.npy','w_1.npy','u_2.npy','w_2.npy','u_3.npy','w_3.npy'],
        'wall_seconds':900,'rss_kib':524288,'output_bytes':134217728},'exact registered contract differs')
    go=None
    if not manufactured:
        require(go_path is not None and go_sha is not None,'separate immutable public byte GO required')
        goraw=Path(go_path).read_bytes()
        require(re.fullmatch('[0-9a-f]{64}',go_sha) and sha(goraw)==go_sha,'GO receipt SHA mismatch')
        go=load(goraw)
        require(go['status']=='PASS_READBACK_PUBLIC_BYTE_GO' and
                re.fullmatch('[0-9a-f]{40}',go['freeze_commit']) and go['registration_sha256']==registration_sha,
                'public byte GO status/commit/registration differs')
        require(go['scope']=='LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE' and
                go['source_go'] is True and go['decode_go'] is True and
                go['public_bytes_verified'] is True and go['independent_review_pass'] is True,
                'complete explicit reviewed source/decode GO required')
        require(go['registered_file_pins']==files,'public GO exact byte pins differ')
    return {'root':root,'registration':registration,'registration_path':Path(registration_path),
            'registration_sha256':registration_sha,'manufactured':manufactured,'go':go,'go_sha256':go_sha}

def reauthenticate(v):
    require(sha(v['registration_path'].read_bytes())==v['registration_sha256'],'registration changed')
    require(source_tree(v['root'],v['registration']['files'])==set(v['registration']['files']),'source closure changed')
    for relative,pin in v['registration']['files'].items():
        data=safe_file(v['root'],relative).read_bytes()
        require(len(data)==pin['bytes'] and sha(data)==pin['sha256'],'frozen source changed')

def resource_guard(started,output):
    elapsed=time.monotonic()-started;rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    require(elapsed<=900 and rss<=524288,'registered wall/RSS budget exceeded')
    size=sum(p.stat().st_size for p in Path(output).rglob('*') if p.is_file())
    require(size<=134217728,'registered output budget exceeded')
    return {'wall_seconds':elapsed,'peak_rss_kib':rss,'output_bytes':size}
