"""Manufactured semantic-closure/GO adversaries; no numeric imports or inputs."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import tempfile
from unittest.mock import patch
sys.dont_write_bytecode=True
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'execution'))
from registration_guard import authenticate,source_tree

CONTRACT={'arithmetic_bits':1024,'momentum_degree':2048,'source_degree':24,'source_coefficient_bits':512,
    'source_cells':64,'export_bits':96,'export_gate':'1/1000000000000000000','anchor':'-9/2','endpoints':['-4','-7/2'],
    'selected_decode_members':['k.npy','momentum_weights.npy','observation_eta.npy','u_1.npy','w_1.npy',
                               'u_2.npy','w_2.npy','u_3.npy','w_3.npy'],
    'wall_seconds':900,'rss_kib':524288,'output_bytes':134217728}

def document(candidate):
    files={}
    for relative in source_tree(candidate):
        raw=(candidate/relative).read_bytes();files[relative]={'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
    return {'schema_version':1,'scope':'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE','files':files,'contract':CONTRACT}

def execute(candidate,d,registration,manufactured=True,go=None):
    raw=(json.dumps(d,sort_keys=True,separators=(',',':'))+'\n').encode();registration.write_bytes(raw)
    digest=hashlib.sha256(raw).hexdigest()
    if go is None:return authenticate(candidate,registration,digest,manufactured)
    go['registration_sha256']=digest;go['registered_file_pins']=d['files']
    gp=registration.with_suffix('.go.json');gr=(json.dumps(go,sort_keys=True)+'\n').encode();gp.write_bytes(gr)
    return authenticate(candidate,registration,digest,manufactured,gp,hashlib.sha256(gr).hexdigest())

def run():
    passed=[]
    with tempfile.TemporaryDirectory(prefix='endpoint-guard-',dir='/tmp') as temporary:
        temp=Path(temporary);candidate=temp/'candidate'
        shutil.copytree(root,candidate,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        reg=temp/'registration.json';d=document(candidate)
        execute(candidate,d,reg);passed.append('complete manufactured source closure')
        def reject(label,action):
            try:action()
            except (ValueError,KeyError):passed.append('reject '+label)
            else:raise RuntimeError('mutation accepted: '+label)
        omitted=json.loads(json.dumps(d));omitted['files'].pop('provenance_decoder/INPUT_SPEC.json')
        reject('omitted semantic InputSpec',lambda:execute(candidate,omitted,reg))
        helper=candidate/'unregistered.txt';helper.write_text('hidden semantic helper')
        reject('unregistered non-Python source artifact',lambda:execute(candidate,d,reg));helper.unlink()
        helper=candidate/'unregistered.py';helper.write_text('raise RuntimeError("must not import")')
        reject('unregistered Python helper',lambda:execute(candidate,d,reg));helper.unlink()
        directory=candidate/'empty_unregistered';directory.mkdir()
        reject('unregistered empty directory roster',lambda:execute(candidate,d,reg));directory.rmdir()
        directory=candidate/'source/later_source';directory.mkdir()
        (directory/'__init__.py').write_text('raise RuntimeError("UNREGISTERED_PACKAGE")')
        directory.chmod(0o111)
        try:
            reject('unreadable package directory',lambda:execute(candidate,d,reg))
            reject('unreadable standalone source traversal',lambda:source_tree(candidate))
        finally:directory.chmod(0o755);shutil.rmtree(directory)
        directory=candidate/'engine';old_mode=directory.stat().st_mode & 0o777
        directory.chmod(0o111)
        try:reject('unreadable registered directory',lambda:execute(candidate,d,reg))
        finally:directory.chmod(old_mode)
        with patch('registration_guard.os.scandir',side_effect=PermissionError('injected enumeration failure')):
            reject('propagated directory enumeration failure',lambda:source_tree(candidate))
        original_open=__import__('os').open;stashed=temp/'raced_engine'
        def replaced_open(path,flags,*args,**kwargs):
            if path=='engine' and 'dir_fd' in kwargs:
                (candidate/'engine').rename(stashed)
                (candidate/'engine').symlink_to(stashed,target_is_directory=True)
            return original_open(path,flags,*args,**kwargs)
        try:
            with patch('registration_guard.os.open',side_effect=replaced_open):
                reject('directory replaced by symlink between stat and open',lambda:execute(candidate,d,reg))
        finally:
            if (candidate/'engine').is_symlink():(candidate/'engine').unlink()
            if stashed.exists():stashed.rename(candidate/'engine')
        malformed=json.loads(json.dumps(d))
        malformed['files']['engine//endpoint_engine.py']=malformed['files'].pop('engine/endpoint_engine.py')
        reject('nonnormal registered file path',lambda:execute(candidate,malformed,reg))
        old=(candidate/'engine/endpoint_engine.py').read_bytes();(candidate/'engine/endpoint_engine.py').write_bytes(old+b'\n')
        reject('registered source byte changed',lambda:execute(candidate,d,reg));(candidate/'engine/endpoint_engine.py').write_bytes(old)
        cache=candidate/'engine/__pycache__';cache.mkdir();(cache/'endpoint_engine.cpython-312.pyc').write_bytes(b'fabricated cached executable')
        reject('unregistered cached executable',lambda:execute(candidate,d,reg));shutil.rmtree(cache)
        native=candidate/'execution/registration_guard.so';native.write_bytes(b'fabricated native alias')
        reject('native local executable alias',lambda:execute(candidate,d,reg));native.unlink()
        package=candidate/'execution/registration_guard';package.mkdir();(package/'__init__.py').write_text('raise RuntimeError("package alias")')
        reject('pure Python package alias',lambda:execute(candidate,d,reg));shutil.rmtree(package)
        alias=candidate/'execution/ctypes.py';alias.write_text('raise RuntimeError("stdlib alias")')
        reject('pure Python standard-library alias',lambda:execute(candidate,d,reg));alias.unlink()
        empty=temp/'empty';empty.mkdir();(empty/'irrelevant.py').write_text('# irrelevant')
        reject('fabricated subset root',lambda:execute(empty,document(empty),reg))
        reject('actual branch without immutable GO',lambda:execute(candidate,d,reg,False))
        go={'status':'PASS_READBACK_PUBLIC_BYTE_GO','freeze_commit':'a'*40,
            'scope':'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE','source_go':False,'decode_go':True,
            'public_bytes_verified':True,'independent_review_pass':True}
        reject('source authorization absent',lambda:execute(candidate,d,reg,False,go))
        go['source_go']=True;go['decode_go']=False
        reject('decode authorization absent',lambda:execute(candidate,d,reg,False,go))
        go['decode_go']=True;go['scope']='NO_RETAINED_ARRAYS_SOURCE_OPERATOR_ONLY'
        reject('legacy GO scope',lambda:execute(candidate,d,reg,False,go))
        go['scope']='LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE';go['independent_review_pass']=False
        reject('independent review absent',lambda:execute(candidate,d,reg,False,go))
        path=candidate/'PROTOCOL.md';old=path.read_bytes();path.unlink();path.symlink_to(root/'PROTOCOL.md')
        reject('registered payload symlink',lambda:execute(candidate,d,reg));path.unlink();path.write_bytes(old)
    if 'flint' in sys.modules:raise RuntimeError('guard imported numeric dependencies')
    print(json.dumps({'status':'PASS_PRE_NUMERIC_SEMANTIC_AND_GO_GUARDS','checks':len(passed),'cases':passed,
                      'retained_arrays_opened':0,'physical_callbacks':0,'flint_imported':False},sort_keys=True))

if __name__=='__main__':run()
