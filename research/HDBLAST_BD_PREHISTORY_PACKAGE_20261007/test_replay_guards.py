"""Manufactured file/JSON guard checks only; never launches a source route."""
from pathlib import Path
import hashlib
import json
import sys
import tempfile
import copy
import replay_checkpoint as r

checks=[];controls=[]
def check(value,label):
    if not value: raise RuntimeError(label)
    checks.append(label)
def reject(fn,label):
    try: fn()
    except (ValueError,TypeError): controls.append(label)
    else: raise RuntimeError('Unsafe fixture accepted: '+label)
with tempfile.TemporaryDirectory(prefix='portable-replay-manufactured-') as tmp:
    tmp=Path(tmp);root=tmp/'checkpoint';root.mkdir()
    (root/'FULL_REGISTRATION.json').write_text('{"schema_version":1,"files":{}}\n')
    (root/'note.txt').write_text('manufactured immutable text\n')
    files={p.name:{'bytes':p.stat().st_size,'sha256':r.sha_file(p)} for p in root.iterdir()}
    manifest={'schema_version':1,'files':files}
    check(r.verify_payload(root,manifest)==2,'complete manufactured inventory')
    raw=(json.dumps(manifest,sort_keys=True)+'\n').encode();path=tmp/'manifest.json';path.write_bytes(raw)
    check(r.pinned_json(path,r.sha_bytes(raw))[0]==manifest,'external expected JSON SHA binds bytes')
    reject(lambda:r.pinned_json(path,'0'*64),'incorrect external expected SHA')
    reject(lambda:r.pinned_json(path,r.sha_bytes(raw).upper()),'noncanonical SHA encoding')
    reject(lambda:r.load(b'{"x":0,"x":1}'),'duplicate JSON keys')
    reject(lambda:r.load(b'{"x":NaN}'),'nonfinite JSON')
    reject(lambda:r.safe_file(root,'../manifest.json'),'relative traversal')
    reject(lambda:r.safe_file(root,'/tmp/manifest.json'),'absolute path')
    reject(lambda:r.safe_file(root,'note.txt/../note.txt'),'interior traversal')
    reject(lambda:r.safe_file(root,'./note.txt'),'dot component')
    (root/'extra.py').write_text('raise RuntimeError("must never be imported")\n')
    reject(lambda:r.verify_payload(root,manifest),'unlisted executable source')
    (root/'extra.py').unlink()
    (root/'note.txt').write_text('changed')
    reject(lambda:r.verify_payload(root,manifest),'changed payload bytes')
    (root/'note.txt').write_text('manufactured immutable text\n')
    link=root/'link';link.symlink_to(root/'note.txt')
    reject(lambda:r.safe_file(root,'link'),'symlink payload')
    reject(lambda:r.verify_payload(root,manifest),'extra symlink leaf')
    link.unlink()
    jsonlink=tmp/'jsonlink';jsonlink.symlink_to(path)
    reject(lambda:r.pinned_json(jsonlink,r.sha_bytes(raw)),'symlink pinned JSON')
    bad={'schema_version':True,'files':files}
    reject(lambda:r.verify_payload(root,bad),'boolean schema version')
    check(r.python_command('/usr/bin/python',False,root/'script.py',['x'])==['/usr/bin/python','-B',str(root/'script.py'),'x'],
          'normal command exact flags')
    check(r.python_command('/usr/bin/python',True,root/'script.py',['x'])==['/usr/bin/python','-B','-O',str(root/'script.py'),'x'],
          'optimized command exact flags')
    reject(lambda:r.compare_bytes(root/'note.txt',path,'fixture'),'exact scientific mismatch')
    check(r.compare_bytes(root/'note.txt',root/'note.txt','fixture')==r.sha_file(root/'note.txt'),'exact scientific identity')
    package=tmp/'package';package.mkdir()
    check(r.fresh_output_directory(tmp/'fresh',(root,package))==tmp/'fresh','fresh output outside immutable sibling directories')
    reject(lambda:r.fresh_output_directory(root/'fresh',(root,package)),'output inside checkpoint')
    reject(lambda:r.fresh_output_directory(package/'fresh',(root,package)),'output inside delivery package')
    reject(lambda:r.fresh_output_directory(root,(root,package)),'existing output directory')
    outlink=tmp/'outlink';outlink.symlink_to(package,target_is_directory=True)
    reject(lambda:r.fresh_output_directory(outlink/'fresh',(root,package)),'output symlink ancestor')
    for route in ('primary','independent'):
        resource={'wall_seconds':1.5,'peak_rss_kib':1024}
        first={'source_degree':1,'source_model_rows':[{'coefficients':['1/1','2/1']}],**(resource if route=='primary' else {'resource':resource})}
        second=copy.deepcopy(first)
        measures=second if route=='primary' else second['resource']
        measures.update(wall_seconds=2.5,peak_rss_kib=2048)
        left=tmp/(route+'-budget-left.json');right=tmp/(route+'-budget-right.json')
        left.write_text(json.dumps(first));right.write_text(json.dumps(second))
        check(bool(r.compare_budget_cores(left,right,route)),'scientific budget ignores only resource measurements '+route)
        second['source_model_rows'][0]['coefficients'][1]='3/1';right.write_text(json.dumps(second))
        reject(lambda:r.compare_budget_cores(left,right,route),'changed scientific budget coefficient '+route)
        for changed,label in ((1.0,'float'),(True,'boolean')):
            second=copy.deepcopy(first);second['source_degree']=changed;right.write_text(json.dumps(second))
            reject(lambda:r.compare_budget_cores(left,right,route),'canonical JSON scientific type change '+route+' '+label)
        if route=='independent':
            second=copy.deepcopy(first);second['resource']['unregistered']=0;right.write_text(json.dumps(second))
            reject(lambda:r.scientific_budget_core(second,route),'extra variable resource field')
result={'status':'PASS_PORTABLE_REPLAY_MANUFACTURED_GUARDS','checks_passed':len(checks),
        'rejected_controls_count':len(controls),'checks':checks,'rejected_controls':controls,
        'physical_source_callbacks':0,'route_launches':0,'saved_array_decodes':0,
        'helper_sha256':hashlib.sha256(Path(r.__file__).read_bytes()).hexdigest(),
        'test_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('MANUFACTURED_GUARDS.json')
with out.open('x') as f:f.write(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({key:result[key] for key in ['status','checks_passed','rejected_controls_count','physical_source_callbacks','route_launches']}))
