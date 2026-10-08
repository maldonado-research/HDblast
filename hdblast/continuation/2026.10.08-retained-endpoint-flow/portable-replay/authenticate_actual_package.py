"""Independent data-only authentication; original dependencies are opaque bytes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import stat

HERE=Path(__file__).resolve().parent
EXTERNAL={'replay_checkpoint.py':'5de25511aea5579b7e0ebe94208ee1584d0869ca248fa8939cdbb29ed25ff6fd',
'PAYLOAD_MANIFEST.json':'976e92a7a14ecb13ed04cbb9f7ed15fff3a285bd40175a28929dff53677d0752',
'REPLAY_PINS.json':'1407dc4bcf0a312d1bc23a29e93e3c1149071f4ce2705ccbbeef430e82446d68'}
REG='33eded898f47a96142f85002b896edd4c294180653f50f2b38e83f22d5a75796'
GO='edded4f67840bc54e95673c5f43b71f6804e373d13211c90ba48ca99f72a8389'
BOOT='438a0478fa185f46ca06372036caa7aee1f9e9659a2778b9d1a9e98651122f0a'
ACTUAL_REVIEW='097fe778797dc78f1dc800885bd064e7406c6eb80f9642b3c7cc8577d05754ca'
PREFIX='research/HDBLAST_CHECKPOINT_20261008_RETAINED_ENDPOINT_FLOW'

def require(ok,message):
 if not ok:raise ValueError(message)
def pairs(items):
 result={}
 for key,value in items:
  require(key not in result,'duplicate JSON key');result[key]=value
 return result
def parse(raw):return json.loads(raw,object_pairs_hook=pairs,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('nonfinite JSON')))
def signature(s):return (s.st_dev,s.st_ino,s.st_mode,s.st_size,s.st_mtime_ns,s.st_ctime_ns,s.st_nlink)
def relative(name):
 p=Path(name);require(type(name) is str and name and not p.is_absolute() and '..' not in p.parts and '\\' not in name and p.as_posix()==name and name!='.','safe normalized relative path')
 return p

def read_member(root,name,capture=False):
 parts=relative(name).parts;root=Path(root).absolute()
 for p in (root,*root.parents):require(stat.S_ISDIR(p.lstat().st_mode),'real root ancestry')
 flags=os.O_RDONLY|os.O_DIRECTORY|os.O_NOFOLLOW|os.O_CLOEXEC
 descriptors=[os.open(root,flags)]
 try:
  for part in parts[:-1]:descriptors.append(os.open(part,flags,dir_fd=descriptors[-1]))
  fd=os.open(parts[-1],os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK|os.O_CLOEXEC,dir_fd=descriptors[-1])
  try:
   before=os.fstat(fd);require(stat.S_ISREG(before.st_mode) and before.st_nlink==1 and before.st_size<=256*(1<<20),'bounded unaliased regular leaf')
   h=hashlib.sha256();count=0;blocks=[]
   while chunk:=os.read(fd,1<<20):
    count+=len(chunk);require(count<=256*(1<<20),'file grew');h.update(chunk)
    if capture:blocks.append(chunk)
   require(signature(before)==signature(os.fstat(fd))==signature(os.stat(parts[-1],dir_fd=descriptors[-1],follow_symlinks=False)) and count==before.st_size,'leaf changed during hashing')
   return {'bytes':count,'sha256':h.hexdigest()},b''.join(blocks) if capture else None
  finally:os.close(fd)
 finally:
  for descriptor in reversed(descriptors):os.close(descriptor)

def inventory(root):
 found=set();dirs=set()
 def error(exc):raise exc
 for parent,names,leaves in os.walk(root,followlinks=False,onerror=error):
  for name in names+leaves:
   path=Path(parent)/name;s=path.lstat();key=path.relative_to(root).as_posix();relative(key)
   require(not stat.S_ISLNK(s.st_mode),'package link')
   if stat.S_ISDIR(s.st_mode):dirs.add(key)
   else:
    require(stat.S_ISREG(s.st_mode) and s.st_nlink==1,'package special or hardlink')
    require('__pycache__' not in path.parts and path.suffix not in ('.pyc','.pyo','.so','.pyd','.dll','.dylib'),'package executable alias')
    found.add(key)
 return found,dirs

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--package',type=Path,required=True);parser.add_argument('--repository-root',type=Path,required=True);parser.add_argument('--receipt',type=Path,required=True);args=parser.parse_args()
 root=args.package.absolute();controls={};decoded={}
 for name,digest in EXTERNAL.items():
  item,raw=read_member(root,name,True);require(item['sha256']==digest,'external control pin:'+name);controls[name]=item
  if name.endswith('.json'):decoded[name]=parse(raw)
 manifest=decoded['PAYLOAD_MANIFEST.json'];pins=decoded['REPLAY_PINS.json'];files=manifest['files']
 require(manifest['schema_version']==1 and type(files) is dict and not set(files)&set(EXTERNAL),'payload manifest schema')
 names,dirs=inventory(root);require(names==set(files)|set(EXTERNAL),'complete package file roster')
 expected_dirs={p.as_posix() for name in names for p in relative(name).parents if p.as_posix()!='.'}
 require(dirs==expected_dirs,'complete package directory roster')
 actual={}
 for name,expected in files.items():
  item,_=read_member(root,name);require(item==expected,'payload file pin:'+name);actual[name]=item
 require(len(actual)==100 and len(names)==103,'sealed exact package count')
 require(pins['schema_version']==1 and pins['sealed'] is True and pins['execution_scope']=='ACTUAL_RETAINED_ENDPOINT_FLOW','sealed actual scope')
 require(pins['helper_sha256']==EXTERNAL['replay_checkpoint.py'] and pins['payload_manifest_sha256']==EXTERNAL['PAYLOAD_MANIFEST.json'],'linked external controls')
 require(pins['core_path']==PREFIX+'/core' and pins['checkpoint_path']==PREFIX and pins['registration_path']==PREFIX+'/FULL_REGISTRATION.json' and pins['public_go_path']=='PUBLIC_GO.json','fixed source layout')
 def document(name):return parse(read_member(root,name,True)[1])
 reg=document(pins['registration_path']);go=document(pins['public_go_path'])
 require(actual[pins['registration_path']]['sha256']==pins['registration_sha256']==REG and actual['PUBLIC_GO.json']['sha256']==pins['public_go_sha256']==GO,'reg/GO fixed pins')
 require(go['status']=='PASS_READBACK_PUBLIC_BYTE_GO' and go['freeze_commit']=='3fc078553d1c9b10f11099870195c7a317951a71' and go['source_go'] is True and go['decode_go'] is True and go['registered_file_pins']==reg['files'],'genuine accepted source GO')
 core={name.removeprefix(PREFIX+'/core/'):item for name,item in actual.items() if name.startswith(PREFIX+'/core/')}
 require(core==reg['files'] and len(core)==41,'accepted complete core005')
 require(actual['bootstrap_replay.py']['sha256']==BOOT,'reviewed bootstrap source')
 limits={'address_space_bytes':536870912,'aggregate_output_bytes':134217728,'child_processes_allowed':False,'core':0,'cpu_seconds':900,'file_size_bytes':134217728,'network_allowed':False,'nonroot':True,'nproc':0,'single_threads':True,'wall_seconds':900}
 require(json.dumps(pins['resource_limits'],sort_keys=True)==json.dumps(limits,sort_keys=True),'unchanged caps')
 require(pins['runtime']=={'packages':{'mpmath':'1.3.0','python-flint':'0.9.0','sympy':'1.14.0'},'platform':'linux','python_version':'3.12.14'},'unchanged arithmetic runtime')
 review=document('reviews/actual/FINAL_ACTUAL_REVIEW_RECEIPT.json')
 require(actual['reviews/actual/FINAL_ACTUAL_REVIEW_RECEIPT.json']['sha256']==ACTUAL_REVIEW and review['status']=='PASS_INDEPENDENT_ACTUAL_ENDPOINT_EXPORTED_CERTIFICATE_REVIEW','actual scientific acceptance')
 for mode in ('normal','optimized'):
  prefix=pins['saved_runs'][mode]+'/'
  saved={name[len(prefix):]:item for name,item in actual.items() if name.startswith(prefix)}
  require(saved==review['bounded_actual_runs'][mode]['all_output_file_pins'],'saved genuine actual mode bytes:'+mode)
 spec=document(PREFIX+'/core/provenance_decoder/INPUT_SPEC.json');expected={};known_sizes={}
 for kind in ('producer','manifest'):expected[spec[kind+'_repository_path']]=spec[kind+'_sha256']
 audit=spec['prior_static_audit'];expected[audit['repository_path']]=audit['sha256'];known_sizes[audit['repository_path']]=audit['bytes']
 for capsule in spec['capsules']:
  expected[capsule['repository_path']]=capsule['reader_inputspec']['file_sha256'];known_sizes[capsule['repository_path']]=capsule['reader_inputspec']['file_bytes']
  original=capsule['original'];expected[original['repository_path']]=original['sha256'];known_sizes[original['repository_path']]=original['bytes']
 dependencies=pins['retained_input_dependencies'];require(len(expected)==11 and set(expected)==set(dependencies) and not set(expected)&set(names),'exact eleven external dependencies')
 require(pins['retained_inputs_commit']=='13a30ef46f5c90c1b01ae83b609db65fd0f8a709','immutable original origin')
 declaration=document('ORIGINAL_INPUT_DEPENDENCIES.json');require(declaration['files']==dependencies and declaration['origin_commit']==pins['retained_inputs_commit'] and declaration['original_input_bytes_embedded_in_this_bundle']==0,'opaque dependency declaration')
 for name,item in dependencies.items():
  require(item['sha256']==expected[name] and (name not in known_sizes or item['bytes']==known_sizes[name]),'dependency metadata from frozen InputSpec')
  measured,_=read_member(args.repository_root,name);require(measured==item,'original opaque dependency bytes:'+name)
 require(sum(item['bytes'] for item in dependencies.values())==213114777,'complete opaque dependency byte count')
 result={'status':'PASS_INDEPENDENT_ACTUAL_PACKAGE_AND_OPAQUE_DEPENDENCIES','external_controls':controls,'bootstrap_sha256':BOOT,
 'package':str(root),'payload_files':100,'total_files':103,'payload_bytes':sum(item['bytes'] for item in actual.values()),
 'source_files':41,'source_registration_sha256':REG,'public_GO_sha256':GO,'freeze_commit':go['freeze_commit'],
 'original_dependency_root':str(args.repository_root.absolute()),'opaque_dependency_files':11,'opaque_dependency_bytes':sum(item['bytes'] for item in dependencies.values()),
 'opaque_dependency_pins':dependencies,'saved_actual_modes_match_independent_acceptance':True,'actual_review_sha256':ACTUAL_REVIEW,
 'original_archive_array_decodes':0,'physical_source_callbacks':0,'producer_executions':0,'package_code_imported':0,
 'next_step':'Captured-bootstrap verify-only saved readback; fresh pair only after that passes under existing root authorization.'}
 with args.receipt.open('x') as out:json.dump(result,out,sort_keys=True,indent=2);out.write('\n')
 print(json.dumps({key:value for key,value in result.items() if key!='opaque_dependency_pins'},sort_keys=True))
if __name__=='__main__':main()
