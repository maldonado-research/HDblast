"""Independent fresh portable custody and byte-identity audit; no study imports."""
from pathlib import Path
import hashlib
import json
import math
import os
import stat
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
PACKAGE=BASE/'delivery/actual-package'
FRESH=HERE/'fresh-pair-001'
PREFIX='research/HDBLAST_CHECKPOINT_20261008_RETAINED_ENDPOINT_FLOW'
REG='33eded898f47a96142f85002b896edd4c294180653f50f2b38e83f22d5a75796'
GO='edded4f67840bc54e95673c5f43b71f6804e373d13211c90ba48ca99f72a8389'
LIMITS={'address_space_bytes':536870912,'aggregate_output_bytes':134217728,'core_bytes':0,'cpu_seconds':900,
'custodian_receipt_reserve_bytes':65536,'file_size_bytes':134217728,'process_limit':0,'wall_seconds':900}
DENIED={'socket','socketpair','connect','bind','listen','accept','accept4','fork','vfork','clone','clone3','execve','execveat','io_uring_setup','io_uring_enter','io_uring_register'}

def require(ok,message):
 if not ok:raise ValueError(message)
def unique(items):
 d={}
 for key,value in items:require(key not in d,'duplicate key');d[key]=value
 return d
def load(path):return json.loads(path.read_bytes(),object_pairs_hook=unique,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('nonfinite JSON')))
def pin(path):
 metadata=path.lstat();require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink==1,'unaliased regular artifact')
 h=hashlib.sha256();count=0
 with path.open('rb') as stream:
  while raw:=stream.read(1<<20):h.update(raw);count+=len(raw)
 require(count==metadata.st_size and (metadata.st_ino,metadata.st_size,metadata.st_mtime_ns)==(path.stat().st_ino,path.stat().st_size,path.stat().st_mtime_ns),'artifact changed')
 return {'bytes':count,'sha256':h.hexdigest()}

def main():
 acceptance=load(BASE/'actual-independent-review/FINAL_ACTUAL_REVIEW_RECEIPT.json')
 require(pin(BASE/'actual-independent-review/FINAL_ACTUAL_REVIEW_RECEIPT.json')['sha256']=='097fe778797dc78f1dc800885bd064e7406c6eb80f9642b3c7cc8577d05754ca','actual acceptance changed')
 science=acceptance['scientific_file_pins'];require(len(science)==9,'nine fixed scientific files')
 helper=load(FRESH/'PORTABLE_REPLAY_RECEIPT.json')
 require(helper['status']=='PASS_COMPLETE_RETAINED_ENDPOINT_SOURCE_REPLAY' and helper['manufactured_control'] is False and
         helper['new_physical_source_callbacks']==256 and helper['new_retained_array_decodes']==72,'complete authorized actual pair')
 runs={}
 for mode in ('normal','optimized'):
  root=FRESH/mode;require(root.is_dir() and not root.is_symlink(),'fresh real mode directory')
  names={p.name for p in root.iterdir()};require(names==set(science)|{'RESOURCE.json','ENTRY_RECEIPT.json','EXECUTION.json','child.log'},'complete actual output roster')
  files={name:pin(root/name) for name in sorted(names)}
  for name,expected in science.items():
   require(files[name]==expected==pin(PACKAGE/'evidence/actual'/mode/name),'fresh/package/root scientific bytes differ:'+mode+'/'+name)
  entry=load(root/'ENTRY_RECEIPT.json');execution=load(root/'EXECUTION.json');resource=load(root/'RESOURCE.json');readback=load(root/'OUTPUT_READBACK.json')
  require(entry['status']=='PASS_REGISTERED_ENDPOINT_ENTRY' and entry['manufactured'] is False and entry['scope']=='LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE','entry scope')
  require(entry['registration_sha256']==REG and entry['go_sha256']==GO,'entry exact authorization')
  require(entry['outputs']=={name:files[name] for name in (*science,'RESOURCE.json')},'entry complete output pins')
  require(execution['status']=='PASS_BOUNDED_REGISTERED_EXECUTION' and execution['manufactured_only'] is False and execution['optimized'] is (mode=='optimized'),'genuine actual mode status')
  require(execution['limits']==LIMITS and all(type(v) is int for v in execution['limits'].values()),'unchanged hard limits')
  for field in ('exit_code','wait_status','wrapper_exit_code'):
   require(type(execution[field]) is int and execution[field]==0,'authoritative wait4/exit failed')
  require(type(execution['uid']) is int and execution['uid']>0,'nonroot worker')
  require(all(execution[key] is None for key in ('stop_reason','monitor_error','preserved_child_log')),'custody failure')
  require(execution['entry_receipt']==files['ENTRY_RECEIPT.json'],'custodian entry binding')
  require(execution['captured_child_log_sha256']==execution['observed_child_log_sha256']==files['child.log']['sha256'],'log custody')
  require(execution['aggregate_limit_enforcement']=='POLL_AND_FINAL_ACCEPTANCE_CHECK_WITH_RECEIPT_RESERVE','aggregate cap semantics')
  for field,cap in [('wall_seconds',900),('cpu_seconds_wait4',900),('peak_rss_kib_wait4',524288)]:
   require(type(execution[field]) in (int,float) and math.isfinite(execution[field]) and 0<=execution[field]<=cap,'measured resource cap')
  command=execution['command'];flags=['-I','-B']+(['-O'] if mode=='optimized' else []);index=1+len(flags)
  require(command[1:index]==flags and command[index]==str(PACKAGE/PREFIX/'core/execution/run_endpoint.py'),'exact isolated worker command')
  arguments=command[index+1:];require(len(arguments)==14 and len(set(arguments[::2]))==7,'exact unique operation arguments')
  require(dict(zip(arguments[::2],arguments[1::2]))=={'--root':str(PACKAGE/PREFIX/'core'),'--output':str(root),
   '--registration':str(PACKAGE/PREFIX/'FULL_REGISTRATION.json'),'--registration-sha256':REG,
   '--repository-root':'/workspace/HDblast','--go':str(PACKAGE/'PUBLIC_GO.json'),'--go-sha256':GO},'fresh command source/output/original input/GO bindings')
  require(execution['custodian_output_directory']==str(root),'fresh held directory')
  require(len(entry['kernel_denied_syscalls'])==len(DENIED) and set(entry['kernel_denied_syscalls'])==DENIED,'complete kernel denial list')
  h=hashlib.sha256();size=0
  for name,item in files.items():
   if name=='EXECUTION.json':continue
   h.update(name.encode()+b'\0'+str(item['bytes']).encode()+b'\0'+bytes.fromhex(item['sha256']));size+=item['bytes']
  require(execution['worker_output']=={'bytes':size,'file_count':12,'directory_count':0,'tree_sha256':h.hexdigest()},'fresh complete custody tree')
  require(size<=134217728-65536 and sum(item['bytes'] for item in files.values())<=134217728,'aggregate limits')
  require(resource['output_bytes']==sum(item['bytes'] for name,item in files.items() if name not in ('RESOURCE.json','ENTRY_RECEIPT.json','EXECUTION.json')),'worker resource snapshot bytes')
  require(type(resource['peak_rss_kib']) is int and 0<resource['peak_rss_kib']<=execution['peak_rss_kib_wait4'] and 0<=resource['wall_seconds']<=900,'worker resource counters')
  for evidence in (entry,readback):require(type(evidence['nodes']) is int and evidence['nodes']==49152 and type(evidence['endpoint_comparisons']) is int and evidence['endpoint_comparisons']==98304,'complete endpoint coverage')
  require(readback['manufactured'] is False and readback['prefix_cases']==24 and readback['source_evaluations']==readback['retained_arrays_opened']==0,'readback scope')
  for name,count in [('SOURCE_ATTEMPTS.jsonl',128),('DECODE_ATTEMPTS.jsonl',36)]:require(len((root/name).read_bytes().splitlines())==count,'complete actual attempt journal')
  require(helper['fresh_custody'][mode]['execution']==files['EXECUTION.json'] and helper['fresh_custody'][mode]['entry']==files['ENTRY_RECEIPT.json'],'helper own fresh custody links')
  require(helper['fresh_checks'][mode]==readback,'helper final readback bytes')
  runs[mode]={'file_pins':files,'wall_seconds':execution['wall_seconds'],'cpu_seconds_wait4':execution['cpu_seconds_wait4'],'peak_rss_kib_wait4':execution['peak_rss_kib_wait4'],'worker_output':execution['worker_output']}
 require(helper['fresh_comparison']['files']==science==helper['saved_comparison']['files'],'helper scientific identity authority')
 result={'status':'PASS_INDEPENDENT_FRESH_ACTUAL_PORTABLE_PAIR_CUSTODY_AND_IDENTITY','fresh_producer_executions':2,
 'normal_optimized_modes_genuine':True,'science_files_identical_to_both_saved_root_modes_and_package':9,
 'scientific_file_pins':science,'runs':runs,'portable_replay_receipt':pin(FRESH/'PORTABLE_REPLAY_RECEIPT.json'),
 'authorized_new_physical_source_callbacks':256,'authorized_new_retained_array_decodes':72,'source_free_custody_auditor_new_calls':0,
 'parameters_adapted':False,'new_GO_issued':False,'source_and_actual_review_remain_separately_pinned':True}
 with (HERE/'FRESH_PAIR_INDEPENDENT_AUDIT.json').open('x') as out:json.dump(result,out,sort_keys=True,indent=2);out.write('\n')
 print(json.dumps({key:value for key,value in result.items() if key not in ('runs','scientific_file_pins')},sort_keys=True))
if __name__=='__main__':main()
