"""Independent source and manufactured custodian byte/counter verification.

No production module, decoder, numerical library or physical source is imported.
"""
from pathlib import Path
import hashlib
import json
import os
import stat
import sys
BASE=Path(__file__).resolve().parent
STABLE=('SOURCE_CERTIFICATE.json','TARGET_BUDGETS.json','NODE_CERTIFICATES.jsonl.gz',
'DATA.json','SCIENCE_SUMMARY.json','SOURCE_ATTEMPTS.jsonl','DECODE_ATTEMPTS.jsonl',
'INPUT_AUTHENTICATION.json','OUTPUT_READBACK.json')
LIMITS={'address_space_bytes':536870912,'aggregate_output_bytes':134217728,
'core_bytes':0,'cpu_seconds':900,'custodian_receipt_reserve_bytes':65536,
'file_size_bytes':134217728,'process_limit':0,'wall_seconds':900}
def require(ok,message):
 if not ok:raise ValueError(message)
def pin(path):
 data=path.read_bytes();return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def inventory(root):
 files={};dirs=set()
 def inaccessible(error):raise error
 for parent,names,leaves in os.walk(root,followlinks=False,onerror=inaccessible):
  for name in names+leaves:
   p=Path(parent)/name;s=p.lstat();relative=str(p.relative_to(root))
   require(not stat.S_ISLNK(s.st_mode),'symlink:'+relative)
   if stat.S_ISDIR(s.st_mode):dirs.add(relative)
   else:
    require(stat.S_ISREG(s.st_mode) and s.st_nlink==1,'special/hardlink:'+relative)
    files[relative]=pin(p)
 return files,dirs
def run(source,registration,normal,optimized,label):
 source,registration,normal,optimized=map(Path,(source,registration,normal,optimized))
 reg=json.loads(registration.read_bytes());regpin=pin(registration)
 sourcepins,dirs=inventory(source);expected=set()
 for name in reg['files']:
  expected.update(str(p) for p in Path(name).parents if str(p)!='.')
 require(sourcepins==reg['files'] and expected<=dirs,'complete source file roster differs')
 runs={};science=[]
 for mode,path in [('normal',normal),('optimized',optimized)]:
  files,outdirs=inventory(path);entry=json.loads((path/'ENTRY_RECEIPT.json').read_bytes())
  execution=json.loads((path/'EXECUTION.json').read_bytes());readback=json.loads((path/'OUTPUT_READBACK.json').read_bytes())
  require(execution['status']=='PASS_BOUNDED_MANUFACTURED_EXECUTION','bounded status')
  require(execution['limits']==LIMITS and execution['manufactured_only'] is True and execution['optimized'] is (mode=='optimized'),'policy/mode')
  for field in ['exit_code','wait_status','wrapper_exit_code']:
   require(type(execution[field]) is int and execution[field]==0,'nonzero wait/exit')
  require(type(execution['uid']) is int and execution['uid']!=0,'root run')
  require(execution['monitor_error'] is None and execution['stop_reason'] is None,'custodian rejected')
  require(execution['wall_seconds']<900 and execution['cpu_seconds_wait4']<900 and execution['peak_rss_kib_wait4']<524288,'resource bound')
  require(execution['entry_receipt']==files['ENTRY_RECEIPT.json'],'entry byte binding')
  require(execution['captured_child_log_sha256']==execution['observed_child_log_sha256']==files['child.log']['sha256'],'log custody')
  require(entry['status']=='PASS_MANUFACTURED_ENDPOINT_ENTRY' and entry['manufactured'] is True and entry['go_sha256'] is None,'entry scope')
  require(entry['registration_sha256']==regpin['sha256'] and entry['scope']==reg['scope'],'registration binding')
  require(entry['outputs']=={name:files[name] for name in (*STABLE,'RESOURCE.json')},'all entry output hashes')
  for receipt in (entry,readback):
   require(type(receipt['nodes']) is int and receipt['nodes']==49152 and type(receipt['endpoint_comparisons']) is int and receipt['endpoint_comparisons']==98304,'coverage')
  require(readback['status']=='PASS_STANDALONE_ENDPOINT_READBACK' and readback['prefix_cases']==24 and readback['source_evaluations']==readback['retained_arrays_opened']==0,'readback status')
  require({'socket','socketpair','connect','bind','listen','accept','accept4','fork','vfork','clone','clone3','execve','execveat','io_uring_setup','io_uring_enter','io_uring_register'}<=set(entry['kernel_denied_syscalls']),'kernel policy')
  tree=hashlib.sha256();count=0;size=0
  for name,p in sorted(files.items()):
   if name=='EXECUTION.json':continue
   tree.update(name.encode('utf-8',errors='surrogateescape')+b'\0'+str(p['bytes']).encode()+b'\0'+bytes.fromhex(p['sha256']));count+=1;size+=p['bytes']
  require(execution['worker_output']=={'bytes':size,'file_count':count,'directory_count':len(outdirs),'tree_sha256':tree.hexdigest()},'custodian complete output inventory')
  require(size+files['EXECUTION.json']['bytes']<134217728,'aggregate cap')
  command=execution['command']
  require('-I' in command and '-B' in command and ('-O' in command)==(mode=='optimized') and '--manufactured-only' in command and '--go' not in command,'supported command')
  require(command[command.index('--registration-sha256')+1]==regpin['sha256'],'command registration pin')
  require(files['SOURCE_ATTEMPTS.jsonl']['bytes']==files['DECODE_ATTEMPTS.jsonl']['bytes']==0,'real attempt journals')
  science.append({name:files[name] for name in STABLE})
  runs[mode]={'directory':str(path),'entry_receipt':files['ENTRY_RECEIPT.json'],'execution_receipt':files['EXECUTION.json'],
  'resource_receipt':files['RESOURCE.json'],'child_log':files['child.log'],'worker_output':execution['worker_output'],
  'wall_seconds':execution['wall_seconds'],'cpu_seconds_wait4':execution['cpu_seconds_wait4'],'peak_rss_kib_wait4':execution['peak_rss_kib_wait4'],'file_pins':files}
 require(science[0]==science[1],'nine science outputs differ')
 result={'status':'PASS_INDEPENDENT_BOUNDED_PAIR_AND_SOURCE_BYTES','registration':regpin,'source_files':sourcepins,
 'source_directories':sorted(dirs),'extra_empty_source_directories':sorted(dirs-expected),'stable_scientific_files':science[0],'runs':runs,'retained_decodes':0,
 'physical_source_calls':0,'production_modules_imported':0,'actual_GO_issued':False}
 (BASE/label).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
 print(json.dumps({'status':result['status'],'source_count':len(sourcepins),'stable_file_count':len(science[0]),'receipt':label},sort_keys=True))
if __name__=='__main__':run(*sys.argv[1:])
