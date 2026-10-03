#!/usr/bin/env python3
"""Stdlib-only authenticated entry; no numerical/physical import before frozen verification."""
from __future__ import annotations
import argparse,hashlib,importlib.util,importlib.metadata,json,re,sys,time,resource
from pathlib import Path


def require(ok,message):
 if not ok:raise RuntimeError(message)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 require(spec is not None and spec.loader is not None,'Cannot load exact registered module')
 module=importlib.util.module_from_spec(spec);sys.modules[name]=module
 spec.loader.exec_module(module)
 return module


def authenticate(root,registration_sha256,freeze_commit):
 root=Path(root).absolute()
 require(re.fullmatch('[0-9a-f]{64}',registration_sha256) is not None,'Explicit registration SHA256 required')
 require(re.fullmatch('[0-9a-f]{40}',freeze_commit) is not None,'Explicit publicfreeze40hex required')
 for path in (root,*root.parents):require(not path.is_symlink(),'Checkpoint symlink ancestry rejected')
 registration=root/'FULL_REGISTRATION.json';helper=root/'code'/'active_integrity.py'
 entry=root/'independent'/'diagnostic_independent.py';core=root/'independent'/'run_independent.py'
 for path in (registration,helper,entry,core,root/'code',root/'independent'):
  require(not path.is_symlink(),'Bootstrap source symlink rejected')
 require(Path(__file__).absolute()==entry,'Execute the exact registered independent entry')
 require(registration.is_file() and sha(registration)==registration_sha256,'Bootstrap registration hash differs')
 files=json.loads(registration.read_text())['files']
 for relative,path in [('code/active_integrity.py',helper),('independent/diagnostic_independent.py',entry),('independent/run_independent.py',core)]:
  require(path.is_file() and sha(path)==files.get(relative),'Bootstrap registered source mismatch '+relative)
 module=load(helper,'active_integrity')
 verified=module.verify_frozen(root,registration_sha256,freeze_commit)
 return module,verified


def main():
 started=time.monotonic()
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--checkpoint-root',type=Path,required=True)
 parser.add_argument('--freeze-commit',required=True)
 parser.add_argument('--registration-sha256',required=True)
 parser.add_argument('--output-dir',type=Path,required=True)
 args=parser.parse_args()
 helper,verified=authenticate(args.checkpoint_root,args.registration_sha256,args.freeze_commit)
 require(sys.flags.optimize==0 and sys.version_info[:3]==(3,12,14),'Pinned Python3.12.14 normal runtime required')
 for name,version in helper.VERSIONS.items():require(importlib.metadata.version(name)==version,'Pinned dependency mismatch '+name)
 helper.fresh_external(args.checkpoint_root,args.output_dir)
 core=load(args.checkpoint_root/'independent'/'run_independent.py','active_source_independent_core')
 core.AUTHENTICATED_BY_WRAPPER={'checkpoint_root':str(args.checkpoint_root.absolute()),'registration_sha256':args.registration_sha256,'freeze_commit':args.freeze_commit}
 core.ROUTE_STARTED=started
 core.main()
 helper.verify_frozen(args.checkpoint_root,args.registration_sha256,args.freeze_commit)
 elapsed=time.monotonic()-started;peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 require(elapsed<=900 and peak<=262144,'Complete wrapper/numerical route resource cap exceeded')
 output=args.output_dir/'diagnostic.json';report=json.loads(output.read_text())
 report['resources'].update(elapsed_seconds=elapsed,peak_rss_kib=peak,scope='Authenticated wrapper/import/provenance through finalserialization; all12cases threecontrols bothMPcontexts',limit_seconds=900,limit_peak_rss_kib=262144)
 report['provenance']['full_frozen_verification']=verified
 report['provenance']['post_run_frozen_verification']='PASS'
 report['producer_sha256']=sha(__file__)
 output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
 require(time.monotonic()-started<=900,'Post-verification serialization exceedsroutebudget')
 print(json.dumps({'status':report['status'],'resources':report['resources'],'post_run_frozen_verification':'PASS'}),flush=True)
if __name__=='__main__':main()
