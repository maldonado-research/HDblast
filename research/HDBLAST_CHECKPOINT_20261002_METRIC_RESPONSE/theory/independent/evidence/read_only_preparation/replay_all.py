from pathlib import Path
import hashlib,json,os,platform,subprocess,sys
root=Path('/workspace/hdblast-research-work/metric-independent-theory')
output=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
state=lambda:{str(p.relative_to(root)):sha(p) for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
before=state();runs=[]
for script,prefix in [('verify_metric_response.py','EXACT_PROOF'),('verify_cse_specialization.py','CSE_AGREEMENT'),('verify_primary_inventory_agreement.py','PRIMARY_INVENTORY_AGREEMENT')]:
 for optimized in (False,True):
  name=prefix+('_OPTIMIZED' if optimized else '')
  destination=output/(name+'.json')
  command=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(output/'run_readonly_proof.py'),'--proof',str(root/script),'--output',str(destination)]
  with (output/(name+'.log')).open('wb') as log:
   run=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},check=False)
  (output/(name+'.exit')).write_text(str(run.returncode)+'\n')
  if run.returncode:raise RuntimeError('Read-only proof failed: '+name)
  receipt=json.loads(destination.read_text())
  runs.append({'script':script,'optimized':optimized,'command':command,'exit_code':run.returncode,'status':receipt['status'],'identity_count':receipt['identity_count'],'receipt_sha256':sha(destination),'log_sha256':sha(output/(name+'.log'))})
# Explicit export is allowed only to the external requested file.
export=output/'CONTACT_INVENTORY_EXPLICIT_EXPORT.json'
command=[sys.executable,'-B',str(output/'run_readonly_proof.py'),'--proof',str(root/'explore_contacts.py'),'--output',str(export)]
with (output/'CONTACT_EXPORT.log').open('wb') as log:run=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},check=False)
(output/'CONTACT_EXPORT.exit').write_text(str(run.returncode)+'\n')
if run.returncode:raise RuntimeError('Explicit contact export failed')
if export.read_bytes()!=(root/'FINITE_K_CONTACT_INVENTORY.json').read_bytes():raise RuntimeError('Explicit export changed exact contact inventory')
after=state()
if after!=before:raise RuntimeError('Source tree bytes changed during read-only proofs')
receipt={'status':'PASS','scope':'Pure symbolic proofs only; zero physical evaluations','python':platform.python_version(),'runs':runs,'source_tree_identical_before_after':True,'protected_source_files':len(before),'source_write_audit_guard_active':True,'bytecode_writes_disabled':True,'explicit_external_inventory_export_matches_legacy_bytes':True,'inventory_sha256':sha(export),'runner_sha256':sha(output/'run_readonly_proof.py'),'orchestrator_sha256':sha(Path(__file__))}
(output/'READ_ONLY_REPLAY.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'PASS','proof_runs':len(runs),'protected_source_files':len(before),'identities_per_round':[x['identity_count'] for x in runs],'explicit_export_bytes_identical':True}))
