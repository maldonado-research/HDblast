"""Execute the final guarded driver entirely on fabricated inputs/sources."""
import hashlib,json,sys,time,resource
from pathlib import Path
import numpy as np
import engine
import run_independent as driver
from kernel import LD,CD,forcing

here=Path(__file__).resolve().parent
base=here/'fabricated-final-driver'
if base.exists():raise RuntimeError('Do not replace a fabricated run receipt')
checkpoint=base/'checkpoint';inputs=checkpoint/'inputs';inputs.mkdir(parents=True)


def fakejet(t,power):
 import math
 z=np.asarray(t,dtype=LD)+LD(4)
 return np.stack([LD(math.factorial(power)//math.factorial(power-n))*z**(power-n) if n<=power else np.zeros_like(z) for n in range(6)])*LD('0.0000001')


def fake_sources(prep,source):
 power={'positive_B':4,'signed_uB':5}[source]
 grid=prep['grid'];dt=prep['dt'];t=prep['t'];s=prep['s']
 nested=grid[:-1,None,None]+t[None,:,None]*s[None,None,:]
 end=grid[:-1,None]+dt*s[None,:]
 return {'grid':fakejet(grid,power),'dense':fakejet(prep['dense'],power),
         'nested_g':forcing(nested,fakejet(nested,power)),
         'end_g':forcing(end,fakejet(end,power))}


def forbidden_source(*args,**kwargs):raise RuntimeError('Actual study source must never run in this fabricated preflight')
engine.source_jet=forbidden_source
driver.physical_sources=fake_sources
for setting,nodes,perunit in [('coarse',8192,128),('fine',16384,256)]:
 k=np.linspace(LD('0.0009'),LD('255.999'),nodes,dtype=LD)
 weights=np.full(nodes,LD(256)/nodes,dtype=LD)
 eta=-LD(6)+np.arange(int(LD('4.5')*perunit)+1,dtype=LD)/perunit
 for source in ('positive_B','signed_uB'):
  state={'k':k,'momentum_weights':weights,'history_eta':eta,'history_values':np.zeros((len(eta),3,9),dtype=LD),
         'history_ledger_integrand':np.zeros((len(eta),3),dtype=LD),
         'observation_eta':np.array(['-5.5','-4.5','-4','-3.5','-2.5','-1.5'],dtype=LD),
         'history_baseline_contact':np.zeros((len(eta),3),dtype=LD),'history_source_jet':np.zeros((len(eta),6),dtype=LD),
         'history_forcing_jet':np.zeros((len(eta),4),dtype=LD),'history_geometry':np.zeros((len(eta),10),dtype=LD),
         'history_cutoffs':np.array([64,128,256],dtype=np.int64),
         'history_quantity_names':np.array(['q','q_prime','q_second','rho','p','Q0','rho0','p0','current']),
         'history_geometry_names':np.array(['a0','L','L_prime','L_second','L_third','M','M_prime','M_second','M_third','M_fourth'])}
  for index in (1,2,3):
   state['u_'+str(index)]=np.full(nodes,CD('0.000000000001')+CD('0.0000000000007j'),dtype=CD)+index*CD('0.00000000000001')
   state['w_'+str(index)]=np.full(nodes,CD('0.0000000000002')-CD('0.0000000000003j'),dtype=CD)+index*CD('0.000000000000005j')
  np.savez(inputs/f'metric_modes_{source}_{setting}.npz',**state)
manifest={'status':'FABRICATED_ONLY','actual_physical_arrays':False,'fake_mode_nodes':[8192,16384],'fake_histories':[577,1153]}
(inputs/'INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
files={str(p.relative_to(checkpoint)):hashlib.sha256(p.read_bytes()).hexdigest() for p in checkpoint.rglob('*') if p.is_file()}
reg={'schema_version':1,'physical_diagnostic_evaluations_before_freeze':0,'scientific_outcome_before_freeze':'UNCOMPUTED_FABRICATED_ONLY','files':files}
regpath=checkpoint/'FULL_REGISTRATION.json';regpath.write_text(json.dumps(reg,indent=2)+'\n')
started=time.monotonic()
sys.argv=['run_independent.py','--checkpoint-root',str(checkpoint),'--freeze-commit','0'*40,
          '--registration-sha256',hashlib.sha256(regpath.read_bytes()).hexdigest(),'--output-dir',str(base/'output')]
driver.AUTHENTICATED_BY_WRAPPER={'checkpoint_root':str(checkpoint.absolute()),'registration_sha256':hashlib.sha256(regpath.read_bytes()).hexdigest(),'freeze_commit':'0'*40}
driver.ROUTE_STARTED=started
driver.main()
result=json.loads((base/'output'/'diagnostic.json').read_text())
receipt={'status':'PASS_FINAL_GUARDED_DRIVER_FABRICATED_FULL_ROUTE','physical_evaluations':0,'actual_study_source_blocked':True,
         'actual_input_capsules_read':False,'rows':len(result['rows']),'controls':result['controls'],'precisions':result['reduction_decimal_digits'],
         'elapsed_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
         'implementation_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in here.glob('*.py')},
         'result_sha256':hashlib.sha256((base/'output'/'diagnostic.json').read_bytes()).hexdigest(),
         'serialization_max_error_by_dps':result['serialization_max_error_by_dps'],'native_profile_serialization_max_error':result['native_profile_serialization_max_error']}
(here/'SYNTHETIC_FINAL_DRIVER_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2),flush=True)
