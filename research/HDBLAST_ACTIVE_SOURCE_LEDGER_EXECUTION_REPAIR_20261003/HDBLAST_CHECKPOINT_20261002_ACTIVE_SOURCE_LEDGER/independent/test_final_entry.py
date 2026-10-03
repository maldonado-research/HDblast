"""Final entry/array/serialization guards on fabricated values only."""
import argparse,hashlib,json,sys,tempfile,time,resource
from pathlib import Path
import numpy as np
import mpmath as mp
import diagnostic_independent as wrapper
import run_independent as core
from kernel import LD,CD


def need(test,message):
 if not test:raise RuntimeError(message)


def main():
 started=time.monotonic();checks=[]
 # Exact fullshape fabricated schemas; no source or phase/integral evaluation.
 for setting,n,count in [('coarse',8192,577),('fine',16384,1153)]:
  a={'k':np.linspace(LD('0.0009'),LD('255.99'),n,dtype=LD),'momentum_weights':np.full(n,LD(1),dtype=LD),
     'observation_eta':np.array(['-5.5','-4.5','-4','-3.5','-2.5','-1.5'],dtype=LD),
     'history_eta':-LD(6)+np.arange(count,dtype=LD)/(128 if setting=='coarse' else 256),'history_values':np.zeros((count,3,9),dtype=LD),'history_ledger_integrand':np.zeros((count,3),dtype=LD),
     'history_baseline_contact':np.zeros((count,3),dtype=LD),'history_source_jet':np.zeros((count,6),dtype=LD),
     'history_forcing_jet':np.zeros((count,4),dtype=LD),'history_geometry':np.zeros((count,10),dtype=LD),'history_cutoffs':np.array([64,128,256],dtype=np.int64),
     'history_quantity_names':np.array(['q','q_prime','q_second','rho','p','Q0','rho0','p0','current']),
     'history_geometry_names':np.array(['a0','L','L_prime','L_second','L_third','M','M_prime','M_second','M_third','M_fourth'])}
  a.update({key:np.zeros(n,dtype=CD) for key in ('u_1','w_1','u_2','w_2','u_3','w_3')})
  core.validate_arrays(a,setting);checks.append(setting+'_exact19schema')
  bad=dict(a);bad['u_1']=np.zeros(n,dtype=np.complex128)
  try:core.validate_arrays(bad,setting)
  except RuntimeError:checks.append(setting+'_badmodedtype_rejected')
  else:raise RuntimeError('Wrong complex precision accepted')
  bad=dict(a);bad['history_values']=a['history_values'].copy();bad['history_values'][0,0,0]=np.nan
  try:core.validate_arrays(bad,setting)
  except RuntimeError:checks.append(setting+'_nonfinitehistory_rejected')
  else:raise RuntimeError('Nonfinite native history accepted')
  bad=dict(a);bad['observation_eta']=a['observation_eta'].copy();bad['observation_eta'][1]+=LD('0.125')
  try:core.validate_arrays(bad,setting)
  except RuntimeError:checks.append(setting+'_observation_epoch_mutation_rejected')
  else:raise RuntimeError('Changed six-epoch observation inventory accepted')
  bad=dict(a);bad['history_eta']=a['history_eta'].copy();bad['history_eta'][7]+=LD('0.00001')
  try:core.validate_arrays(bad,setting)
  except RuntimeError:checks.append(setting+'_full_history_grid_mutation_rejected')
  else:raise RuntimeError('Changed full historical grid accepted')
  bad=dict(a);bad['k']=a['k'].copy();bad['k'][n//4-1]=LD('64.0001')
  try:core.validate_arrays(bad,setting)
  except RuntimeError:checks.append(setting+'_K_prefix_count_mutation_rejected')
  else:raise RuntimeError('Changed K-prefix count accepted')

 for dps in (80,100):
  ctx=mp.mp.clone();ctx.dps=dps
  for x in ('-0.00003271','2.031','1.1834e-17'):
   core.serial(ctx,ctx.mpf(x))
  checks.append('exactFraction_serialization_'+str(dps))
 for x in ('-0.00003271','2.031','1.1834e-17'):core.serial_native(LD(x))
 checks.append('native_binary80_profile_roundtrip')
 # A symlink root must be rejected before helper or numerical imports.
 with tempfile.TemporaryDirectory(dir='/tmp') as td:
  t=Path(td);(t/'real').mkdir();(t/'alias').symlink_to(t/'real',target_is_directory=True)
  try:wrapper.authenticate(t/'alias','0'*64,'0'*40)
  except RuntimeError as e:need('symlink' in str(e),'Wrong bootstrap failure order');checks.append('symlink_before_resolution')
  else:raise RuntimeError('Symlink checkpoint accepted')
 # Unauthorized direct core must reject before reaching np.load/source.
 old=sys.argv
 sys.argv=['core','--checkpoint-root','/tmp/not-an-input','--freeze-commit','0'*40,'--registration-sha256','0'*64,'--output-dir','/tmp/never-physical-output']
 try:core.main()
 except RuntimeError as e:need('authenticated' in str(e),'Wrong core failure order');checks.append('standalone_core_rejected_before_reader')
 else:raise RuntimeError('Unauthorized core accepted')
 finally:sys.argv=old
 here=Path(__file__).resolve().parent
 result={'status':'PASS_FINAL_ENTRY_SCHEMA_AND_SERIALIZATION_GUARDS','physical_evaluations':0,'actual_arrays_read':False,'study_source_calls':0,
         'checks':checks,'elapsed_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
         'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in here.glob('*.py')},
         'scope':'Final standalone-auth wrapper and array/serialization additions over the unchanged numericalcore benchmark.'}
 name='FINAL_ENTRY_GUARDS_OPTIMIZED.json' if not __debug__ else 'FINAL_ENTRY_GUARDS_NORMAL.json'
 (here/name).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='files'},indent=2))
if __name__=='__main__':main()
