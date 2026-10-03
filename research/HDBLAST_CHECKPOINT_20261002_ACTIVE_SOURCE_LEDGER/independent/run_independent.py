#!/usr/bin/env python3
"""Prospective independent active-source direct-F diagnostic.
Physical inputs are inaccessible until explicit public-freeze/provenance arguments pass.
"""
from __future__ import annotations
import argparse,hashlib,json,re,time,resource
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import mpmath as mp
import numpy as np
from kernel import LD,CD,DenseKernel
from engine import prepare,physical_sources,trajectory,cut_prefix

CONTROLS=((16,16),(16,24),(24,24))
CUTOFFS=(64,128,256)
EPS_RATIO=(3777893186295716171,37778931862957161709568)
PI_RATIO=(14488038916154245685,4611686018427387904)
SERIAL_MAX={'80':Fraction(0),'100':Fraction(0)}
NATIVE_SERIAL_MAX=Fraction(0)
AUTHENTICATED_BY_WRAPPER=None
ROUTE_STARTED=None


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def require(test,message):
 if not test:raise RuntimeError(message)
def ratio(ctx,x):
 n,d=x.as_integer_ratio();return ctx.mpf(n)/ctx.mpf(d)
def serial(ctx,x):
 require(ctx.isfinite(x),'Nonfinite MP scientific field')
 text=ctx.nstr(x,55)
 sign,mantissa,exponent,_=x._mpf_
 reference=Fraction((-1 if sign else 1)*mantissa)*(Fraction(2)**exponent)
 SERIAL_MAX[str(ctx.dps)]=max(SERIAL_MAX[str(ctx.dps)],abs(Fraction(Decimal(text))-reference))
 return text
def serial_native(x):
 global NATIVE_SERIAL_MAX
 require(np.isfinite(x),'Nonfinite native scientific field')
 text=str(x)
 require(LD(text).as_integer_ratio()==x.as_integer_ratio(),'Native profile decimal failed binary80 roundtrip')
 NATIVE_SERIAL_MAX=max(NATIVE_SERIAL_MAX,abs(Fraction(Decimal(text))-Fraction(*x.as_integer_ratio())))
 return text


def verify(args):
 root=args.checkpoint_root.resolve()
 regpath=root/'FULL_REGISTRATION.json'
 require(re.fullmatch('[0-9a-f]{40}',args.freeze_commit) is not None,'Missing public 40hex freeze pin')
 require(sha(regpath)==args.registration_sha256,'Registration hash mismatch')
 reg=json.loads(regpath.read_text())
 require(reg.get('physical_diagnostic_evaluations_before_freeze')==0,'Prospective evaluation marker mismatch')
 require(isinstance(reg.get('files'),dict) and reg['files'],'Missing complete registration inventory')
 for rel,pin in reg['files'].items():
  p=(root/rel).resolve()
  require(p.is_relative_to(root) and p.is_file() and not p.is_symlink(),'Invalid registered path')
  require(sha(p)==pin,'Registered file changed: '+rel)
 require(not args.output_dir.resolve().is_relative_to(root),'Output must be external to checkpoint')
 require(not args.output_dir.exists(),'Refusing output overwrite')
 return {'freeze_commit':args.freeze_commit,'registration_sha256':args.registration_sha256,'registered_files_checked':len(reg['files']),'input_manifest_sha256':sha(root/'inputs'/'INPUT_MANIFEST.json')}


def validate_arrays(a,setting):
 names={'k','momentum_weights','u_1','w_1','u_2','w_2','u_3','w_3','observation_eta','history_eta','history_values','history_ledger_integrand','history_baseline_contact','history_source_jet','history_forcing_jet','history_geometry','history_cutoffs','history_quantity_names','history_geometry_names'}
 require(set(a)==names,'Exactly19 frozen capsule members required')
 nodes=8192 if setting=='coarse' else 16384
 history=577 if setting=='coarse' else 1153
 for key in ('k','momentum_weights'):
  require(a[key].shape==(nodes,) and a[key].dtype==np.dtype(LD) and np.all(np.isfinite(a[key])) and np.all(a[key]>0),'Invalid quadraturevalues/header')
 require(np.all(np.diff(a['k'])>0) and a['k'][-1]<256,'Momentum ordering/range changed')
 for key in ('u_1','w_1','u_2','w_2','u_3','w_3'):
  require(a[key].shape==(nodes,) and a[key].dtype==np.dtype(CD) and np.all(np.isfinite(a[key])),'Invalid retainedmodevalues/header')
 shapes={'observation_eta':(6,),'history_eta':(history,),'history_values':(history,3,9),'history_ledger_integrand':(history,3),'history_baseline_contact':(history,3),'history_source_jet':(history,6),'history_forcing_jet':(history,4),'history_geometry':(history,10)}
 for key,shape in shapes.items():
  require(a[key].shape==shape and a[key].dtype==np.dtype(LD) and np.all(np.isfinite(a[key])),'Invalid historyvalues/header '+key)
 require(np.array_equal(a['history_cutoffs'],np.array(CUTOFFS,dtype=np.int64)),'Cutoff inventory changed')
 require(np.array_equal(a['observation_eta'],np.array(['-5.5','-4.5','-4','-3.5','-2.5','-1.5'],dtype=LD)),'Six observation epochs changed')
 dt=LD(1)/(128 if setting=='coarse' else 256)
 require(np.array_equal(a['history_eta'],-LD(6)+np.arange(history,dtype=LD)*dt),'Full historical epoch grid changed')
 require(tuple(int(np.count_nonzero(a['k']<K)) for K in CUTOFFS)==(nodes//4,nodes//2,nodes),'K-prefix node counts changed')
 require(tuple(a['history_quantity_names'])==('q','q_prime','q_second','rho','p','Q0','rho0','p0','current'),'Quantity labels differ')
 require(tuple(a['history_geometry_names'])==('a0','L','L_prime','L_second','L_third','M','M_prime','M_second','M_third','M_fourth'),'Geometry labels differ')


def simpson(F,end,dt):
 require(end%2==0,'Global Simpson endpoint parity')
 return dt/3*(F[0]+F[end]+4*np.sum(F[1:end:2],axis=0,dtype=LD)+2*np.sum(F[2:end:2],axis=0,dtype=LD))


def mode_sums(k,weights,u,w,eta,eps,pi,counts):
 L=-LD(1)/eta;U,W=u/eps,w/eps
 common=-L*W.real-k*W.imag
 r=((2*k*k+3*L*L)*U.real+common)/(2*k)
 p=((2*k*k/3-L*L)*U.real+common)/(2*k)
 f=(3*L*L*L*U.real+L*L*W.real)/k+L*W.imag
 mw=weights*k*k/(2*pi*pi)
 return {key:cut_prefix(mw*value,counts) for key,value in [('R',r),('P',p),('F',f)]}


def flow(ctx,k,weights,stored_u,stored_w,pred_u,pred_w,eta,eps,pi,counts):
 L=-ctx.one/ratio(ctx,eta);epsilon=ratio(ctx,eps);pi_mp=ratio(ctx,pi)
 terms=[];triangle=[]
 for kk,ww,us,ws,up,wp in zip(k,weights,stored_u,stored_w,pred_u,pred_w):
  km=ratio(ctx,kk);weight=ratio(ctx,ww)*km*km/(2*pi_mp*pi_mp)
  du=ctx.mpc(ratio(ctx,us.real)-ratio(ctx,up.real),ratio(ctx,us.imag)-ratio(ctx,up.imag))
  dw=ctx.mpc(ratio(ctx,ws.real)-ratio(ctx,wp.real),ratio(ctx,ws.imag)-ratio(ctx,wp.imag))
  coefficient=(2*km*km+3*L*L)/km
  terms.append(weight*(coefficient*du.real-dw.imag-L*dw.real/km)/(2*epsilon))
  triangle.append(weight*(coefficient*abs(du)+(1+L/km)*abs(dw))/(2*epsilon))
 return [ctx.fsum(terms[:n]) for n in counts],[ctx.fsum(triangle[:n]) for n in counts]


def term_integral(ctx,k,weights,terms,pi,counts):
 p=ratio(ctx,pi)
 converted=[ratio(ctx,v)*ratio(ctx,w)*ratio(ctx,kk)**2/(2*p*p) for v,w,kk in zip(terms,weights,k)]
 return [ctx.fsum(converted[:n]) for n in counts]


def profile_strings(array):return [[str(x) for x in row] for row in array]


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--checkpoint-root',type=Path,required=True)
 parser.add_argument('--freeze-commit',required=True)
 parser.add_argument('--registration-sha256',required=True)
 parser.add_argument('--output-dir',type=Path,required=True)
 args=parser.parse_args()
 require(AUTHENTICATED_BY_WRAPPER=={'checkpoint_root':str(args.checkpoint_root.absolute()),'registration_sha256':args.registration_sha256,'freeze_commit':args.freeze_commit},'Physical core requires authenticated registered wrapper')
 started=ROUTE_STARTED if ROUTE_STARTED is not None else time.monotonic()
 def check():
  require(time.monotonic()-started<=900,'900s complete-route limit exceeded')
  require(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<=262144,'262144KiB complete-route limit exceeded')
 provenance=verify(args);check()
 np.seterr(invalid='raise',divide='raise',over='raise',under='ignore')
 require(np.finfo(LD).nmant==63 and np.dtype(LD).itemsize==16 and np.little_endian,'Pinned native binary80 runtime mismatch')
 eps=LD(EPS_RATIO[0])/LD(EPS_RATIO[1]);pi=LD(PI_RATIO[0])/LD(PI_RATIO[1])
 require(eps.as_integer_ratio()==EPS_RATIO and pi.as_integer_ratio()==PI_RATIO,'Inherited scalar-ratio mismatch')
 args.output_dir.mkdir(parents=True)
 report={'schema_version':1,'old_metric_status':'FAIL','fixed_cases':12,'input_manifest_sha256':provenance['input_manifest_sha256'],'resource_limits':{'wall_seconds':900,'peak_rss_kib':262144},'route':'independent_direct_F_nested_Duhamel','status':'RUNNING','provenance':provenance,
         'interval':['-4.5','-3.5'],'midpoint':'-4','controls':[[a,b] for a,b in CONTROLS],
         'canonical_control':[24,24],'reduction_decimal_digits':[80,100],
         'source_phase_arithmetic':'native x87binary80 longdouble','quadrature_constants':'binary64 leggauss nodes/weights cast to nativeLD',
         'high_precision_scope':'Final exact-ratio transfer and weighted accumulation only; source/phases/localgeometry/timeaccumulation nativeLD.',
         'units':{'profiles_R_P':'a0^4 delta_rho,p/epsilon','F':'dR/deta','baseline_r0_p0':'a0^4 times physical renormalized baseline rho0,p0','contact_F':'L*(contact_R-3contact_P)-3hprime*(r0+p0)','baseline_work_integral':'integral of signed -3hprime*(r0+p0)'},'momentum_target':'Inherited discrete nodes and weights; no continuum-momentum error certificate','rows':[]}
 path=args.output_dir/'diagnostic.json'
 try:
  for setting,steps in [('coarse',128),('fine',256)]:
   archives={}
   for source in ('positive_B','signed_uB'):
    p=args.checkpoint_root/'inputs'/f'metric_modes_{source}_{setting}.npz'
    with np.load(p,allow_pickle=False) as z:archives[source]={name:z[name].copy() for name in z.files}
    validate_arrays(archives[source],setting)
   first=archives['positive_B'];k=first['k'];weights=first['momentum_weights']
   require(np.all(k>0) and np.all(weights>0),'Nonpositive inherited quadrature')
   require(all(np.array_equal(archives[s]['k'],k) and np.array_equal(archives[s]['momentum_weights'],weights) for s in archives),'Source momentum rules differ')
   counts=tuple(int(np.count_nonzero(k<K)) for K in CUTOFFS)
   require(counts[-1]==len(k),'Unexpected modes beyondK256')
   a,b=LD('-4.5'),LD('-3.5');dt=(b-a)/steps
   cache={};outputs={}
   for source_q,ledger_q in CONTROLS:
    if ledger_q not in cache:cache[ledger_q]=prepare(k,weights,steps,ledger_q,counts,pi,a,b,check)
    prep=dict(cache[ledger_q]);local=DenseKernel(k[:1],dt,source_q,ledger_q)
    prep.update(source_q=source_q,s=local.s,sw=local.sw)
    for source,archive in archives.items():
     outputs[(source,source_q,ledger_q)]=trajectory(k,weights,archive['u_1'],archive['w_1'],prep,physical_sources(prep,source),eps,counts,pi,check)
     print(json.dumps({'source':source,'setting':setting,'source_order':source_q,'ledger_order':ledger_q,'stage':'trajectory_done','elapsed_seconds':time.monotonic()-started}),flush=True)
   grid=cache[24]['grid']
   for source,archive in archives.items():
    ids=[int(np.flatnonzero(archive['history_eta']==t)[0]) for t in grid]
    start,end=ids[0],ids[-1]
    require(end-start==steps and all(y-x==1 for x,y in zip(ids,ids[1:])),'Inherited activegrid differs')
    storedR=archive['history_values'][ids,:,3];storedP=archive['history_values'][ids,:,4];storedF=archive['history_ledger_integrand'][ids]
    Sab=simpson(archive['history_ledger_integrand'],end,dt)-simpson(archive['history_ledger_integrand'],start,dt)
    Sreset=simpson(storedF,steps,dt)
    if setting=='fine':
     full_double=archive['history_ledger_integrand'][::2]
     Sdouble=simpson(full_double,end//2,2*dt)-simpson(full_double,start//2,2*dt)
    canonical=outputs[(source,24,24)]
    recomputed={}
    for knot,ix,index in [('initial',0,1),('midpoint',steps//2,2),('endpoint',steps,3)]:
     mode=mode_sums(k,weights,archive[f'u_{index}'],archive[f'w_{index}'],grid[ix],eps,pi,counts)
     recomputed[knot]={key:mode[key]+canonical['contact_'+key][ix] for key in ('R','P','F')}
    precision_results={}
    for dps in (80,100):
     ctx=mp.mp.clone();ctx.dps=dps
     C={}
     M_discrete=term_integral(ctx,k,weights,1/(2*k),pi,counts)
     for sq,lq in CONTROLS:
      output=outputs[(source,sq,lq)]
      modeI=term_integral(ctx,k,weights,output['modal_I'],pi,counts)
      final_flow,final_triangle=flow(ctx,k,weights,archive['u_3'],archive['w_3'],output['u_end'],output['w_end'],b,eps,pi,counts)
      mid_flow,mid_triangle=flow(ctx,k,weights,archive['u_2'],archive['w_2'],output['u_mid'],output['w_mid'],grid[steps//2],eps,pi,counts)
      C[(sq,lq)]={'I_mode':modeI,'I_contacts':[ratio(ctx,x) for x in output['contact_I']],
                  'I_ab':[modeI[j]+ratio(ctx,x) for j,x in enumerate(output['contact_I'])],
                  'baseline_work_integral':[ratio(ctx,x) for x in output['work_I']],
                  'J_Lg':[ratio(ctx,output['J_Lg'])]*3,'M_discrete':M_discrete,
                  'Lg_projection':[m*ratio(ctx,output['J_Lg']) for m in M_discrete],
                  'E_flow':final_flow,'triangle_bound':final_triangle,'midpoint_E_flow':mid_flow,'midpoint_triangle_bound':mid_triangle}
     precision_results[str(dps)]=[]
     for j,K in enumerate(CUTOFFS):
      fields={key:value[j] for key,value in C[(24,24)].items()}
      fields.update(S_ab=ratio(ctx,Sab[j]),S_direct_reset=ratio(ctx,Sreset[j]),
                    DeltaR=ratio(ctx,storedR[-1,j])-ratio(ctx,storedR[0,j]))
      fields['D_S']=fields['DeltaR']-fields['S_ab'];fields['D_cont']=fields['DeltaR']-fields['I_ab'];fields['E_Q']=fields['I_ab']-fields['S_ab']
      fields['decomposition_error']=fields['D_S']-fields['D_cont']-fields['E_Q']
      fields['E_operator']=(ratio(ctx,storedR[-1,j])-ratio(ctx,recomputed['endpoint']['R'][j]))-(ratio(ctx,storedR[0,j])-ratio(ctx,recomputed['initial']['R'][j]))
      fields['E_evolution_ledger']=fields['D_cont']-fields['E_flow']-fields['E_operator']
      fields['midpoint_E_operator']=ratio(ctx,storedR[steps//2,j])-ratio(ctx,recomputed['midpoint']['R'][j])
      fields['midpoint_density_gap']=ratio(ctx,storedR[steps//2,j])-ratio(ctx,canonical['R'][steps//2,j])
      fields['midpoint_projection_closure']=fields['midpoint_density_gap']-fields['midpoint_E_flow']-fields['midpoint_E_operator']
      fields.update(initial_R_reconstruction_error=ratio(ctx,storedR[0,j])-ratio(ctx,recomputed['initial']['R'][j]),endpoint_R_reconstruction_error=ratio(ctx,storedR[-1,j])-ratio(ctx,recomputed['endpoint']['R'][j]))
      if setting=='fine':fields.update(S_ab_double_global=ratio(ctx,Sdouble[j]),S_native_minus_double=fields['S_ab']-ratio(ctx,Sdouble[j]))
      fields['forcing_control_gap']=C[(24,24)]['I_ab'][j]-C[(16,24)]['I_ab'][j]
      fields['ledger_control_gap']=C[(16,24)]['I_ab'][j]-C[(16,16)]['I_ab'][j]
      fields['joint_control_gap']=C[(24,24)]['I_ab'][j]-C[(16,16)]['I_ab'][j]
      fields.update({key+'_profile_max_error':ratio(ctx,np.max(np.abs(canonical[key][:,j]-stored))) for key,stored in [('R',storedR[:,j]),('P',storedP[:,j]),('F',storedF[:,j])]})
      precision_results[str(dps)].append({'K':K,'fields':{key:serial(ctx,value) for key,value in fields.items()},
                                         'controls':{str(sq)+'_'+str(lq):{key:serial(ctx,value[j]) for key,value in c.items()} for (sq,lq),c in C.items()}})
     check()
    for j,K in enumerate(CUTOFFS):
     row={'source':source,'setting':setting,'K':K,'times':[serial_native(x) for x in grid],
          'profiles':{key:[serial_native(x) for x in canonical[key][:,j]] for key in ('R','P','F')},
          'contact_profiles':{key:[serial_native(x) for x in canonical['contact_'+key][:,j]] for key in ('R','P','F')},
          'baseline_profiles':{key:[serial_native(x) for x in cache[24]['grid_geometry'][key][:,j]] for key in ('r0','p0')},
          'stored_profiles':{key:[serial_native(x) for x in values[:,j]] for key,values in [('R',storedR),('P',storedP),('F',storedF)]},
          'control_profiles':{str(sq)+'_'+str(lq):{key:[serial_native(x) for x in outputs[(source,sq,lq)][key][:,j]] for key in ('R','P','F')} for sq,lq in CONTROLS},
          'anchor_recomputed':{knot:{key:serial_native(value[j]) for key,value in values.items()} for knot,values in recomputed.items()},
          'precisions':{key:value[j] for key,value in precision_results.items()}}
     report['rows'].append(row)
    path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');check()
  require(len(report['rows'])==12,'Incomplete case universe')
  report['status']='COMPLETED_DIAGNOSTIC_WITHOUT_RECLASSIFYING_OLD_FAIL'
  report['serialization_max_error_by_dps']={key:str(value) for key,value in SERIAL_MAX.items()}
  report['native_profile_serialization_max_error']=str(NATIVE_SERIAL_MAX)
  report['native_profile_serialization_scope']='Every native scientific/profile scalar preserves an exact binary80 roundtrip; error relative to original binary80 rational reported separately.'
  require(max(SERIAL_MAX.values())<=Fraction('0.000000000001') and NATIVE_SERIAL_MAX<=Fraction('0.000000000001'),'Serialization exceeds1e-12')
  report['limits']=['Empirical local forcing/ledger controls only; no rigorous continuum or inheritedinitialstate error bound.',
                   'Root validator classifies numerical consistency without altering the historicalscientificFAIL.',
                   'All source and phase arithmetic remainsnativeLD; MP80/100 checks finalreductiononly.']
 except Exception as exc:
  report['status']='FATAL_EXECUTION_FAILURE';report['exception']={'type':type(exc).__name__,'message':str(exc)}
  raise
 finally:
  report['resources']={'elapsed_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
  path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
  check()
 print(json.dumps({'status':report['status'],'resources':report['resources'],'rows':len(report['rows'])}),flush=True)
if __name__=='__main__':main()
