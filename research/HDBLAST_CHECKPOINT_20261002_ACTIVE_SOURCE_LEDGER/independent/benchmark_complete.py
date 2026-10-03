"""Complete fixed-size fabricated pipeline, including reduction and serialization.
No capsule, actual B source, actual initial state, or registered array is read.
"""
import hashlib,json,time,resource,math
from pathlib import Path
import numpy as np
import mpmath as mp
from kernel import LD,CD,forcing,DenseKernel
from engine import prepare,trajectory
here=Path(__file__).resolve().parent


def ratio(ctx,x):
 n,d=x.as_integer_ratio();return ctx.mpf(n)/ctx.mpf(d)


def fake_jet(t,power):
 z=np.asarray(t,dtype=LD)+LD(8)
 return np.stack([LD(math.factorial(power)//math.factorial(power-n))*z**(power-n) if n<=power else np.zeros_like(z) for n in range(6)])*LD('0.00017')


def sources(prep,power):
 grid=prep['grid'];dt=prep['dt'];t=prep['t'];s=prep['s']
 nested=grid[:-1,None,None]+t[None,:,None]*s[None,None,:]
 end=grid[:-1,None]+dt*s[None,:]
 return {'grid':fake_jet(grid,power),'dense':fake_jet(prep['dense'],power),
         'nested_g':forcing(nested,fake_jet(nested,power)),
         'end_g':forcing(end,fake_jet(end,power))}


def main():
 started=time.monotonic();rows=[]
 def check():
  if time.monotonic()-started>900:raise RuntimeError('Complete fabricated route exceeds 900s')
  if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>262144:raise RuntimeError('Complete fabricated route exceeds memory')
 for name,nodes,steps in [('coarse',8192,128),('fine',16384,256)]:
  k=np.linspace(LD('0.0009'),LD('255.999'),nodes,dtype=LD)
  weights=np.full(nodes,LD(256)/nodes,dtype=LD)
  counts=tuple(int(np.count_nonzero(k<K)) for K in (64,128,256))
  # Fabricated constants and interval differ from every physical study input.
  pi=LD('3.2');eps=LD('0.003');a=LD('-8.5');b=LD('-7.5')
  u0=np.full(nodes,CD('0.00012')+CD('0.00031j'),dtype=CD)
  w0=np.full(nodes,CD('0.000017')-CD('0.000013j'),dtype=CD)
  caches={}
  for source_q,q in ((16,16),(16,24),(24,24)):
   ts=time.monotonic()
   if q not in caches:caches[q]=prepare(k,weights,steps,q,counts,pi,a,b,check)
   prep=dict(caches[q])
   sk=DenseKernel(k[:1],prep['dt'],source_q,q)
   prep.update(source_q=source_q,s=sk.s,sw=sk.sw)
   print(json.dumps({'setting':name,'source_q':source_q,'q':q,'stage':'geometry_ready','elapsed_seconds':time.monotonic()-started}),flush=True)
   for power in (4,5):
    result=trajectory(k,weights,u0,w0,prep,sources(prep,power),eps,counts,pi,check)
    contexts={}
    for dps in (80,100):
     ctx=mp.mp.clone();ctx.dps=dps
     rp=ratio(ctx,pi);terms=[ratio(ctx,r)*ratio(ctx,w)*ratio(ctx,kk)**2/(2*rp*rp) for r,w,kk in zip(result['modal_I'],weights,k)]
     contexts[str(dps)]={str(K):ctx.nstr(ctx.fsum(terms[:n])+ratio(ctx,result['contact_I'][j]),55) for j,(K,n) in enumerate(zip((64,128,256),counts))}
    row={'setting':name,'source_q':source_q,'q':q,'fabricated_power':power,'precisions':contexts,
         'profiles':{key:[[str(float(x)) for x in r] for r in result[key]] for key in ('R','P','F')},
         'finite':all(np.all(np.isfinite(result[key])) for key in ('R','P','F','modal_I','contact_I'))}
    rows.append(row);check()
    print(json.dumps({'setting':name,'source_q':source_q,'q':q,'fabricated_power':power,'stage':'trajectory_reduced','elapsed_seconds':time.monotonic()-started}),flush=True)
 report={'status':'PASS_COMPLETE_FABRICATED_PIPELINE','physical_evaluations':0,'actual_capsule_payloads_read':False,
         'actual_B_source_evaluated':False,'all_four_histories_three_controls_both_reduction_contexts':True,
         'node_counts':{'coarse':8192,'fine':16384},'step_counts':{'coarse':128,'fine':256},
         'source_and_ledger_controls':[[16,16],[16,24],[24,24]],'quadrature_constants':'binary64 leggauss nodes/weights cast to nativeLD',
         'scope':'Geometry cache + all12 trajectories + profiles + twoMPcontexts + serialization; inputdecoding/hashverification cost is represented only by synthetic-byte hashing below.',
         'rows':rows,'elapsed_seconds_before_write':time.monotonic()-started,
         'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 # Bound fake decoding/hash workload to the real input order of magnitude.
 data=np.arange(1024*1024,dtype=np.uint8).tobytes()
 for j in range(16):hashlib.sha256(data).hexdigest()
 serialized=json.dumps(report,indent=2,allow_nan=False)+'\n'
 (here/'SYNTHETIC_COMPLETE_PIPELINE.json').write_text(serialized)
 final={'status':report['status'],'physical_evaluations':0,'elapsed_seconds_including_write':time.monotonic()-started,
        'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'output_bytes':len(serialized.encode()),'output_sha256':hashlib.sha256(serialized.encode()).hexdigest(),'all_rows_finite':all(r['finite'] for r in rows)}
 (here/'SYNTHETIC_COMPLETE_PIPELINE_RECEIPT.json').write_text(json.dumps(final,indent=2)+'\n')
 check();print(json.dumps(final,indent=2),flush=True)
if __name__=='__main__':main()
