"""Fabricated algebra/operator/flow tests; never calls the study source."""
import json,time,resource
from pathlib import Path
import numpy as np
import mpmath as mp
from kernel import LD,CD,full_subtractions,baseline,DenseKernel
from specialized_contacts import directional_subtractions
from moment_contacts import weighted_subtractions
from engine import geometry
from run_independent import flow


def need(test,message):
 if not test:raise RuntimeError(message)


def main():
 started=time.monotonic()
 k=np.array(['0.11','0.7','2.1','17.9','251'],dtype=LD)
 L=np.array(['0.21','0.29'],dtype=LD)[None,:]
 jet=np.stack([np.full((1,2),LD('0.0037')/(n+1),dtype=LD) for n in range(6)])
 generic=full_subtractions(k[:,None],L,jet)
 special=directional_subtractions(k[:,None],L,jet)
 max_gap=max(np.max(np.abs(a-b)) for a,b in zip(generic,special))
 need(max_gap<LD('0.000000000000001'),'IndependentPair specialization synthetic mismatch')
 weights=np.ones_like(k)/LD(3);pi=LD('3.2')
 times=-LD(1)/L.ravel()
 counts=(2,4,5)
 geo=geometry(k,weights,times,counts,pi)
 dr,dp=weighted_subtractions(L.T,geo['moments'],jet[:,0,:,None])
 measure=weights*k*k/(2*pi*pi)
 moment_gap=max(np.max(np.abs(dr[:,j]-np.sum((generic[0]*measure[:,None])[:n],axis=0,dtype=LD))) for j,n in enumerate(counts))
 moment_gap=max(moment_gap,max(np.max(np.abs(dp[:,j]-np.sum((generic[1]*measure[:,None])[:n],axis=0,dtype=LD))) for j,n in enumerate(counts)))
 need(moment_gap<LD('0.000000000000001'),'Discrete moment regrouping synthetic mismatch')
 # Full-size final projections are intentionally synthetic and count toward cost.
 projection_rows=[]
 for setting,n in [('coarse',8192),('fine',16384)]:
  kk=np.linspace(LD('0.0009'),LD('255.999'),n,dtype=LD);ww=np.full(n,LD(256)/n,dtype=LD)
  count=tuple(int(np.count_nonzero(kk<K)) for K in (64,128,256))
  pred_u=np.full(n,CD('0.00012')+CD('0.00031j'),dtype=CD);pred_w=np.full(n,CD('0.000017')-CD('0.000013j'),dtype=CD)
  stored_u=pred_u+CD('0.00000000000003')-CD('0.00000000000001j');stored_w=pred_w+CD('0.00000000000007')+CD('0.00000000000004j')
  for source in (0,1):
   for dps in (80,100):
    ctx=mp.mp.clone();ctx.dps=dps
    for q in (16,24):
     for eta in (LD('-8'),LD('-7.5')):
      signed,triangle=flow(ctx,kk,ww,stored_u,stored_w,pred_u,pred_w,eta,LD('0.003'),pi,count)
      need(all(abs(a)<=b for a,b in zip(signed,triangle)),'Signedflow projection exceeds triangle bound')
      projection_rows.append({'setting':setting,'fabricated_source':source,'dps':dps,'order':q,'eta':str(eta),'signed':[ctx.nstr(x,55) for x in signed],'triangle':[ctx.nstr(x,55) for x in triangle]})
 report={'status':'PASS_FABRICATED_OPERATOR_AND_FULLSIZE_PROJECTION_CHECKS','physical_evaluations':0,'study_source_calls':0,
         'generic_vs_specialized_max_gap':str(max_gap),'moment_regrouping_max_gap':str(moment_gap),
         'fullsize_projection_rows':projection_rows,'elapsed_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 name='SYNTHETIC_OPTIMIZED.json' if __debug__ is False else 'SYNTHETIC_NORMAL.json'
 Path(__file__).with_name(name).write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k!='fullsize_projection_rows'},indent=2))
if __name__=='__main__':main()
