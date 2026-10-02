#!/usr/bin/env python3
"""Render saved metric response evidence; no source or mode evaluations."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--summary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise RuntimeError('Fresh figure output required')
    if matplotlib.__version__!='3.10.1':raise RuntimeError('Pinned Matplotlib3.10.1 required')
    data=json.loads(args.summary.read_text())
    if data['status']!='registered_metric_calibration_passed' or not data['postrun_only']:raise RuntimeError('Saved passing metric summary required')
    args.output.mkdir(parents=True)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180,'svg.hashsalt':'registered-metric-response'})
    sources=('positive_B','signed_uB');colors=('#235789','#d95f02');figures={}
    fig,axes=plt.subplots(3,2,figsize=(9,8),sharex=True,constrained_layout=True)
    for col,source in enumerate(sources):
        rows=sorted([r for r in data['registered_rows'] if r['source']==source],key=lambda r:r['eta'])
        for j,q in enumerate(('q','rho','p')):
            ax=axes[j,col];x=[r['eta'] for r in rows]
            ax.axvspan(-5,-3,color='0.9',label='source support' if j==0 else None)
            ax.plot(x,[r['continuum_values'][q] for r in rows],'o--',color=colors[col],label='primary continuum formula')
            ax.plot(x,[next(f for f in r['finite_k'] if f['K']==256)['direct_mode_values'][q] for r in rows],'x',color='black',label='direct K=256')
            ax.axhline(0,color='0.5',lw=.5);ax.set_ylabel(q+' / registered scaling')
            if j==0:ax.set_title(source)
            if j==2:ax.set_xlabel('conformal time eta')
            if j==0 and col==0:ax.legend(fontsize=7)
    fig.suptitle('Homogeneous linear metric response: fixed x=r=2 and incoming BD state')
    figures['metric_response']=fig
    fig,axes=plt.subplots(1,3,figsize=(10,3.4),constrained_layout=True)
    points=[f for r in data['registered_rows'] for f in r['finite_k']]
    for ax,q in zip(axes,('q','rho','p')):
        ks=(64,128,256)
        gaps=[max(f['direct_to_continuum_differences'][q] for f in points if f['K']==K) for K in ks]
        ax.loglog(ks,np.maximum(gaps,1e-20),'o-',label='max observed difference')
        tails=data.get('continuum_tail_comparisons',[])
        bounds=[]
        for K in ks:
            matching=[t for t in tails if t['K']==K and t['quantity']==q]
            values=[t['analytic_tail_upper_bound'] for t in matching]
            values=[v for v in values if isinstance(v,(int,float))]
            bounds.append(max(values) if values else None)
        if all(v is not None for v in bounds):ax.loglog(ks,np.maximum(bounds,1e-20),'s--',label='omitted-band envelope')
        ax.set_title(q);ax.set_xlabel('fixed comoving K');ax.set_ylabel('absolute scaled difference / envelope');ax.grid(alpha=.2)
    axes[0].legend(fontsize=7)
    fig.suptitle('Cutoff evidence; envelopes exclude finite integration and arithmetic error')
    figures['metric_cutoff']=fig
    fig,axes=plt.subplots(1,3,figsize=(10,3.4),constrained_layout=True)
    names=('q','q_prime','q_second','rho','p')
    metrics=data['metrics']
    axes[0].bar(range(5),np.maximum([metrics['maximum_cross_route_difference'][q] for q in names],1e-20),color=colors[0]);axes[0].set_xticks(range(5),names,rotation=35);axes[0].set_yscale('log');axes[0].set_title('Primary / direct agreement')
    ref=metrics['maximum_refinement_difference'];axes[1].bar(range(5),np.maximum([ref[q] for q in names],1e-20),color=colors[1]);axes[1].set_xticks(range(5),names,rotation=35);axes[1].set_yscale('log');axes[1].set_title('Direct coarse / fine refinement')
    axes[2].bar([0,1],np.maximum([metrics['maximum_fine_ward_endpoint_difference'],metrics['maximum_ward_refinement_difference']],1e-20),color=['#4c956c','#808080']);axes[2].set_xticks([0,1],['endpoint','refinement']);axes[2].set_yscale('log');axes[2].set_title('Independent Ward ledger')
    for ax in axes:ax.set_ylabel('maximum absolute residual');ax.grid(axis='y',alpha=.2)
    fig.suptitle('Empirical checks of finite-K calibration; displayed zero floor is1e-20')
    figures['metric_checks']=fig
    outputs={}
    for name,fig in figures.items():
        for ext in ('png','svg','pdf'):
            target=args.output/(name+'.'+ext);fig.savefig(target);outputs[target.name]=sha(target)
        plt.close(fig)
    provenance={'summary_sha256':sha(args.summary),'renderer_sha256':sha(__file__),'matplotlib':matplotlib.__version__,'postrun_only':True,'new_physical_evaluations':0,'figures_sha256':outputs}
    (args.output/'FIGURE_PROVENANCE.json').write_text(json.dumps(provenance,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'figure_artifacts':len(outputs),'new_physical_evaluations':0}))

if __name__=='__main__':main()
