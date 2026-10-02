#!/usr/bin/env python3
"""Render saved registered stress points; no new physical source evaluation."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR','/tmp/hdblast-stress-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/hdblast-stress-cache')
Path(os.environ['XDG_CACHE_HOME']).mkdir(parents=True,exist_ok=True)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator, ScalarFormatter

COLORS={'positive_B':'#176b91','signed_uB':'#b44a2a'}
LABELS={'positive_B':'Positive pulse B','signed_uB':'Signed pulse uB'}
QUANTITIES=('rho','p','current')
TITLES={'rho':'Density response','p':'Pressure response','current':'Scalar current response'}
UNITS={'rho':r'$R=a^4\delta\rho/\epsilon$','p':r'$P=a^4\delta p/\epsilon$',
       'current':r'$J=a^2\delta j/\epsilon$'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(figure,output,name):
    for extension in ('svg','png','pdf'):
        figure.savefig(output/(name+'.'+extension),dpi=200,bbox_inches='tight',
                       metadata={'Creator':'Post-run registered stress data renderer'})
    plt.close(figure)


def log_points(axis,x,y,*args,**kwargs):
    points=[(a,b) for a,b in zip(x,y) if b>0]
    if points:axis.semilogy([p[0] for p in points],[p[1] for p in points],*args,**kwargs)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise RuntimeError('Refusing to overwrite stress figures')
    data=json.loads(args.summary.read_text())
    if data['status']!='registered_matched_stress_calibration_passed' or data['new_physical_source_response_or_mode_evaluations']!=0:
        raise RuntimeError('Unexpected summary status or presentation scope')
    groups={source:[r for r in data['registered_rows'] if r['source']==source] for source in COLORS}
    for source,rows in groups.items():
        if [r['eta'] for r in rows]!=[-5.5,-4.5,-4,-3.5,-2.5,-1.5]:raise RuntimeError('Unexpected registered grid')
    args.output.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.size':10,'axes.labelsize':11,'axes.titlesize':12,
                         'figure.titlesize':16,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.hashsalt':'matched-stress-registration-4a5dad8'})
    figure,axes=plt.subplots(2,3,figsize=(14,7.1),sharex=True)
    for row_index,(source,rows) in enumerate(groups.items()):
        times=[r['eta'] for r in rows]
        color=COLORS[source]
        for column,name in enumerate(QUANTITIES):
            axis=axes[row_index,column]
            axis.axvspan(-5,-3,color='#e9edf0',zorder=0)
            axis.axhline(0,color='#75818a',linewidth=.8)
            axis.plot(times,[r['continuum_values'][name] for r in rows],'o--',color=color,
                      linewidth=1.2,markersize=5.5,zorder=3)
            axis.plot(times,[r['finite_k'][-1]['direct_mode_values'][name] for r in rows],
                      'D',markerfacecolor='none',markeredgecolor='#20262b',markersize=7.5,zorder=4)
            axis.set_title(LABELS[source]+' · '+TITLES[name])
            axis.set_ylabel(UNITS[name])
            axis.set_xticks(times)
            axis.grid(axis='y',alpha=.2)
            if row_index==1:axis.set_xlabel(r'Conformal time $\eta$ ($H=1$)')
    legend=[Line2D([0],[0],color='#58646c',marker='o',linestyle='--',label='Removed-cutoff logarithmic / matched stress route'),
            Line2D([0],[0],color='#20262b',marker='D',markerfacecolor='none',linestyle='none',label='Independent direct modes, K=256')]
    figure.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.965),ncol=2,frameon=False,fontsize=10)
    figure.suptitle('Matched stress and current at the twelve registered points',y=1.02)
    figure.text(.5,-.015,'Shading: exact source support −5 < η < −3. Lines guide the eye; open and filled markers nearly coincide.\n'
                'First-order perturbations about the fixed BD baseline: δρ=εR/a⁴, δp=εP/a⁴, δj=εJ/a²; ε=10⁻⁴. No heating or stability inference.',
                ha='center',fontsize=9,color='#4f5962')
    figure.tight_layout(rect=(0,.025,1,.925),h_pad=2.0,w_pad=2.0)
    save(figure,args.output,'stress_response')

    figure,axes=plt.subplots(2,3,figsize=(14,7.1),sharex=True)
    for row_index,(source,rows) in enumerate(groups.items()):
        rows=[r for r in rows if r['eta']>-5]
        times=[r['eta'] for r in rows]
        for column,name in enumerate(QUANTITIES):
            axis=axes[row_index,column];points=[r['finite_k'][-1] for r in rows]
            log_points(axis,times,[p['analytic_omitted_band_bounds'][name] for p in points],'D--',color='#56626b',markersize=5,linewidth=1)
            log_points(axis,times,[p['primary_to_continuum_differences'][name] for p in points],'o',color=COLORS[source],markersize=6)
            log_points(axis,times,[p['direct_to_continuum_differences'][name] for p in points],'^',color='#20262b',markerfacecolor='none',markersize=8)
            axis.set_title(LABELS[source]+' · '+name)
            axis.set_ylabel('Absolute gap or bound in '+{'rho':'R','p':'P','current':'J'}[name])
            axis.set_xticks(times)
            axis.grid(axis='y',which='major',alpha=.22)
            if row_index==1:axis.set_xlabel(r'Conformal time $\eta$')
    handles=[Line2D([0],[0],color='#56626b',marker='D',linestyle='--',label='Analytic omitted-band bound'),
             Line2D([0],[0],color='#58646c',marker='o',linestyle='none',label='Primary finite-K / continuum gap'),
             Line2D([0],[0],color='#20262b',marker='^',markerfacecolor='none',linestyle='none',label='Direct-mode K / continuum gap')]
    figure.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.965),ncol=3,frameon=False,fontsize=9.5)
    figure.suptitle('K=256: observed discrepancies and constructive UV bounds',y=1.02)
    figure.text(.5,-.015,'Only the omitted momentum band is analytically enclosed. Quadrature, refinement and arithmetic remain separate numerical uncertainties.\n'
                'Exact zeros are omitted on logarithmic axes. Lines connect registered points only.',ha='center',fontsize=9,color='#4f5962')
    figure.tight_layout(rect=(0,.025,1,.925),h_pad=2.0,w_pad=2.0)
    save(figure,args.output,'stress_uv_bounds')

    figure,axes=plt.subplots(1,2,figsize=(11.4,4.8))
    cutoffs=[64,128,256]
    styles={'rho':('#176b91','o'),'p':('#b44a2a','s')}
    for name,(color,marker) in styles.items():
        cross=[max(r['finite_k'][i]['direct_to_primary_differences'][name] for r in data['registered_rows']) for i in range(3)]
        refine=[max(r['finite_k'][i]['mode_refinement_differences'][name] for r in data['registered_rows']) for i in range(3)]
        log_points(axes[0],cutoffs,cross,color=color,marker=marker,label=name+': direct / primary')
        log_points(axes[0],cutoffs,refine,color=color,marker=marker,markerfacecolor='none',linestyle='--',label=name+': coarse / fine')
    ward=[max(r['finite_k'][i]['ward']['endpoint_difference'] for r in data['registered_rows']) for i in range(3)]
    ward_ref=[max(r['finite_k'][i]['ward']['refinement_difference'] for r in data['registered_rows']) for i in range(3)]
    trace=[max(r['gap'] for r in data['trace_comparisons'] if r['K']==K) for K in cutoffs]
    log_points(axes[1],cutoffs,ward,'o-',color='#176b91',label='Independent Ward ledger / direct density')
    log_points(axes[1],cutoffs,ward_ref,'o--',color='#176b91',markerfacecolor='none',label='Ward ledger coarse / fine')
    log_points(axes[1],cutoffs,trace,'s-',color='#b44a2a',label='Direct finite-band trace residual')
    for axis in axes:
        axis.set_xticks(cutoffs);axis.xaxis.set_minor_locator(NullLocator())
        axis.xaxis.set_major_formatter(ScalarFormatter());axis.set_xlabel('Registered comoving cutoff K')
        axis.grid(axis='y',alpha=.23);axis.legend(frameon=False,fontsize=8.5)
        axis.set_ylabel('Maximum absolute difference over registered observations')
    axes[0].set_title('Direct stress agreement and refinement')
    axes[1].set_title('Independent Ward ledger and matched trace')
    figure.suptitle('Finite-band numerical checks',y=1.01)
    figure.text(.5,-.02,'Residuals use the registered normalization. These checks remain numerical evidence, separate from the analytic UV tail theorem.\n'
                'Ward uses physical Q0,K and independently saved direct-stress histories; the density was not defined by the ledger.',ha='center',fontsize=9,color='#4f5962')
    figure.tight_layout(rect=(0,.025,1,1),w_pad=3)
    save(figure,args.output,'stress_checks')
    files={p.name:sha(p) for p in sorted(args.output.iterdir()) if p.suffix in ('.svg','.png','.pdf')}
    (args.output/'FIGURE_PROVENANCE.json').write_text(json.dumps({'summary_sha256':sha(args.summary),'renderer_sha256':sha(__file__),
        'matplotlib':matplotlib.__version__,'postrun_only':True,'new_physical_evaluations':0,'figures_sha256':files},indent=2,sort_keys=True)+'\n')
    print('Rendered',len(files),'figures from',len(data['registered_rows']),'saved registered rows')


if __name__=='__main__':main()
