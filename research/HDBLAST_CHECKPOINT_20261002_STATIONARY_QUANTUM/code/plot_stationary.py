#!/usr/bin/env python3
"""Plot recorded stationary points; no sources, integration or new acceptance tests."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--summary',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    data=json.loads(args.summary.read_text())
    args.output.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','savefig.bbox':'tight'})
    fig,axes=plt.subplots(1,2,figsize=(10.3,4.4),layout='constrained')
    styles={-1:('#2667a8','s'),0:('#373737','o'),1:('#c25425','^')}
    reference={'H2':5.924014794328956e-5,'eta':8.405268308670538e-5}
    for ax,field in zip(axes,['H2','eta']):
        ax.set_xscale('log')
        ax.set_yscale('log' if field=='H2' else 'symlog',**({} if field=='H2' else {'linthresh':1e-10}))
        ax.set_xlabel(r'Quantum source amplitude $\gamma$')
        ax.set_ylabel(r'$\Delta H_b^2/H_0^2$' if field=='H2' else r'$\Delta\phi_b/|\phi_0-1|$')
        ax.grid(True,alpha=.2)
        for b,(color,marker) in styles.items():
            rows=sorted([s for s in data['shifts'] if s['b']==b],key=lambda s:s['gamma'])
            x=[r['gamma'] for r in rows]
            y=[r['observables'][field]['delta_from_own_solved_control']/reference[field] for r in rows]
            ax.plot(x,y,color=color,alpha=.45,linewidth=.8)
            for i,row in enumerate(rows):
                resolved=row['observables'][field]['shift_resolved']
                ax.scatter([x[i]],[y[i]],color=color,marker=marker,s=36,
                           facecolors=color if resolved else 'none',linewidths=1.1,
                           label=f'$b={b:+d}$' if i==0 else None,zorder=3)
        ax.axvspan(1,1.7e6,color='#dfd9ce',alpha=.3,zorder=0)
        ax.set_xlim(.006,1.7e6)
    axes[0].legend(frameon=False,loc='upper left')
    axes[1].axhline(0,color='#555555',linewidth=.7)
    fig.suptitle('Registered stationary quantum shell responses',fontsize=13)
    fig.supxlabel('Shaded sampled amplitudes fail the declared hierarchy screen. Lines only guide the eye; hollow marker: unresolved shift.',fontsize=9)
    for suffix in ['svg','png','pdf']:
        fig.savefig(args.output/f'stationary_responses.{suffix}',dpi=180)
    print(args.output)

if __name__=='__main__':
    main()
