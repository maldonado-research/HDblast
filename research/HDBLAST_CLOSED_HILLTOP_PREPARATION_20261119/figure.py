"""Compare smooth-field limits with the waiting time of the homogeneous solutions."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parent


def main(output):
    limits=json.loads((ROOT/'minimax/MINIMAX.json').read_text())['branches']
    modes=json.loads((ROOT/'harmonics/HARMONICS.json').read_text())['branches']
    blue, amber, dark, grey = '#146a8a', '#b55b22', '#303e49', '#bec5ca'
    fig,axes=plt.subplots(1,2,figsize=(10,4.4),sharey=True)
    fig.subplots_adjust(left=.1,right=.965,top=.82,bottom=.27,wspace=.16)
    for ax, (key, data) in zip(axes,limits.items()):
        position=[2,1,0]
        ax.set_title('Minimum at '+key.split('=')[1],fontsize=13,fontweight='bold',pad=10)
        for y,(det,entry) in zip(position,data.items()):
            horizon=entry['band_results']['n2']['rms_horizons']['0.01']
            needed=entry['time_required']
            force=modes[key]['detunings'][det]['x_n2_optimal_spatial_rms_at_force_radius']
            ax.hlines(y,0,needed,color=grey,lw=3)
            ax.hlines(y,0,horizon,color=blue,lw=7)
            ax.plot(horizon,y,'o',color=blue,ms=6)
            ax.plot(force,y,'|',color=amber,ms=15,mew=2)
            ax.plot(needed,y,'D',color=dark,ms=6)
            ax.text(needed+.37,y-.015,f'{needed:.2f}',va='center',fontsize=9,color=dark)
            ax.text(horizon+.04,y+.14,f'{horizon:.2f}',va='bottom',fontsize=9,color=blue)
        ax.set(xlim=(-.3,18.6),ylim=(-.45,2.55),xticks=[0,4,8,12,16],
               xlabel=r'Elapsed time from bounce ($H_{\rm top}^{-1}$)')
        ax.grid(axis='x',lw=.6,color='#e2e5e7')
        ax.spines[['right','top','left']].set_visible(False)
        ax.tick_params(axis='y',length=0)
    axes[0].set(yticks=[2,1,0],yticklabels=[r'$10^{-3}$',r'$10^{-10}$',r'$10^{-20}$'],ylabel='Detuning')
    fig.suptitle('How long can the field stay smooth?',fontsize=17,fontweight='bold',y=.98)
    handles=[Line2D([0],[0],color=blue,lw=5,label='Latest time below 0.01 field RMS'),
             Line2D([0],[0],marker='|',color=amber,ms=12,mew=2,linestyle='none',label='Using the 10% force tolerance'),
             Line2D([0],[0],marker='D',color=dark,ms=6,linestyle='none',label='Time needed by homogeneous launch')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.04),ncol=3,frameon=False,fontsize=9)
    fig.savefig(output,dpi=170,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True,type=Path)
    main(p.parse_args().output)
