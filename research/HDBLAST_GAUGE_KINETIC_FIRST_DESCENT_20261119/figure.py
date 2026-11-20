"""Recreate the outgoing photon and dissipative radiation comparisons."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT=Path(__file__).resolve().parent


def draw():
    finite=json.loads((ROOT/'feedback/FEEDBACK.json').read_bytes())['branches']
    fluid=json.loads((ROOT/'fluid/FLUID.json').read_bytes())['branches']
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.7),sharey=True)
    for (key,branch),color in zip(finite.items(),['#1f6981','#b05227']):
        cutoffs=[16,24,32,48,64,96,128,192]
        values=[branch['runs'][f'{n}/3']['end']['fraction_of_matter'] for n in cutoffs]
        label=f"Minimum at {key.split('=')[1]}"
        axes[0].plot(cutoffs,values,'o-',color=color,lw=1.8,ms=4.3,label=label)
        rows=fluid[key]['scan']
        rates=sorted(map(float,rows))
        vals=[rows[str(int(r)) if r.is_integer() and str(int(r)) in rows else str(r)]['maximum_before_turnaround']['fraction'] for r in rates]
        axes[1].plot(rates,vals,'o-',color=color,lw=1.8,ms=3)
    for ax in axes:
        ax.axhline(.5,color='#434343',ls='--',lw=1)
        ax.set_ylim(0,.82)
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.grid(axis='y',alpha=.18)
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_title('Photon share at first field minimum',loc='left',fontsize=11)
    axes[1].set_title('Highest fluid radiation share before turnaround',loc='left',fontsize=11)
    axes[0].set_xlabel('Highest photon harmonic retained, $K$')
    axes[1].set_xlabel('Assumed peak conversion rate, $\\Gamma/H_0$')
    axes[0].set_ylabel('Share of total matter energy')
    fig.suptitle('When does radiation carry half the energy?',x=.075,ha='left',fontsize=16)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False)
    fig.subplots_adjust(left=.08,right=.98,top=.80,bottom=.24,wspace=.20)
    fig.savefig(ROOT/'radiation-share.png',dpi=180,bbox_inches='tight',facecolor='white')
    plt.close(fig)

if __name__=='__main__':
    draw()
