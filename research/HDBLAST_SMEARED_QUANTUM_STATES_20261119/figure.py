"""Plot three photon launch fractions on each of the two calibrated potentials."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

HERE = Path(__file__).resolve().parent


def draw():
    data = json.loads((HERE/'radiation/PHOTONS.json').read_bytes())['branches']
    colors = {'0.25':'#236482','0.75':'#a45120','0.95':'#61498b'}
    fig, axes = plt.subplots(1,2,figsize=(10,4.8),sharex=True,sharey=True)
    for ax,(key,b) in zip(axes,data.items()):
        for load,color in colors.items():
            result=b['runs'][load]
            initial=result['initial']
            curve=result['curve']
            ax.plot([p['x'] for p in curve],[p['a_over_ai'] for p in curve],color=color,lw=2.3)
            turn=result['first_maximum']
            ax.plot(turn['x'],turn['a']/initial['a'],'s',ms=5.8,color=color)
            passage=result['first_crossing']
            if passage is not None:
                ax.plot(passage['x'],passage['a']/initial['a'],'o',ms=8,mfc='white',mec=color,mew=1.7)
        ax.set_title('Potential minimum at '+key.split('=')[1],fontsize=12)
        ax.set_xlabel('Time since bounce, $H_0 t$')
        ax.xaxis.set_major_locator(MultipleLocator(.25))
        ax.yaxis.set_major_locator(MultipleLocator(.05))
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y',alpha=.23)
    axes[0].set_ylabel('Radius relative to its bounce value')
    handles=[Line2D([0],[0],color=c,lw=2.3,label=f'{100*data["phi_m=0.5"]["runs"][load]["initial"]["matter_fraction"]:.2f}% photons')
             for load,c in colors.items()]
    handles.extend([Line2D([0],[0],marker='o',ls='',mfc='white',mec='#393939',mew=1.5,label='Field reaches minimum'),
                    Line2D([0],[0],marker='s',ls='',color='#393939',label='Expansion stops')])
    fig.suptitle('More light at launch can leave less time for the field to roll',x=.065,ha='left',fontsize=16)
    fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,fontsize=10)
    fig.subplots_adjust(left=.085,right=.98,top=.82,bottom=.28,wspace=.15)
    fig.savefig(HERE/'photon-turnaround.png',dpi=180,bbox_inches='tight',facecolor='white')
    plt.close(fig)


if __name__=='__main__':
    draw()
