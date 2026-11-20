"""Plot the finite-distance, dissipation-free, coupled Einstein–Maxwell example."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
NAVY='#243744'
BLUE='#226977'
GOLD='#ad562b'
MUTED='#506674'


def main():
    data=json.loads((ROOT/'loading/LOADING.json').read_bytes())['forward']
    row=data['runs']['maxwell_primary']
    curve=dict(zip(data['trace_columns'],np.array(row['time_series']).T,strict=True))
    t,a=curve['t'],curve['a']
    after=t>=row['bounce_tau']
    fraction=np.divide(curve['sheet'],curve['e'],out=np.zeros_like(t),where=curve['e']>0)
    start=np.interp(row['first_em90_radius'],a[(t>row['bounce_tau'])&(t<row['maximum_tau'])],
                    t[(t>row['bounce_tau'])&(t<row['maximum_tau'])])
    plt.rcParams.update({'font.size':11,'text.color':NAVY,'axes.labelcolor':NAVY,'axes.titlecolor':NAVY,
                         'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#b3c3ca',
                         'xtick.color':MUTED,'ytick.color':MUTED})
    fig,(left,right)=plt.subplots(1,2,figsize=(11.3,5.4))
    fig.subplots_adjust(left=.084,right=.97,bottom=.23,top=.74,wspace=.35)
    fig.text(.084,.955,'From an incoming wave to an expanding wall',weight='bold',fontsize=20,va='top')
    fig.text(.084,.86,'Energy returns to the same electromagnetic field in a closed five-dimensional model',color=MUTED)
    for ax in (left,right):
        ax.axvspan(start,row['maximum_tau'],color='#d9e9e4',zorder=0)
        ax.set_xlabel('Time on the wall (model units)')
        ax.grid(axis='y',alpha=.16)
    left.plot(t,a,color=BLUE,lw=2.6)
    left.set_title('Wall radius',loc='left',weight='bold')
    left.set_ylabel('Radius (bulk curvature lengths)')
    left.scatter([row['bounce_tau'],row['maximum_tau']],[row['bounce_radius'],row['maximum_radius']],
                 color=BLUE,zorder=5,s=30)
    left.annotate('Bounce',(row['bounce_tau'],row['bounce_radius']),xytext=(6,10),textcoords='offset points',
                  fontsize=9,color=MUTED)
    left.annotate('Expansion stops',(row['maximum_tau'],row['maximum_radius']),xytext=(-7,-18),ha='right',
                  textcoords='offset points',fontsize=9,color=MUTED)
    left.set_ylim(1,10.65);left.set_xlim(-.15,max(t)+.1)
    right.plot(t[after],100*fraction[after],color=GOLD,lw=2.6)
    right.axhline(90,ls='--',lw=1.0,color='#63727a')
    right.set_title('Electromagnetic share after the bounce',loc='left',weight='bold')
    right.set_ylabel('Percent of added wall energy')
    right.set_xlim(.7,max(t)+.1);right.set_ylim(-3,102)
    fig.text(.084,.072,'Shading marks expansion with more than 90% of the added wall energy in coherent electromagnetic fields.',color=MUTED,fontsize=9.4)
    fig.text(.084,.032,'The radiation-fluid density and damping rate are zero throughout. The model stops shortly after contraction resumes.',color=MUTED,fontsize=9.4)
    fig.savefig(ROOT/'curved-wall.png',dpi=170)
    plt.close(fig)


if __name__=='__main__':
    main()
