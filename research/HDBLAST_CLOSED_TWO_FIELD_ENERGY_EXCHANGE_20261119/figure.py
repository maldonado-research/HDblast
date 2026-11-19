"""Plot the integrated work and net energy change to each first contact."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def make(output):
    values=json.loads((HERE/'coupled/COUPLED.json').read_bytes())['summary']
    cases=['0.001/velocity','0.001/displacement','-0.001/velocity','-0.001/displacement']
    rows=[values[c+'/0.01'] for c in cases]
    direct=np.array([r['transfer_to_crossing_as_percent_of_daughter_seed'] for r in rows])
    total=np.array([r['energy_change_to_crossing_as_percent_of_daughter_seed'] for r in rows])
    fig,ax=plt.subplots(figsize=(8.0,4.0),layout='constrained')
    positions=np.arange(4)[::-1]
    for numbers,offset,color,label in [(direct,.15,'#006178','From the rolling field'),
                                      (total,-.15,'#a14526','Net change including expansion')]:
        ax.barh(positions+offset,numbers,height=.25,color=color,label=label,zorder=2)
        for y,val in zip(positions+offset,numbers):
            ax.text(val+(.3 if val>=0 else -.3),y,f'{val:+.2f}',ha='left' if val>=0 else 'right',
                    va='center',fontsize=9)
    ax.axvline(0,lw=.9,color='#354448',zorder=1)
    ax.set_yticks(positions,['Poles · velocity','Poles · displacement','Equator · velocity','Equator · displacement'])
    ax.set_xlabel('Energy change by first crossing (% of the second field’s initial energy)')
    ax.set_xlim(-17.5,23.5)
    ax.grid(axis='x',color='#d7e0e0',lw=.7)
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y',length=0)
    ax.legend(loc='lower center',bbox_to_anchor=(.5,1.00),ncols=2,frameon=False)
    fig.savefig(output,dpi=160)
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    make(parser.parse_args().output)
