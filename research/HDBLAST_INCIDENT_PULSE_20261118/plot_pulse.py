"""Draw the measured shell response and signed energy transfer using the finest grid."""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load(path):
    with path.open(newline='') as f:
        rows=list(csv.DictReader(f))
    return {k:np.array([float(v[k]) for v in rows]) for k in rows[0]}


def transfer(s,h):
    weights=np.array([77,-214,234,-122,25])/12
    grad=sum(w*(s[f'phi_minus_{j}']-s['phi_minus_1']) for j,w in enumerate(weights,1))/h
    rate=-2*np.exp(3*(s['A_wall']-s['A_wall'][0])-s['B_wall'])*s['phi_T_wall']*grad
    return np.r_[0,np.cumsum(np.diff(s['T'])*(rate[:-1]+rate[1:])/2)]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('runs',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    grid=args.runs/'g8'
    quiet=load(grid/'quiet.csv')
    t=quiet['T']
    q=transfer(quiet,0.000125)
    fig,axes=plt.subplots(1,2,figsize=(10.7,4.35),constrained_layout=True)
    for name,color,label in [('delayed','#155995','+0.001'),('delayed_negative','#bd4925','−0.001')]:
        s=load(grid/(name+'.csv'))
        axes[0].plot(t,1000*(s['phi_wall']-quiet['phi_wall']),lw=2,color=color,label=label)
        axes[1].plot(t,1000*(transfer(s,.000125)-q),lw=2,color=color,label=label)
    for ax,title,units in zip(axes,('Scalar field at the shell','Net energy transferred into the shell'),
                              ('Change from quiet (× 10⁻³)','Cumulative change from quiet (× 10⁻³ model units)')):
        ax.axhline(0,color='#555555',lw=.7,zorder=0)
        ax.axvline(.10,color='#777777',ls=':',lw=1)
        ax.axvline(.15,color='#555555',ls='--',lw=1)
        ax.set(xlim=(.087,.158),xlabel='Conformal time T',ylabel=units,title=title)
        ax.grid(axis='y',alpha=.18)
        ax.spines[['top','right']].set_visible(False)
        ax.set_xticks([.09,.10,.11,.12,.13,.14,.15])
    axes[0].legend(title='Initial scalar amplitude',fontsize=9,loc='lower left',frameon=False)
    fig.suptitle('A matched incoming pulse reaches the model shell',fontsize=15,fontweight='bold')
    fig.savefig(args.output,dpi=180,facecolor='white')
    plt.close(fig)


if __name__=='__main__':
    main()
