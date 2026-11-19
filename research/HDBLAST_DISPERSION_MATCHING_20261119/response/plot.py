"""Scientific plots of the predeclared bath response experiment."""
import argparse
import os
from pathlib import Path
import tempfile

os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'hdblast-dispersion-mpl'))
os.environ.setdefault('XDG_CACHE_HOME',str(Path(tempfile.gettempdir())/'hdblast-dispersion-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def plot(data,output):
    with np.load(data) as z:
        q,cutoffs,ratios,orders = (z[key] for key in ('q','cutoffs','ratios','orders'))
        error,quad,weighted = (z[key] for key in ('cutoff_relative_error','quadrature_relative_error','weighted_quadrature_relative_l1'))
    fig,axes = plt.subplots(1,3,figsize=(13.8,4.1),layout='constrained')
    colors = plt.get_cmap('viridis')(np.linspace(.12,.9,len(cutoffs)))
    for ci,cutoff in enumerate(cutoffs):
        axes[0].plot(ratios,100*np.max(error[ci],axis=0),'o-',color=colors[ci],label=f'Λ={cutoff}')
        axes[1].loglog(q,np.maximum(np.max(abs(quad[0,ci]),axis=-1),1e-15),'o-',ms=3,color=colors[ci],label=f'Λ={cutoff}')
        axes[2].loglog(orders,np.max(abs(quad[:,ci]),axis=(1,2)),'o-',ms=4,color=colors[ci],label=f'Λ={cutoff}')
    axes[0].set(xlabel='Held-out ν/q',ylabel='Maximum missing continuum response (%)',title='Finite cutoff after matching value and slope')
    axes[0].legend(frameon=False)
    axes[1].set(xlabel='Native comoving q / H₀',ylabel='N=64 maximum relative quadrature error',title='Threshold resolution in original finite bath')
    axes[2].set(xlabel='Bath Gauss order N',ylabel='Maximum relative quadrature error',title='All 24 modes, all four frequency probes')
    for ax in axes:ax.grid(alpha=.2)
    fig.savefig(output,dpi=190)
    plt.close(fig)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=Path(__file__).resolve().parent/'TABLES.npz')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'response.png')
    args = parser.parse_args()
    plot(args.data,args.output)
