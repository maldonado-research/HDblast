"""Plot the executed follow-up arrays, preserving the original classification."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path)
    args=parser.parse_args()
    summary=json.loads((args.run/'summary.json').read_text())
    if summary['status'] not in ('PASS','FAIL'):
        raise ValueError('Follow-up must be complete')
    output=args.run/'figures'
    output.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True,constrained_layout=True)
    with np.load(args.run/summary['selected']/'arrays.npz',allow_pickle=False) as a:
        t=a['t']
        for ax,key in ((axes[0,0],'p'),(axes[1,0],'rho'),(axes[1,1],'Q')):
            for K,color in ((96,'#d1872c'),(192,'#709eb5'),(384,'#16697a')):
                ax.plot(t,a[f'K{K}_{key}'],color=color,lw=1.5 if K==384 else .8,label=f'K={K}')
            ax.set_ylabel({'p':r'$p$','rho':r'$\rho$','Q':r'$Q$'}[key])
            ax.axvline(.5,color='gray',lw=.7,ls=':')
            ax.axvline(1,color='gray',lw=.7,ls='--')
            ax.legend(fontsize=8)
        ax=axes[0,1]
        ax.plot(t,np.abs(a['K384_p']-a['K192_p']),color='#16697a',label='quadrature-refined cutoff change')
        allowance=summary['cutoff_K192_to_K384']['p']['limit']
        ax.axhline(allowance,color='#bc4749',ls='--',label='registered maximum allowance')
        ax.set_ylabel(r'$|p_{384}-p_{192}|$')
        ax.set_title(f"pressure cutoff gate: {'PASS' if summary['cutoff_K192_to_K384']['p']['passed'] else 'FAIL'}")
        ax.legend(fontsize=7)
        for ax in axes[1]:
            ax.set_xlabel(r'conformal time $\eta$')
    fig.suptitle(f"A=0.2 prescribed smooth control: K384 follow-up {summary['status']}\nOriginal K192 matrix remains FAIL")
    fig.savefig(output/'followup_cutoff.png',dpi=200)
    fig.savefig(output/'followup_cutoff.pdf')
    plt.close(fig)
    print(output)


if __name__=='__main__':
    main()
