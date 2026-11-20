"""Render the threshold, heating and closed-wall continuation from the frozen products."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def plot(decay, fate, target):
    model = decay['runs']['primary']
    d = dict(zip(decay['trace_columns'], np.asarray(model['trace']).T, strict=True))
    anchor = fate['scalar_threshold']['anchor']
    late = (d['t'] >= 9.5) & (d['t'] <= model['maximum_tau']+.20)
    t = d['t'][late]

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.3), layout='constrained')
    ax[0].plot(t, 100*d['f'][late], color='#19756a', lw=2.3, label='Neutral parents decayed')
    ax[0].plot(t, 100*d['intrinsic_random'][late]/d['e'][late], color='#823d70',
               lw=2.6, label='Particle internal share of wall energy')
    ax[0].axvline(model['threshold_tau'], c='#686873', ls=':', lw=1.4)
    ax[0].axvline(model['maximum_tau'], c='#343541', ls='--', lw=1.3)
    ax[0].text(model['threshold_tau']-.035, 47, 'Decay opens', rotation=90, ha='right', va='center', size=9)
    ax[0].text(model['maximum_tau']-.035, 47, 'Expansion stops', rotation=90, ha='right', va='center', size=9)
    ax[0].set(title='Charged particles appear as the wall expands', xlabel='Wall proper time (model units)',
              ylabel='Fraction (%)', ylim=(-3,104))
    ax[0].legend(loc='lower right', frameon=True, fontsize=8)
    ratio = np.linspace(1,1.36,500)
    for share,col,name in [(0,'#474c5b','No vacuum'),(.1,'#2479a4','10% vacuum'),
                           (fate['scalar_threshold']['radiation']['critical_vacuum_share'],'#b56521','41.34% vacuum (critical)')]:
        a = anchor['a']*ratio
        e = anchor['e']*(share+(1-share)/ratio**4)
        rhs = e*(2+e)+anchor['mass']/a**4-1/a**2
        ax[1].plot(ratio, 1000*rhs, color=col, lw=2.2, label=name)
    ax[1].axhline(0,c='#242535',ls='--',lw=1)
    ax[1].set(title='Can radiation keep the closed wall expanding?',
              xlabel='Radius / radius at the 90% heating sample',
              ylabel='Friedmann right-hand side × 1,000')
    ax[1].set_ylim(-2.2,.7)
    ax[1].legend(loc='lower left', fontsize=8)
    for panel in ax:
        panel.grid(alpha=.19)
        panel.spines[['right','top']].set_visible(False)
    fig.savefig(target, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=HERE/'neutral-parent-decay.png')
    args = p.parse_args()
    plot(json.loads((HERE/'coupled/DECAY.json').read_bytes()),
         json.loads((HERE/'expansion/FATE.json').read_bytes()), args.output)
