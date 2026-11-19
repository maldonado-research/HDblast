"""Plot the curvature threshold and the number of quanta below a given momentum."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import zeta

ROOT = Path(__file__).resolve().parent
BLUE, RED, MUTED = "#155a82", "#b34d29", "#69737c"


def main(output):
    transition = json.loads((ROOT / 'curvature/CURVATURE.json').read_text())['universal']
    scan = transition['scanned']
    strength = [v['g'] for v in scan]
    frac = np.array([v['ratio_at_turn'] for v in scan])
    frac = frac / (1 + frac)
    case = json.loads((ROOT / 'spectra/SPECTRA.json').read_text())['branches']['0.5']['0.001']['epochs']['equality']
    with (ROOT / 'spectra/SPECTRA.csv').open() as stream:
        profile = [r for r in csv.DictReader(stream) if r['phi_min']=='0.5' and r['detuning']=='0.001' and r['epoch']=='equality']
    ratio = case['bose_count_fraction']
    temp_over_half_mass = case['mean_p_over_half_mass'] * ratio / (np.pi**4 / (30 * zeta(3)))
    p = np.geomspace(0.00001, 1.5, 3500)
    x = p / temp_over_half_mass
    eq_count = cumulative_trapezoid(x**2 / np.expm1(x), x, initial=0) / (2 * zeta(3))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), layout='constrained')
    left, right = axes
    left.plot(strength, frac, color=BLUE, marker='.', markersize=5)
    left.axhline(.5, color=MUTED, linestyle='--', linewidth=1)
    left.axvline(transition['g_critical'],color=MUTED,linestyle=':',linewidth=1)
    left.scatter([transition['g_critical']], [.5], color=RED, zorder=4)
    left.annotate('Equality at first maximum\n$g=0.8351$', xy=(.8351,.5),xytext=(1.1,.31),fontsize=10)
    left.set(xscale='log', xlim=(.09,17), ylim=(0,1), xlabel=r'Decay strength $g=\gamma M/K^{3/2}$',
             ylabel='Radiation share at first maximum', title='When does radiation overtake the parent?')
    right.plot(p, eq_count, color=MUTED, label='Equilibrium gas at the same energy')
    px = np.array([float(r['p_over_half_mass']) for r in profile])
    counts = [float(r['cumulative_number_fraction'])*ratio for r in profile]
    right.plot(np.r_[p[0],px,1.5], [0,*counts,ratio], color=BLUE, linewidth=2.4, label='Quanta from parent decays')
    right.annotate('8.74% of equilibrium count',xy=(.7,.082),xytext=(.009,.19),fontsize=10,color=BLUE)
    right.set(xscale='log',xlim=(.001,1.5),ylim=(0,1.07),xlabel='Daughter momentum / (half the parent mass)',
              ylabel='Count below momentum / full equilibrium count', title='How many particles carry the energy?')
    right.legend(loc='center right',fontsize=8.8,frameon=False)
    for ax in axes:
        ax.grid(axis='y', color='#d6dce0', linewidth=.6)
        ax.spines[['right','top']].set_visible(False)
    fig.savefig(output, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output', type=Path, required=True)
    main(cli.parse_args().output)
