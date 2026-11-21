"""Draw the change in early transverse pressure and momentum counts with elastic collisions."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def draw(study, output):
    trials = json.loads((study / 'live/LIVE.json').read_bytes())['runs']
    cases = ('small_pulse', 'pressure_fit')
    selected = [r for r in trials if r['birth_pairs'] == 16384 and r['steps'] == 512]
    baseline = [next(r for r in selected if r['source'] == s and r['law'] == 'off') for s in cases]
    collisions = [[r for r in selected if r['source'] == s and r['law'] == 'vector'] for s in cases]
    if any(len(group) != 3 for group in collisions):
        raise ValueError('Expected three independent collision seeds per source')
    fig, (pressure, counts) = plt.subplots(1, 2, figsize=(10.3, 4.4))
    labels = ['Small extra waves', 'Pressure-fitted waves']
    centers = np.arange(2)
    grey = '#71818d'
    teal = '#136c70'
    for ax in (pressure, counts):
        ax.set_xticks(centers, labels)
        ax.grid(axis='y', color='#dce0e2', lw=.7)
        ax.set_axisbelow(True)
        ax.spines[['right','top']].set_visible(False)
    for ax, measure, scale in ((pressure, lambda r: r['sideways_pressure_over_target'][0], 1.),
                               (counts, lambda r: r['tv20'], 100.)):
        old = np.array([measure(row) for row in baseline]) * scale
        samples = np.array([[measure(row) for row in group] for group in collisions]) * scale
        mean = samples.mean(axis=1)
        lo = mean - samples.min(axis=1)
        hi = samples.max(axis=1) - mean
        ax.bar(centers-.17, old, .31, color=grey, label='Collisions off', zorder=2)
        ax.bar(centers+.17, mean, .31, color=teal, label='Vector-mediated collisions', zorder=2)
        ax.errorbar(centers+.17, mean, yerr=[lo,hi], fmt='none', capsize=3,
                    ecolor='#1d3538', lw=1.2, zorder=3)
        for i, (base, ave, high) in enumerate(zip(old, mean, hi)):
            if ax is pressure:
                a, b = ('<0.001' if base < .001 else f'{base:.2f}'), f'{ave:.2f}'
                offset = .048
            else:
                a, b = f'{base:.1f}%', f'{ave:.1f}%'
                offset = 1.1
            ax.text(i-.17, base+offset, a, ha='center', fontsize=9, color='#263a47')
            ax.text(i+.17, ave+high+offset, b, ha='center', fontsize=9, color='#12494c')
    pressure.axhline(1., ls='--', color='#28383c', lw=1)
    pressure.text(.02, 1.018, 'Thermal pressure', transform=pressure.get_yaxis_transform(),
                  ha='left', va='bottom', fontsize=9)
    pressure.set(title='Pressure in the weaker sideways direction', ylabel='Pressure / thermal pressure',
                 ylim=(0,1.54), yticks=[0,.25,.5,.75,1,1.25,1.5])
    counts.set(title='Distance from the thermal particle counts', ylabel='Particles needing to change bin (%)',
               ylim=(0,72), yticks=range(0,71,10))
    handles, names = counts.get_legend_handles_labels()
    fig.legend(handles, names, loc='lower center', bbox_to_anchor=(.50,.06), ncol=2, frameon=False)
    fig.text(.5, .028, 'At 1% of reference electrical work. Lines above bars span three collision seeds; 32,768 particles per run.',
             ha='center', va='bottom', fontsize=9, color='#35434a')
    fig.subplots_adjust(top=.89, bottom=.23, left=.075, right=.99, wspace=.28)
    fig.savefig(output, dpi=175, bbox_inches='tight', facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--output', type=Path, default=Path('scattering-results.png'))
    args = parser.parse_args()
    draw(args.study, args.output)
