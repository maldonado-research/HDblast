"""Plot the late axial field and the directional excess in three collision models."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

HERE = Path(__file__).resolve().parent
LABEL = 't6.67000-p14-n4096-r91283-particle-pole'


def draw(snapshots, triad, destination):
    with np.load(snapshots / (LABEL + '.npz')) as record:
        time = record['t']
        field = record['E'][:, 2]
        fluid = record['reference_E'][:, 2]
    data = json.loads(triad.read_bytes())
    samples = [next(r for r in data['runs'] if r['snapshot'] == LABEL)]
    for mode in ('elastic', 'convert'):
        samples.append(next(r for r in data['shared_runs']
                            if r['snapshot'] == LABEL and r['flavor_mode'] == mode))
    gaps = [100*r['number_weighted']['any_boost_gap'] for r in samples]

    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10,
                         'axes.spines.top':False, 'axes.spines.right':False})
    blue, ochre = '#174E72', '#985210'
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.4, 3.7), layout='constrained')
    ax.set_title('The axial field reverses at different times', loc='left', fontsize=12, weight='bold')
    ax.axhline(0, color='#777777', linewidth=.7)
    ax.plot(time, fluid, '--', color=ochre, lw=1.8, label='Fluid reference')
    ax.plot(time, field, color=blue, lw=2, label='Simulated particle current')
    for series, color, offset in ((field, blue, 7), (fluid, ochre, -11)):
        ax.plot(time[-1], series[-1], 'o', ms=4, color=color)
        ax.annotate(f'{series[-1]:+.3f}', (time[-1], series[-1]),
                    xytext=(-2, offset), textcoords='offset points', ha='right', color=color)
    ax.set(xlabel='Shell proper time (model units)', ylabel='Axial electric field (model units)',
           xlim=(6.445, 6.674), ylim=(-.052, .130))
    ax.xaxis.set_major_locator(MultipleLocator(.05))
    ax.legend(frameon=False, loc='upper left')
    ax.grid(axis='y', alpha=.15)
    bx.set_title('Particles still favor particular directions', loc='left', fontsize=12, weight='bold')
    y = np.arange(3)
    bars = bx.barh(y, gaps, height=.5, color=['#174E72','#448183','#448183'])
    bx.set_yticks(y, ['Separate\ninteraction', 'Shared, stay\nin flavor', 'Shared, allow\nflavor changes'])
    bx.invert_yaxis()
    bx.set(xlim=(0, 64), xlabel='Minimum angular mismatch (percentage points)')
    bx.xaxis.set_major_locator(MultipleLocator(15))
    bx.bar_label(bars, labels=[f'{value:.1f}' for value in gaps], padding=4)
    bx.grid(axis='x', alpha=.15)
    bx.set_axisbelow(True)
    fig.savefig(destination, dpi=220, bbox_inches='tight', facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshots', required=True, type=Path,
                        help='Regenerated live/snapshots directory')
    parser.add_argument('--triad', type=Path, default=HERE/'triad/TRIAD.json')
    parser.add_argument('--output', type=Path, default=HERE/'axial-response.png')
    args = parser.parse_args()
    draw(args.snapshots, args.triad, args.output)
