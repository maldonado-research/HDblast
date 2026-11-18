"""Render the past first-variation coefficient from the capsule's interval data."""

import csv
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


HERE = Path(__file__).resolve().parent
HORIZON = Fraction(25, 16)
WINDOW = (Fraction(5, 4), Fraction(11, 8))


def plot_cells(inputs):
    future = json.loads((inputs / 'FUTURE_ENCLOSURE.json').read_text())
    transfer = json.loads((inputs / 'PAST_TRANSFER.json').read_text())
    cells = inputs / 'FUTURE_CELLS.csv'
    grid = future['cell_data']
    if hashlib.sha256(cells.read_bytes()).hexdigest() != grid['sha256']:
        raise ValueError('Interval file differs from the enclosure report')
    ell = Fraction(future['state']['ell'])
    error = Fraction(transfer['difference']['safe_strict_global_upper'])
    # The CSV bounds ell*S in integer units. The transfer error is already in S.
    scale = ell * grid['value_denominator']
    with cells.open(newline='') as handle:
        for row in csv.DictReader(handle):
            left = Fraction(int(row['left_grid']), grid['time_grid_denominator'])
            if left >= HORIZON:
                break
            right = Fraction(int(row['right_grid']), grid['time_grid_denominator'])
            low = int(row['ell_S_lo']) / scale - error
            high = int(row['ell_S_hi']) / scale + error
            yield left, min(right, HORIZON), low, high


def band(ax, cells, interval, color):
    # Repeated endpoints draw each enclosure as a rectangle, including jumps.
    edges, lows, highs = [], [], []
    for left, right, low, high in cells:
        left, right = max(left, interval[0]), min(right, interval[1])
        if left < right:
            edges.extend((float(left), float(right)))
            lows.extend([float(low / 1000)] * 2)
            highs.extend([float(high / 1000)] * 2)
    ax.fill_between(edges, lows, highs, facecolor=color, edgecolor=color, linewidth=0.65, zorder=3)


def render(inputs, destination):
    cells = list(plot_cells(inputs))
    if not cells or cells[0][0] != 0 or cells[-1][1] != HORIZON:
        raise ValueError('Incomplete plotting interval')
    teal = '#17616a'
    gold = '#f8e8b6'
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.9), gridspec_kw={'width_ratios': [1.65, 1]})
    fig.subplots_adjust(left=0.10, right=0.97, top=0.78, bottom=0.23, wspace=0.30)
    fig.suptitle('Past shell scalar response at first order', x=0.10, y=0.96,
                 ha='left', fontsize=16, weight='bold')
    full, detail = axes
    full.set_title('Full observer interval', loc='left', fontsize=12)
    detail.set_title('Signed interval, expanded', loc='left', fontsize=12)
    full.axvspan(*map(float, WINDOW), color=gold, zorder=0)
    full.axhline(0, color='#66747b', ls=':', lw=0.9, zorder=2)
    band(full, cells, (0, HORIZON), teal)
    band(detail, cells, WINDOW, teal)
    full.set_xlim(0, float(HORIZON))
    full.set_ylim(-51, 26)
    full.set_yticks([-40, -20, 0, 20])
    full.set_xticks([0, 0.5, 1, 1.25, float(HORIZON)], ['0', '1/2', '1', '5/4', '25/16'])
    detail.set_xlim(*map(float, WINDOW))
    detail.set_ylim(-48, -33)
    detail.set_yticks([-46, -42, -38, -34])
    detail.set_xticks([1.25, 1.3125, 1.375], ['5/4', '21/16', '11/8'])
    full.set_ylabel(r'$S(-q)$  (thousands of model units)', labelpad=7)
    for ax in axes:
        ax.set_xlabel(r'Past observer time  $q/\ell$')
        ax.tick_params(length=3, color='#66747b')
        ax.grid(axis='y', color='#dce2e5', zorder=1)
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['bottom', 'left']].set_color('#83929a')
    fig.legend(handles=[Patch(facecolor=teal, label='Cellwise coefficient enclosure'),
                        Patch(facecolor=gold, label='Window [5/4, 11/8]')],
               loc='lower left', bbox_to_anchor=(0.085, 0.012), frameon=False, ncol=2)
    fig.savefig(destination, dpi=190, facecolor='white', metadata={'Software': 'Matplotlib'})
    plt.close(fig)


if __name__ == '__main__':
    render(HERE / 'inputs', HERE / 'response.png')
