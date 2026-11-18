"""Plot the two weak tests from the adjacent shape.py replay and SHAPE.json."""
from fractions import Fraction as F
import json
from pathlib import Path
import runpy

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
H = F(1, 32)
STARTS = (F(5, 4), F(11, 8))
COLORS = ('#126677', '#9c4924')


def spline(t):
    if t < 0 or t > 3:
        return 0
    # Three unit boxes convolved together, divided by the peak (3/4).
    return 2 * sum(c * max(t - j, 0)**2 for j, c in enumerate((1, -3, 3, -1))) / 3


def render():
    api = runpy.run_path(str(ROOT / 'shape.py'))
    rows = api['load_cells']()
    result = json.loads((ROOT / 'SHAPE.json').read_text())['waveform']
    bands = [[F(x) for x in pair] for pair in result['mass_squared_intervals']]
    for start, band in zip(STARTS, bands):
        if tuple(band) != api['mass_test'](rows, start, H, F(0)):
            raise ValueError('Weak test bounds disagree with SHAPE.json')

    fig = plt.figure(figsize=(10.4, 5.4), layout='constrained')
    layout = fig.add_gridspec(2, 2, width_ratios=(1.55, 1), height_ratios=(3.4, 1))
    trace = fig.add_subplot(layout[0, 0])
    weight = fig.add_subplot(layout[1, 0], sharex=trace)
    mass = fig.add_subplot(layout[:, 1])
    fig.suptitle('Two clock intervals require different masses for one free scalar',
                 x=0.02, ha='left', fontsize=15, weight='bold')
    trace.set_title('Normalized shell waveform', loc='left', fontsize=12)
    mass.set_title('Required scaled mass squared', loc='left', fontsize=12)

    lower, upper = F(19, 16), F(3, 2)
    xx, lo, hi = [], [], []
    for left, right, low, high in api['overlaps'](rows, lower, upper):
        xx.extend((float(left), float(right)))
        lo.extend([float(low)] * 2)
        hi.extend([float(high)] * 2)
    trace.fill_between(xx, lo, hi, color='#263e4a', linewidth=1)
    trace.set_ylim(-0.12, 4.75)
    trace.set_ylabel(r'Normalized amplitude $U$')
    trace.tick_params(labelbottom=False)
    for i, (start, color) in enumerate(zip(STARTS, COLORS)):
        end = start + 3 * H
        trace.axvspan(float(start), float(end), color=color, alpha=0.13, lw=0)
        trace.text(float(start + 3*H/2), 0.19, f'Probe {i+1}', ha='center', color=color)
        u = [F(k, 100) for k in range(301)]
        weight.fill_between([float(start + H*t) for t in u], [float(spline(t)) for t in u],
                            color=color, alpha=0.85, lw=0)
    weight.set_xlim(float(lower), float(upper))
    weight.set_ylim(0, 1.08)
    weight.set_yticks([0, 1])
    weight.set_xticks([1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5])
    weight.set_ylabel('Probe weight\n(unit peak)')
    weight.set_xlabel(r'Past observer clock $x=q/\ell$')

    for i, ((left, right), color) in enumerate(zip(bands, COLORS)):
        height = 1-i
        mass.hlines(height, float(left), float(right), color=color, lw=8)
        mass.vlines([float(left), float(right)], height-.07, height+.07, color=color, lw=1.7)
        # Text labels rounded outward to two decimal places.
        low = (left * 100).__floor__() / 100
        high = (right * 100).__ceil__() / 100
        mass.text(float((left+right)/2), height+.14, f'{low:.2f} to {high:.2f}',
                  color=color, ha='center', weight='bold')
    mass.set_xlim(40, 90)
    mass.set_ylim(-0.45, 1.6)
    mass.set_yticks([1, 0], ['Probe 1', 'Probe 2'])
    mass.set_xlabel('Dimensionless mass squared in the scaled clock')
    mass.grid(axis='x', color='#dce2e4', zorder=0)
    for ax in (trace, weight, mass):
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['bottom', 'left']].set_color('#a4aeb2')
        ax.tick_params(color='#808b90')
    fig.savefig(ROOT/'mass_test.png', dpi=180, facecolor='white', metadata={'Software':'Matplotlib'})
    plt.close(fig)


if __name__ == '__main__':
    render()
