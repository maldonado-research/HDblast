#!/usr/bin/env python3
"""Plot certified cell enclosures; floating conversion is display only."""
import csv
from datetime import datetime, timezone
from fractions import Fraction
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
report = json.loads((ROOT/'computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json').read_text())
rows = list(csv.DictReader((ROOT/'computed/COORDINATE_SHELL_SOURCE_CELLS.csv').open()))
D, G = report['cell_data']['value_denominator'], report['cell_data']['time_grid_denominator']
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.titleweight': 'bold', 'figure.facecolor': '#f7f8fc',
                     'axes.facecolor': '#ffffff', 'savefig.facecolor': '#f7f8fc'})
fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.2), layout='constrained')
colors = ['#166b87', '#a64736', '#7852a0']
xs = np.array([int(row[key])/G for row in rows for key in ('left_grid', 'right_grid')])
for ax, field, factor, title, ylabel in (
        (axes[0, 0], 'u', 1e8, 'Actual first metric boundary trace', r'$10^8 u_b$'),
        (axes[0, 1], 'utau', 1e4, 'Proper-time derivative enclosure', r'$10^4\,\partial_\tau u_b$')):
    lower = np.repeat([int(row[field+'_lo'])/D*factor for row in rows], 2)
    upper = np.repeat([int(row[field+'_hi'])/D*factor for row in rows], 2)
    ax.fill_between(xs, lower, upper, color=colors[0], alpha=.22)
    ax.plot(xs, (lower+upper)/2, color=colors[0], linewidth=.8)
    ax.set(title=title, xlabel=r'$x=t/\ell$', ylabel=ylabel, xlim=(0, 3))
    ax.grid(alpha=.18)
    ax.axvspan(2, 3, color='#dae1ed', alpha=.25)
ax = axes[1, 0]
for i, name in enumerate(('Sp', 'SK', 'EN')):
    lo, hi = [float(Fraction(value)) for value in report['exact']['norms'][name]['squared_L2']]
    ax.plot([lo, hi], [i, i], color=colors[i], linewidth=5, solid_capstyle='round')
    ax.plot([lo, hi], [i, i], '|', color=colors[i], markersize=12)
ax.set_xscale('log')
ax.set_yticks([0, 1, 2], [r'$S_pY^2$', r'$S_KY^2$', r'$E_{N,4}H_4$'])
ax.set(title='Squared source norms, registered model units', xlabel=r'$\|\mathrm{source}\|_{L^2(d\mu_0)}^2$', ylim=(-.5, 2.5))
ax.grid(axis='x', alpha=.18)
ax.text(.03, .03, 'Boundary source amplitudes; the response is not evaluated.',
        transform=ax.transAxes, fontsize=9, color='#44516b')
ax = axes[1, 1]
ax.axis('off')
ax.set_title('Coordinate source signs in the central bracket', loc='left')
data = [['At left endpoint A', 'negative', 'positive'],
        [r'At a scalar zero $x_*\in(A,B)$', 'positive', 'negative'],
        ['At right endpoint B', 'negative', 'positive'],
        [r'Throughout tail $2\leq x\leq3$', 'negative', 'positive']]
table = ax.table(cellText=data, colLabels=['Certified location', r'$S_p$', r'$S_K$'],
                 cellLoc='left', colLoc='left', bbox=[0, .36, 1, .55], colWidths=[.59, .205, .205])
table.auto_set_font_size(False)
table.set_fontsize(9)
for (r, c), cell in table.get_celld().items():
    cell.set_edgecolor('#dce2ee')
    cell.set_facecolor('#e8eef6' if r == 0 else '#ffffff')
    if r == 0:
        cell.set_text_props(weight='bold')
ax.text(0, .24, r'$A=1.49998,\quad B=1.49999$', transform=ax.transAxes)
ax.text(0, .13, 'Continuity gives at least two distinct zeros of each coefficient in (A,B).',
        transform=ax.transAxes, fontsize=9)
ax.text(0, .02, 'The signs depend on the retained coordinate representative.\nThey do not determine a force, energy or second-order response.',
        transform=ax.transAxes, fontsize=9, color='#44516b')
fig.suptitle('HDBLAST · continuous coordinate shell source checkpoint', fontsize=16, weight='bold')
target = ROOT/'computed/COORDINATE_SHELL_SOURCES'
fig.savefig(target.with_suffix('.png'), dpi=180, metadata={'Software': 'HDBLAST coordinate shell sources'})
fixed = datetime(2026, 10, 15, tzinfo=timezone.utc)
fig.savefig(target.with_suffix('.pdf'), metadata={'Title': 'HDBLAST continuous coordinate shell sources',
            'Creator': 'HDBLAST reproducible figure', 'CreationDate': fixed, 'ModDate': fixed})
print('Wrote computed/COORDINATE_SHELL_SOURCES.png and .pdf; display only.')
