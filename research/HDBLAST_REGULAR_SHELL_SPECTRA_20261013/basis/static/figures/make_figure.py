#!/usr/bin/env python3
"""Visualize exact analytic bounds; plotting floats never decide a theorem gate."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from fractions import Fraction
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    certificate_path = root / 'bounds/RESCALED_BOUND_CERTIFICATE.json'
    baseline_path = root / 'bounds/public_baseline/EXACT_UNIFORM_CONSTANTS.json'
    certificate = json.loads(certificate_path.read_text())
    baseline = json.loads(baseline_path.read_text())
    primary = certificate['primary_uniform_certificate']
    delta_star = Fraction(primary['delta_star'])
    old_star = Fraction(baseline['models'][0]['delta_star'])
    coefficient = Fraction(primary['exact_constants']['suppression_q'])
    remainder = Fraction(primary['K_H_best'])
    margin = coefficient - remainder * delta_star
    if not margin > Fraction('0.000423'):
        raise RuntimeError('The certificate does not support the plotted lower bound')
    args.output.mkdir(parents=True, exist_ok=True)
    points = [(delta_star * i / 200) for i in range(201)]
    rows = [(d, coefficient - remainder*d, coefficient + remainder*d) for d in points]
    with (args.output / 'controlled-response.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['delta_exact', 'normalized_lower_exact', 'normalized_upper_exact',
                         'delta_display', 'lower_display', 'upper_display'])
        for d, lo, hi in rows:
            writer.writerow([str(d), str(lo), str(hi), float(d), float(lo), float(hi)])
    metadata = {
        'scope': 'Analytic enclosure for the specified static local branch; not data or numerical trajectories.',
        'certificate_sha256': hashlib.sha256(certificate_path.read_bytes()).hexdigest(),
        'baseline_sha256': hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
        'old_delta_star': str(old_star), 'new_delta_star': str(delta_star),
        'endpoint_ratio_exact': str(delta_star / old_star),
        'endpoint_ratio_display': float(delta_star / old_star),
        'q_exact': str(coefficient), 'K_H_exact': str(remainder),
        'minimum_lower_bound_exact': str(margin),
        'zero_endpoint': 'A limiting coefficient only; the theorem domain is delta > 0.',
    }
    (args.output / 'figure-data.json').write_text(json.dumps(metadata, indent=2) + '\n')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'figure.dpi': 180})
    fig, (left, right) = plt.subplots(1, 2, figsize=(11.5, 4.5),
                                     gridspec_kw={'width_ratios': [1, 1.35]})
    left.scatter([float(old_star), float(delta_star)], [0, 1], s=65,
                 c=['#68717d', '#1768ac'], zorder=3)
    left.set_xscale('log')
    left.set_xlim(2e-8, 1e-2)
    left.set_ylim(-.5, 1.6)
    left.set_yticks([0, 1], ['8 Oct\nbaseline', '11 Oct\nimproved'])
    left.text(float(old_star)*1.35, 0, f'{float(old_star):.3g}', va='center')
    left.text(float(delta_star), 1.22, '0.002', ha='center', color='#1768ac')
    left.annotate('', xy=(float(delta_star), .6), xytext=(float(old_star), .6),
                  arrowprops={'arrowstyle': '<->', 'color': '#1768ac'})
    left.text((float(old_star)*float(delta_star))**.5, .78,
              f'{float(delta_star/old_star):,.0f}× larger endpoint', ha='center', fontsize=9)
    left.set_xlabel(r'Upper detuning endpoint $\bar\delta$ (log scale)')
    left.set_title('Proved sufficient ranges', loc='left', fontweight='bold')
    left.grid(axis='x', alpha=.2)

    xx = [float(d) for d, _, _ in rows]
    lower = [float(lo) for _, lo, _ in rows]
    upper = [float(hi) for _, _, hi in rows]
    right.fill_between(xx, lower, upper, color='#a8d4eb', alpha=.65,
                       label=r'Exact enclosure $q\pm K_H\delta$')
    right.plot(xx, lower, color='#1768ac', lw=1.8)
    right.plot(xx, upper, color='#1768ac', lw=1.8)
    right.axhline(float(coefficient), color='#68717d', linestyle='--', lw=1,
                  label=r'Quadratic coefficient $q=c^2/384$')
    right.axhline(0, color='#333333', lw=1)
    right.set_xlim(0, float(delta_star))
    right.set_ylim(-.00007, .00165)
    right.set_xlabel(r'Detuning $\delta$ (zero shown as a limit)')
    right.set_ylabel(r'$(H_{\rm metric}^2-H^2)/\delta^2$')
    right.set_title('Curvature reduction stays positive', loc='left', fontweight='bold')
    right.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    right.xaxis.set_major_locator(ticker.MultipleLocator(.0005))
    right.text(.00082, .00019, r'Lower bound $>0.000423$ throughout', fontsize=9)
    right.legend(loc='upper left', fontsize=8, frameon=False)
    fig.suptitle('Controlled scalar boundary response', x=.04, ha='left', fontsize=16,
                 fontweight='bold')
    fig.text(.04, .018, 'Conditional static theorem • specified local tube • no stability or cosmological inference',
             fontsize=9, color='#515962')
    fig.tight_layout(rect=(0, .055, 1, .93), w_pad=2.4)
    fig.savefig(args.output / 'controlled-response.png')
    fig.savefig(args.output / 'controlled-response.pdf', metadata={
        'Title': 'Controlled scalar boundary response', 'Author': 'Ricardo Maldonado',
        'CreationDate': None, 'ModDate': None})
    plt.close(fig)


if __name__ == '__main__':
    main()
