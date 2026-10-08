"""Plot already certified JSON summaries; no source or array evaluation."""
from fractions import Fraction
from pathlib import Path
import csv
import argparse
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def exact(value):
    return Fraction(int(value['numerator']), int(value['denominator']))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    DATA = args.data
    OUT = args.output
    OUT.mkdir(exist_ok=True)
    raw = DATA.read_bytes()
    data = json.loads(raw)
    radii = {key: Fraction(value) for key, value in data['max_complete_target_export_L1_radius'].items()}
    fields = [('suprema', 'delta_U_L1_upper', r'Incoming $\delta U$: maximum L1 norm upper bound'),
              ('suprema', 'delta_W_L1_upper', r'Incoming $\delta W$: maximum L1 norm upper bound'),
              ('sums', 'R_uniform_abs_upper', 'Finite weighted density: uniform upper bound'),
              ('sums', 'P_uniform_abs_upper', 'Finite weighted pressure: uniform upper bound')]
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.0), constrained_layout=True)
    colors = {'positive_B': '#2166ac', 'signed_uB': '#d95f02'}
    labels = {'positive_B': r'$h=B$', 'signed_uB': r'$h=zB$'}
    for ax, (group, key, title) in zip(axes.flat, fields):
        for capsule in ('positive_B/coarse', 'positive_B/fine', 'signed_uB/coarse', 'signed_uB/fine'):
            source, grid = capsule.split('/')
            rows = [row for row in data['rows'] if row['capsule'] == capsule]
            x = [float(exact(row['K'])) for row in rows]
            y = [float(exact(row[group][key])) for row in rows]
            ax.plot(x, y, marker='o', markersize=4, linewidth=1.5,
                    linestyle='-' if grid == 'coarse' else '--', color=colors[source],
                    label=labels[source] + ', ' + grid)
        ax.set_yscale('log')
        ax.set_xticks([64, 128, 256])
        ax.set_xlabel('Fixed momentum cutoff K')
        ax.set_ylabel('Certified upper bound')
        ax.set_title(title, fontsize=10)
        ax.grid(True, which='both', alpha=.2)
    axes[0, 0].legend(fontsize=8, loc='best')
    fig.suptitle('Exact retained-state comparison in the prescribed scalar model\n'
                 '49,152 capsule-node occurrences; 12 fixed weighted prefixes', fontsize=12)
    fig.savefig(OUT / 'stored_state_bounds.png', dpi=180)
    fig.savefig(OUT / 'stored_state_bounds.pdf', metadata={'Title': 'HDBLAST exact retained-state bounds',
                                                        'CreationDate': None, 'ModDate': None})
    plt.close(fig)
    columns = ['capsule', 'K', 'nodes', 'U_L1_sup_lower', 'U_L1_sup_upper',
               'W_L1_sup_lower', 'W_L1_sup_upper', 'c_sup', 'first_order_Wronskian_sup',
               'R_uniform_upper', 'P_uniform_upper', 'work_upper', 'integral_pressure_upper']
    canon = lambda q: f'{q.numerator}/{q.denominator}'
    with (OUT / 'EXACT_TABLE.csv').open('x', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        for row in data['rows']:
            s = row['suprema']; sums = row['sums']
            u = exact(s['delta_U_L1_upper']); w = exact(s['delta_W_L1_upper'])
            writer.writerow({'capsule': row['capsule'], 'K': canon(exact(row['K'])), 'nodes': row['node_count'],
                             'U_L1_sup_lower': canon(max(Fraction(0), u - 2 * radii['U'])), 'U_L1_sup_upper': canon(u),
                             'W_L1_sup_lower': canon(max(Fraction(0), w - 2 * radii['W'])), 'W_L1_sup_upper': canon(w),
                             'c_sup': canon(exact(s['c_abs'])),
                             'first_order_Wronskian_sup': canon(exact(s['first_order_Wronskian_defect_abs'])),
                             'R_uniform_upper': canon(exact(sums['R_uniform_abs_upper'])),
                             'P_uniform_upper': canon(exact(sums['P_uniform_abs_upper'])),
                             'work_upper': canon(exact(sums['work_abs_upper'])),
                             'integral_pressure_upper': canon(exact(sums['int_pressure_abs_upper']))})
    receipt = {'status': 'PASS_PLOT_OF_CERTIFIED_JSON_SUMMARIES',
               'input_DATA_sha256': hashlib.sha256(raw).hexdigest(),
               'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'figure_values': 'Binary64 display only, from exact rational certified upper bounds',
               'exact_table': 'Exact rational strings; max(0,S-2R) lower rule uses global target L1 radius',
               'units': 'U=u/epsilon, W=w/epsilon; inherited a0^4/epsilon scaled linear R/P response',
               'time_interval': '[-9/2,-7/2] for uniform finite weighted state contribution',
               'source_callbacks': 0, 'original_array_decodes': 0,
               'continuous_momentum_or_full_stress_certificate': 'UNRESOLVED'}
    (OUT / 'FIGURE_RECEIPT.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(receipt['status'])

if __name__ == '__main__':
    main()
