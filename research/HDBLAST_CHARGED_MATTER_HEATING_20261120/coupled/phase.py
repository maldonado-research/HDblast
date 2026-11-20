"""Test the electric-energy predictor against six post-exit switch times."""
import argparse
from bisect import bisect_left
import json
from pathlib import Path

NAMES = ('activation_4p7', 'activation_7', 'activation_9p5', 'coarse',
         'activation_11', 'activation_11p5')


def phase_table(result):
    runs = result['runs']
    params = result['parameters']
    columns = result['trace_columns']
    neutral = [dict(zip(columns, row, strict=True)) for row in runs['neutral_precursor']['trace']]
    times = [row['t'] for row in neutral]
    rows = []
    for name in NAMES:
        model = params[name]
        if (model['cells'], model['D'], model['Q'], model['nu']) != (24, 1e-4, 1e8, 6):
            raise ValueError('Cases must share the same seed, mesh, charge and drag.')
        t = model['activation']
        j = bisect_left(times, t)
        lo, hi = neutral[j-1], neutral[j]
        z = (t-lo['t']) / (hi['t']-lo['t'])
        electric = lo['X'] + z*(hi['X']-lo['X'])
        magnetic = lo['B'] + z*(hi['B']-lo['B'])
        prediction = electric**2 / (electric**2 + magnetic**2)
        charged = runs[name]
        onset_energy = charged['activation_radius']**4 * charged['activation_sheet_density']
        measured = charged['postexit_peak_work_a4'] / onset_energy
        delta = measured - prediction
        if abs(delta) > .012:
            raise ValueError(f'Phase {t} exceeds the reported 1.2-percentage-point test tolerance')
        rows.append(dict(case=name, activation=t, electric_fraction_at_activation=prediction,
                         comoving_work_fraction_at_peak=measured, difference=delta))
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', nargs='?', type=Path,
                        default=Path(__file__).resolve().with_name('CHARGED.json'))
    args = parser.parse_args()
    print(json.dumps(phase_table(json.loads(args.results.read_bytes())), indent=2))
