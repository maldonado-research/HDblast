"""Compare the entire pole-driven wall history after halving both thermal table spacings."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
METRICS = (
    'hot_radius_factor', 'peak_hot_fraction', 'max_nu', 'integral_nu',
    'end_rest_drag_work', 'cool_or_fast_absolute_electrical_work_fraction',
    'first_hot_time', 'maximum_radius', 'doubling_theta',
)


def compare(coarse, fine):
    if any(coarse[key] != fine[key] for key in ('action', 'daughter_mass_L', 'pole_mass_over_daughter_mass', 'width_over_mass')):
        raise ValueError('Different physical parameters in grid comparison')
    output = {}
    for name in ('pole_cap', 'pole_cap_half', 'pole_vacancy', 'pole_vacancy_half'):
        base = coarse['runs'][name]
        refined = fine['runs'][name]
        output[name] = {
            metric: {'coarse': base[metric], 'fine': refined[metric],
                     'relative_change': refined[metric] / base[metric] - 1}
            for metric in METRICS
        }
        if abs(output[name]['hot_radius_factor']['relative_change']) > .002:
            raise ValueError(f'{name} hot interval is sensitive to rate-table resolution')
    return {
        'coarse': {'log10_temperature_step': .125, 'scaled_rapidity_step': .05, 'quadrature_order': 16},
        'fine': {'log10_temperature_step': .0625, 'scaled_rapidity_step': .025, 'quadrature_order': 20},
        'runs': output,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=HERE.parent / 'prior')
    parser.add_argument('--feedback', type=Path, default=HERE.parent / 'feedback')
    parser.add_argument('--coarse', type=Path, default=HERE / 'EVOLUTION.json')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1')
    with tempfile.TemporaryDirectory(prefix='hdblast-grid-check-') as work:
        table = Path(work) / 'fine.dat'
        history = Path(work) / 'fine.json'
        subprocess.run([sys.executable, '-I', '-B', str(HERE / 'thermal.py'),
                        '--output', str(table), '--log-step', '.0625', '--dy', '.025',
                        '--order', '20', '--jobs', str(args.jobs)], env=env, check=True)
        subprocess.run([sys.executable, '-I', '-B', str(HERE / 'evolve.py'),
                        '--base', str(args.base), '--feedback', str(args.feedback), '--table', str(table),
                        '--output', str(history), '--jobs', str(args.jobs)], env=env, check=True)
        report = compare(json.loads(args.coarse.read_text()), json.loads(history.read_text()))
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
