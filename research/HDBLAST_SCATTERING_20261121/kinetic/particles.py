"""Rebuild charge-tagged particle samples at the early marker of HDblast study 198."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np


def rebuild(study, fields, out, jobs, powers):
    feedback = runpy.run_path(str(study / 'feedback/currents.py'))
    orbits = feedback['pairs']
    inputs = orbits['read_inputs'](study / 'parent/parent/history', fields)
    out.mkdir(parents=True, exist_ok=True)
    labels = ('dedicated_half', 'shared_half')
    jobspec = [(label, case, power) for power in powers for label in labels
               for case in ('pressure_fit', 'small_pulse')]
    with tempfile.TemporaryDirectory(prefix='collision-input-') as tmpname:
        tmp = Path(tmpname)
        exe = tmp / 'homogeneous'
        subprocess.run(['g++', '-O3', '-std=c++17', str(study / 'feedback/homogeneous.cpp'),
                        '-o', str(exe)], check=True)

        def one(job):
            label, case, power = job
            seed = orbits['SCRAMBLE'] + {12: 1, 13: 71, 14: 121}[power]
            if case == 'pressure_fit':
                spec = orbits['uniform_pair'](orbits['FIT'][label][1])
            else:
                spec = (2, 2, 'early', 'late_birth', .04/np.sqrt(2), .04/np.sqrt(2), np.pi/2)
            stats, p = feedback['integrate'](inputs, label, case, spec, power, 512,
                                             seed, True, exe, tmp)
            sign = np.tile([1., -1.], 2**power)
            h, _, m = inputs[label]
            drift = np.interp(m['time'], h['t'], h['u'])
            file = out / f'{label}_{case}_{power}.npz'
            births = orbits['seeds'](h, m, power, seed)[:, 0]
            birth_time = np.repeat(np.sort(births, kind='stable'), 2)
            if np.any(birth_time[1:] < birth_time[:-1]):
                raise ArithmeticError('Particle records must be sorted by birth.')
            np.savez(file, p=p, sign=sign, born=birth_time, theta=m['theta'], drift=drift, a=m['a'],
                     time=m['time'], D=m['D'], seed=seed, steps=512,
                     reference_tv=stats['coarse_cell_total_variation'])
            print(file.name, stats['coarse_cell_total_variation'], flush=True)

        with ThreadPoolExecutor(jobs) as pool:
            list(pool.map(one, jobspec))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--fields', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--powers', type=int, nargs='+', default=[12,13])
    args = parser.parse_args()
    rebuild(args.study.resolve(), args.fields.resolve(), args.out.resolve(), args.jobs, args.powers)
