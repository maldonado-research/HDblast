"""Regenerate the prebirth trajectory and directly record the incident radial characteristic."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np

COLUMNS = ('t', 'a', 'H', 'K', 'E', 'B', 'current', 'incoming', 'phi', 'D')


def capture(prior, out, jobs):
    prior = prior.resolve()
    out.mkdir(parents=True, exist_ok=True)
    library = runpy.run_path(str(prior / 'bulk/incoming.py'))
    solver_inputs = prior / 'parent/prior'
    meta = json.loads((prior / 'history/markers.json').read_text())
    config = json.loads((solver_inputs / 'prior/envelope/ENVELOPE.json').read_text())['parameters']['mass_4e12']
    run = runpy.run_path(str(solver_inputs / 'prior/coupled/joined.py'))['run']
    with tempfile.TemporaryDirectory() as build:
        tmp = Path(build)
        exe = library['capture_solver'](solver_inputs, tmp)
        def run_one(label):
            spec = meta['runs'][label]
            with np.load(prior / 'history' / (label + '.npz')) as tab:
                end = float(tab['t'][-1])
            target = (out / (label + '.bin')).resolve()
            args = ['2', repr(meta['action']), str(solver_inputs/'boost/BOOST_TABLE.dat'), '1',
                    repr(spec['channel_fraction']), str(target)]
            wrapper = tmp / label
            wrapper.write_text('#!/bin/sh\nexec ' + shlex.quote(str(exe)) + ' "$@" ' + shlex.join(args) + '\n')
            wrapper.chmod(0o755)
            params = dict(config, cfl=config['cfl'] * (.5 if label.endswith('half') else 1),
                          sample=0, end=end+.0002)
            try:
                run(wrapper, params)
            except subprocess.CalledProcessError as err:
                raise RuntimeError(err.stderr) from err
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            list(pool.map(run_one, meta['runs']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--jobs', default=2, type=int)
    args = parser.parse_args()
    capture(args.prior, args.output_dir, args.jobs)
