"""Recompute both canonical-memory bounds from the saved directional replay slices."""
import argparse
import json
from pathlib import Path

import numpy as np

from memory_bound import measure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--slices',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True,help='Pinned upstream kinetics/SWEEP.json')
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    references = json.loads(args.reference.read_text())['runs']
    result = {}
    for case in ('mass_4e10','mass_4e12'):
        ref = references[case]
        result[case] = {}
        for grid in ('native','half_step'):
            traces = np.fromfile(args.slices/f'{case}-{grid}.bin',dtype='=f8').reshape(-1,12)
            for label,scale in (('witness',1.),('floor',ref['reference_drag_floor']/ref['action'])):
                result[case][f'{grid}_{label}'] = measure(traces,ref['mu'],scale)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
