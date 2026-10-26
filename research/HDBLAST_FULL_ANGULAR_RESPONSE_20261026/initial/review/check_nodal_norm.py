#!/usr/bin/env python3
"""Compare full23 Gram-normalized initial radial data with the retained cap."""
import argparse
import importlib.util
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_REPO = Path('/workspace/hdblast-research-work/shell-observable-20261025/repo')


def sqrt_lower(value):
    grid = 1 << 112
    return F(math.isqrt((value*grid*grid).__floor__()), grid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=DEFAULT_REPO)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('nodal_norm_full23', ROOT/'full23_initial.py')
    producer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(producer)
    context = producer.Context(args.repo_root)
    profiles = json.loads((ROOT/'results/INITIAL_PROFILES.json').read_text())
    report = json.loads((ROOT/'results/INITIAL_REPORT.json').read_text())
    weights = context.weights
    squares = [sum((weight*F(entry)**2 for weight, entry in zip(weights, row)), F(0))
               for row in profiles['profiles']['du']]
    index = max(range(len(squares)), key=squares.__getitem__)
    r = F(profiles['nodes'][index])
    nodal_nominal = producer.sqrt_upper(squares[index])
    # P1 error globally bounds every nodal coefficient roundoff error too.
    nodal_actual_upper = (nodal_nominal+F(report['P1_H_error_caps'][1])
                          +F(report['total_H_caps'][1]))
    actual = [item.a.val() for item in context.paired_profiles(context.I(r))[1]]
    upper_square = sum((weight*box.mag()**2 for weight, box in zip(weights, actual)), F(0))
    lower_square = sum((weight*(0 if box.lo <= 0 <= box.hi else min(abs(box.lo), abs(box.hi)))**2
                        for weight, box in zip(weights, actual)), F(0))
    point_lower, point_upper = sqrt_lower(lower_square), producer.sqrt_upper(upper_square)
    producer.need(point_lower <= point_upper < F(103*10**31), 'POINT_NORM_BELOW_RETAINED_CAP')
    producer.need(nodal_actual_upper < F(103*10**31), 'ALL_NODAL_NORMS_BELOW_RETAINED_CAP')
    context.reauthenticate()
    result = {
        'status': 'PASS_FULL23_ACTUAL_NODAL_RADIAL_NORM_COMPARISON',
        'r_at_max_nominal_nodal_H_radial_norm': str(r),
        'max_nominal_nodal_H_radial_upper': str(nodal_nominal),
        'actual_point_H_radial_norm': [str(point_lower), str(point_upper)],
        'actual_point_H_radial_norm_display': [float(point_lower), float(point_upper)],
        'all_actual_nodes_H_radial_upper': str(nodal_actual_upper),
        'all_actual_nodes_H_radial_upper_display': float(nodal_actual_upper),
        'retained_response_X_cap': str(F(103*10**31)),
        'same_normalization': 'Authenticated full23 probability sphere Gram; U_r=ell*dU/dz=-ell*rho*dU/dy; positive frame Frobenius metric plus scalar norm.',
        'scope': 'Nodal check and one actual point enclosure, not a new global response norm proof. The1.248e33 whole-cell nominal radial upper can be loose.',
        'complete_roots_authenticated_before_and_after': True,
        'producer_sha256': producer.sha(ROOT/'full23_initial.py'),
    }
    (HERE/'NODAL_NORM_REVIEW.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
