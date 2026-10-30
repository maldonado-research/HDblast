#!/usr/bin/env python3
"""Cross-check independent active-shell enclosures and the source report."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def interval(row):
    need(len(row) == 2, 'interval shape')
    lo, hi = map(F, row)
    need(lo <= hi, 'interval order')
    return lo, hi


def overlap(a, b):
    return max(a[0], b[0]) <= min(a[1], b[1])


def compare(detailed, coarse, source):
    need(all(d['status'] == 'PASS' and d['point'] == ['9/8', '0']
             for d in (detailed, coarse, source)), 'same successful shell point')
    drows, crows, grows = (detailed['full23_point_table'],
                           coarse['full23_point_table'], source['rows'])
    need(len(drows) == len(crows) == len(grows) == 23, 'complete modal tables')
    need(F(detailed['C1_error_exact_recomputed_upper']) ==
         F(coarse['C1_error_exact_recomputed_upper']), 'same exact C1 transfer')
    checks = 0
    for i, (d, c, g) in enumerate(zip(drows, crows, grows)):
        need(d['label'] == c['label'] == g['label'] and
             d['parity'] == c['parity'], 'row identity')
        need(interval(d['affine_shell_source']) == interval(g['interval']),
             'same authenticated source interval')
        for key in ('value', 'dx', 'dr'):
            dd, cc = interval(d['candidate'][key]), interval(c['candidate'][key])
            need(cc[0] <= dd[0] <= dd[1] <= cc[1],
                 'candidate outside independent cone box: ' + d['label'] + ' ' + key)
            checks += 1
            need(overlap(interval(d['conditional_actual_formal_response'][key]),
                         interval(c['conditional_actual_formal_response'][key])),
                 'disjoint conditional actual boxes: ' + d['label'] + ' ' + key)
            checks += 1
        need(overlap(interval(d['homogeneous_boundary_defect']),
                     interval(d['direct_homogeneous_boundary_defect'])),
             'disjoint cancellation-free and direct boundary identities')
        if d['parity'] > 0:
            need(overlap(interval(d['conditional_actual_formal_Robin_radial']),
                         interval(c['shell_Robin_actual_radial_derivative'])),
                 'disjoint actual-law Robin radial boxes: ' + d['label'])
            checks += 1
            direct = interval(d['conditional_actual_formal_response']['dr'])
            robin = interval(d['conditional_actual_formal_Robin_radial'])
            need(interval(d['conditional_actual_formal_dr_intersection']) ==
                 (max(direct[0], robin[0]), min(direct[1], robin[1])),
                 'radial intersection')
        else:
            for key in ('value', 'dx'):
                need(interval(d['candidate'][key]) == (0, 0) and
                     interval(d['conditional_actual_formal_response'][key]) == (0, 0),
                     'exact homogeneous Dirichlet trace')
            need(interval(g['interval']) == (0, 0), 'Dirichlet source zero')
    d, c = drows[12], crows[12]
    need(d['label'] == c['label'] == 'w:0', 'primary scalar mean label')
    for key in ('value', 'dx'):
        need(interval(d['conditional_actual_formal_response'][key])[0] > 0 and
             interval(c['conditional_actual_formal_response'][key])[0] > 0,
             'independent positive actual scalar mean ' + key)
    need(interval(d['conditional_actual_formal_dr_intersection'])[0] > 0 and
         interval(c['shell_Robin_actual_radial_derivative'])[0] > 0,
         'independent positive actual inward radial derivative')
    lo, hi = interval(d['conditional_actual_formal_response']['dr'])
    need(lo < 0 < hi, 'plain C1 radial interval still undecided')
    need(interval(d['full_affine_boundary_defect'])[0] > 0,
         'nonzero actual-law candidate boundary defect')
    # Readable, deliberately outward-rounded consequences of the detailed boxes.
    for key, lower, upper in [('value', F(8287*10**27), F(8294*10**27)),
                               ('dx', F(19744*10**28), F(19749*10**28))]:
        lo, hi = interval(d['conditional_actual_formal_response'][key])
        need(lower < lo <= hi < upper, 'readable conditional ' + key + ' caps')
    lo, hi = interval(d['conditional_actual_formal_dr_intersection'])
    need(F(2153*10**24) < lo <= hi < F(2155*10**24), 'readable Robin radial caps')
    lo, hi = interval(d['full_affine_boundary_defect'])
    need(F(584*10**20) < lo <= hi < F(667*10**20), 'readable candidate defect caps')
    return {'status': 'PASS', 'modal_rows': 23,
            'cross_method_interval_checks': checks,
            'independent_actual_formal_signs': 'scalar mean value, dx, inward dr positive',
            'candidate_boundary_defect': 'actual-law scalar mean defect strictly positive',
            'scope': 'one registered shell sample; conditional inherited certificates; no C2 error, full boundary solve or physical observable'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    for name in ('detailed', 'coarse', 'source'):
        ap.add_argument('--' + name, required=True, type=Path)
    args = ap.parse_args()
    print(json.dumps(compare(*(json.loads(getattr(args, name).read_text())
                               for name in ('detailed', 'coarse', 'source'))), sort_keys=True))
