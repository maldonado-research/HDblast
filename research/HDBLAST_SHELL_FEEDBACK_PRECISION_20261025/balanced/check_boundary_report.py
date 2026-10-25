#!/usr/bin/env python3
"""Independent exact-rational replay of the stored boundary experiment.

This reads fixed certificates and P1 data. It does not evaluate profile jets,
repeat quadrature, or execute any frozen producer.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

SOURCE_PIN = 'c01a8f5f95737dcacebdc26379ae0ee5863af9c1977890b3901ad6113103e619'
REPORT_PIN = '36cdd588603c2c15226a6793430d8ed3ff3f2101e94afca1f8744a73fce9e276'
SEGMENTS = ((F(65, 64), F(66, 64), 16),
            (F(66, 64), F(9, 8), 112),
            (F(9, 8), F(5, 4), 128))


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interval(row):
    lo, hi = map(F, row)
    need(lo <= hi, 'INTERVAL_ORDER')
    return lo, hi


def p1_integrals(nodes, values, x):
    """Direct segment integrals: integral v and integral (x-r)*v.

    This deliberately does not import the frozen hinge/primitive kernel.
    """
    first, second = [F(0)] * 5, [F(0)] * 5
    for a, b, va, vb in zip(nodes, nodes[1:], values, values[1:]):
        if x <= a:
            break
        width = min(b, x)-a
        for j in range(5):
            slope = (vb[j]-va[j])/(b-a)
            contribution = va[j]*width+slope*width**2/2
            first[j] += contribution
            second[j] += ((x-a)*contribution
                          - va[j]*width**2/2-slope*width**3/3)
    return first, second


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--source', type=Path, default=Path(__file__).with_name('exact_boundary.py'))
    parser.add_argument('--report', type=Path, default=Path(__file__).with_name('EXACT_BOUNDARY_NORMAL.json'))
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('INDEPENDENT_BOUNDARY_REVIEW.json'))
    args = parser.parse_args()
    need(not args.output.exists(), 'FRESH_REVIEW_OUTPUT_REQUIRED')
    need(sha(args.source) == SOURCE_PIN, 'REVIEWED_SOURCE_CHANGED')
    need(sha(args.report) == REPORT_PIN, 'REVIEWED_REPORT_CHANGED')
    data = json.loads(args.report.read_text())
    current = args.repo_root.resolve()/'research/HDBLAST_COUPLED_MEAN_PILOT_20261025'
    profile_path = current/'results/PROFILES.json'
    evaluated_path = current/'results/EVALUATED_CANDIDATE.json'
    need(current not in args.output.resolve().parents, 'OUTPUT_INSIDE_FROZEN_CAPSULE')
    manifest = json.loads((current/'MANIFEST.json').read_text())
    for path in (profile_path, evaluated_path):
        rel = str(path.relative_to(current))
        record = next(row for row in manifest['files'] if row['path'] == rel)
        need(sha(path) == record['sha256'] and path.stat().st_size == record['bytes'],
             'FROZEN_DATA_MEMBER_HASH:'+rel)
    need(sha(current/'MANIFEST.json') == data['authentication']['current']['manifest_sha256'],
         'FROZEN_MANIFEST_HASH')
    profile = json.loads(profile_path.read_text())
    evaluated = json.loads(evaluated_path.read_text())
    need(data['schema'] == 'HDBLAST_EXACT_BOUNDARY_CANDIDATE_QUADRATURE_V1', 'SCHEMA')
    need(data['status'] == 'EVALUATED_NEW_EXACT_BOUNDARY_CANDIDATE_ONLY', 'SCOPE_STATUS')
    need(data['complete_closures_authenticated_before_and_after'] is True, 'PRODUCER_AUTHENTICATION_DECLARATION')
    need(data['arithmetic']['fixed_subintervals'] == 256, 'FIXED_SUBINTERVALS')
    need(data['arithmetic']['distinct_point_evaluations'] == 257, 'FIXED_POINT_COUNT')
    rows = data['panel_certificates']
    expected = []
    for a, b, count in SEGMENTS:
        step = (b-a)/count
        expected += [(a+k*step, a+(k+2)*step) for k in range(0, count, 2)]
    need(len(rows) == len(expected) == 128, 'PANEL_COUNT')
    for row, (left, right) in zip(rows, expected):
        a, b = map(F, row['r'])
        need((a, b) == (left, right), 'FIXED_PANEL_COVERAGE')
        need(F(row['midpoint']) == (a+b)/2, 'SIMPSON_MIDPOINT')
        derivative = interval(row['fourth_derivative_interval'])
        need(F(row['remainder_upper']) == (b-a)**5*max(map(abs, derivative))/2880,
             'SIMPSON_REMAINDER_FORMULA')
        interval(row['simpson_interval'])

    nodes = list(map(F, profile['nodes']))
    uu = [list(map(F, row)) for row in profile['profiles']['u']]
    vv = [list(map(F, row)) for row in profile['profiles']['v']]
    M = [list(map(F, row)) for row in profile['matrices']['M']]
    need(all(M[3][j] == 0 for j in (1, 2, 3)), 'J_ROBIN_ROW_STRUCTURE')
    need(M[3][4] == F(data['definitions']['m']), 'M_JW')
    need(M[3][0] == F(data['definitions']['M_JK']), 'M_JK')
    need(all(row[0] == 0 for row in uu), 'OLD_U_K_ZERO')
    results = []
    for sample in data['samples']:
        x = F(sample['x'])
        need(x in (F(9, 8), F(5, 4)) and F(sample['r']) == 0, 'SAMPLE_LOCATION')
        used = [row for row in rows if F(row['r'][1]) <= x]
        error = sum((F(row['remainder_upper']) for row in used), F(0))
        need(error == F(sample['simpson_remainder_upper']), 'CUMULATIVE_REMAINDER')
        lo = sum((F(row['simpson_interval'][0]) for row in used), F(0))
        hi = sum((F(row['simpson_interval'][1]) for row in used), F(0))
        sl, sh = interval(sample['simpson_sum_interval'])
        need(sl <= lo <= hi <= sh, 'OUTWARD_SIMPSON_SUM')
        il, ih = interval(sample['m_Cprime_f_integral_interval'])
        need(il <= sl-error and ih >= sh+error, 'OUTWARD_INTEGRAL')

        pu, unused = p1_integrals(nodes, uu, x)
        pv, qv = p1_integrals(nodes, vv, x)
        old = next(row for row in evaluated['samples'] if F(row['x']) == x and F(row['r']) == 0)
        old_feedback = F(old['feedback_exact']['value'][3])
        old_bulk = old_feedback + sum((M[3][j]*(pu[j]+qv[j]) for j in range(5)), F(0))
        boundary_velocity = -sum((M[3][j]*qv[j] for j in range(5)), F(0))
        need(old_bulk == F(sample['old_bulk_J_exact']), 'INDEPENDENT_OLD_BULK_RECONSTRUCTION')
        need(pv[3] == F(sample['old_velocity_primitive_J_exact']), 'INDEPENDENT_VELOCITY_PRIMITIVE')
        need(boundary_velocity == F(sample['old_boundary_velocity_J_exact']),
             'INDEPENDENT_BOUNDARY_VELOCITY')
        velocity_and_bulk = pv[3]+old_bulk+boundary_velocity
        need(velocity_and_bulk == F(sample['old_velocity_and_bulk_exact']), 'EXACT_VELOCITY_BULK_ASSEMBLY')
        el, eh = interval(sample['balanced_endpoint_interval'])
        cl, ch = interval(sample['new_J_candidate_interval'])
        need(cl <= el+il+velocity_and_bulk and ch >= eh+ih+velocity_and_bulk,
             'OUTWARD_FINAL_CANDIDATE_ASSEMBLY')
        ol, oh = interval(sample['old_J_candidate_interval'])
        original_lo, original_hi = interval(old['complete_candidate_interval_evaluation']['value'][3])
        need(ol <= original_lo <= original_hi <= oh, 'OLD_CANDIDATE_ENCLOSED')
        dl, dh = interval(sample['candidate_change_interval'])
        need(dl <= cl-oh and dh >= ch-ol, 'CANDIDATE_CHANGE_ENCLOSED')
        results.append({'x': str(x), 'r': '0', 'panels': len(used),
                        'new_candidate_interval': [str(cl), str(ch)],
                        'candidate_interval_strictly_negative': ch < 0,
                        'quadrature_remainder_upper': str(error),
                        'all_exact_arithmetic_checks_passed': True})
    need(len(results) == 2 and {row['x'] for row in results} == {'9/8', '5/4'}, 'TWO_SAMPLES')
    result = {'schema': 'HDBLAST_INDEPENDENT_EXACT_BOUNDARY_REVIEW_V1',
              'status': 'PASS_INDEPENDENT_BOUNDARY_CERTIFICATE_REPLAY',
              'source_sha256': SOURCE_PIN, 'report_sha256': REPORT_PIN,
              'checker_sha256': sha(Path(__file__)),
              'frozen_profiles_sha256': sha(profile_path),
              'frozen_evaluated_candidate_sha256': sha(evaluated_path),
              'arithmetic': 'Python Fraction; no floating-point decisions',
              'checks': ['Fixed 256 subintervals / 128 Simpson panels and full ordered coverage.',
                         'Each stored remainder equals width^5 times the stored fourth-derivative magnitude divided by 2880.',
                         'Outward panel sums, cumulative remainder, integral enclosure and final candidate assembly.',
                         'Independent direct P1 segment integration reconstructs frozen velocity primitives and boundary velocity.',
                         'Frozen bulk reconstructed from old feedback plus old boundary primitive.',
                         'Old candidate and candidate-change intervals enclosed.'],
              'manual_source_review': ['AST extraction preserves the nominal f expressions and zero collar.',
                                       'Point jet order 2 leaves f order 0; panel jet order 6 leaves f order 4.',
                                       'The factor 4! ell^4 converts the fourth z jet coefficient to the fourth r derivative.',
                                       'The balanced identity has the correct sign; nominal u_K is identically zero.',
                                       'Splitting at the cutoff permits separate one-sided fourth-derivative bounds.'],
              'samples': results,
              'limits': ['Candidate evaluation review only; no actual-response or intrinsic-observable sign is certified.',
                         'This replay checks stored derivative bounds algebraically; it does not recompute profile jets or quadrature.',
                         'No frozen producer or old numerical stage was executed.']}
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({'status': result['status'], 'output': str(args.output),
                      'sha256': sha(args.output)}, sort_keys=True))


if __name__ == '__main__':
    main()
