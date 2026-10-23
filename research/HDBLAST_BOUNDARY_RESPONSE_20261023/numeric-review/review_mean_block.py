#!/usr/bin/env python3
"""Independent exact review of the closed five-mean comparison certificate.

No new source evaluations. Full retained source authentication precedes imports.
Connection/curvature mean matrices are reconstructed directly on covariant 5x5
tensors, and the positive resolvent is solved by determinants, independently
of the candidate's 23-amplitude extraction and Gaussian elimination.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
from itertools import permutations
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('boundary_review_auth', HERE / 'review_boundary_response.py')
auth = importlib.util.module_from_spec(spec)
spec.loader.exec_module(auth)
need = auth.need


def add(*matrices):
    return [[sum((m[i][j] for m in matrices), Q(0)) for j in range(5)] for i in range(5)]


def scale(c, m): return [[c * v for v in row] for row in m]
def absolute(m): return [[abs(v) for v in row] for row in m]
def mv(m, v): return [sum((a * b for a, b in zip(row, v)), Q(0)) for row in m]
def zero(): return [[Q(0)] * 5 for _ in range(5)]
def identity(): return [[Q(i == j) for j in range(5)] for i in range(5)]


def scalar_matrix(entries):
    m = zero()
    for i, j, value in entries: m[i][j] = Q(value)
    return m


def tensor_basis(j):
    t = zero()
    if j == 0: t[0][0] = Q(1)
    elif j == 1: t[0][1] = t[1][0] = Q(1)
    elif j == 2: t[1][1] = Q(1)
    elif j == 3:
        for k in range(2, 5): t[k][k] = Q(1)
    return t


def tensor_connection(t, direction, a, h):
    gamma = zero()
    if direction == 1:
        gamma[0][1] = gamma[1][0] = a
    else:
        gamma[0][direction] = a; gamma[1][direction] = h
        gamma[direction][0] = -a; gamma[direction][1] = h
    return [[-sum((gamma[i][k] * t[k][j] + gamma[j][k] * t[i][k]
                   for k in range(5)), Q(0)) for j in range(5)] for i in range(5)]


def tensor_curvature(t, b, d):
    sig = [1, -1, 1, 1, 1]
    def R(i, j, k, l):
        c = b if 0 in (i, j, k, l) else d
        return c * sig[i] * sig[j] * ((i == k and j == l) - (i == l and j == k))
    ric = [[sum((sig[k] * R(k, i, k, j) for k in range(5)), Q(0))
            for j in range(5)] for i in range(5)]
    return [[sum((ric[i][k] * sig[k] * t[k][j] + ric[j][k] * sig[k] * t[i][k]
                  for k in range(5)), Q(0)) / 2
             - sum((R(i, k, j, l) * sig[k] * sig[l] * t[k][l]
                    for k in range(5) for l in range(5)), Q(0))
             for j in range(5)] for i in range(5)]


def tensor_operator(fn):
    result = zero()
    for col in range(4):
        t = fn(tensor_basis(col))
        values = [t[0][0], t[0][1], t[1][1], sum(t[k][k] for k in range(2, 5)) / 3]
        reconstructed = add(*(scale(v, tensor_basis(j)) for j, v in enumerate(values)))
        need(reconstructed == t, 'DIRECT_TENSOR_MEAN_RECONSTRUCTION')
        for row, value in enumerate(values): result[row][col] = value
    return result


def direct_mean_primitives():
    def connection(t, k, a, h): return tensor_connection(t, k, Q(a), Q(h))
    def square(t, k, a, h): return connection(connection(t, k, a, h), k, a, h)
    def summed(t, kind):
        terms = []
        for k in range(2, 5):
            if kind == 'Sa': terms.append(square(t, k, 1, 0))
            elif kind == 'Sh': terms.append(square(t, k, 0, 1))
            else:
                terms.extend([connection(connection(t, k, 1, 0), k, 0, 1),
                              connection(connection(t, k, 0, 1), k, 1, 0)])
        return add(*terms)
    out = {'I': identity(), 'O': tensor_operator(lambda t: connection(t, 1, 1, 0)),
           'O2': tensor_operator(lambda t: square(t, 1, 1, 0)),
           'Rb': tensor_operator(lambda t: tensor_curvature(t, Q(1), Q(0))),
           'Rd': tensor_operator(lambda t: tensor_curvature(t, Q(0), Q(1)))}
    for kind in ('Sa', 'Sh', 'Sah'): out[kind] = tensor_operator(lambda t: summed(t, kind))
    for name in ('Delta', 'Ca', 'Ch', 'Ephi'): out[name] = zero()
    out['Aphi'] = scalar_matrix([(1, 4, 2), (4, 1, 1)])
    out['Bphi'] = scalar_matrix([(0, 4, -4), (4, 0, Q(1, 2)), (4, 2, Q(1, 2)), (4, 3, Q(-3, 2))])
    out['Pmet'] = scalar_matrix([(i, i, 1) for i in range(4)])
    out['Pg'] = scalar_matrix([(0, 4, 1), (2, 4, -1), (3, 4, 1)])
    out['Pww'] = scalar_matrix([(4, 4, 1)])
    out['W_K'] = scalar_matrix([(4, 0, 1)])
    out['W_L'] = scalar_matrix([(4, 1, 1)])
    out['By_a'] = scalar_matrix([(2, 0, -1), (3, 0, 1), (0, 0, -4), (0, 2, -2), (0, 3, 6)])
    out['By_p'] = scalar_matrix([(4, 0, Q(1, 2))])
    out['By_sigma1'] = scalar_matrix([(2, 4, Q(-1, 3)), (3, 4, Q(1, 3)), (0, 4, Q(4, 3))])
    out['By_sigma2'] = scalar_matrix([(4, 4, Q(-1, 2))])
    out['PD'] = scalar_matrix([(1, 1, 1)])
    out['PN'] = add(identity(), scale(-1, out['PD']))
    return out


def determinant(matrix):
    result = Q(0)
    for p in permutations(range(5)):
        sign = (-1) ** sum(p[i] > p[j] for i in range(5) for j in range(i + 1, 5))
        term = Q(sign)
        for i, j in enumerate(p): term *= matrix[i][j]
        result += term
    return result


def cramer(matrix, rhs):
    den = determinant(matrix); need(den != 0, 'NONZERO_RESOLVENT_DETERMINANT')
    out = []
    for col in range(5):
        numerator = [row[:] for row in matrix]
        for i in range(5): numerator[i][col] = rhs[i]
        out.append(determinant(numerator) / den)
    return out


def check_global_caps(source):
    g = source['global_cover']; division = g['division_per_unit']
    need(division == 2, 'INHERITED_GLOBAL_COVER_DIVISION')
    expected = {(i, j) for i in range(6) for j in range(16) if i + j < 16}
    rows = g['cell_enclosures']; keys = {tuple(row['index']) for row in rows}
    need(len(rows) == len(keys) == g['cells'] and keys == expected, 'WHOLE_D_SOURCE_COVER')
    for row in rows:
        i, j = row['index']
        need(tuple(map(Q, row['x'])) == (Q(i, 2), Q(i + 1, 2))
             and tuple(map(Q, row['r'])) == (Q(j, 2), Q(j + 1, 2)), 'GLOBAL_SOURCE_BOX')
    caps = g['summary']['absolute_coefficient_caps']
    names = ['metric_forcing/yy/0', 'metric_forcing/y0/0', 'metric_forcing/00/0',
             'metric_forcing/ab/0', 'scalar_forcing/0', 'QC/y/0']
    for name in names:
        def value(row):
            node = row['sources']
            for part in name.split('/'): node = node[int(part)] if isinstance(node, list) else node[part]
            need(Q(node[0]) <= Q(node[1]), 'GLOBAL_SOURCE_INTERVAL_ORDER')
            return max(abs(Q(v)) for v in node)
        need(Q(caps[name]) == max(map(value, rows)), 'GLOBAL_COMPONENT_CAP_FROM_ALL_CELLS')
    shell = source['shell_whole_time_cover']['cell_enclosures']
    intervals = [tuple(map(Q, row['x'])) for row in shell]
    need(intervals == [(Q(i, 2), Q(i + 1, 2)) for i in range(6)], 'COMPLETE_SHELL_SOURCE_COVER')
    def cap(name):
        values = [tuple(map(Q, row['sources']['one_copy_outward_hessian'][name][0])) for row in shell]
        need(all(lo <= hi for lo, hi in values), 'SHELL_SOURCE_INTERVAL_ORDER')
        return max(abs(v) for row in values for v in row)
    return [Q(caps[k]) for k in names[:5]], [cap('00'), cap('ab'), cap('scalar'), Q(caps['QC/y/0'])]


def evaluate(args):
    before = [auth.authenticate(args.source_root, auth.OLD_SHA, 169),
              auth.authenticate(args.previous_source_root, auth.PREVIOUS_SHA, 49)]
    own_bytes = Path(__file__).read_bytes(); candidate_raw = args.candidate_report.read_bytes()
    source_raw = args.source_report.read_bytes()
    # Numeric integration review is joined by exact report bytes, not by its
    # rounded display strings. Its independent geometry/source replay is separate.
    source = auth.decode(source_raw); candidate = auth.decode(candidate_raw)
    need(auth.sha(source_raw) == '6caf16c304fb6db940899d58d09b648066eb2356385189ccb79eb2b0a6879675', 'SELECTED_N16_SOURCE_REPORT_PIN')
    need(candidate['status'] == 'PASS_CONDITIONAL_MEAN_BLOCK_POSITIVE_POINT_CERTIFICATE', 'CANDIDATE_STATUS')
    op_path = args.source_root / 'projection/inputs/linear_operator.py'
    spec = importlib.util.spec_from_file_location('review_retained_linear_operator', op_path)
    op = importlib.util.module_from_spec(spec); spec.loader.exec_module(op)
    indices = [op.IDX[name] for name in ('K:0', 'L:0', 'N:0', 'J:0', 'w:0')]
    complete = op.matrices(); direct = direct_mean_primitives()
    need(set(complete) == set(direct), 'ALL_PRIMITIVES_ACCOUNTED')
    for name, matrix in complete.items():
        need(all(matrix[i][j] == matrix[j][i] == 0 for i in indices for j in range(23) if j not in indices), 'TWO_WAY_EXACT_PRIMITIVE_CLOSURE_' + name)
        need(direct[name] == [[matrix[i][j] for j in indices] for i in indices], 'DIRECT_TENSOR_PRIMITIVE_' + name)
    need([op.gram()[i] for i in indices] == [1, 2, 1, 3, 1], 'MEAN_GRAM_AND_FRAME_MULTIPLICITIES')
    ell, R, alpha, P, Ht = Q(1, 10**6), Q(131), Q(14, 125), Q(7, 50000), Q(3, 10**6)
    h = Ht / 129; bcurv = alpha**2 + Q(1, 10000); dcurv = alpha**2 + Q(1, 129**2)
    def term(c, name): return scale(c, absolute(direct[name]))
    aa = add(term(3*ell*Ht, 'I'), term(2*ell*R*alpha, 'O'), term(ell*R*P, 'Aphi'))
    bb = add(term(3*ell*R*alpha, 'I'), term(ell*R*P, 'Bphi'))
    dd = scale(ell**2*R**2, add(
        scale(alpha**2, absolute(add(direct['Sa'], scale(-1, direct['O2'])))),
        scale(alpha*h, absolute(add(direct['Sah'], scale(-3, direct['O'])))),
        term(h*h, 'Sh'), term(2*bcurv, 'Rb'), term(2*dcurv, 'Rd'),
        term(Q(4, 3), 'Pmet'), term(Q(4, 3000), 'Pg'), term(4, 'Pww'),
        term(Q(1, 500)+4*alpha*P, 'W_K'), term(3*P*h, 'W_L')))
    mm = scale(ell*R, add(term(alpha, 'By_a'), term(P, 'By_p'), term(2*P, 'By_sigma1'), term(4, 'By_sigma2')))
    sigma = add(aa, bb, scale(3, dd)); cc = add(scale(Q(3, 2), sigma), scale(3, mm))
    matrices = {'A': aa, 'B': bb, 'D0': dd, 'M': mm, 'X_contraction': cc}
    need(candidate['component_operator_majorants'] == {k: [[str(v) for v in row] for row in m] for k, m in matrices.items()}, 'INDEPENDENT_COMPONENT_MAJORANTS')
    q = max(sum(row) for row in cc); need(q < 1, 'NONNEGATIVE_MATRIX_CONTRACTION')
    actual = auth.decode((args.source_root / 'numeric/ACTUAL_SOURCES_NORMAL.json').read_bytes())
    fc, rawgc = check_global_caps(actual)
    q00, qab, qp, qc = rawgc
    gc = [2*ell*R*(q00+3*qab+qc), Q(0), 2*ell*R*q00, 2*ell*R*qab, ell*R*qp]
    need(candidate['component_source_caps'] == {'F': list(map(str, fc)), 'G': list(map(str, gc))}, 'ALL_HARMONIC_QUADRATIC_MEAN_SOURCE_CAPS')
    rhs = [Q(3, 2)*f + g for f, g in zip(fc, gc)]
    bound = cramer(add(identity(), scale(-1, cc)), rhs)
    need(all(v > 0 for v in bound) and bound == [s+c for s, c in zip(rhs, mv(cc, bound))], 'EXACT_POSITIVE_COMPONENT_RESOLVENT')
    need(candidate['component_X_bounds'] == list(map(str, bound)), 'DETERMINANT_RESOLVENT_AGREEMENT')
    eb = Q(5, 8)*mv(sigma, bound)[0]; er = Q(1, 32)*mv(mm, bound)[0]; ee = eb+er
    lo, hi = map(Q, source['signed_integral_ledger']['S_equals_KF_plus_LG_mean'])
    need(candidate['total_feedback_value_bound'] == str(ee)
         and candidate['bulk_feedback_value_bound'] == str(eb)
         and candidate['Robin_feedback_value_bound'] == str(er), 'CAUSAL_POINT_FEEDBACK_ALLOWANCE')
    need(candidate['actual_Wyy_probability_mean_interval'] == [str(lo-ee), str(hi+ee)]
         and candidate['strict_positive_margin'] == str(lo-ee) and lo-ee > 0, 'POSITIVE_CONDITIONAL_AUXILIARY_MEAN_INTERVAL')
    # Six scientifically specific in-process algebra faults. These are bounded
    # comparisons, not an exhaustive theorem about all possible bad proofs.
    faults = {
        'drop_Robin_point_allowance': ee != eb,
        'omit_Lorentzian_QK00_cap_in_GK': gc[0] != 2*ell*R*(3*qab+qc),
        'omit_full_metric_lowered_QC_cap_in_GK': gc[0] != 2*ell*R*(q00+3*qab),
        'double_one_copy_outward_QK_caps': gc[0] != 2*ell*R*(2*q00+6*qab+qc),
        'omit_value_norm_factor_three_in_bulk_feedback': sigma != add(aa, bb, dd),
        'couple_nonmean_mode_into_mean': all(complete['Ca'][indices[0]][j] == 0 for j in range(23) if j not in indices),
    }
    need(all(faults.values()), 'MEAN_SPECIFIC_FALSE_FORMULA_DISTINCTIONS')
    after = [auth.authenticate(args.source_root, auth.OLD_SHA, 169),
             auth.authenticate(args.previous_source_root, auth.PREVIOUS_SHA, 49)]
    need(before == after and own_bytes == Path(__file__).read_bytes()
         and candidate_raw == args.candidate_report.read_bytes() and source_raw == args.source_report.read_bytes(), 'READ_ONLY_INPUTS_AND_PROGRAM_STABLE')
    return {'status': 'PASS_INDEPENDENT_CONDITIONAL_FIVE_MEAN_BLOCK_REVIEW',
            'source_evaluations': 0, 'all_retained_inputs_authenticated_before_import_and_after': [v[0] for v in before],
            'mean_primitive_direct_tensor_reconstructions': sorted(direct),
            'two_way_primitive_closure': True, 'mean_gram': [1, 2, 1, 3, 1],
            'comparison_resolvent_method': 'Exact 5x5 determinants and Cramer rule; nonnegative maximum row sum <1',
            'maximum_row_sum': str(q), 'bulk_feedback': str(eb), 'Robin_feedback': str(er), 'total_feedback': str(ee),
            'strict_positive_auxiliary_mean_margin': str(lo-ee),
            'auxiliary_mean_interval': [str(lo-ee), str(hi+ee)],
            'specific_false_formula_distinctions': faults,
            'candidate_report_sha256': auth.sha(candidate_raw), 'source_integral_sha256': auth.sha(source_raw),
            'reviewer_source_sha256': auth.sha(own_bytes),
            'analytic_premises_reviewed': [
                'SO(4) averaging commutes with all primitives and boundary projectors; first-field quadratic sources are projected after their complete construction.',
                'Component norm Xi=max(sup|Wi|/3,sup|d_x Wi|,sup|d_r Wi|); this is not an H-norm with omitted Gram factors.',
                'Retarded maps on the backward-closed domain obey Xi(Kf)<=3/2 sup|fi| and Xi(Lg)<=sup|gi|.',
                'Zero source wedge and flat corner preserve Wi=0 when x+r<1 under Picard iteration; inherited actual smooth C1 solution supplies physical-side boundary traces.',
                'Boundary trace |Wi(s,0)|<=(s-1)bi for 1<=s<=5/4 gives the Robin integral factor1/32; weighted active bulk area is5/8.',
                'All whole-collar background bounds and inherited source enclosures remain conditional retained premises; this review does not reprove those source PDE certificates.'],
            'limits': ['Internal independent review, not external specialist review or a formal proof assistant certificate.',
                       'One conditional auxiliary forced-increment mean at one point; no claim about total physical U, force, nonlinear integrability, a waveform, cosmological origin or novelty.']}


if __name__ == '__main__':
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode, 'USE_PYTHON_I_S_B')
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ('source-root', 'previous-source-root', 'candidate-report', 'source-report', 'output'):
        ap.add_argument('--'+name, type=Path, required=True)
    args = ap.parse_args()
    need(args.output.is_absolute() and not args.output.exists(), 'NEW_REVIEW_OUTPUT')
    result = evaluate(args)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(result['status'])
