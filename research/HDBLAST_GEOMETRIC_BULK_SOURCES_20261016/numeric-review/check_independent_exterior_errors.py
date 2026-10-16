#!/usr/bin/env python3
"""Independent exact adjugate replay of backward-closed exterior errors."""
import argparse
from fractions import Fraction as F
from math import comb
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NUMERIC = ROOT.parent/'numeric'


def need(ok, text):
    if not ok:
        raise ValueError(text)


def mv(matrix, vector):
    return [sum((a*b for a, b in zip(row, vector)), F(0)) for row in matrix]


def plus(*vectors):
    return [sum(row, F(0)) for row in zip(*vectors)]


def times(c, vector):
    return [c*v for v in vector]


def adjugate(matrix):
    det = sum((matrix[0][j]*(matrix[1][(j+1)%3]*matrix[2][(j+2)%3]-matrix[1][(j+2)%3]*matrix[2][(j+1)%3]) for j in range(3)), F(0))
    need(det != 0, 'nonzero independent adjugate determinant')
    return [[(matrix[(j+1)%3][(i+1)%3]*matrix[(j+2)%3][(i+2)%3]-matrix[(j+1)%3][(i+2)%3]*matrix[(j+2)%3][(i+1)%3])/det for j in range(3)] for i in range(3)]


def derive(data, lift):
    master, co, bump = (data[k] for k in ('master', 'coefficients', 'bump'))
    ell, x, theta = F(1, 10**6), F(3), F(3, 2)
    m = master['exact']['interior']
    Q = list(map(F, bump['dimensionless_Q_derivative_sup']))
    velocity = [x*ell, ell, ell*(Q[1]+ell), ell*(Q[2]+2*ell*Q[1]+ell**2)]
    kernel = [x*x*F(m['forcing']['f'])/2, F(m['first_error']), max(map(F, m['second_error'].values())), max(map(F, m['third_error'].values()))]
    ew = [sum(row, F(0)) for row in zip(velocity, lift, kernel)]
    c = 6*ell/130
    comparison = [5*c]+[c*(5*Q[k]+k*Q[k-1]) for k in range(1, 4)]
    actual = list(map(F, master['source_bound_adapter']['scaled_order_bounds']))
    indices = tuple((i, n-i) for n in range(3) for i in range(n+1))
    mag = lambda item: max(abs(F(v)) for v in item['exact'])
    matrices = {name: {ij: [[mag(entry[f'{ij[0]},{ij[1]}']) for entry in row] for row in co['matrix_mixed_derivative_intervals'][name]] for ij in indices} for name in ('A', 'B', 'D')}
    ps, p0 = (co['source_profile_mixed_derivative_intervals'][name] for name in ('P_s', 'P_0'))
    residual, defect = {}, {}
    for i, j in indices:
        value, source = [F(0)]*3, F(0)
        for a in range(i+1):
            for b in range(j+1):
                weight, n, key = comb(i, a)*comb(j, b), i+j-a-b, f'{a},{b}'
                delta = max(abs(F(v)-(c if a+b == 0 else 0)) for v in ps[key]['exact'])
                source += weight*(delta*Q[n+1]+mag(ps[key])*ew[n+1]+mag(p0[key])*actual[n])
                for channel in range(3):
                    value[channel] += weight*((matrices['A'][a, b][channel][0]+matrices['B'][a, b][channel][0])*comparison[n+1]+matrices['D'][a, b][channel][0]*comparison[n])
        value[0] += source
        residual[i, j], defect[i, j] = value, source
    A, B, D = (matrices[k] for k in ('A', 'B', 'D'))
    linear = lambda key: [[A[key][i][j]+B[key][i][j]+theta*D[key][i][j] for j in range(3)] for i in range(3)]
    M = [times(x, row) for row in linear((0, 0))]
    need(max(map(sum, M)) < 1, 'exterior resolvent is a contraction')
    inv = adjugate([[F(i == j)-M[i][j] for j in range(3)] for i in range(3)])
    need(all(v >= 0 for row in inv for v in row), 'positive independent exterior resolvent')
    r0, rx, rr, rxx = (residual[k] for k in ((0, 0), (1, 0), (0, 1), (2, 0)))
    u = mv(inv, times(x, r0))
    v = mv(inv, plus(times(2, r0), times(x, plus(mv(linear((1, 0)), u), rx))))
    acceleration = plus(mv(A[0, 0], r0), rx)
    free = [max(r0[i]/theta+2*acceleration[i], rr[i]+acceleration[i]) for i in range(3)]
    w = mv(inv, plus(free, times(x, plus(times(2, mv(linear((1, 0)), v)), mv(linear((2, 0)), u), rxx))))
    result = {'0,0': times(theta, u), '1,0': u, '0,1': u, '2,0': [min(v[i], theta*w[i]) for i in range(3)], '1,1': v, '3,0': w, '2,1': w}
    result['0,2'] = plus(result['2,0'], mv(A[0, 0], u), mv(B[0, 0], u), mv(D[0, 0], result['0,0']), r0)
    result['1,2'] = plus(w, mv(A[1, 0], u), mv(A[0, 0], result['2,0']), mv(B[1, 0], u), mv(B[0, 0], v), mv(D[1, 0], result['0,0']), mv(D[0, 0], u), rx)
    result['0,3'] = plus(w, mv(A[0, 1], u), mv(A[0, 0], v), mv(B[0, 1], u), mv(B[0, 0], result['0,2']), mv(D[0, 1], result['0,0']), mv(D[0, 0], u), rr)
    return ew, result, {f'{i},{j}': v for (i, j), v in defect.items()}, {f'{i},{j}': v for (i, j), v in residual.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse overwrite')
    pins = json.loads((NUMERIC/'INPUT_PINS.json').read_text())
    data = {}
    for name, pin in pins.items():
        raw = (NUMERIC/pin['path']).read_bytes()
        need(len(raw) == pin['bytes'] and hashlib.sha256(raw).hexdigest() == pin['sha256'], 'independent input '+name)
        data[name] = json.loads(raw)
    spec = importlib.util.spec_from_file_location('reviewed_exterior', NUMERIC/'exterior_error_caps.py')
    producer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(producer)
    report = producer.run()
    ew, eu, defect, residual = derive(data, [F(0)]*4)
    need(ew == list(map(F, report['master_actual_minus_comparison_by_total_order'])), 'all four independent exterior master errors')
    need(eu == {k: list(map(F, row)) for k, row in report['gauge_actual_minus_comparison_mixed'].items()}, 'all thirty independent exterior gauge errors')
    need(defect == {k: F(v) for k, v in report['source_defect_mixed'].items()}, 'all six independent exterior source defects')
    need(residual == {k: list(map(F, row)) for k, row in report['complete_gauge_residual_mixed'].items()}, 'all eighteen independent exterior gauge residuals')
    inherited = json.loads((NUMERIC/'inputs/APPROXIMATION_ERRORS_NORMAL.json').read_text())
    m, b = data['master']['exact']['interior'], data['master']['exact']['boundary']
    full = derive(data, list(map(F, (b['y'], b['yprime'], m['ysecond'], m['ythird']))))
    need(full[0] == [F(v['exact_upper']) for v in inherited['master_error_budget_by_scaled_total_order']['actual_minus_comparison']], 'independent full-domain master baseline')
    need(full[1] == {k: [F(v['exact_upper']) for v in row] for k, row in inherited['actual_minus_gauge_comparison_mixed_scaled_derivative_caps'].items()}, 'independent full-domain gauge baseline')
    need(all(ew[k] <= full[0][k] for k in range(4)) and all(eu[k][j] <= full[1][k][j] for k in eu for j in range(3)), 'monotone exterior restriction')
    out = {'status': 'PASS_INDEPENDENT_EXTERIOR_C3_ERRORS', 'exact_rational_comparisons': 92,
           'independent_matrix_method': 'cyclic adjugate inverse',
           'all_point_domain_argument': 'For r>=x every backward characteristic point (s,q) has q>=r-x+s>=s. It meets the physical shell only at the flat corner. The exact causal dynamic lift and its jets through3 vanish. Other inherited global positive majorants remain upper bounds; no C0 cap is differentiated.',
           'master_exterior_errors': list(map(str, ew)),
           'review_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'producer_source_sha256': hashlib.sha256((NUMERIC/'exterior_error_caps.py').read_bytes()).hexdigest(),
           'input_sha256': {k: v['sha256'] for k, v in pins.items()}}
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
    print(out['status'])


if __name__ == '__main__':
    main()
