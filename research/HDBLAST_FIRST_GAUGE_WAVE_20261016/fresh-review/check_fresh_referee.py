#!/usr/bin/env python3
"""Bounded independent rational checks supporting the fresh scientific review.

No producer module is imported. The PDE/regularity proof is assessed in the
review note; this program authenticates its files and checks selected exact
quantitative implications by separate coarse estimates.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path


def need(condition, message):
    if not condition:
        raise ValueError(message)


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def mul(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def scale(c, a):
    return mul((F(c), F(c)), a)


def sq(a):
    lower = F(0) if a[0] <= 0 <= a[1] else min(a[0]**2, a[1]**2)
    return lower, max(a[0]**2, a[1]**2)


def recip(a):
    need(a[0] > 0 or a[1] < 0, 'nonzero interval denominator')
    return 1/a[1], 1/a[0]


def read(root, path):
    return json.loads((root/path).read_text())


def independent_checks(root):
    bg = read(root, 'coefficients/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json')
    shell = read(root, 'coefficients/inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json')
    coef = read(root, 'coefficients/GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json')
    master = read(root, 'coefficients/inputs/HIGHER_JETS_NORMAL.json')
    sol = read(root, 'solver/ACTUAL_GAUGE_SOLUTION_NORMAL.json')
    ell = F(1, 10**6)
    need(bg['domain']['Z'] == str(8*ell), 'actual full collar domain')
    need(master['scope']['conditional_on_retained_smooth_actual_solution'] is True,
         'retained actual-master premise is explicit')
    need(master['potential_report_sha256'] == coef['input_sha256']['ACTUAL_POTENTIAL_JET_BOUNDS.json'],
         'same potential in actual master and gauge coefficients')

    boxes = {k: tuple(map(F, bg['background_intervals'][k]['exact']))
             for k in ('rho', 'eta', 'p', 'a')}
    r, eta, p, a = [boxes[k] for k in ('rho', 'eta', 'p', 'a')]
    deriv = [F(0), F(28, 9), F(50, 9), F(-2, 3), F(-20, 9), F(-4, 9)]
    up = (F(0), F(0))
    for c in reversed(deriv):
        up = add(mul(up, eta), (c, c))
    rhs = {'rho': scale(-1, mul(sq(r), a)), 'eta': scale(-1, mul(r, p)),
           'p': scale(-1, mul(r, add(up, scale(-4, mul(a, p))))),
           'a': add(recip(r), scale(F(1, 3), mul(r, sq(p))))}
    margins = {}
    for name in boxes:
        need(rhs[name] == tuple(map(F, bg['background_intervals'][name+'_z']['exact'])),
             'independent ODE image: '+name)
        origin = tuple(map(F, shell['boxes'][name+'_b']['exact']))
        tube = add(origin, mul((F(0), 8*ell), rhs[name]))
        margin = tube[0]-boxes[name][0], boxes[name][1]-tube[1]
        need(min(margin) > 0, 'strict 8ell first-exit margins: '+name)
        need(list(map(str, margin)) == coef['first_exit'][name]['strict_face_margins'],
             'coefficient report first-exit agreement: '+name)
        margins[name] = list(map(str, margin))

    # Coarse matrices are derived directly from the retained outer box.
    # They are intentionally independent of the producer's tighter intervals.
    R, amin, pmax, H = F(131), F(14, 125), F(7, 50000), 3*ell
    ay = F(1, 129**2)+pmax*pmax/3
    A0 = ell*(3*H+2*amin*R*R)
    B0 = 5*ell*amin*R
    Drows = [ell**2*(8+R*R*(4*ay+8*amin**2)+6*amin*R*R*H+16*amin*R*R),
             ell**2*(3+3*H*H+8+16*H), ell**2*(2*amin+2*H+4)]
    D0 = max(Drows)
    q3 = 3*(A0+B0+F(3, 2)*D0)
    need(q3 < F(3, 250), 'independent coarse C1 contraction below 0.012')

    # Check the exact tighter matrix certificate by multiplication, not by
    # using the producer's matrix-inversion routine.
    M = [[F(x) for x in row] for row in coef['weighted_C1_matrix']]
    inverse = [[F(x) for x in row] for row in coef['weighted_C1_resolvent']]
    matrices = {k: [[F(x) for x in row] for row in coef['matrix_elementwise_absolute_caps'][k]]
                for k in ('A', 'B', 'D')}
    for i in range(3):
        for j in range(3):
            expected = 3*(matrices['A'][i][j]+matrices['B'][i][j]+F(3, 2)*matrices['D'][i][j])
            need(M[i][j] == expected and M[i][j] >= 0, 'exact C1 matrix entry')
            product = sum((F(i == k)-M[i][k])*inverse[k][j] for k in range(3))
            need(product == F(i == j), 'exact positive resolvent identity')
            need(inverse[i][j] >= 0, 'positive resolvent')
    exact_q = max(map(sum, M))
    need(exact_q == F(coef['weighted_C1_contraction_upper']) and exact_q < F('0.011412'),
         'reported contraction factor')

    # Recheck every actual point sign using a coarser scalar contraction and
    # independently bounded first-column couplings, not the producer's
    # componentwise inverse or its point intervals.
    psbox = tuple(map(F, coef['source_profile_mixed_derivative_intervals']['P_s']['0,0']['exact']))
    need(psbox[0] > 0, 'positive actual derivative source multiplier')
    ps = max(map(abs, psbox))
    pzero = max(abs(F(v)) for v in coef['source_profile_mixed_derivative_intervals']['P_0']['0,0']['exact'])
    cstar = 6*ell/130
    deviation = max(abs(v-cstar) for v in psbox)
    W0 = F(master['exact']['interior']['U'])
    yp = F(master['exact']['boundary']['yprime'])
    potential = F(master['exact']['background']['W'])
    need(potential == 102000, 'actual master potential bound')
    first_A = [3*ell*H, 2*ell*amin, F(0)]
    first_B = [3*ell*amin*R, F(0), F(0)]
    first_D = [ell**2*(8+R*R*(4*ay+8*amin**2)), F(0), 2*ell**2*amin]

    def bump_value(s):
        v = abs(s)
        if v <= 1 or v >= 2:
            return F(0)
        need(v == F(3, 2), 'point evaluation uses exact bump peak only')
        return F(1)

    point_results = []
    for point in sol['signed_actual_points']:
        x, depth = F(point['scaled_x']), F(point['scaled_r'])
        need(0 < x <= 3 and depth >= 0 and x+depth <= 8, 'evaluated point in causal triangle')
        plus, minus = depth+x, depth-x
        sgn = lambda v: F((v > 0)-(v < 0))
        comparator = cstar/4*(plus*bump_value(plus)+minus*bump_value(minus)
                      -depth*(sgn(plus)*bump_value(plus)+sgn(minus)*bump_value(minus)))
        need(comparator == F(point['comparator_radial_amplitude']), 'exact signed comparator value')
        boundary = F(0) if x <= 1 else yp
        master_derivative_error = ell+boundary+x*ell**2*potential*W0
        source_error = 16*deviation+ps*master_derivative_error+pzero*W0
        pvalue, pfirst = cstar*x/2, cstar*(F(1, 2)+8*x)
        residuals = [(source_error if i == 0 else F(0))
                     +(first_A[i]+first_B[i])*pfirst+first_D[i]*pvalue for i in range(3)]
        q = x*(A0+B0+x*D0/2)
        error = x*x*max(residuals)/(2*(1-q))
        lower, upper = comparator-error, comparator+error
        sign = 'positive' if lower > 0 else 'negative' if upper < 0 else 'unresolved'
        need(sign == point['actual_sign'] and sign != 'unresolved', 'independent actual gauge-amplitude sign')
        point_results.append({'x': str(x), 'r': str(depth), 'actual_sign': sign,
                              'independent_coarse_actual_interval': [str(lower), str(upper)],
                              'coarse_scalar_contraction': str(q)})
    need(len(point_results) == 5, 'all five reported signs independently covered')
    return {'status': 'PASS_FRESH_REFEREE_EXACT_QUANTITATIVE_CHECKS',
            'whole_8ell_first_exit_margins': margins,
            'independent_coarse_contraction': str(q3),
            'verified_producer_contraction': str(exact_q),
            'independent_actual_point_signs': point_results,
            'scope': 'Selected exact consequences of the analytically reviewed conditional gauge-wave theorem; no nonlinear or second-order evolution claim.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse to overwrite output')
    pins = json.loads(Path(__file__).with_name('FRESH_REVIEW_PINS.json').read_text())
    for path, entry in pins['reviewed_files'].items():
        need(not path.startswith('fresh-review/') and path not in ('MANIFEST.json', 'REPLAY_ROSTER.json', 'reproduce.py'),
             'nonrecursive review pin set')
        raw = (args.root/path).read_bytes()
        need(len(raw) == entry['bytes'] and sha256(raw).hexdigest() == entry['sha256'],
             'reviewed scientific file changed: '+path)
    result = independent_checks(args.root)
    result['reviewed_files'] = pins['reviewed_files']
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')


if __name__ == '__main__':
    main()
