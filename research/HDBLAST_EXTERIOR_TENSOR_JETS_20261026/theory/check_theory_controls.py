#!/usr/bin/env python3
"""Independent exact-rational identities for the exterior TT two-jet ledger.

No project module is imported.  Actual-map checks are structural read-only
audits, not substitutes for the evaluator's complete-tree authentication.
"""
import argparse
import ast
from fractions import Fraction as Q
import hashlib
import json
from math import comb, factorial
from pathlib import Path

PROTECTED_CAPSULES = (
    'HDBLAST_COUPLED_MEAN_PILOT_20261025',
    'HDBLAST_INPUT_PRECISION_20261023',
    'HDBLAST_COMMON_GAUGE_FORCING_20261016',
    'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016',
    'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
    'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023',
    'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023',
    'HDBLAST_FULL_ANGULAR_RESPONSE_20261026',
    'HDBLAST_SHELL_FEEDBACK_PRECISION_20261025',
)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def guard_output(repo, output):
    """Replay writes only fresh outputs outside all old and current capsules."""
    need(not output.exists() and not output.is_symlink(), 'FRESH_OUTPUT_REQUIRED')
    target = output.resolve()
    protected = [(repo/'research'/name).resolve() for name in PROTECTED_CAPSULES]
    # The script is packaged at <capsule>/theory/check_theory_controls.py.
    # Protect its actual source tree even when --repo-root points elsewhere.
    protected.append(Path(__file__).resolve().parent.parent)
    need(all(target != root and root not in target.parents for root in protected),
         'OUTPUT_INSIDE_PROTECTED_SOURCE_TREE')


def clean(p):
    return {k: Q(v) for k, v in p.items() if v}


def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, Q(0)) + v
    return clean(out)


def scale(p, c):
    return clean({k: c*v for k, v in p.items()})


def mul(p, q):
    out = {}
    for (i, j), a in p.items():
        for (k, l), b in q.items():
            key = (i+k, j+l)
            out[key] = out.get(key, Q(0)) + a*b
    return clean(out)


def diff(p, axis):
    out = {}
    for (i, j), a in p.items():
        exponent = (i, j)[axis]
        if exponent:
            key = (i-1, j) if axis == 0 else (i, j-1)
            out[key] = exponent*a
    return clean(out)


def evaluate(p, x, r):
    return sum((a*x**i*r**j for (i, j), a in p.items()), Q(0))


def at_x_zero(p):
    return clean({(i, j): a for (i, j), a in p.items() if i == 0})


def K(p):
    """Exact double cone integral from the monomial beta-integral formula."""
    out = {}
    for (i, j), a in p.items():
        for k in range(0, j+1, 2):
            key = (i+k+2, j-k)
            c = a*comb(j, k)*Q(factorial(i)*factorial(k), factorial(i+k+2))
            out[key] = out.get(key, Q(0)) + c
    return clean(out)


def characteristic_integral(p, sign):
    """Exact integral_0^x p(s,r+sign*(x-s)) ds, without retarded 1/2."""
    out = {}
    for (i, j), a in p.items():
        for k in range(j+1):
            key = (i+k+1, j-k)
            c = a*comb(j, k)*sign**k*Q(factorial(i)*factorial(k), factorial(i+k+1))
            out[key] = out.get(key, Q(0)) + c
    return clean(out)


def operator(p, A, B, D):
    return add(mul(A, diff(p, 0)), mul(B, diff(p, 1)), mul(D, p))


def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def theory_controls():
    tests = []
    # Distinct degrees and both parities prevent accidental cancellation of
    # the characteristic signs, endpoint term, or retarded normalization.
    forcing = {(0, 0): Q(7, 13), (1, 2): Q(-3, 5), (2, 3): Q(11, 17), (3, 1): Q(2, 7)}
    z = K(forcing)
    plus = characteristic_integral(forcing, 1)
    minus = characteristic_integral(forcing, -1)
    zx = scale(add(plus, minus), Q(1, 2))
    zr = scale(add(plus, scale(minus, -1)), Q(1, 2))
    fr = diff(forcing, 1)
    zxr = scale(add(characteristic_integral(fr, 1), characteristic_integral(fr, -1)), Q(1, 2))
    zrr = scale(add(characteristic_integral(fr, 1), scale(characteristic_integral(fr, -1), -1)), Q(1, 2))
    for label, left, right in [
        ('K_x_boundary_sum', diff(z, 0), zx),
        ('K_r_boundary_difference', diff(z, 1), zr),
        ('K_xr_boundary_Fr_sum', diff(diff(z, 0), 1), zxr),
        ('K_rr_boundary_Fr_difference', diff(diff(z, 1), 1), zrr),
        ('K_xx_endpoint_term', diff(diff(z, 0), 0), add(forcing, zrr)),
        ('K_commutes_spatial_derivative', diff(z, 1), K(fr)),
    ]:
        need(left == right, label)
        tests.append(label)
    need(not at_x_zero(z) and not at_x_zero(diff(z, 0)), 'K_zero_initial_data')
    tests.append('K_zero_initial_data')

    u = {(2, 3): Q(1), (3, 2): Q(3, 5), (4, 1): Q(-2, 7), (5, 0): Q(5, 11)}
    A = {(0, 0): Q(1, 13), (1, 1): Q(1, 17)}
    B = {(0, 0): Q(-1, 19), (0, 2): Q(1, 23)}
    D = {(0, 0): Q(1, 29), (1, 0): Q(1, 31), (0, 1): Q(-1, 37)}
    wave_u = add(diff(diff(u, 0), 0), scale(diff(diff(u, 1), 1), -1))
    F = add(wave_u, scale(operator(u, A, B, D), -1))
    v = diff(u, 1)
    Ar, Br, Dr = [diff(coefficient, 1) for coefficient in (A, B, D)]
    Q_u = operator(u, Ar, Br, Dr)
    v_force = add(diff(F, 1), operator(v, A, B, D), Q_u)
    need(u == K(add(F, operator(u, A, B, D))), 'manufactured_u_integral_equation')
    need(v == K(v_force), 'manufactured_v_differentiated_integral_equation')
    need(not at_x_zero(v) and not at_x_zero(diff(v, 0)), 'v_zero_initial_data')
    acceleration = at_x_zero(diff(diff(u, 0), 0))
    need(bool(acceleration) and acceleration == at_x_zero(F), 'nonzero_initial_acceleration')
    tests.extend(['manufactured_u_integral_equation', 'manufactured_v_differentiated_integral_equation',
                  'v_zero_initial_data', 'nonzero_initial_acceleration'])

    q, s = Q(3, 17), Q(19, 11)  # s>1: no q+s<1 assumption is allowed.
    block = [[q, Q(0)], [s, q]]
    power = [[Q(1), Q(0)], [Q(0), Q(1)]]
    for n in range(1, 9):
        power = mm(power, block)
        need(power == [[q**n, Q(0)], [n*s*q**(n-1), q**n]], 'triangular_power')
    inverse = [[1/(1-q), Q(0)], [s/(1-q)**2, 1/(1-q)]]
    need(mm([[1-q, Q(0)], [-s, 1-q]], inverse) == [[1, 0], [0, 1]], 'triangular_resolvent')
    T, f, g = Q(1, 4), Q(7, 3), Q(11, 5)
    Mu = T*f/(1-q)
    Mv = (T*g+s*Mu)/(1-q)
    need(q*Mv+s*Mu == q*T*g/(1-q)+s*T*f/(1-q)**2, 'complete_second_jet_tail')
    tests.extend(['triangular_power', 'triangular_resolvent_with_s_above_one', 'complete_second_jet_tail'])

    # rho=2/(1+kappa*r), so -1/(ell*rho) is polynomial and its vector-field
    # action can be differentiated independently from the conversion formula.
    ell, kappa = Q(1, 7), Q(2, 5)
    invrho = {(0, 0): Q(1, 2), (0, 1): kappa/2}
    a = kappa/(2*ell)
    yy_direct = mul(scale(invrho, -1/ell), diff(mul(scale(invrho, -1/ell), diff(u, 1)), 1))
    yy_formula = add(scale(mul(mul(invrho, invrho), diff(diff(u, 1), 1)), 1/ell**2),
                     scale(mul(invrho, diff(u, 1)), a/ell))
    need(yy_direct == yy_formula, 'natural_yy_connection_correction')
    ty_direct = scale(diff(mul(scale(invrho, -1/ell), diff(u, 1)), 0), 1/ell)
    ty_formula = scale(mul(invrho, diff(diff(u, 0), 1)), -1/ell**2)
    need(ty_direct == ty_formula, 'natural_ty_orientation_sign')
    tests.extend(['natural_yy_connection_correction', 'natural_ty_orientation_sign'])

    # A quadratic inverse rho makes a and a_y nonconstant, so the third-M
    # identities also test the a_y sign rather than only constant-a cases.
    ir = {(0, 0): Q(2, 3), (0, 1): Q(3, 11), (0, 2): Q(5, 17)}
    aa = scale(diff(ir, 1), 1/ell)
    aay = scale(mul(ir, diff(aa, 1)), -1/ell)
    M = scale(u, 2)
    My = scale(mul(ir, diff(M, 1)), -1/ell)
    Myy = scale(mul(ir, diff(My, 1)), -1/ell)
    Myt = scale(diff(My, 0), 1/ell)
    Mtt = scale(diff(diff(M, 0), 0), 1/ell**2)
    Myyr = add(scale(mul(mul(ir, ir), diff(diff(diff(u, 1), 1), 1)), 2/ell**2),
               scale(mul(mul(aa, ir), diff(diff(u, 1), 1)), 6/ell),
               scale(mul(add(mul(aa, aa), scale(aay, -1)), diff(u, 1)), 2))
    Mytr = add(scale(mul(ir, diff(diff(diff(u, 0), 1), 1)), -2/ell**2),
               scale(mul(aa, diff(diff(u, 0), 1)), -2/ell))
    Mttr = scale(diff(diff(diff(u, 0), 0), 1), 2/ell**2)
    for label, direct, formula in [('M_yy_r_C3_conversion', diff(Myy, 1), Myyr),
                                   ('M_yt_r_C3_conversion', diff(Myt, 1), Mytr),
                                   ('M_tt_r_C3_conversion', diff(Mtt, 1), Mttr)]:
        need(direct == formula, label)
        tests.append(label)

    # Compare unsimplified and compact D_r identities at exact synthetic
    # background jets satisfying a_y=-rho^-2-p^2/3.  This is separate from
    # the numerical producer's Taylor-jet differentiation implementation.
    for rhoj, aj, pj, pyj, Vj, Vpj, Hj in [
        (Q(7, 3), Q(2, 11), Q(-5, 13), Q(3, 17), Q(-7, 19), Q(2, 23), Q(5, 29)),
        (Q(11, 7), Q(-3, 5), Q(7, 17), Q(-5, 19), Q(13, 23), Q(-7, 29), Q(-3, 31)),
    ]:
        ayj = -1/rhoj**2-pj**2/3
        rhor, ar = -ell*rhoj**2*aj, -ell*rhoj*ayj
        pr, Vr = -ell*rhoj*pyj, -ell*rhoj*pj*Vpj
        ayr = 2*rhor/rhoj**3-Q(2, 3)*pj*pr
        bj, dj, hj = -(ayj+aj**2), 1/rhoj**2-aj**2, Hj/rhoj
        br, dr, hr = -ayr-2*aj*ar, -2*rhor/rhoj**3-2*aj*ar, -Hj*rhor/rhoj**2
        inner = -2*aj**2+2*hj**2-2*bj-8*dj+Q(4, 3)*Vj
        inner_r = -4*aj*ar+4*hj*hr-2*br-8*dr+Q(4, 3)*Vr
        direct = ell**2*(2*rhoj*rhor*inner+rhoj**2*inner_r)
        Sj = 8*aj**2-Q(2, 3)*pj**2+Q(4, 3)*Vj
        compact = -ell**3*rhoj**3*(2*aj*Sj+16*aj*ayj-Q(4, 3)*pj*pyj+Q(4, 3)*pj*Vpj)
        need(direct == compact, 'D_r_unsimplified_equals_compact')
        need(-3*ell*(rhor*aj+rhoj*ar) == 3*ell**2*rhoj**2*(aj**2+ayj), 'B_r_identity')
    tests.extend(['D_r_unsimplified_equals_compact', 'B_r_identity'])

    x0, r0, H, potential = Q(1, 4), Q(5, 4), Q(2, 13), Q(-3, 23)
    rho = 1/evaluate(invrho, x0, r0)
    value, ux, ur, uxx, urr = [evaluate(p, x0, r0) for p in
                              (u, diff(u, 0), diff(u, 1), diff(diff(u, 0), 0), diff(diff(u, 1), 1))]
    AA, BB = -3*ell*H, -3*ell*rho*a
    FF = uxx-urr-AA*ux-BB*ur-potential*value
    ut, uy = ux/ell, -ur/(ell*rho)
    utt = uxx/ell**2
    uyy = evaluate(yy_direct, x0, r0)
    need(utt/rho**2-uyy == -3*H*ut/rho**2+4*a*uy+(potential*value+FF)/(ell**2*rho**2),
         'natural_coordinate_PDE_factors')
    tests.append('natural_coordinate_PDE_factors')

    # Deliberate mutations must violate independent polynomial identities.
    mutations = {
        'dropped_Br_ur': K(add(v_force, scale(mul(Br, diff(u, 1)), -1))) != v,
        'radial_boundary_plus_instead_of_minus': zx != diff(z, 1),
        'omitted_retarded_half': add(plus, minus) != diff(z, 0),
        'dropped_natural_yy_correction': scale(mul(mul(invrho, invrho), diff(diff(u, 1), 1)), 1/ell**2) != yy_direct,
        'omitted_second_jet_coefficient_tail': q*Mv != q*Mv+s*Mu,
    }
    for name, rejected in mutations.items():
        need(rejected, 'MUTATION_ACCEPTED:'+name)
    return {'passed_identities': tests, 'mutations_rejected': mutations,
            'triangular_control_q': str(q), 'triangular_control_s': str(s),
            'initial_acceleration_polynomial': [{'powers': list(k), 'coefficient': str(a)} for k, a in sorted(acceleration.items())]}


def source_dependency_audit(repo):
    base = repo/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
    mapfile = base/'projection/GEOMETRIC_COEFFICIENT_MAPS.json'
    jetfile = base/'numeric/actual_first_jets.py'
    maps = json.loads(mapfile.read_text())
    need(maps['status'] == 'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION', 'MAP_STATUS')
    row = maps['maps']['metric_forcing']['ab'][5]
    need(row == maps['maps']['metric_forcing']['ab'][6] and len(row) == 78, 'TT_ROW_IDENTITY')
    symbols = sorted({name for term in row for name in term['powers']})
    second = sorted(name for name in symbols if name.endswith(('yy', 'yt', 'tt')))
    need(second == ['Mtt', 'Myt', 'Myy'], 'SECOND_DERIVATIVE_SYMBOL_SET')
    positions = [i for i, term in enumerate(row) if set(term['powers']) & set(second)]
    need(positions == [7, 22, 50], 'SECOND_DERIVATIVE_POSITIONS')
    expected = [
        {'coefficient': '-4/3', 'powers': {'A': 1, 'Myy': 1, 'ell': 2, 'rho': 2}},
        {'coefficient': '8/3', 'powers': {'B': 1, 'Myt': 1, 'ell': 2, 'rho': 1}},
        {'coefficient': '-4/3', 'powers': {'D': 1, 'Mtt': 1, 'ell': 2}},
    ]
    need([row[i] for i in positions] == expected, 'SECOND_DERIVATIVE_EXACT_TERMS')
    function = next(node for node in ast.parse(jetfile.read_text()).body
                    if isinstance(node, ast.FunctionDef) and node.name == 'frame_amplitudes')
    amps = next(node.value for node in function.body if isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == 'amps' for target in node.targets))
    M = next(value for key, value in zip(amps.keys, amps.values) if isinstance(key, ast.Constant) and key.value == 'M')
    need(ast.dump(M) == ast.dump(ast.parse('2*d', mode='eval').body), 'REGISTERED_M_EQ_2D')
    unpack = next(node for node in function.body if isinstance(node, ast.Assign)
                  and any(isinstance(target, ast.Tuple) and [name.id for name in target.elts] == ['f', 'b', 'd']
                          for target in node.targets))
    need(isinstance(unpack.value, ast.Name) and unpack.value.id == 'U', 'D_IS_GAUGE_AMPLITUDE')
    return {'status': 'PASS_READ_ONLY_STRUCTURAL_DEPENDENCY_AUDIT',
            'scope': 'No producer import; complete frozen-tree authentication remains a separate evaluator prerequisite.',
            'source_map': {'logical_path': mapfile.relative_to(repo).as_posix(), 'sha256': hashlib.sha256(mapfile.read_bytes()).hexdigest()},
            'first_jet_source': {'logical_path': jetfile.relative_to(repo).as_posix(), 'sha256': hashlib.sha256(jetfile.read_bytes()).hexdigest()},
            'total_terms': len(row), 'symbols': symbols, 'second_derivative_symbols': second,
            'second_derivative_term_indices_zero_based': positions, 'second_derivative_terms': expected,
            'registered_M': '2*d, with f,b,d=U (gauge amplitudes)',
            'highest_required_scaled_actual_master_or_gauge_order': 3,
            'C4_required_for_this_TT_source_derivative': False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    guard_output(args.repo_root.resolve(), args.output)
    report = {'schema': 'HDBLAST_TT_TWOJET_INDEPENDENT_THEORY_CONTROLS_V1',
              'status': 'PASS_INDEPENDENT_TT_TWOJET_THEORY_CONTROLS',
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'exact_rational_theory_controls': theory_controls(),
              'actual_TT_dependency_audit': source_dependency_audit(args.repo_root.resolve()),
              'limits': ['No actual forcing value or response computed.',
                         'Manufactured polynomial controls validate formulas, not actual input error certificates.',
                         'Ordinary registered coefficient derivatives are not gauge-invariant curvature.']}
    with args.output.open('x') as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({'status': report['status'], 'identities': len(report['exact_rational_theory_controls']['passed_identities']),
                      'mutations_rejected': len(report['exact_rational_theory_controls']['mutations_rejected'])}))


if __name__ == '__main__':
    main()
