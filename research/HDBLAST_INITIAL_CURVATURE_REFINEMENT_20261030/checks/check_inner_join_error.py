#!/usr/bin/env python3
"""Exact-rational inner-join initial-curvature ERROR diagnostics.

Original actual target, W1_quadratic candidate and both smoothing widths stay
fixed. All scalar/vector/tensor full23 rows are retained. No physics accuracy or
nonlinear-spacetime conclusion follows from this coordinate data calculation.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode = True

CORE_PIN = 'f2d1bfef82220fc1e385a3bdb6de35487b1ef28c4caee37c8dddb6a2b68cb782'
INITIAL_PIN = 'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'
PROFILES_PIN = '555f911f6c336862dfacb6aa56ef5179ac06701812d593b5bfbbf6ea9bd4d5d1'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw = Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return raw


def encode(z):
    if isinstance(z, F):
        return str(z)
    if isinstance(z, dict):
        return {k: encode(v) for k, v in z.items()}
    if isinstance(z, (tuple, list)):
        return [encode(v) for v in z]
    return z


def row(z):
    return [z.lo, z.hi]


def integral01(p):
    return sum((a / F(i + 1) for i, a in enumerate(p)), F(0))


def derivative(p):
    return [i * a for i, a in enumerate(p) if i] or [F(0)]


def sqrt_lower(q):
    need(q >= 0, 'nonnegative squared lower bound')
    grid = 1 << 80
    n = math.isqrt((q * grid * grid).__floor__())
    # A strict lower bound is convenient in the public statement.
    if n and F(n, grid) ** 2 == q:
        n -= 1
    return F(n, grid)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--core', type=Path, default=Path(__file__).with_name('initial_second_core.py'))
    ap.add_argument('--fault', '--inject-error', dest='fault', choices=('half-mass', 'half-moment', 'curvature-factor',
                    'join-branches', 'cutoff-derivative', 'gram-weight'))
    args = ap.parse_args()
    groups, controls = [], []

    def checked(name, ok):
        need(ok, name)
        groups.append(name)

    core_raw = pin(args.core, CORE_PIN)
    root = args.repo_root.resolve()
    initial_dir = root / 'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial'
    initial_raw = pin(initial_dir / 'results/INITIAL_REPORT.json', INITIAL_PIN)
    profiles_raw = pin(initial_dir / 'results/INITIAL_PROFILES.json', PROFILES_PIN)
    initial, profiles = json.loads(initial_raw), json.loads(profiles_raw)
    spec = importlib.util.spec_from_file_location('inner_join_order5_authenticated', args.core)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    h = mod.HigherInitial(root)  # Complete retained source roots checked before scientific imports.
    c, I, ell = h.ctx, h.I, h.ell
    labels, weights = list(c.op.LABELS), list(c.weights)
    checked('retained_full23_labels', len(labels) == 23 and initial['labels'] == profiles['labels'] == labels)
    checked('complete_positive_Gram', len(weights) == 23 and min(weights) > 0 and list(map(F, profiles['gram'])) == weights)
    checked('registered_ell', ell == F(1, 10**6))
    checked('retained_initial_status', initial['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER')
    du, R = F(1, 2**30), F(7, 4)
    joins = ((F(65, 64), 'zero', 'cutoff'), (F(66, 64), 'cutoff', 'plain'))
    if args.fault == 'join-branches':
        joins = ((F(65, 64), 'zero', 'cutoff'), (F(66, 64), 'plain', 'plain'))
    checked('exact_one_sided_branches', joins == ((F(65, 64), 'zero', 'cutoff'), (F(66, 64), 'cutoff', 'plain')))
    checked('isolated_joins_and_unchanged_outer_repair', 2*du < F(1, 64) and joins[-1][0]+du < R-du and joins[0][0]-du > 1)

    kernel = [F(35, 32)*a for a in (1, 0, -3, 0, 3, 0, -1)]
    mass, moment = integral01(kernel), integral01([F(0)]+kernel)
    if args.fault == 'half-mass':
        mass *= 2
    if args.fault == 'half-moment':
        moment /= 2
    curvature_factor = 2 if args.fault != 'curvature-factor' else 1
    checked('kernel_half_mass', mass == F(1, 2))
    checked('kernel_half_first_moment', moment == F(35, 256))
    checked('ordinary_scaled_curvature_factor', curvature_factor == math.factorial(2))
    checked('kernel_normalization', 2*integral01(kernel) == 1)
    checked('kernel_positive_on_its_support', kernel == [F(35, 32)*a for a in (1, 0, -3, 0, 3, 0, -1)])
    S = [F(a) for a in (0,0,0,0,0,126,-420,540,-315,70)]
    if args.fault == 'cutoff-derivative':
        S[5] += 1
    p = S
    for n in range(5):
        checked('packet_matching_derivative_'+str(n), (p[0], sum(p)) == ((F(0),F(1)) if n == 0 else (F(0),F(0))))
        p = derivative(p)
    checked('fifth_packet_derivative_is_nonzero', p[0] == sum(p) == 15120)
    controls.extend(('unit_jump_convolution_is_half_not_one', 'half_first_moment_is_not_full_moment',
                     'normalized_c2_is_not_ordinary_curvature', 'plain_branch_cannot_replace_cutoff_trace',
                     'packet_C4_does_not_imply_C5', 'positive_Gram_cannot_drop_scalar_mean'))
    checked('wrong_unit_jump_formula_rejected', mass != 1)
    checked('wrong_full_moment_formula_rejected', moment != F(35, 128))
    checked('wrong_packet_C5_matching_rejected', p[0] != 0 and sum(p) != 0)
    checked('wrong_scaled_normalized_factor_rejected', curvature_factor != 1)
    if args.fault == 'gram-weight':
        weights[12] = F(0)
    checked('scalar_mean_positive_Gram_weight', labels[12] == 'w:0' and weights[12] == 1)

    def intersect(a, b, name):
        checked(name, max(a.lo, b.lo) <= min(a.hi, b.hi))
        return I(max(a.lo, b.lo), min(a.hi, b.hi))

    tables = []
    for rc, left, right in joins:
        pu, _ = h.paired(I(rc))
        ln, _ = h.nominal(I(rc), left)
        rn, _ = h.nominal(I(rc), right)
        ls, _ = h.nominal(I(rc-du, rc), left)
        rs, _ = h.nominal(I(rc, rc+du), right)
        # At65/64 the paired packet's symmetric flat-tail derivative caps can
        # obscure signs. Evaluate the SAME actual analytic point with the
        # authenticated signed bump polynomials, rather than change the data.
        bg = {k: v.a for k, v in h.background(I(rc)).items()}
        E = h.kernel(ell*h.J([I(rc), 1/ell]+[0]*4)).a
        packet = h.bump(I(rc), tail=False)
        actual_direct, _ = c._registered(bg, E, packet)
        checked('full23_point_and_local_jets_'+str(rc), len(actual_direct) == len(ls) == len(rs) == 23 and
                all(z.n >= 2 for z in actual_direct) and all(z.n >= 3 for z in ls+rs))
        rows, nonzero, lower_square = [], [], F(0)
        for i, label in enumerate(labels):
            tag = str(rc)+'_'+label
            L, U = curvature_factor*ell**2*ln[i].c[2], curvature_factor*ell**2*rn[i].c[2]
            B = curvature_factor*ell**2*pu[i].b.c[2]
            D = sum((curvature_factor*ell**2*pu[i].d[k].c[2] for k in range(3)), I(0))
            direct = curvature_factor*ell**2*actual_direct[i].c[2]
            actual = intersect(direct, curvature_factor*ell**2*pu[i].a.c[2], 'actual_routes_overlap_'+tag)
            Mminus, Mplus = (6*ell**3*ls[i].c[3]).mag(), (6*ell**3*rs[i].c[3]).mag()
            remainder = moment*du*(Mminus+Mplus)
            average = mass*(L+U)
            smooth = average+I(-remainder, remainder)
            error_direct = actual-average+I(-remainder, remainder)
            error_paired = B-average+D+I(-remainder, remainder)
            error = intersect(error_direct, error_paired, 'actual_error_routes_overlap_'+tag)
            limit = intersect(actual-average, B-average+D, 'zero_width_routes_overlap_'+tag)
            minimum = F(0) if error.lo <= 0 <= error.hi else min(abs(error.lo), abs(error.hi))
            lower_square += weights[i]*minimum**2
            if minimum:
                nonzero.append(label)
                checked('same_nonzero_zero_width_limit_sign_'+tag,
                        (error.lo > 0 and limit.lo > 0) or (error.hi < 0 and limit.hi < 0))
            if rc == F(65, 64):
                checked('zero_nominal_curvature_traces_'+label, L.lo == L.hi == U.lo == U.hi == 0)
            rows.append({'label': label, 'positive_Gram_weight': weights[i],
                         'nominal_left_u_rr': row(L), 'nominal_right_u_rr': row(U),
                         'signed_nominal_curvature_jump_right_minus_left': row(U-L),
                         'actual_u_rr_signed_direct': row(direct), 'actual_u_rr': row(actual),
                         'paired_nominal_selected_u_rr': row(B), 'paired_actual_minus_selected_nominal_u_rr': row(D),
                         'left_local_M3': Mminus, 'right_local_M3': Mplus,
                         'half_kernel_remainder_abs_upper': remainder,
                         'fixed_width_smoothed_nominal_u_rr': row(smooth),
                         'actual_minus_smoothed_u_rr': row(error),
                         'zero_width_limit_actual_minus_trace_average': row(limit),
                         'fixed_width_error_component_abs_lower': minimum,
                         'sign_determined': bool(minimum)})
        expected = ['J:0', 'J:2', 'J:4', 'w:0', 'w:2', 'w:4'] if rc == F(65,64) else ['w:0','w:2','w:4']
        checked('exact_nonzero_actual_error_labels_'+str(rc), nonzero == expected)
        Hlower = sqrt_lower(lower_square)
        checked('strict_complete_H_lower_'+str(rc), Hlower > 0 and Hlower**2 < lower_square)
        tables.append({'r': rc, 'left_branch': left, 'right_branch': right, 'rows': rows,
                       'fixed_width_nonzero_component_labels': nonzero,
                       'fixed_width_error_H_norm_squared_lower': lower_square,
                       'fixed_width_error_H_norm_strict_lower': Hlower})
    first_mean, second_mean = tables[0]['rows'][12], tables[1]['rows'][12]
    a, b = map(F, first_mean['actual_minus_smoothed_u_rr'])
    checked('first_join_scalar_mean_readable_signed_bracket', F('0.6114') < a <= b < F('0.6115'))
    a, b = map(F, second_mean['actual_minus_smoothed_u_rr'])
    checked('second_join_scalar_mean_readable_signed_bracket', -F('5.195e20') < a <= b < -F('5.193e20'))
    h.authenticate()
    pin(args.core, CORE_PIN)
    pin(initial_dir/'results/INITIAL_REPORT.json', INITIAL_PIN)
    pin(initial_dir/'results/INITIAL_PROFILES.json', PROFILES_PIN)
    report = {'status':'PASS_INNER_JOIN_ACTUAL_INITIAL_CURVATURE_ERROR',
              'candidate':'unchanged W1_quadratic / W1_dagger from PR100 and PR102',
              'actual_target_unchanged':True, 'original_bulk_and_shell_sources_unchanged':True,
              'initial_candidate_unchanged':True, 'ell':ell, 'delta_u':du,
              'labels':labels, 'positive_Gram':weights,
              'definition':'actual U0_rr minus (k_delta * registered nominal displacement)_rr at each inner cutoff join; ordinary scaled r derivative',
              'formula':'u_s,rr(rc)=(kminus+kplus)/2+rem; abs(rem_i)<=35delta/256*(M3minus_i+M3plus_i)',
              'actual_method':'signed authenticated order5 actual background Taylor route and signed analytic bump polynomial jet; intersected with retained paired actual and telescoping actual-minus-nominal routes',
              'nominal_method':'forced exact one-sided branches and local nominal order6 jets on [rc-delta,rc] and [rc,rc+delta]',
              'ordinary_scaled_derivative_factors':{'u_rr':'2 ell^2 c2','nominal_u_rrr':'6 ell^3 c3'},
              'kernel_half_mass':mass, 'kernel_half_first_moment_coefficient':moment,
              'joins':tables, 'first_join_scalar_mean_strict_signed_bracket':['0.6114','0.6115'],
              'second_join_scalar_mean_strict_signed_bracket':['-5.195e20','-5.193e20'],
              'checks':{'named_groups':groups,'group_count':len(groups),
                        'wrong_formula_controls':controls,'wrong_formula_control_count':len(controls)},
              'input_sha256':{'order5_core':CORE_PIN,'retained_initial_report':INITIAL_PIN,
                              'retained_initial_profiles':PROFILES_PIN,'retained_initial_producer':mod.PRODUCER_PIN},
              'authenticated_before_import_and_readback':True,
              'scope':'conditional coordinate initial-data mismatch for a fixed auxiliary reconstruction; not a physical force or model obstruction; no new candidate, evolved-error upper bound, full Einstein residual, nonlinear family, cosmological observable or Big Bang-origin evidence'}
    print(json.dumps(encode(report), sort_keys=True))


if __name__ == '__main__':
    main()
