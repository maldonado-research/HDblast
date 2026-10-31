#!/usr/bin/env python3
"""Portable actual combined TT-jet enclosure from authenticated old box inputs.

No scientific producer is imported and no source oracle is run.  The retained
whole-box F/G ranges remain inherited conditional inputs.  All operations in
this checker use integers/Fraction; displayed decimals are not proof data.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True

T, R, ELL = Q(5, 4), Q(1, 8), Q(1, 10**6)
L = T-R
GRID = 1 << 192
OLD_GRID = 1 << 112
CAPS = {
 'HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031': 'fb3046fa483cd781060e167324c8002a241898f82627172bcde8187836a3234e',
 'HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031': 'e05bf3bb3ca31fe19321568886da29e917b4f1b784239daf69e05c733b71a672',
 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016': '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2',
 'HDBLAST_FULL_ANGULAR_RESPONSE_20261026': 'a35f84119e84bf2ea8a712ded86ff79977d3ee91ac25f91ccbd202d379c859ac',
 'HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026': '204e647b9f68f2af169b052deabee7c5d7f4bd71195fae31cefbfc2d10265777',
 'HDBLAST_EXACT_PACKET_INITIAL_20261031': '8201452fa5bc745a9ccfd64f67f5f2a15c564a3364af619873825c3a7a6015aa',
}
# Literal formula switches are exercised by real corrupted checker subprocesses.
NEUMANN_COMBO_SIGN = -1
REFLECTED_COMBO_SIGN = 1
POSITIVE_WEIGHT_IS_EXP_MINUS = True
RADIAL_REFLECTED_SIGN = -1
FEEDBACK_GAIN = 1
SECOND_DERIVATIVE_NORMALIZATION = 1


def need(ok, name):
    if not ok:
        raise ValueError(name)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_bytes())


def encode(v):
    if isinstance(v, I):
        return [str(v.lo), str(v.hi)]
    if isinstance(v, Q):
        return str(v)
    if isinstance(v, dict):
        return {k: encode(x) for k, x in v.items()}
    if isinstance(v, (tuple, list)):
        return [encode(x) for x in v]
    return v


class I:
    """Exact rational interval. Quantization is explicit, never implicit."""
    def __init__(self, lo=0, hi=None):
        self.lo, self.hi = Q(lo), Q(lo if hi is None else hi)
        need(self.lo <= self.hi, 'INTERVAL_ORDER')

    @staticmethod
    def coerce(v):
        return v if isinstance(v, I) else I(v)

    def __add__(self, v):
        v = self.coerce(v)
        return I(self.lo+v.lo, self.hi+v.hi)
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, v):
        return self+-self.coerce(v)

    def __rsub__(self, v):
        return self.coerce(v)+-self

    def __mul__(self, v):
        v = self.coerce(v)
        p = [a*b for a in (self.lo, self.hi) for b in (v.lo, v.hi)]
        return I(min(p), max(p))
    __rmul__ = __mul__

    def inv(self):
        need(not self.lo <= 0 <= self.hi, 'NONZERO_DENOMINATOR')
        return I(1/self.hi, 1/self.lo)

    def __truediv__(self, v):
        return self*self.coerce(v).inv()

    def __rtruediv__(self, v):
        return self.coerce(v)*self.inv()

    def __pow__(self, n):
        need(n in (0, 1, 2), 'SUPPORTED_POWER')
        if n == 0:
            return I(1)
        if n == 1:
            return self
        return I(0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2),
                 max(self.lo**2, self.hi**2))

    def mag(self):
        return max(abs(self.lo), abs(self.hi))

    def lattice(self, grid=GRID):
        return I(Q((self.lo*grid).__floor__(), grid),
                 Q((self.hi*grid).__ceil__(), grid))


def box(v):
    need(isinstance(v, list) and len(v) == 2, 'TWO_ENDPOINT_BOX')
    return I(*v)


def same(v, encoded, label):
    need(encode(v) == encoded, label)


def exp_point(x):
    """Taylor at |x|/2^k <= 1/2 with geometric positive-tail majorant.

    Terms after degree 64 have ratio <= y/66.  The entire discarded tail
    is <= (y^65/65!)/(1-y/66).  Negative arguments use exact reciprocals.
    Each squaring is rounded outward on the 192-bit lattice.
    """
    x = Q(x)
    y, k = abs(x), 0
    need(y <= 16, 'EXP_ARGUMENT_REGISTERED_RANGE')
    while y > Q(1, 2):
        y /= 2
        k += 1
    term, partial = Q(1), Q(1)
    for j in range(1, 65):
        term *= y/j
        partial += term
    omitted = term*y/65
    remainder = omitted/(1-y/66)
    out = I(partial, partial+remainder).lattice()
    for _ in range(k):
        out = (out**2).lattice()
    if x < 0:
        out = out.inv().lattice()
    return out


def exp_box(v):
    return I(exp_point(v.lo).lo, exp_point(v.hi).hi)


def checked_member(root, row, prefix):
    rel = Path(row['path'])
    need(not rel.is_absolute() and '..' not in rel.parts, prefix+':SAFE_PATH')
    path = root/rel
    need(path.is_file() and not path.is_symlink(), prefix+':MEMBER_TYPE:'+str(rel))
    need(all(not p.is_symlink() for p in path.parents), prefix+':ANCESTOR_SYMLINK')
    need(path.stat().st_size == row['bytes'] and sha(path) == row['sha256'],
         prefix+':MEMBER_BYTES:'+str(rel))
    return str(rel)


def authenticate_capsule(root, pin):
    need(root.is_dir() and not root.is_symlink(), 'CAPSULE_ROOT:'+root.name)
    need(sha(root/'MANIFEST.json') == pin, 'CAPSULE_MANIFEST:'+root.name)
    seen = {'MANIFEST.json'}
    for row in read(root/'MANIFEST.json')['files']:
        rel = checked_member(root, row, root.name)
        need(rel not in seen, 'DUPLICATE_MEMBER:'+rel)
        seen.add(rel)
    for path in root.rglob('*'):
        need(not path.is_symlink(), 'CAPSULE_TREE_SYMLINK:'+root.name)
        if path.is_file():
            need(str(path.relative_to(root)) in seen, 'UNLISTED_MEMBER:'+root.name)
    return {'manifest_sha256': pin, 'members_including_manifest': len(seen)}


def authenticate(repo):
    result = {}
    for name, pin in CAPS.items():
        result[name] = authenticate_capsule(repo/'research'/name, pin)
    transfer = repo/'research/HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031'
    pins = read(transfer/'binding/SOURCE_PINS.json')
    for row in pins['files']:
        checked_member(repo/'research', row, 'PR112_SELECTED_SOURCE')
    for row in pins['capsules']:
        need(row['manifest'] == 'MANIFEST.json', 'STANDARD_CAPSULE_MANIFEST')
        result[row['path']] = authenticate_capsule(repo/'research'/row['path'], row['manifest_sha256'])
    result['PR112_selected_source_count'] = len(pins['files'])
    return result


def combined(p, m, a, b, g):
    wp = m if POSITIVE_WEIGHT_IS_EXP_MINUS else p
    return NEUMANN_COMBO_SIGN*p*g + wp*a + REFLECTED_COMBO_SIGN*p*b


def radial(a, b, g):
    return g+a+RADIAL_REFLECTED_SIGN*b


def formula_witnesses():
    # Independent cosh/sinh contraction at e^chi=2; Fplus=3, Fminus=5,
    # constant G=7, T=5/4.  This tests a nonzero boost and both endpoint signs.
    p, m, g = Q(2), Q(1, 2), Q(7)
    a, b = T*3/2, T*5/2
    ux, ur = -g+a+b, g+a-b
    direct = (p+m)*ux/2-(p-m)*ur/2
    need(direct == Q(-109, 16), 'INDEPENDENT_BOOST_WITNESS_VALUE')
    need(combined(p, m, a, b, g) == direct, 'COMBINED_KERNEL_BOOST_IDENTITY')
    need(radial(a, b, g) == ur, 'RADIAL_KERNEL_ENDPOINT_SIGN')
    # Complete derivative tail is a max norm bound; boost l1 gain is e^|chi|.
    need(FEEDBACK_GAIN*Q(2)*Q(1, 7) == Q(2, 7), 'COMPLETE_FEEDBACK_GAIN')
    need(SECOND_DERIVATIVE_NORMALIZATION*Q(3, 2)/Q(1, 1000000) == 1500000,
         'NO_EXTRA_HALF_SECOND_DERIVATIVE')
    # f(v)=|v| integrates to these first derivatives at x>r. A sign on the
    # folded time endpoint gives a different x derivative.
    ip = (T+R)*T-T*T/2
    im = (L*L+(T-L)*(T-L))/2
    need((ip+im)/2 == Q(101, 128) and (ip-im)/2 == Q(19, 128),
         'REFLECTED_PHYSICAL_ENDPOINT_DERIVATIVES')
    need(exp_point(0).lo <= 1 <= exp_point(0).hi, 'EXP_ZERO')
    ep, em = exp_point(Q(1, 4)), exp_point(Q(-1, 4))
    product = ep*em
    need(product.lo <= 1 <= product.hi, 'EXP_RECIPROCAL_WITNESS')
    return {'nonzero_boost_contraction': str(direct),
            'time_fold_polynomial_derivative': '101/128',
            'radial_fold_polynomial_derivative': '19/128',
            'feedback_l1_gain': '2', 'second_derivative_half_absent': True}


def audit_old(d):
    """Reconstruct all relevant old report arithmetic on its original lattice."""
    need(d['point'] == ['5/4', '1/8'], 'FIXED_P_DAGGER')
    need(d['all_sources_authenticated_before_imports_and_after'] is True, 'OLD_SOURCE_AUTHENTICATION_RECEIPT')
    need(d['calls'] == {'actual_first_jets': 140, 'bulk_source_cells': 130, 'shell': 9, 'apex': 1}, 'OLD_CAMPAIGN_EXACT_ROSTER')
    need(d['exact_source_role_joins']['bulk_terms'] == 78 and d['exact_source_role_joins']['Neumann_source_generally_nonzero'], 'COMPLETE_PHYSICAL_SOURCE_ROLES')
    need(d['screen_error_charges']['actual_first_errors_included_in_every_source_and_quadratic_interval'] is True, 'OLD_ACTUAL_FIRST_ERRORS_INCLUDED')
    bulk, shell, curve = d['bulk_cap_cover'], d['shell_history_cover'], d['characteristic_segments']
    need(len(bulk) == 64 and len(shell) == 8 and len(curve) == 33, 'OLD_DECLARED_BOX_ROSTER')
    for k, row in enumerate(bulk):
        i, j = divmod(k, 8)
        need(list(map(Q, row['x'])) == [T*i/8, T*(i+1)/8] and
             list(map(Q, row['r'])) == [(T+R)*j/8, (T+R)*(j+1)/8], 'BULK_COVER_GEOMETRY_%d'%k)
        box(row['F_T'])
    f = max(box(row['F_T']).mag() for row in bulk)
    need(f == Q(d['bulk_source_cap']), 'BULK_SOURCE_CAP_EXACT_MAX')
    for k, row in enumerate(shell):
        need(list(map(Q, row['x'])) == [L*k/8, L*(k+1)/8], 'FULL_NEUMANN_HISTORY_GEOMETRY_%d'%k)
    g = max(box(row['G_T']).mag() for row in shell)
    need(g == Q(d['shell_source_cap']), 'NEUMANN_SOURCE_CAP_EXACT_MAX')
    need(T <= 3 and T+(T+R) < 8, 'SOURCE_OVERCOVER_INSIDE_ACTUAL_FIRST_DOMAIN')
    knots = sorted({T*i/32 for i in range(33)} | {L})
    need(len(knots) == 34 and L in knots, 'REFLECTION_JOIN_EXPLICIT_SPLIT')
    old_integral, a, b, length = I(0), I(0), I(0), Q(0)
    for k, (lo, hi) in enumerate(zip(knots, knots[1:])):
        row = curve[k]
        need(list(map(Q, row['x'])) == [lo, hi] and Q(row['length']) == hi-lo, 'CHARACTERISTIC_CELL_TIME_%d'%k)
        need(list(map(Q, row['positive_r'])) == [T+R-hi, T+R-lo], 'POSITIVE_CHARACTERISTIC_GEOMETRY_%d'%k)
        need(list(map(Q, row['reflected_r'])) == sorted([abs(lo-L), abs(hi-L)]), 'FOLDED_CHARACTERISTIC_GEOMETRY_%d'%k)
        need(hi <= L or lo >= L, 'NO_CELL_CROSSES_REFLECTION_JOIN')
        fp, fm, dt = box(row['positive_F']), box(row['reflected_F']), hi-lo
        # This reproduces every old implicit rounding, including division
        # through the reciprocal of the singleton 2.
        pair = (fp+fm).lattice(OLD_GRID)
        weighted = (dt*pair).lattice(OLD_GRID)
        term = (weighted*Q(1, 2)).lattice(OLD_GRID)
        old_integral = (old_integral+term).lattice(OLD_GRID)
        a += dt*fp/2
        b += dt*fm/2
        length += dt
    need(length == T, 'FULL_CHARACTERISTIC_LENGTH')
    same(old_integral, d['bulk_source_derivative_integral'], 'OLD_TIME_INTEGRAL_DYADIC_RECONSTRUCTION')
    endpoint = box(d['shell_endpoint_G'])
    need(endpoint.mag() <= g, 'NEUMANN_ENDPOINT_WITHIN_HISTORY_CAP')
    old_seed = (old_integral-endpoint).lattice(OLD_GRID)
    same(old_seed, d['zero_bulk_feedback_seed_x'], 'NEUMANN_TIME_LIFT_MINUS_G')
    op = d['operator']
    caps = {k: Q(v) for k, v in op['caps'].items()}
    for k in ('A', 'B', 'D'):
        need(caps[k] == box(op['coefficient_intervals'][k]).mag(), 'OPERATOR_%s_CAP'%k)
    sigma = caps['A']+caps['B']+T*caps['D']/2
    q = T*sigma
    need(sigma == Q(op['sigma']) and q == Q(op['q']) and 0 < q < 1, 'COMPLETE_FEEDBACK_Q')
    norm = (2*g+T*f)/(1-q)
    tail = q*norm
    need(norm == Q(d['full_target_X_upper']) and tail == Q(d['complete_feedback_Ux_error']), 'FULL_C1_NORM_AND_ALL_FEEDBACK_TAIL')
    old_ux = (old_seed+I(-tail, tail).lattice(OLD_GRID)).lattice(OLD_GRID)
    same(old_ux, d['actual_full_U_Tx_interval'], 'EXISTING_ACTUAL_UX_VERIFIED')
    bg = d['apex_background']
    need(Q(bg['ell']) == ELL, 'PHYSICAL_ELL')
    rho, c = box(bg['rho']), box(bg['c'])
    rho_sq = (rho**2).lattice(OLD_GRID)
    denominator = (rho_sq*I(ELL).lattice(OLD_GRID)).lattice(OLD_GRID)
    denominator = (denominator*c).lattice(OLD_GRID)
    old_scale = (Q(3, 2)*denominator.inv().lattice(OLD_GRID)).lattice(OLD_GRID)
    same(old_scale, d['screen_scale'], 'EXISTING_POSITIVE_TT_COEFFICIENT_FACTOR_VERIFIED')
    need(old_scale.lo > 0, 'STRICTLY_POSITIVE_TT_FACTOR')
    return {'A_half_integral_Fplus': a, 'B_half_integral_Fminus': b,
            'G_endpoint': endpoint, 'source_F_cap': f, 'source_G_cap': g,
            'operator_caps': caps, 'sigma': sigma, 'q': q, 'full_C1_norm': norm,
            'feedback_derivative_error': tail, 'old_Ux_verified': old_ux,
            'tt_factor': old_scale, 'rho_at_P': rho, 'cosh_t_at_P': c,
            'old_exact_arithmetic_geometry_checks': 'PASS',
            'time_fold_sign': '+Fminus', 'radial_fold_sign': '-Fminus',
            'Neumann_time_sign': '-G(L)', 'Neumann_radial_sign': '+G(L)',
            'source_proof_status': 'AUTHENTICATED_BOX_INPUTS; source evaluations are inherited, not reproved'}


def enclose(old, chi):
    p, m = exp_box(chi), exp_box(-chi)
    gain = exp_point(chi.mag()).hi
    a, b, g = (old[k] for k in ('A_half_integral_Fplus', 'B_half_integral_Fminus', 'G_endpoint'))
    seed_combo = combined(p, m, a, b, g)
    seed_r = radial(a, b, g)
    error = FEEDBACK_GAIN*gain*old['feedback_derivative_error']
    ur = seed_r+I(-old['feedback_derivative_error'], old['feedback_derivative_error'])
    actual = seed_combo+I(-error, error)
    tt = SECOND_DERIVATIVE_NORMALIZATION*old['tt_factor']*actual
    need(p.lo > 0 and m.lo > 0, 'EXP_WEIGHTS_STRICTLY_POSITIVE')
    need((a+b-g).lo >= old['old_Ux_verified'].lo and
         (a+b-g).hi <= old['old_Ux_verified'].hi, 'EXACT_SEED_INSIDE_VERIFIED_ACTUAL_UX')
    return {'chi_uniform_all_background_roots': chi, 'exp_chi': p,
            'exp_minus_chi': m, 'feedback_boost_gain_upper': gain,
            'combined_seed': seed_combo.lattice(), 'radial_seed': seed_r.lattice(),
            'actual_U_T_r': ur.lattice(), 'actual_combined_TT_jet': actual.lattice(),
            'complete_combined_feedback_error': error,
            'positive_tt_coefficient_factor': old['tt_factor'],
            'metric_linear_second_TT_coefficient': tt.lattice(),
            'tt_term_sign': 'negative' if tt.hi < 0 else ('positive' if tt.lo > 0 else 'inconclusive'),
            'formula': '-exp(chi)G(L)+1/2 integral_0^T[exp(-chi)Fplus+exp(chi)Fminus]ds',
            'actual_target': 'U=LG+K_N(F+A U_x+B U_r+D U); zero TT Cauchy data',
            'chi_dependency_policy': 'Single interval contains rapidities of every enclosed background shooting root; positive exponential weights act on whole-cell source ranges; interval dependency is overenclosed.',
            'no_full_observer_sign_claim': True}


def background_input(path, expected_sha, repo, old):
    need(path.is_file() and not path.is_symlink(), 'BACKGROUND_REPORT_TYPE')
    need(sha(path) == expected_sha, 'BACKGROUND_REPORT_PIN')
    d = read(path)
    need(d.get('schema') == 'proper-clock-background-certificate-v1', 'BACKGROUND_SCHEMA')
    need(d.get('status') == 'PASS_ACTUAL_BACKGROUND_PROPER_CLOCK_SHOOT_AND_GAUSSIAN_TUBE', 'BACKGROUND_STATUS')
    target = d['target']
    need(target['same_background_event_as_retained_Pdagger'] is True and
         Q(target['scaled_x']) == T and Q(target['scaled_r']) == R and
         Q(target['t']) == ELL*T and Q(target['z']) == ELL*R, 'BACKGROUND_SAME_P_DAGGER')
    boost = d['boost']
    need(boost['chi_uniform_all_shooting_roots'] is True, 'BACKGROUND_UNIFORM_ROOT_COVERAGE')
    cbg = box(boost['cosh_t_dagger'])
    cold = old['cosh_t_at_P']
    own_cosh = (exp_point(ELL*T)+exp_point(-ELL*T))/2
    need(encode(old['rho_at_P']) == boost['rho_at_Pdagger'] and
         cold.lo <= cbg.lo <= own_cosh.lo <= own_cosh.hi <= cbg.hi <= cold.hi,
         'SAME_AUTHENTICATED_APEX_BACKGROUND_BOXES')
    pins = d['input_sha256']
    need(bool(pins), 'BACKGROUND_INPUT_PIN_MAP_PRESENT')
    for rel, pin in pins.items():
        p = Path(rel)
        need(not p.is_absolute() and '..' not in p.parts and
             (repo/p).is_file() and not (repo/p).is_symlink() and
             all(not a.is_symlink() for a in (repo/p).parents), 'BACKGROUND_INPUT_SAFE_PATH')
        need(sha(repo/p) == pin, 'BACKGROUND_SHARED_INPUT_PIN:'+rel)
    shared = 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/MANIFEST.json'
    need(pins.get(shared) == CAPS['HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'], 'BACKGROUND_SAME_ACTUAL_SOURCE_CAPSULE')
    chi = box(boost['chi_dagger'])
    need(encode(chi) == d['shoot']['chi_dagger'], 'BACKGROUND_SHOOT_BOOST_SAME_CHI')
    # The background's exponential bounds are separately checked for coherence;
    # this checker uses its own sharper rational Taylor bounds in the TT sum.
    for key, calculated in [('exp_chi_dagger', exp_box(chi)),
                            ('exp_minus_chi_dagger', exp_box(-chi))]:
        inherited = box(boost[key])
        need(inherited.lo <= calculated.lo <= calculated.hi <= inherited.hi,
             'BACKGROUND_EXP_COHERENCE:'+key)
    return chi, {'report_sha256': expected_sha, 'schema': d['schema'],
                 'point': [T, R], 'chi_uniform_all_roots': True,
                 'shared_authenticated_source_pin_count': len(pins),
                 'same_apex_rho_box_and_nested_cosh_bounds': True,
                 'background_source_pin_map': pins}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', required=True, type=Path)
    ap.add_argument('--background-report', type=Path)
    ap.add_argument('--background-sha256')
    ap.add_argument('--chi-preview', nargs=2)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(), 'FRESH_OUTPUT')
    repo = args.repo_root.resolve()
    need(repo not in args.output.resolve().parents and repo != args.output.resolve(), 'OUTPUT_OUTSIDE_FROZEN_SOURCE')
    need(all(not p.is_symlink() for p in args.output.parents), 'OUTPUT_PARENT_SYMLINK')
    own = sha(__file__)
    before = authenticate(repo)
    witnesses = formula_witnesses()
    source = repo/'research/HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031/main/POSTCONTACT_NORMAL.json'
    normal, optimized = read(source), read(source.with_name('POSTCONTACT_OPTIMIZED.json'))
    need(normal == optimized, 'FROZEN_NORMAL_OPTIMIZED_CANONICAL_EQUAL')
    need(normal['producer_sha256'] == sha(source.with_name('producer.py')) and
         normal['freeze_sha256'] == sha(source.with_name('FREEZE.json')),
         'OLD_REPORT_PRODUCER_AND_FREEZE_PIN_COHERENCE')
    old = audit_old(normal)
    binding_path = repo/'research/HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/binding/BINDING_NORMAL.json'
    binding = read(binding_path)
    need(binding['count'] == 40 and binding['inherited_contraction'] == normal['operator']['q'] and
         binding['inherited_screen_interval'] == normal['actual_full_screen_h'],
         'EXACT_FAMILY_BINDING_SAME_RETAINED_TT_TARGET')
    if args.chi_preview:
        need(args.background_report is None and args.background_sha256 is None, 'PREVIEW_WITHOUT_BACKGROUND')
        chi, bg = I(*args.chi_preview), {'status': 'PREVIEW_ONLY_NO_CERTIFIED_BACKGROUND'}
        status = 'PREVIEW_CONDITIONAL_COMBINED_TT_JET_ONLY'
    else:
        need(args.background_report is not None and args.background_sha256 is not None, 'PINNED_BACKGROUND_REQUIRED')
        chi, bg = background_input(args.background_report, args.background_sha256, repo, old)
        status = 'CERTIFIED_CONDITIONAL_ACTUAL_COMBINED_TT_JET_ONLY'
    result = {'schema': 'HDBLAST_PROPER_CLOCK_COMBINED_ACTUAL_TT_V1', 'status': status,
              'point': [T, R], 'ell': ELL, 'L': L, 'source_policy': 'AUTHENTICATED_BOX_INPUTS',
              'background_authentication': bg, 'old_inputs_audit': old,
              'enclosure': enclose(old, chi), 'formula_witnesses': witnesses,
              'source_authentication': before, 'checker_sha256': own,
              'new_source_oracle_calls': 0,
              'exact_family_TT_binding': 'Inherited conditional PR112 selected exact-family TT identification; same retained actual target and source data.',
              'limits': ['This is the metric-linear exact second TT term only; full proper-clock quadratic and angular-tilt terms must be added separately.',
                         'Retained exact-family, regularity, actual-source and constraint premises remain conditional.',
                         'Source hashes authenticate inherited box data; this checker does not repeat or reprove the physical source evaluations.',
                         'The coefficient is a second parameter derivative; no extra one-half is inserted.']}
    need(authenticate(repo) == before and sha(__file__) == own, 'SOURCES_OR_CHECKER_CHANGED_DURING_RUN')
    if args.background_report:
        need(sha(args.background_report) == args.background_sha256, 'BACKGROUND_CHANGED_DURING_RUN')
        after_chi, after_bg = background_input(args.background_report, args.background_sha256, repo, old)
        need(encode(after_chi) == encode(chi) and after_bg == bg, 'BACKGROUND_SHARED_INPUTS_CHANGED_DURING_RUN')
    result['all_inputs_authenticated_before_and_after'] = True
    args.output.write_text(json.dumps(encode(result), sort_keys=True, indent=2)+'\n')
    print(json.dumps({'status': status, 'report_sha256': sha(args.output),
                      'tt_term_sign': result['enclosure']['tt_term_sign']}))


if __name__ == '__main__':
    main()
