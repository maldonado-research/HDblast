#!/usr/bin/env python3
"""Portable independent witness exclusion; Python standard library only.

Run with Python -I -S -B, optionally also -O.  This program imports no
producer code and makes no background, bump, or evolved-field oracle query.
"""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import factorial
import json
import os
from pathlib import Path
import stat
import sys


RAW_REL = 'research/HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/numeric/references/RAW_WITNESS.json'
STATIC_REL = 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json'
RAW_SHA = '412c4840503166a1a1893b8f16b9a3be35a8af98643c461dd55d820c1dc79d76'
STATIC_SHA = '3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22'


class Rejected(Exception):
    pass


def require(condition, reason):
    if not condition:
        raise Rejected(reason)


class Parser(argparse.ArgumentParser):
    def error(self, message):
        raise Rejected('invalid_arguments: ' + message)


def parse_arguments(argv):
    parser = Parser(add_help=False, allow_abbrev=False)
    for flag in ('--repo-root', '--candidate', '--output'):
        parser.add_argument(flag, required=True)
    for flag in ('--repo-root', '--candidate', '--output'):
        occurrences = sum(arg == flag or arg.startswith(flag + '=') for arg in argv)
        require(occurrences <= 1, 'duplicate_argument: ' + flag)
    return parser.parse_args(argv)


def normalized_absolute(value, label):
    require(isinstance(value, str) and bool(value), label + '_path_empty')
    require('\x00' not in value, label + '_path_contains_nul')
    require(value.startswith('/') and not value.startswith('//'), label + '_path_not_normalized_absolute')
    require(os.path.normpath(value) == value, label + '_path_not_normalized_absolute')
    return Path(value)


def checked_existing_path(path, kind, label):
    current = Path('/')
    for part in path.parts[1:]:
        current = current / part
        try:
            mode = current.lstat().st_mode
        except OSError:
            raise Rejected(label + '_path_unavailable') from None
        require(not stat.S_ISLNK(mode), label + '_path_has_symlink')
    if kind == 'directory':
        require(path.is_dir(), label + '_path_not_directory')
    elif kind == 'file':
        require(path.is_file(), label + '_path_not_regular_file')
    else:
        raise Rejected('internal_path_kind')
    return path


def checked_output_path(value):
    path = normalized_absolute(value, 'output')
    checked_existing_path(path.parent, 'directory', 'output_parent')
    require(not os.path.lexists(path), 'output_already_exists')
    return path


def load_pinned(path, expected, label):
    checked_existing_path(path, 'file', label)
    try:
        data = path.read_bytes()
    except OSError:
        raise Rejected(label + '_read_failed') from None
    require(sha256(data).hexdigest() == expected, label + '_byte_pin_mismatch')
    try:
        return json.loads(data)
    except (ValueError, UnicodeError):
        raise Rejected(label + '_invalid_json') from None


class I:
    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = self.lo if hi is None else F(hi)
        require(self.lo <= self.hi, 'reversed_interval')

    def __add__(self, x):
        x = interval(x)
        return I(self.lo + x.lo, self.hi + x.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, x):
        return self + -interval(x)

    def __rsub__(self, x):
        return interval(x) + -self

    def __mul__(self, x):
        x = interval(x)
        q = [self.lo * x.lo, self.lo * x.hi, self.hi * x.lo, self.hi * x.hi]
        return I(min(q), max(q))

    __rmul__ = __mul__

    def inv(self):
        require(self.lo * self.hi > 0, 'division_interval_meets_zero')
        return I(1 / self.hi, 1 / self.lo)

    def __truediv__(self, x):
        return self * interval(x).inv()

    def __rtruediv__(self, x):
        return interval(x) * self.inv()

    def __pow__(self, n):
        require(isinstance(n, int) and n >= 0, 'invalid_interval_power')
        answer = I(1)
        for unused in range(n):
            answer = answer * self
        return answer

    def report(self):
        return {'exact': [str(self.lo), str(self.hi)],
                'display_only': [float(self.lo), float(self.hi)]}


def interval(x):
    return x if isinstance(x, I) else I(x)


def read(x):
    values = x['exact']
    require(isinstance(values, list) and len(values) == 2, 'invalid_exact_interval_shape')
    return I(*values)


def contains(outer, inner):
    return outer.lo <= inner.lo and inner.hi <= outer.hi


def fresh_Q_target():
    """Q(T)=exp(-4/3), by positive Taylor sum and geometric remainder."""
    x = F(4, 3)
    n = 80
    lo = sum((x**k / F(factorial(k)) for k in range(n + 1)), F(0))
    tail = (x**(n + 1) / F(factorial(n + 1))) / (1 - x / F(n + 2))
    require(tail > 0, 'nonpositive_exponential_remainder')
    Q = I(1 / (lo + tail), 1 / lo)
    t = F(1, 4)
    d = t * (1 - t)
    return [Q,
            Q * ((1 - 2*t) / d**2),
            Q * ((1 - 6*d + 6*d**2) / d**4),
            Q * ((1 - 2*t) * (1 - 10*d + 24*d**2 - 12*d**3) / d**6)]


def build_report(raw, static):
    require(raw['schema'] == 'HDBLAST_P1_PREPARATION_RAW_WITNESS_V1', 'raw_schema_differs')
    require(static['domain']['ell'] == '1/1000000', 'static_ell_differs')
    require(static['domain']['z_min'] == '1/1000000' and
            static['domain']['z_max'] == '1/500000', 'static_domain_differs')
    jets = static['y_jets']
    caps = {'rho': I(129, 131), 'a': I('111/1000', '14/125'),
            'p': I('-1/7000', '-3/25000')}
    for name, cap in caps.items():
        require(contains(cap, read(jets[name][0])), 'coarse_' + name + '_cap_failed')
    rb_global = read(static['annulus_refinement']['rho']['initial_shell'])
    require(contains(I(129, 131), rb_global), 'coarse_rho_b_cap_failed')
    require(contains(I('-1/10000', '1/10000'), read(jets['a'][1])), 'coarse_a_y_cap_failed')
    require(contains(I('-1/4000', '1/4000'), read(jets['p'][1])), 'coarse_p_y_cap_failed')

    target = raw['target']
    rho = read(target['background']['rho'])
    a = read(target['background']['a'])
    p = read(target['background']['p'])
    rb = read(raw['shared_inputs']['rho_b'])
    py = read(jets['p'][1])
    for name, val in [('rho', rho), ('a', a), ('p', p)]:
        require(contains(read(jets[name][0]), val), 'target_' + name + '_outside_shared_annulus_box')
    require(rb.lo == rb_global.lo and rb.hi == rb_global.hi, 'rho_b_shared_inputs_differ')
    # Broad independently checked cap; no producer auxiliary derivative is used.
    ay = I('-1/10000', '1/10000')
    fresh_jets = fresh_Q_target()
    old_jets = list(map(read, target['Q_r_jets']))
    require(len(old_jets) == len(fresh_jets), 'retained_Q_jet_count_differs')
    for old, fresh in zip(old_jets, fresh_jets):
        require(contains(old, fresh), 'retained_old_Q_jet_does_not_contain_fresh_enclosure')
    q, q1, q2, q3 = fresh_jets
    ell = F(1, 10**6)
    b = ell * rho * a
    c = -3 * rb**3 / (2 * ell * rho**4 * p)
    d = ell * rho * py / p
    R = (rb / rho)**3
    h = -R * q1
    hr = -R * (q2 + 3*b*q1)
    S = q2 + 4*b*q1 + 3*ell**2*q
    s = c * S
    sr = c * (q3 + 4*b*q2 + 3*ell**2*q1 + (4*b+d)*S
              - 4*ell**2*rho**2*(a*a+ay)*q1)
    ar = -ell * rho * ay
    pr = -ell * rho * py
    D = F(3, 2)*h - 3*a*s/p
    Dr = F(3, 2)*hr - 3*(ar*s/p + a*sr/p - a*s*pr/(p*p))
    L = 10**21
    radius = 2 * 10**33
    require(D.lo*L > radius, 'D_does_not_exclude_scaled_witness')
    require(-Dr.hi*L > radius, 'D_r_does_not_exclude_scaled_witness')

    # Separate fresh elementary proof: Q(T)>1/4, Q(T)<1, all positive
    # terms of s retained as lower bounds, and |h| bounded above.
    x = F(4, 3)
    exp_hi = sum((x**k / F(factorial(k)) for k in range(11)), F(0))
    exp_hi += (x**11 / F(factorial(11))) / (1 - x/F(12))
    require(exp_hi < 4, 'coarse_exponential_upper_bound_failed')
    cmin = F(3) * 129**3 / (2*ell*131**4*F(1, 7000))
    h_abs = F(131, 129)**3 * F(128, 9)
    coarse = 3*F(111, 1000)*7000*cmin*F(1408, 81) - F(3, 2)*h_abs
    require(coarse > 3*10**12, 'coarse_shortened_lower_bound_failed')
    require(coarse*L > radius, 'coarse_D_does_not_exclude_scaled_witness')

    return {
        'schema': 'HDBLAST_INDEPENDENT_GEOMETRIC_CLASS_WITNESS_EXCLUSION_V1',
        'status': 'PASS_INDEPENDENT_WITNESS_EXCLUSION',
        'method': 'stdlib Fraction; shared old background enclosures; fresh analytic rational target bump replay',
        'shared_old_input_pins': {RAW_REL: RAW_SHA, STATIC_REL: STATIC_SHA},
        'candidate_logical_path': RAW_REL,
        'formula': {
            'D': '3*h/2 - 3*a*s/p',
            'D_r': '3*h_r/2 -3*(a_r*s/p+a*s_r/p-a*s*p_r/p^2)',
            'a_r': '-ell*rho*a_y',
            'p_r': '-ell*rho*p_y'},
        'target': '5/4',
        'direction_D': D.report(),
        'direction_D_r': Dr.report(),
        'lambda': str(-L),
        'class_radius': str(radius),
        'witness_absolute_D_lower': str(D.lo*L),
        'witness_absolute_D_r_lower': str(-Dr.hi*L),
        'coarse_D_lower': str(coarse),
        'coarse_D_lower_display_only': float(coarse),
        'coarse_exp_4over3_upper': str(exp_hi),
        'coarse_caps_checked_from_pinned_static_source': True,
        'retained_old_Q_jets_contain_fresh_analytic_enclosures': True,
        'fresh_Q_jet_construction': 'Positive exp(4/3) Taylor sum through n=80; geometric tail ratio (4/3)/82; invert and multiply exact target derivative factors.',
        'scope': {
            'bump_evaluations': 'fresh analytic rational replay; not new oracle or bump queries',
            'fresh_analytic_target_bump_evaluations': 1,
            'new_background_oracle_queries': 0,
            'new_bump_queries': 0,
            'new_evolved_field_oracle_queries': 0,
            'producer_code_imports': 0,
            'shared_old_enclosure_premises': True,
            'excluded_preparation': 'the one fixed lambda=-10^21 witness',
            'physical_class_selection_rule': False,
            'uniform_nonlinear_epsilon_radius': False}}


def canonical_bytes(report):
    return (json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + '\n').encode('utf-8')


def write_fresh(path, data):
    # Recheck after all science checks; O_EXCL and O_NOFOLLOW prevent a final
    # existing-file or final-symlink overwrite if the path changes meanwhile.
    checked_existing_path(path.parent, 'directory', 'output_parent')
    require(not os.path.lexists(path), 'output_already_exists')
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW
    try:
        fd = os.open(path, flags, 0o600)
    except OSError:
        raise Rejected('output_create_failed') from None
    try:
        with os.fdopen(fd, 'wb') as stream:
            stream.write(data)
            stream.flush()
    except Exception:
        try:
            path.unlink()
        except OSError:
            pass
        raise Rejected('output_write_failed') from None


def run(argv):
    require(sys.flags.isolated and sys.flags.no_site and sys.dont_write_bytecode, 'required_interpreter_flags:-I_-S_-B')
    args = parse_arguments(argv)
    root = normalized_absolute(args.repo_root, 'repo_root')
    checked_existing_path(root, 'directory', 'repo_root')
    candidate = normalized_absolute(args.candidate, 'candidate')
    output = checked_output_path(args.output)
    static_path = root / STATIC_REL
    raw = load_pinned(candidate, RAW_SHA, 'candidate')
    static = load_pinned(static_path, STATIC_SHA, 'static_source')
    report = build_report(raw, static)
    data = canonical_bytes(report)
    digest = sha256(data).hexdigest()
    receipt = {'status': report['status'], 'report_sha256': digest,
               'new_oracle_queries': 0}
    receipt_text = json.dumps(receipt, sort_keys=True, separators=(',', ':'), allow_nan=False)
    write_fresh(output, data)
    return receipt_text


def main():
    try:
        receipt = run(sys.argv[1:])
    except Rejected as exc:
        print(json.dumps({'status': 'REJECTED', 'reason': str(exc)},
                         sort_keys=True, separators=(',', ':')))
        return 2
    except Exception as exc:
        print(json.dumps({'status': 'REJECTED', 'reason': 'validation_failed:' + type(exc).__name__},
                         sort_keys=True, separators=(',', ':')))
        return 2
    print(receipt)
    return 0


if __name__ == '__main__':
    sys.exit(main())
