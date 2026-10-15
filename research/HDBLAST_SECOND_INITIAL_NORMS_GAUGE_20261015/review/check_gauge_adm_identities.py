#!/usr/bin/env python3
"""All-point ADM epsilon identities, independent of the direct 5D checker.

Only the generic Laurent arithmetic engine from check_reduced_sources is used;
no producer or gauge-result module is imported. Coefficients are ordinary
epsilon Taylor coefficients and second derivatives include the factor two.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ENGINE = Path(__file__).with_name('check_reduced_sources.py')
spec = importlib.util.spec_from_file_location('formal_arithmetic', ENGINE)
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)
L, var, zero = engine.L, engine.var, engine.zero


class S:
    def __init__(self, values=0):
        if isinstance(values, (tuple, list)):
            self.c = tuple(engine.cast(x) for x in values)
        else:
            self.c = (engine.cast(values), L(0), L(0))

    def __add__(self, other):
        other = series(other)
        return S([a+b for a, b in zip(self.c, other.c)])

    __radd__ = __add__

    def __neg__(self):
        return S([-x for x in self.c])

    def __sub__(self, other):
        return self + -series(other)

    def __rsub__(self, other):
        return series(other) + -self

    def __mul__(self, other):
        other = series(other)
        return S([sum((self.c[j]*other.c[k-j] for j in range(k+1)), L(0))
                  for k in range(3)])

    __rmul__ = __mul__

    def __truediv__(self, other):
        a, b, c = series(other).c
        inverse = S([L(1)/a, -b/a**2, b*b/a**3-c/a**2])
        return self*inverse

    def __rtruediv__(self, other):
        return series(other)/self

    def __pow__(self, exponent):
        out = S(1)
        for _ in range(exponent):
            out = out*self
        return out


def series(x):
    return x if isinstance(x, S) else S(x)


def check(injection):
    # The arithmetic engine's indeterminates have no scientific meaning here.
    U, Uy, T, a, fy, f, unused, Ui, fi, q, rho = map(var, engine.NAMES)
    ep = S([0, 1, 0])
    J = 1+ep*U
    L2 = 1-4*ep*U
    M = 1+2*ep*U+ep*ep*f/2
    Jy, L2y, My = ep*Uy, -4*ep*Uy, 2*ep*Uy+ep*ep*fy/2
    Ji, L2i, Mi = ep*Ui, -4*ep*Ui, 2*ep*Ui+ep*ep*fi/2

    # Exact initial C=0 lapse/shift formulas, specialized only after forming
    # their rational ADM expressions. Reference shift is zero, N0=rho,
    # reference slice time-symmetric and Gamma_bar^y_tt=rho^2*a.
    Kyy = -2*ep*T/rho+ep*ep*q/2
    Kang = ep*T/rho+ep*ep*q/2
    dotN = J*J*rho*rho*(Kyy/L2+3*Kang/M)
    dotBy = J*J*rho*rho*(L2y/(2*L2*L2)-S(4*a)/L2
                          -3*My/(2*L2*M)-Jy/(L2*J)+S(3*a)/M)+rho*rho*a
    dotBi = J*J*(-L2i/(2*M*L2)-Mi/(2*M*M)-Ji/(M*J))

    n1, n2 = rho*T, 4*rho*rho*q-24*rho*U*T
    by1 = -6*rho*rho*(Uy+4*a*U)
    by2 = rho*rho*(-74*U*Uy-194*a*U*U-3*fy/2-3*a*f)
    bi2 = 22*U*Ui-fi/2
    if injection == 'halve_lapse_cross':
        n2 += 12*rho*U*T
    if injection == 'omit_radial_metric_cross':
        by2 += 74*rho*rho*U*Uy
    if injection == 'wrong_angular_coefficient':
        bi2 -= 11*U*Ui

    tests = [('LAPSE_CONSTANT', dotN.c[0]), ('LAPSE_FIRST', dotN.c[1]-n1),
             ('LAPSE_SECOND', 2*dotN.c[2]-n2),
             ('RADIAL_SHIFT_CONSTANT', dotBy.c[0]),
             ('RADIAL_SHIFT_FIRST', dotBy.c[1]-by1),
             ('RADIAL_SHIFT_SECOND', 2*dotBy.c[2]-by2),
             ('ANGULAR_SHIFT_CONSTANT', dotBi.c[0]),
             ('ANGULAR_SHIFT_FIRST', dotBi.c[1]),
             ('ANGULAR_SHIFT_SECOND', 2*dotBi.c[2]-bi2)]
    for label, residual in tests:
        zero(residual, label)

    # g_ti=gamma_ij beta^j: its time derivative contains the metric/shift
    # cross term even though beta itself vanishes initially at every order.
    dotgty = L2*dotBy
    gty2 = by2-8*U*by1
    if injection == 'omit_mixed_metric_conversion':
        gty2 = by2
    zero(2*dotgty.c[2]-gty2, 'MIXED_METRIC_CONVERSION')
    zero((by2-8*U*by1)/(rho*rho)+2*a*U*U
         -(-26*U*Uy-3*fy/2-3*a*f), 'NORMAL_MIXED_FRAME_CANCELLATION')
    zero(-2*n2/(rho*rho)-4*U*T/rho
         -(-8*q+44*U*T/rho), 'NORMAL_LAPSE_FRAME_CANCELLATION')
    return {'status': 'PASS_ALL_POINT_INITIAL_ADM_GAUGE_IDENTITIES',
            'all_point_series_and_frame_identities': 12,
            'second_coefficient_convention': 'epsilon^2 X2/2',
            'hypotheses': ['zero initial shift at every parameter order',
                           'N=rho(1+epsilon U), N2=0',
                           'gamma_yy=1-4epsilon U',
                           'gamma_ij=(1+2epsilon U+epsilon^2 f/2) gamma0_ij',
                           'K1_yy=-2T/rho, K1_ij=T gamma0_ij/rho',
                           'K2=q gamma0', 'reference slice time-symmetric'],
            'scope': 'Initial gauge coefficient identities only; no gauge-wave, reduced acceleration, boundary or nonlinear-family solve',
            'engine_sha256': hashlib.sha256(ENGINE.read_bytes()).hexdigest()}


CONTROLS = {'halve_lapse_cross': 'LAPSE_SECOND',
            'omit_radial_metric_cross': 'RADIAL_SHIFT_SECOND',
            'wrong_angular_coefficient': 'ANGULAR_SHIFT_SECOND',
            'omit_mixed_metric_conversion': 'MIXED_METRIC_CONVERSION'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--inject', choices=CONTROLS)
    ap.add_argument('--run-controls', action='store_true')
    args = ap.parse_args()
    if args.output.exists():
        raise ValueError('Refusing to overwrite output')
    report = check(args.inject)
    if args.run_controls:
        report['negative_controls'] = []
        with tempfile.TemporaryDirectory() as tmp:
            for name, reason in CONTROLS.items():
                target = Path(tmp)/(name+'.json')
                cmd = [sys.executable]+(['-O'] if not __debug__ else [])
                cmd += ['-I', '-S', '-B', str(Path(__file__).resolve()),
                        '--inject', name, '--output', str(target)]
                r = subprocess.run(cmd, capture_output=True, text=True)
                if r.returncode == 0 or reason not in r.stderr or target.exists():
                    raise ValueError('NEGATIVE_CONTROL_FAILED: '+name)
                report['negative_controls'].append({'name': name, 'expected_reason': reason,
                                                    'returncode': r.returncode, 'no_output': True})
    report['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, sort_keys=True, indent=2)+'\n')


if __name__ == '__main__':
    main()
