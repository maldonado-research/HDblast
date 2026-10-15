#!/usr/bin/env python3
"""Independent all-point Laurent identities for the new initial-data bounds.

No producer or scientific module is imported. Variables are formal; p is
invertible. This supplements, rather than replaces, the interval certificate.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

NAMES = ('u', 'v', 'vy', 'a', 'p', 'h', 'w', 'E', 'U2', 'S', 'V')
N = len(NAMES)


class L:
    def __init__(self, terms=0):
        if isinstance(terms, dict):
            self.terms = {k: Q(v) for k, v in terms.items() if v}
        else:
            self.terms = {(0,) * N: Q(terms)} if terms else {}

    def __add__(self, other):
        other = cast(other)
        d = self.terms.copy()
        for k, v in other.terms.items():
            d[k] = d.get(k, Q(0)) + v
        return L(d)

    __radd__ = __add__

    def __neg__(self):
        return L({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + -cast(other)

    def __rsub__(self, other):
        return cast(other) + -self

    def __mul__(self, other):
        other = cast(other)
        d = {}
        for a, x in self.terms.items():
            for b, y in other.terms.items():
                k = tuple(i + j for i, j in zip(a, b))
                d[k] = d.get(k, Q(0)) + x * y
        return L(d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = cast(other)
        if len(other.terms) != 1:
            raise ValueError('Only monomial Laurent division is supported')
        k, v = next(iter(other.terms.items()))
        return self * L({tuple(-j for j in k): 1 / v})

    def __pow__(self, power):
        if power < 0:
            return cast(1) / (self ** -power)
        out = cast(1)
        for _ in range(power):
            out = out * self
        return out

    def diff(self, name):
        j = NAMES.index(name)
        out = {}
        for k, v in self.terms.items():
            if k[j]:
                q = list(k)
                q[j] -= 1
                out[tuple(q)] = v * k[j]
        return L(out)


def cast(value):
    return value if isinstance(value, L) else L(value)


def var(name):
    k = [0] * N
    k[NAMES.index(name)] = 1
    return L({tuple(k): 1})


def zero(expression, label):
    if expression.terms:
        raise ValueError(label + ': nonzero all-point Laurent residual')


def check(injection):
    u, v, j, a, p, h, w, E, U2, S, V = map(var, NAMES)
    uy = -p * v / 3 - 2 * a * u
    ay = -w**2 - p**2 / 3
    py = p * (h - 4 * a)
    rules = {'u': uy, 'v': j, 'a': ay, 'p': py, 'w': -a*w, 'E': w*E}

    def dy(expression):
        for name in ('vy', 'h', 'U2', 'S', 'V'):
            if expression.diff(name).terms:
                raise ValueError('Unexpected derivative needed: ' + name)
        return sum((expression.diff(name) * rule for name, rule in rules.items()), cast(0))

    uyy = dy(uy)
    d = 3*w/p
    ut, vt = E*u, E*(v-d*u)
    uyt, vyt = dy(ut), dy(vt)

    # Independently differentiate the original published amplitude F.
    F = (18*u*ut-v*vt)*w
    G = (12*u*uyt-48*a*u*ut-6*uy*ut-2*vt*j)*w
    original_P = (G-dy(F))/p
    reduced_P = E*w/p * (d*u*j + (d*(h-w)+10*p)*u*v
                            + (30*a-6*w)*u*u + (2*w-a)*v*v)
    if injection == 'restore_cancelled_v_vy':
        reduced_P += E*w*v*j/p
    if injection == 'omit_momentum_radius':
        reduced_P = reduced_P/w
    zero(original_P-reduced_P, 'MOMENTUM_REDUCTION')

    original_A = (-24*uy*uy-288*a*u*uy-24*u*uyy
                  -(192*ay+384*a*a+32*p*p)*u*u
                  -w*w*(336*u*u+12*ut*ut+2*vt*vt)
                  -2*j*j-16*p*u*j-2*U2*v*v)
    reduced_A = (-2*j*j-8*p*u*j+8*p*(h+2*a)*u*v
                 +(16*p*p-192*w*w)*u*u-(2*U2+Q(8,3)*p*p)*v*v
                 -w*w*(12*ut*ut+2*vt*vt))
    if injection == 'halve_scalar_gradient':
        reduced_A += j*j
    zero(original_A-reduced_A, 'HAMILTONIAN_REDUCTION')
    reduced_q = E*w*(18*u*u-v*v+d*u*v)/3
    zero(F/3-reduced_q, 'TRACE_K_REDUCTION')

    # On S3: S=n1^2+n2^2, V=Y^2, 0<=V<=S^2<=S<=1.
    # |D(Y^2)|^2=16V(S-V)<=4S^2<=4. The first inequality
    # is exactly the nonnegative square below. |DZ2|<=1 follows
    # from 4S(1-S)<=1. Then |DH4|<=|D(Y^2)|+|DS|/2<=5/2.
    cap = 4*S*S
    if injection == 'halve_angular_gradient_cap':
        cap = S*S
    zero(cap-16*V*(S-V)-4*(S-2*V)**2, 'ANGULAR_GRADIENT_BOUND')
    zero(1-4*S*(1-S)-(2*S-1)**2, 'ANGULAR_Z2_BOUND')
    # Polynomial endpoint comparison for H4=V-S/2+1/12:
    # -5/12<=H4<=7/12, since V>=0 and V<=S^2, 0<=S<=1.
    zero(Q(7,12)-(S*S-S/2+Q(1,12))-(1-S)*(S+Q(1,2)),
         'ANGULAR_H4_UPPER_BOUND')
    zero((V-S/2+Q(1,12))+Q(5,12)-V-(1-S)/2,
         'ANGULAR_H4_LOWER_BOUND')
    return {
        'status': 'PASS_ALL_POINT_REDUCTIONS_AND_ANGULAR_MAJORANTS',
        'all_point_source_reductions': 3,
        'angular_nonnegative_identity_certificates': 4,
        'hypotheses': ['p != 0', 'w = 1/rho', 'h = U_prime/p',
                       'p_y = p(h-4a)', 'a_y = -w^2-p^2/3',
                       'u_y = -pv/3-2au', 'E_y=wE',
                       '0<=Y^2<=S^2<=S<=1 on S3'],
        'angular_value_sup': ['1', '1/2', '7/12'],
        'angular_gradient_sup': ['0', '1', '5/2'],
        'gradient_Y_squared_sup': '2',
        'scope': 'Formal identities and angular bounds only; no interval norm, evolution or nonlinear-family theorem',
    }


CONTROLS = {
    'restore_cancelled_v_vy': 'MOMENTUM_REDUCTION',
    'omit_momentum_radius': 'MOMENTUM_REDUCTION',
    'halve_scalar_gradient': 'HAMILTONIAN_REDUCTION',
    'halve_angular_gradient_cap': 'ANGULAR_GRADIENT_BOUND',
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--inject', choices=CONTROLS)
    ap.add_argument('--run-controls', action='store_true')
    args = ap.parse_args()
    if args.output.exists():
        raise ValueError('Refusing to overwrite existing output')
    report = check(args.inject)
    if args.run_controls:
        outcomes = []
        with tempfile.TemporaryDirectory() as tmp:
            for name, reason in CONTROLS.items():
                out = Path(tmp) / (name+'.json')
                cmd = [sys.executable]
                if not __debug__:
                    cmd += ['-O']
                cmd += ['-I', '-S', '-B', str(Path(__file__).resolve()),
                        '--inject', name, '--output', str(out)]
                r = subprocess.run(cmd, text=True, capture_output=True)
                if r.returncode == 0 or reason not in r.stderr or out.exists():
                    raise ValueError('NEGATIVE_CONTROL_FAILED: '+name)
                outcomes.append({'name': name, 'expected_reason': reason,
                                 'returncode': r.returncode, 'no_output': True})
        report['negative_controls'] = outcomes
    report['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, sort_keys=True, indent=2)+'\n')


if __name__ == '__main__':
    main()
