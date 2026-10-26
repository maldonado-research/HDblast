#!/usr/bin/env python3
"""Independent exact-jet review against the frozen full23 initial adapter."""
import argparse
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
REPO = Path('/workspace/hdblast-research-work/shell-observable-20261025/repo')


class Jet:
    def __init__(self, coefficients):
        self.c = list(map(F, coefficients))
    @property
    def n(self): return len(self.c)-1
    def cv(self, value):
        return value if isinstance(value, Jet) else Jet([value]+[0]*self.n)
    def __add__(self, value):
        value = self.cv(value)
        return Jet([a+b for a, b in zip(self.c, value.c)])
    __radd__ = __add__
    def __neg__(self): return Jet([-v for v in self.c])
    def __sub__(self, value): return self+-self.cv(value)
    def __rsub__(self, value): return self.cv(value)+-self
    def __mul__(self, value):
        value = self.cv(value)
        return Jet([sum((self.c[k]*value.c[j-k] for k in range(j+1)), F(0))
                    for j in range(min(self.n, value.n)+1)])
    __rmul__ = __mul__
    def inv(self):
        values = [1/self.c[0]]
        for j in range(1, self.n+1):
            values.append(-values[0]*sum((self.c[k]*values[j-k] for k in range(1, j+1)), F(0)))
        return Jet(values)
    def __truediv__(self, value): return self*self.cv(value).inv()
    def __rtruediv__(self, value): return self.cv(value)*self.inv()
    def __pow__(self, power):
        if power < 0: return self.inv()**(-power)
        value = self.cv(1)
        for _ in range(power): value = value*self
        return value
    def dz(self): return Jet([j*c for j, c in enumerate(self.c) if j])


def constant(value, order=6):
    return Jet([value]+[0]*order)


def equal(a, b):
    if not isinstance(a, Jet): a = constant(a, b.n)
    if not isinstance(b, Jet): b = constant(b, a.n)
    return a.c[:min(a.n, b.n)+1] == b.c[:min(a.n, b.n)+1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=REPO)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    target = args.output.resolve()
    if args.output.exists() or args.output.is_symlink():
        raise ValueError('FRESH_EXTERNAL_OUTPUT_REQUIRED')
    if target == ROOT.parent or ROOT.parent in target.parents:
        raise ValueError('OUTPUT_INSIDE_CAPSULE')
    producer = load('full23_initial_under_review', ROOT/'full23_initial.py')
    context = producer.Context(args.repo_root)
    if any(target == root or root in target.parents for root in context.roots.values()):
        raise ValueError('OUTPUT_INSIDE_FROZEN_CAPSULE')
    # Distinct nonconstant exact polynomial inputs avoid constant-background
    # identities hiding a missing derivative, scalar channel or angular slot.
    bg = {
        'rho': Jet([3, F(1, 7), F(1, 11), 0, F(-1, 13), 0, 0]),
        'eta': Jet([F(2, 5), F(-1, 9), 0, F(1, 17), 0, 0, 0]),
        'p': Jet([-2, F(1, 19), F(-1, 23), 0, 0, 0, 0]),
        'a': Jet([F(1, 3), F(1, 29), 0, F(1, 31), 0, 0, 0]),
        'alpha': Jet([1, F(-1, 37), 0, F(1, 41), 0, 0, 0]),
    }
    packet = Jet([F(1, 2), F(2, 3), F(-3, 5), F(5, 7), F(-7, 11), F(11, 13), F(-13, 17)])
    E = Jet([F(3, 4), F(-1, 4), F(1, 8), F(-1, 24), F(1, 96), 0, 0])
    rho, eta, p, a, alpha = [bg[key] for key in ('rho', 'eta', 'p', 'a', 'alpha')]
    dy = lambda value: -(1/rho)*value.dz()
    h = context.init.poly(context.init.UP[1], eta)/p
    k = rho*(F(5, 2)*a-h)
    u = alpha/rho**2*packet
    v = 3*alpha/(rho**3*p)*(packet.dz()+k*packet)
    uy, vy, ut = dy(u), dy(v), E*u
    d = 3/(rho*p)
    vt = E*(v-d*u)
    U2 = context.init.poly(context.init.UP[2], eta)
    A = (-2*vy**2-8*p*u*vy+8*p*(h+2*a)*u*v
         +(16*p*p-192/rho**2)*u*u-(2*U2+F(8, 3)*p*p)*v*v
         -(12*ut*ut+2*vt*vt)/rho**2)
    D = (36*u*u-2*v*v)/rho**2
    sectors = context.adapter.scalar_project(A, D)
    fs = [-rho**2*source/den for source, den in zip(sectors, (-6, 10, 42))]
    fy = [dy(value) for value in fs]
    K = E/(3*rho)*(18*u*u-v*v+d*u*v)
    Pi = E/(rho*p)*(d*u*vy+(d*(h-1/rho)+10*p)*u*v
                       +(30*a-6/rho)*u*u+(2/rho-a)*v*v)
    Fmom = 3*K
    Fmom_y = dy(Fmom)
    Gmom = Fmom_y+p*Pi
    frozen_data = context.adapter.initial_frame_data(
        rho=rho, a=a, p=p, u=u, uy=uy, ut=ut, v=v, vt=vt,
        Fmom=Fmom, Gmom=Gmom, Fmom_y=Fmom_y, f=fs, fy=fy)
    reference = context.op.pack_initial(frozen_data, context.ell)
    candidate = context._registered(bg, E, packet)
    coefficient_checks = 0
    for kind in range(2):
        for i in range(23):
            if not equal(candidate[kind][i], reference[kind][i]):
                raise ValueError('ADAPTER_IDENTITY:' + context.op.LABELS[i])
            if isinstance(candidate[kind][i], Jet) and isinstance(reference[kind][i], Jet):
                coefficient_checks += min(candidate[kind][i].n, reference[kind][i].n)+1
    nonzero = [[context.op.LABELS[i] for i, value in enumerate(row)
                if isinstance(value, Jet) and any(value.c)] for row in candidate]
    expected_u = [f'{field}:{degree}' for field in ('N', 'J', 'w') for degree in (0, 2, 4)]
    expected_v = [name for name in context.op.LABELS if
                  name.split(':')[0] in ('K', 'L', 'N', 'J', 'w', 'u')]
    if nonzero != [expected_u, expected_v]:
        raise ValueError('EXACT_NONZERO_INITIAL_CHANNELS')
    # Use the actual pair implementation with exact jets and three nonzero
    # primitive channels. This checks algebraic telescoping independently of
    # interval range width or floating point cancellation.
    fake = SimpleNamespace(constant=constant)
    Pair = producer.Context.make_pair_type(fake)
    b = Jet([F(2), F(1, 3), F(1, 5), F(1, 7), 0, 0, 0])
    d1 = Jet([F(1, 13), F(1, 17), 0, 0, 0, 0, 0])
    d2 = Jet([F(-1, 19), 0, F(1, 23), 0, 0, 0, 0])
    d3 = Jet([F(1, 29), 0, 0, F(1, 31), 0, 0, 0])
    x = Pair(b+d1+d2+d3, b, [d1, d2, d3])
    z = Pair(b+2*d1-3*d2+5*d3, b, [2*d1, -3*d2, 5*d3])
    expressions = [x+z, x*z, x/z, (3+x)**3, (1+x*z).inv(), (x*z).dz()]
    for index, expression in enumerate(expressions):
        if not equal(expression.a-expression.b, sum(expression.d, constant(0))):
            raise ValueError('EXACT_TELESCOPING_IDENTITY:'+str(index))
    # A one-sided exact cutoff endpoint control shows the sharp regularity.
    # At the upper cutoff tt=1+z, S and1 agree through derivative4.
    # Their induced scalar displacement and mixed velocity differ at order2.
    simple_bg = {key: constant(value) for key, value in
                 {'rho': 3, 'eta': F(2, 5), 'p': -2, 'a': F(1, 3), 'alpha': 1}.items()}
    tt = Jet([1, 1, 0, 0, 0, 0, 0])
    left_packet = 126*tt**5-420*tt**6+540*tt**7-315*tt**8+70*tt**9
    one = constant(1)
    cutoff_left = context._registered(simple_bg, one, left_packet)
    cutoff_right = context._registered(simple_bg, one, one)
    scalar_jump = cutoff_left[0][context.op.IDX['w:0']]-cutoff_right[0][context.op.IDX['w:0']]
    mixed_jump = cutoff_left[1][context.op.IDX['L:0']]-cutoff_right[1][context.op.IDX['L:0']]
    if scalar_jump.c[:2] != [F(0), F(0)] or mixed_jump.c[:2] != [F(0), F(0)]:
        raise ValueError('EXACT_CUTOFF_C1_MATCH')
    if scalar_jump.c[2] != F(-140, 9) or mixed_jump.c[2] != F(7, 75000):
        raise ValueError('EXACT_CUTOFF_SECOND_DERIVATIVE_JUMP')
    context.reauthenticate()
    report = {
        'status': 'PASS_EXACT_RATIONAL_JET_FULL23_INITIAL_ADAPTER_REVIEW',
        'producer_sha256': hashlib.sha256((ROOT/'full23_initial.py').read_bytes()).hexdigest(),
        'reference': 'Authenticated G/projection/inputs/source_adapters.initial_frame_data and linear_operator.pack_initial',
        'full23_adapter_slots_checked': 46,
        'exact_jet_coefficients_compared': coefficient_checks,
        'nonzero_displacements': nonzero[0], 'nonzero_velocities': nonzero[1],
        'telescoping_expressions_checked': len(expressions),
        'time_scaling': 'pack_initial multiplies every ordinary t velocity by ell; exact identity passed.',
        'regularity_finding': 'Both U and V are generally C1 at66/64. L velocity contains f_y, hence packet third derivatives. Cellwise V_rr bounds remain valid because cutoff branches align with mesh.',
        'exact_cutoff_regularity_control': {
            'packet': 'S(1+z) on the left;1 on the right; constant rational backgrounds',
            'scalar_displacement_left_minus_right_first_coefficients': list(map(str, scalar_jump.c[:3])),
            'mixed_L_velocity_left_minus_right_first_coefficients': list(map(str, mixed_jump.c[:3])),
            'conclusion': 'Value and first derivative match; second derivative may jump.'},
        'support': 'Actual packet data vanish on [0,1], nominal data on [0,65/64]; candidate only registered through7/4.',
        'complete_roots_authenticated_before_and_after': True,
        'limits': ['Algebraic adapter review, not a separate numerical error campaign.',
                   'Exact manufactured polynomial jets validate all slots and derivative couplings; retained physical shell input validity remains a premise.'],
    }
    with args.output.open('x') as handle:
        handle.write(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
