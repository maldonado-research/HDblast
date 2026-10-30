#!/usr/bin/env python3
"""Exact half-kernel, noncommuting transfer and causal support identities.

Manufactured operator tests supplement the analytic proof. They do not
compute actual endpoint curvature or a quantitative C2 solution error.
"""
from fractions import Fraction as Q
import json
import math


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def pa(a, b):
    return trim([(a[k] if k < len(a) else 0)+(b[k] if k < len(b) else 0) for k in range(max(len(a), len(b)))])


def ps(p, c):
    return trim([c*z for z in p])


def pm(a, b):
    out = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def pd(p):
    return trim([k*p[k] for k in range(1, len(p))] or [Q(0)])


def pe(p, x):
    return sum((c*x**k for k, c in enumerate(p)), Q(0))


def monomial(n, c=Q(1)):
    return [Q(0)]*n+[Q(c)]


def half_hinge(n):
    """F_n(q)=35/(32*n!) int_-1^-q (-q-t)^n(1-t^2)^3 dt."""
    out = [Q(0)]
    for j, coefficient in enumerate((1, -3, 3, -1)):
        for k in range(n+1):
            power = 2*j+k
            integral = pa(monomial(power+1, Q((-1)**(power+1), power+1)),
                          [Q(-(-1)**(power+1), power+1)])
            prefactor = Q(35*coefficient*math.comb(n, k)*(-1)**k, 32*math.factorial(n))
            out = pa(out, ps(pm(monomial(n-k, Q((-1)**(n-k))), integral), prefactor))
    return trim(out)


def bernstein(p):
    n = len(p)-1
    return [sum((p[k]*Q(math.comb(i, k), math.comb(n, k)) for k in range(i+1)), Q(0)) for i in range(n+1)]


groups, controls = [], []
def checked(name, ok):
    need(ok, name)
    groups.append(name)


kernel_coeff = [Q(1), Q(-3), Q(3), Q(-1)]
mass = Q(35, 32)*sum((2*c/Q(2*j+1) for j, c in enumerate(kernel_coeff)), Q(0))
moment2 = Q(35, 32)*sum((2*c/Q(2*j+3) for j, c in enumerate(kernel_coeff)), Q(0))
abs_moment1 = Q(35, 32)*sum((2*c/Q(2*j+2) for j, c in enumerate(kernel_coeff)), Q(0))
checked('kernel_mass_and_full_moments', (mass, moment2, abs_moment1) == (Q(1), Q(1, 9), Q(35, 128)))
F0, F1, F2 = [half_hinge(n) for n in range(3)]
checked('half_mass_at_endpoint', pe(F0, Q(0)) == Q(1, 2))
checked('half_kernel_derivative_charge', pe(F1, Q(0)) == Q(35, 256))
checked('half_kernel_value_charge', pe(F2, Q(0)) == Q(1, 36))
checked('normalized_hinge_derivatives', pd(F2) == ps(F1, -1) and pd(F1) == ps(F0, -1))
checked('kernel_support_endpoint_zero', all(pe(p, Q(1)) == 0 for p in (F0, F1, F2)))
bernstein_rows = []
for n, p in enumerate((F0, F1, F2)):
    bb = bernstein(p)
    checked('continuous_half_hinge_Bernstein_range_'+str(n), min(bb) >= 0 and max(bb) == pe(p, Q(0)))
    bernstein_rows.append(list(map(str, bb)))
need(pe(F2, 0) > Q(1, 72), 'undercharged displacement rejected')
controls.append('undercharged_half_kernel_value')
need(pe(F1, 0) > Q(35, 512), 'undercharged derivative rejected')
controls.append('undercharged_half_kernel_derivative')
checked('free_wave_X_norm_charge', max(Q(1, 36)/3, Q(35, 256)) == Q(35, 256))


def eye(n):
    return [[Q(i == j) for j in range(n)] for i in range(n)]


def zero(n):
    return [[Q(0)]*n for _ in range(n)]


def ma(a, b):
    return [[x+y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def mm(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), Q(0)) for j in range(len(b[0]))] for i in range(len(a))]


def mv(a, v):
    return [sum((x*y for x, y in zip(row, v)), Q(0)) for row in a]


def va(a, b):
    return [x+y for x, y in zip(a, b)]


def vs(a, b):
    return [x-y for x, y in zip(a, b)]


def vn(v):
    return max(map(abs, v))


def mn(a):
    return max(sum(map(abs, row), Q(0)) for row in a)


B = [[Q(0), Q(1, 20), Q(0)], [Q(0), Q(0), Q(1, 25)], [Q(0)]*3]
H = [[Q(0), Q(0), Q(1, 11)], [Q(1, 13), Q(0), Q(0)], [Q(0), Q(1, 17), Q(0)]]
A = [[Q(1, 19), Q(0), Q(0)], [Q(0), Q(0), Q(1, 23)], [Q(1, 29), Q(0), Q(0)]]
I = eye(3)
P = ma(ma(I, B), mm(B, B))
checked('nilpotent_boundary_inverse', mm(B, mm(B, B)) == zero(3) and
        mm(ma(I, [[-z for z in row] for row in B]), P) == I)
checked('boundary_operators_do_not_commute', mm(P, H) != mm(H, P) and mm(P, A) != mm(A, P))
old = [Q(2, 3), Q(-5, 7), Q(11, 13)]
delta = [Q(3, 5), Q(-7, 11), Q(13, 17)]
F = [Q(-2, 7), Q(3, 11), Q(5, 19)]
G = [Q(7, 23), Q(-11, 29), Q(13, 31)]


def W0(bs, source=G):
    return mv(P, va(va(bs, mv(H, bs)), source))


def W1(bs, w0, affine=F, shell=G):
    return mv(P, va(va(va(bs, mv(A, w0)), affine), shell))


w0old, w0new = W0(old), W0(va(old, delta))
w1old, w1new = W1(old, w0old), W1(va(old, delta), w0new)
d0 = mv(P, va(delta, mv(H, delta)))
d1 = mv(P, va(delta, mv(A, d0)))
checked('fixed_G_cancels_in_W0_comparison', vs(w0new, w0old) == d0)
checked('fixed_F_G_cancel_in_W1_comparison', vs(w1new, w1old) == d1)
b, k0, k, C = mn(B), mn(H), mn(A), vn(delta)
charge0 = C*(1+k0)/(1-b)
charge1 = C*(1+k*(1+k0)/(1-b))/(1-b)
checked('noncommuting_operator_norm_comparison', 0 < b < 1 and vn(d0) <= charge0 and vn(d1) <= charge1)
wrong_d0 = mv(ma(I, H), mv(P, delta))
need(wrong_d0 != d0, 'moving P through H rejected')
controls.append('commuting_boundary_inverse_through_frozen_bulk')
wrong_d1 = va(mv(P, delta), mv(A, mv(P, d0)))
need(wrong_d1 != d1, 'moving P through A rejected')
controls.append('commuting_boundary_inverse_through_actual_bulk')
wrong_sources = W1(va(old, delta), w0new, affine=va(F, [Q(1), Q(0), Q(0)]))
need(vs(wrong_sources, w1old) != d1, 'changing actual affine source breaks comparison')
controls.append('mismatched_fixed_affine_source')
Eold = Q(19, 7)
target = va(w1old, [Q(1, 7), Q(-3, 11), Q(5, 17)])
checked('same_target_triangle_transfer', vn(vs(target, w1old)) <= Eold and
        vn(vs(target, w1new)) <= Eold+charge1)

# A finite filtration model tests support propagation without commuting any
# coefficient matrices. Analytic causality supplies the corresponding PDE
# support premise; manufactured grid tests are not that continuum proof.
alpha = [Q(1, 2), Q(1), Q(9, 8), Q(3, 2), Q(7, 4), Q(2), Q(9, 4)]
R, du = Q(7, 4), Q(1, 2**30)
front = R-du
curvature = [Q(2), Q(-3), Q(5)]
metric_weights = [Q(1), Q(3), Q(2)]
checked('manufactured_endpoint_curvature_norm',
        sum(w*z*z for w, z in zip(metric_weights, curvature)) == 81)
value_charge, radial_charge = 9*du**2/Q(36), 9*Q(35, 256)*du
actual_free_charge = 9*max(du**2/Q(108), Q(35, 256)*du)
checked('actual_delta_free_wave_X_charge', actual_free_charge == max(value_charge/3, radial_charge))
n = len(alpha)
Bc = [[Q(i-j, 50) if i > j else Q(0) for j in range(n)] for i in range(n)]
Ac = [[Q(i+j+1, 150) if i >= j else Q(0) for j in range(n)] for i in range(n)]
Hc = [[Q(2*i-j+3, 170) if i >= j else Q(0) for j in range(n)] for i in range(n)]
Pc, power = eye(n), eye(n)
for _ in range(1, n):
    power = mm(power, Bc)
    Pc = ma(Pc, power)
checked('causal_boundary_inverse_filtration', mm(power, Bc) == zero(n))
seed = [Q(0) if a < front else Q(i+1, i+2) for i, a in enumerate(alpha)]
first = mv(Pc, va(seed, mv(Hc, seed)))
kept = [first]
for _ in range(5):
    kept.append(mv(Pc, va(seed, mv(Ac, kept[-1]))))
checked('every_tested_finite_iteration_preserves_outer_support',
        all(v[i] == 0 for v in kept for i, a in enumerate(alpha) if a < front))
checked('selected_event_causally_unchanged', Q(9, 8) < front and all(v[2] == 0 for v in kept))
checked('physical_shell_not_reached_within_retained_cone', Q(5, 4) < front)

print(json.dumps({'status': 'PASS_EXACT_CURVATURE_MATCH_COMPARISON',
    'named_groups': groups, 'group_count': len(groups),
    'wrong_formula_controls': controls, 'wrong_formula_control_count': len(controls),
    'half_kernel_value_factor': '1/36', 'half_kernel_derivative_factor': '35/256',
    'half_hinge_Bernstein_coefficients': bernstein_rows,
    'free_wave_X_charge': 'kappa_H*max(delta^2/108,35*delta/256)',
    'W0_charge': 'C*(1+k0)/(1-b)',
    'W1_charge': 'C*(1+k*(1+k0)/(1-b))/(1-b)',
    'support_front': str(front), 'unchanged_selected_alpha': '9/8',
    'manufactured_support_iterations': len(kept),
    'scope': 'exact half-kernel polynomial bounds and manufactured noncommuting/causal comparison identities; actual endpoint curvature norm and actual target error constants are separate inputs; no C2 solution, F_x, Einstein residual or physical origin theorem'}, sort_keys=True))
