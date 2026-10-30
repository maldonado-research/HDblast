#!/usr/bin/env python3
"""Exact polynomial controls for a Robin-preserving actual-bulk update.

No HDBLAST evaluator is imported. A manufactured full variable-coefficient
solution supplies an independent reference. The nilpotent Robin fixture makes
its boundary resolvent finite, while bulk matrices are dense, variable and do
not commute with the mixed projector or Robin matrix.
"""
from fractions import Fraction as Q
from math import comb, factorial
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def poly(items):
    return {ij: Q(z) for ij, z in items.items() if z}


def add(a, b):
    out = dict(a)
    for ij, z in b.items():
        out[ij] = out.get(ij, Q(0)) + z
    return poly(out)


def scale(a, c):
    return poly({ij: c*z for ij, z in a.items()})


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    out = {}
    for (i, j), z in a.items():
        for (k, l), w in b.items():
            ij = i+k, j+l
            out[ij] = out.get(ij, Q(0)) + z*w
    return poly(out)


def deriv(a, axis):
    out = {}
    for ij, z in a.items():
        if ij[axis]:
            kk = list(ij)
            kk[axis] -= 1
            out[tuple(kk)] = ij[axis]*z
    return poly(out)


def dt(a): return deriv(a, 0)
def dr(a): return deriv(a, 1)
def wave(a): return sub(dt(dt(a)), dr(dr(a)))
def trace(a): return poly({(i, 0): z for (i, j), z in a.items() if j == 0})
def initial(a): return poly({(0, j): z for (i, j), z in a.items() if i == 0})
def reflect(a): return poly({(i, j): (-1)**j*z for (i, j), z in a.items()})
def value(a, x, r): return sum((z*x**i*r**j for (i, j), z in a.items()), Q(0))


def J(a):
    need(all(j == 0 for i, j in a), 'time primitive requires trace polynomial')
    return poly({(i+1, 0): z/Q(i+1) for (i, j), z in a.items()})


def char(a):
    need(all(j == 0 for i, j in a), 'characteristic lift requires trace polynomial')
    out = {}
    for (i, _), z in a.items():
        for j in range(i+1):
            ij = i-j, j
            out[ij] = out.get(ij, Q(0)) + z*comb(i, j)*(-1)**j
    return poly(out)


def free(u, v):
    out = {}
    for (_, n), z in u.items():
        for j in range(n//2+1):
            out = add(out, {(2*j, n-2*j): z*Q(comb(n, 2*j))})
    for (_, n), z in v.items():
        for j in range(n//2+1):
            out = add(out, {(2*j+1, n-2*j): z*Q(factorial(n), factorial(n-2*j)*factorial(2*j+1))})
    return out


def K(a):
    """Zero-initial full-line wave response, for a parity-compatible source.

    Exact polynomial time integration solves the recurrence independently of
    the physical characteristic-integral reference below.
    """
    out = {}
    for (i, n), z in a.items():
        for j in range(n//2+1):
            c = z*Q(factorial(n)*factorial(i),
                    factorial(n-2*j)*factorial(i+2*j+2))
            out = add(out, {(i+2*j+2, n-2*j): c})
    return out


def line(a, c, d, lo, hi):
    """Integral of a(s,c+d*s) ds by exact polynomial expansion."""
    out = Q(0)
    for (i, j), z in a.items():
        for n in range(j+1):
            power = i+n+1
            out += z*comb(j, n)*c**(j-n)*d**n*(hi**power-lo**power)/power
    return out


def strip(a, lower, upper, lo, hi):
    anti = {(i, j+1): z/Q(j+1) for (i, j), z in a.items()}
    return line(anti, *upper, lo, hi)-line(anti, *lower, lo, hi)


def physical_K(a, parity, x, r):
    if r >= x:
        return strip(a, (r-x, 1), (r+x, -1), Q(0), x)/2
    t = x-r
    return (strip(a, (0, 0), (r+x, -1), 0, t)
            +parity*strip(a, (0, 0), (t, -1), 0, t)
            +strip(a, (-t, 1), (r+x, -1), t, x))/2


def va(a, b): return [add(x, y) for x, y in zip(a, b)]
def vs(a, b): return [sub(x, y) for x, y in zip(a, b)]
def vc(a, c): return [scale(x, c) for x in a]
def vd(a, axis): return [deriv(x, axis) for x in a]
def vw(a): return [wave(x) for x in a]
def vt(a): return [trace(x) for x in a]
def vJ(a): return [J(x) for x in a]
def vK(a): return [K(x) for x in a]
def PN(a): return [z if s > 0 else {} for z, s in zip(a, P)]
def PD(a): return [z if s < 0 else {} for z, s in zip(a, P)]


def mv(matrix, vector):
    out = []
    for row in matrix:
        p = {}
        for a, z in zip(row, vector):
            p = add(p, mul(a, z))
        out.append(p)
    return out


def matconstant(rows): return [[{(0, 0): Q(z)} if z else {} for z in row] for row in rows]
def at_origin(matrix): return [[trace(initial(z)) for z in row] for row in matrix]
def matvalue(matrix, x, r): return [[value(z, x, r) for z in row] for row in matrix]
def scalar_mv(matrix, vector): return [sum((a*z for a, z in zip(row, vector)), Q(0)) for row in matrix]


P = (1, -1, 1)
ZERO = [{}, {}, {}]
ONE = {(0, 0): Q(1)}
X = {(1, 0): Q(1)}
R = {(0, 1): Q(1)}
RR = mul(R, R)
aa = ((1, 2, -1), (3, -2, 4), (2, 1, 3))
bb = ((-2, 3, 1), (1, 4, -3), (5, -1, 2))
dd = ((2, -1, 4), (-3, 2, 1), (1, 5, -2))
A = [[scale(mul(add(add(ONE, X), RR), R if P[i] != P[j] else ONE), Q(aa[i][j], 10))
      for j in range(3)] for i in range(3)]
B = [[scale(mul(add(add(ONE, mul(X, X)), RR), R if P[i] == P[j] else ONE), Q(bb[i][j], 11))
      for j in range(3)] for i in range(3)]
D = [[scale(mul(add(add(ONE, X), RR), R if P[i] != P[j] else ONE), Q(dd[i][j], 12))
      for j in range(3)] for i in range(3)]
A0, B0, D0 = map(at_origin, (A, B, D))
M = matconstant(((0, 0, Q(1, 7)), (0, 0, 0), (0, 0, 0)))
M0 = matconstant(((0, 0, Q(1, 8)), (0, 0, 0), (0, 0, 0)))


def lower(vector, coefficients=(A, B, D)):
    a, b, d = coefficients
    return va(va(mv(a, vd(vector, 0)), mv(b, vd(vector, 1))), mv(d, vector))


def robin(base, source, post_contact, matrix=M, lift_sign=-1):
    """Exact (I-B)^-1(base+LG) when M^2=0; no inverse of M."""
    h = PN(vt(base))
    # Trace of base+LG is h-JG; (I-B)^-1=I+B because B^2=0.
    load = va(source, mv(matrix, vs(h, vJ(source))))
    if not post_contact:
        return base
    return va(base, vc([char(z) for z in vJ(load)], lift_sign))


def close_homogeneous(base, post_contact):
    return robin(base, ZERO, post_contact)


checks = []
u = [{(0, 4): Q(1, 3)}, {(0, 5): Q(2, 5)}, {(0, 6): Q(-1, 7)}]
v = [{(0, 6): Q(-2, 9)}, {(0, 7): Q(1, 11)}, {(0, 4): Q(3, 13)}]
Bs = [free(a, b) for a, b in zip(u, v)]
need(vw(Bs) == ZERO, 'free-wave principal identity')
need([initial(z) for z in Bs] == u and [initial(dt(z)) for z in Bs] == v,
     'free-wave initial data')
checks.extend(('free_wave_identity', 'free_initial_data'))
f0 = lower(Bs, (A0, B0, D0))
W0 = va(Bs, vK(f0))
G = vc(mv(M, PN(vt(W0))), -1)
need(any(G), 'nonzero actual affine shell source')
need(mv(M, mv(M, PN(vt(W0)))) == ZERO, 'nilpotent Robin fixture')
need(robin(W0, G, True) == W0, 'old frozen-source reconstruction in manufactured fixture')
checks.extend(('nonzero_affine_shell_source', 'nilpotent_Robin_fixture', 'frozen_source_reconstruction'))

# Manufacture an independently known exact target. Its initial data differ
# from the smoothed/free candidate data, so an omitted initial bias fails.
Einitial = [free({(0, 4): Q(1, 17)}, {(0, 6): Q(1, 19)}), {}, {}]
U = va(W0, Einitial)
Fbulk = vs(vw(U), lower(U))
need(any(Fbulk), 'nonzero actual affine bulk source')
need(PD(vt(U)) == ZERO and vs(vs(PN(vt(vd(U, 1))), mv(M, PN(vt(U)))), G) == ZERO,
     'manufactured target actual affine boundary law')
need(vw(U) == va(lower(U), Fbulk), 'manufactured exact variable-coefficient target')
checks.extend(('nonzero_affine_bulk_source', 'manufactured_target_boundary', 'manufactured_target_bulk'))
f1 = va(lower(W0), Fbulk)
r0 = vs(f1, f0)
decomposed = va(va(vs(lower(Bs), f0), lower(vs(W0, Bs))), Fbulk)
need(r0 == decomposed, 'complete residual decomposition')
S1 = va(Bs, vK(f1))

for name, field in (('frozen', f0), ('actual_update', f1), ('bulk_residual', r0)):
    need(all(reflect(z) == scale(z, parity) for z, parity in zip(field, P)),
         'reference polynomial output-parity compatibility')
    for x, r in ((Q(1, 8), Q(3, 8)), (Q(1, 8), Q(1, 8)),
                 (Q(3, 8), Q(1, 8)), (Q(1, 2), Q(0))):
        for i, z in enumerate(field):
            need(value(K(z), x, r) == physical_K(z, P[i], x, r),
                 'independent physical reflected-source integral')
    checks.append('physical_reflection_reference_' + name)

for post in (False, True):
    W1 = robin(S1, G, post)
    correction = close_homogeneous(vK(r0), post)
    need(vs(W1, W0) == correction, 'Robin-preserving residual correction identity')
    need(vw(W1) == f1, 'actual bulk enters once')
    forcing_defect = vs(va(lower(W1), Fbulk), vw(W1))
    need(forcing_defect == lower(vs(W1, W0)), 'actual-RHS-minus-principal forcing-defect sign')
    need(any(forcing_defect), 'one update is not full bulk fixed point')
    need(vs(U, W1) == va(close_homogeneous(Einitial, post),
                         close_homogeneous(vK(lower(vs(U, W0))), post)),
         'error identity includes actual initial-data bias')
    need(any(close_homogeneous(Einitial, post)), 'nonzero initial-data bias retained')
    if not post:
        need([initial(z) for z in W1] == u and [initial(dt(z)) for z in W1] == v,
             'physical initial branch data unchanged')
    checks.extend(('residual_correction_post_' + str(post),
                   'actual_bulk_update_post_' + str(post),
                   'nonzero_remaining_bulk_defect_post_' + str(post),
                   'initial_bias_error_identity_post_' + str(post)))

W1 = robin(S1, G, True)
correction = close_homogeneous(vK(r0), True)
need(PD(vt(W1)) == ZERO and vs(vs(PN(vt(vd(W1, 1))), mv(M, PN(vt(W1)))), G) == ZERO,
     'new actual affine boundary identity')
need(vs(PN(vt(vd(correction, 1))), mv(M, PN(vt(correction)))) == ZERO,
     'correction homogeneous actual Robin identity')
checks.extend(('updated_actual_affine_boundary', 'correction_homogeneous_actual_boundary'))

# The characteristic correction has compatible jets on x=r. Polynomial
# substitution r=x proves equality along the whole interface, not just points.
def contact(poly0):
    out = {}
    for (i, j), z in poly0.items():
        out[(i+j, 0)] = out.get((i+j, 0), Q(0)) + z
    return poly(out)
delta_branch = vs(robin(S1, G, True), robin(S1, G, False))
for order, derivatives in enumerate(((lambda z: z,), (dt, dr),
                                     (lambda z: dt(dt(z)), lambda z: dt(dr(z)), lambda z: dr(dr(z))))):
    need(all(contact(derivative(z)) == {} for derivative in derivatives for z in delta_branch),
         'fixture contact jets through order two')
    checks.append('contact_jet_order_' + str(order))

am = matvalue(A, Q(1, 4), Q(1, 5))
mm = matvalue(M, 0, 0)
need(any(am[i][j] and P[i] != P[j] for i in range(3) for j in range(3)),
     'actual point matrix does not commute with parity')
need(any(sum(mm[i][k]*am[k][j] for k in range(3)) !=
         sum(am[i][k]*mm[k][j] for k in range(3)) for i in range(3) for j in range(3)),
     'actual bulk and Robin matrices do not commute')
checks.append('noncommuting_actual_bulk_fixture')

x0, r0point = Q(1, 4), Q(1, 5)
y = [value(z, x0, r0point) for z in W0]
yx = [value(z, x0, r0point) for z in vd(W0, 0)]
yr = [value(z, x0, r0point) for z in vd(W0, 1)]
affine = [value(z, x0, r0point) for z in Fbulk]
wrong_reflected = scalar_mv(matvalue(A, x0, r0point), [s*z for s, z in zip(P, yx)])
for part in (scalar_mv(matvalue(B, x0, r0point), [-s*z for s, z in zip(P, yr)]),
             scalar_mv(matvalue(D, x0, r0point), [s*z for s, z in zip(P, y)]),
             [s*z for s, z in zip(P, affine)]):
    wrong_reflected = [a+b for a, b in zip(wrong_reflected, part)]
correct_reflected = [s*value(z, x0, r0point) for s, z in zip(P, f1)]
need(wrong_reflected != correct_reflected,
     'reflecting input fields before untransformed physical coefficients must fail')
checks.append('wrong_input_reflection_order_rejected')

wrong_m = robin(S1, G, True, M0)
need(any(vs(vs(PN(vt(vd(wrong_m, 1))), mv(M, PN(vt(wrong_m)))), G)),
     'substituting frozen M0 fails actual boundary')
need(any(vs(vs(PN(vt(vd(robin(S1, ZERO, True), 1))),
                   mv(M, PN(vt(robin(S1, ZERO, True))))), G)),
     'omitting affine G fails actual boundary')
need(any(vs(vs(PN(vt(vd(robin(S1, G, True, lift_sign=1), 1))),
                   mv(M, PN(vt(robin(S1, G, True, lift_sign=1))))), G)),
     'wrong characteristic lift sign fails actual boundary')
need(any(vs(lower(W0), lower(W0, (A0, B0, D0)))), 'freezing actual variable bulk coefficients fails')
need(any(vs(f1, lower(W0))), 'omitting actual affine bulk forcing fails')
need(vs(U, W1) != close_homogeneous(vK(lower(vs(U, W0))), True),
     'dropping initial-data error bias fails')
checks.extend(('wrong_frozen_M_rejected', 'omitted_G_rejected', 'wrong_lift_sign_rejected',
               'frozen_actual_bulk_rejected', 'omitted_F_rejected', 'omitted_initial_bias_rejected'))

# Independent rational checks of the exact majorant recurrence. These are
# algebra identities for upper certificates, not estimates of a true error.
for k, b, initial_error, residual in ((Q(1, 10), Q(1, 20), Q(2), Q(3)),
                                      (Q(3, 10000), Q(1, 10000), Q(7, 3), Q(2, 5))):
    beta = k/(1-b)
    bias = initial_error/(1-b)
    floor = initial_error/(1-k-b)
    delta = residual/(1-b)
    current = (initial_error+residual)/(1-k-b)
    for n in range(6):
        closed = floor+beta**n*delta/(1-beta)
        need(current == closed, 'geometric recurrence and initial-data majorant')
        current = bias+beta*current
    checks.append('biased_majorant_recurrence_' + str(k))

# Exact support geometry controls. The scalar envelope 1_{s+y>1} is only a
# comparison integrand; no jump source is asserted for the physical model.
def affine_at(row, t):
    return row[0]+row[1]*t


def positive_strip_area(upper, lowers, end):
    """Area where max(lowers(s))<y<upper(s), 0<=s<=end."""
    lines = [upper]+list(lowers)
    cuts = {Q(0), end}
    for i, first in enumerate(lines):
        for second in lines[i+1:]:
            if first[1] != second[1]:
                cross = (second[0]-first[0])/(first[1]-second[1])
                if 0 < cross < end:
                    cuts.add(cross)
    out = Q(0)
    for lo, hi in zip(sorted(cuts), sorted(cuts)[1:]):
        mid = (lo+hi)/2
        lower = max(lowers, key=lambda row: affine_at(row, mid))
        if affine_at(upper, mid) > affine_at(lower, mid):
            c, d = upper[0]-lower[0], upper[1]-lower[1]
            out += c*(hi-lo)+d*(hi*hi-lo*lo)/2
    return out


def active_source_kernel(x, r, sign):
    direct = positive_strip_area((r+x, -1), ((r-x, 1), (0, 0), (1, -1)), x)
    reflected = positive_strip_area((x-r, -1), ((0, 0), (1, -1)), x)
    upper_edge = x if x+r > 1 else Q(0)
    negative_length = max(x-r, Q(0))
    lower_edge = sign*negative_length if x-r > 1 else Q(0)
    lower_edge += max(x-max(Q(0), x-r, (1+x-r)/2), Q(0))
    return ((direct+sign*reflected)/2,
            (upper_edge+lower_edge)/2, (upper_edge-lower_edge)/2)


extent = Q(9, 8)
active = extent-1
area = (extent*extent-1)/2
for sign in (-1, 1):
    for i in range(9):
        for j in range(9-i):
            x, r = extent*i/8, extent*j/8
            z, zx, zr = active_source_kernel(x, r, sign)
            need(abs(z) <= area and abs(zx) <= extent and abs(zr) <= extent,
                 'local active-source annular area and characteristic-edge bounds')
    checks.append('local_active_source_geometry_parity_' + str(sign))
z, zx, zr = active_source_kernel(extent, Q(0), 1)
need(z == area and zx == extent and zr == 0,
     'even-shell comparison attains annular area and full outgoing length')
need(z > active*active/2 and zx > active,
     'short-history triangle and derivative bounds must be rejected')
z, zx, zr = active_source_kernel(extent, Q(0), -1)
need(z == zx == 0 and zr == extent, 'odd-shell full normal edge length')
checks.extend(('local_area_and_outgoing_length_sharpness',
               'wrong_short_history_area_derivative_rejected',
               'local_odd_normal_edge_sharpness'))

print(json.dumps({'status': 'PASS', 'named_check_groups': len(checks), 'checks': checks,
    'method': 'manufactured variable-coefficient polynomial target; exact physical reflected-source integrals; finite nilpotent Robin resolvent; independent initial-biased majorant identities and active-source polygon/edge geometry',
    'scope': 'finite algebra controls for one boundary-preserving actual-bulk update; fixtures do not reproduce the full HDBLAST numerical source, prove its zero-wedge assumptions, establish numerical C2 error or solve nonlinear Einstein constraints'}, sort_keys=True))
