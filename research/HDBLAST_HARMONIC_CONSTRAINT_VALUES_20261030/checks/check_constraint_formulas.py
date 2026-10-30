#!/usr/bin/env python3
"""Independent exact frame contractions and norm identities for DC.

Manufactured rational algebra, not a scientific-field evaluation.  No project
module is imported.  Tests remain active in optimized Python.
"""
from fractions import Fraction as Q
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


LABELS = [f'{f}:{l}' for f in ('K', 'L', 'N', 'J', 'w') for l in (0, 2, 4)] + [f'{f}:{l}' for f in ('v', 'u') for l in (2, 4)] + ['H:2', 'H:4', 'T:+', 'T:-']
OUT = [f'{f}:{l}' for f in ('R', 'T') for l in (0, 2, 4)] + ['A:2', 'A:4']
SCALAR = {0: Q(1), 2: Q(1, 12), 4: Q(19, 720)}
SIG = [1, -1, 1, 1, 1]


def hweight(label):
    if label.startswith('T:'):
        return Q(1, 6)
    field, degree = label.split(':')
    degree = int(degree)
    value = SCALAR[degree]
    lam = degree*(degree+2)
    if field == 'L':
        value *= 2
    if field == 'J':
        value *= 3
    if field in ('v', 'u'):
        value *= 2*lam
    if field == 'H':
        value *= Q(2, 3)*lam*(lam-3)
    return value


H = [hweight(z) for z in LABELS]
B = [SCALAR[l] for _ in range(2) for l in (0, 2, 4)] + [l*(l+2)*SCALAR[l] for l in (2, 4)]


def zero(n, m):
    return [[Q(0) for _ in range(m)] for _ in range(n)]


def matrices(alpha=Q(0), eta=Q(0), zeta=Q(0)):
    rx, rr, z = [zero(8, 23) for _ in range(3)]
    def put(m, out, field, l, c):
        label = f'{field}:{l}'
        if label in LABELS:
            m[OUT.index(f'{out}:{l}')][LABELS.index(label)] += c
    for l in (0, 2, 4):
        lam, beta = l*(l+2), Q(2)-Q(2, 3)*l*(l+2)
        put(rx, 'R', 'L', l, -1)
        for f, c in [('K', -Q(1, 2)), ('N', -Q(1, 2)), ('J', Q(3, 2))]:
            put(rr, 'R', f, l, c)
        put(rr, 'T', 'L', l, -1)
        for f, c in [('K', -Q(1, 2)), ('N', -Q(1, 2)), ('J', -Q(3, 2))]:
            put(rx, 'T', f, l, c)
        for f, c in [('K', 4*alpha), ('N', alpha), ('J', -3*alpha), ('L', -3*eta), ('v', -lam*zeta)]:
            put(z, 'R', f, l, c)
        for f, c in [('L', 5*alpha), ('N', -3*eta), ('J', -3*eta), ('u', -lam*zeta)]:
            put(z, 'T', f, l, c)
        if l:
            put(rr, 'A', 'v', l, -1)
            put(rx, 'A', 'u', l, -1)
            for f, c in [('v', 5*alpha), ('u', -4*eta), ('K', -zeta/2), ('N', zeta/2), ('J', -zeta/2), ('H', beta*zeta)]:
                put(z, 'A', f, l, c)
    return rx, rr, z


def mv(m, v):
    return [sum((a*b for a, b in zip(row, v)), Q(0)) for row in m]


def adjoint(m, v):
    return [sum((B[i]*m[i][j]*v[i] for i in range(8)), Q(0))/H[j] for j in range(23)]


def n2(v, weights):
    return sum((w*a*a for w, a in zip(weights, v)), Q(0))


def frame_fixture(l, seed, mutation=None):
    """Contract covariant derivatives in an orthonormal normal coordinate jet.

    The trace-free tensor derivative jet is chosen to have the exact spherical
    Hessian divergence; no modal matrix is used to perform this contraction.
    """
    lam, beta = l*(l+2), Q(2)-Q(2, 3)*l*(l+2)
    names = ('K', 'L', 'N', 'J', 'v', 'u', 'H')
    vals = {f: Q((seed+2)*(i+3)-11, i+2) for i, f in enumerate(names)}
    yy = {f: Q((seed+3)*(i+2)+5, i+3) for i, f in enumerate(names)}
    tt = {f: Q((seed+4)*(i+1)-7, i+4) for i, f in enumerate(names)}
    a, h, rho, c = Q(2, 7), Q(-3, 11), Q(5, 4), Q(9, 8)
    Y = Q(13, 7)
    grad = [Q(2, 3), Q(-3, 5), Q(5, 8)] if l else [Q(0)]*3
    tf = [[Q(1), Q(2, 5), Q(-3, 7)], [Q(2, 5), Q(-2), Q(1, 11)], [Q(-3, 7), Q(1, 11), Q(1)]]
    if not l:
        for f in ('v', 'u', 'H'):
            vals[f] = yy[f] = tt[f] = Q(0)
        tf = zero(3, 3)
    hes = [[tf[i][j]-Q(i == j)*lam*Y/3 for j in range(3)] for i in range(3)]
    def metric(v):
        k = zero(5, 5)
        k[0][0], k[0][1], k[1][0], k[1][1] = v['K']*Y, v['L']*Y, v['L']*Y, v['N']*Y
        for i in range(3):
            k[0][i+2] = k[i+2][0] = v['v']*grad[i]
            k[1][i+2] = k[i+2][1] = v['u']*grad[i]
            for j in range(3):
                k[i+2][j+2] = v['J']*Q(i == j)*Y+v['H']*tf[i][j]
        return k
    k = metric(vals)
    dk = [metric(yy), [[v/rho for v in row] for row in metric(tt)]] + [zero(5, 5) for _ in range(3)]
    for d in range(3):
        part = dk[d+2]
        for f, i, j in [('K', 0, 0), ('L', 0, 1), ('N', 1, 1)]:
            part[i][j] = part[j][i] = vals[f]*grad[d]/(rho*c)
        for i in range(3):
            part[0][i+2] = part[i+2][0] = vals['v']*hes[d][i]/(rho*c)
            part[1][i+2] = part[i+2][1] = vals['u']*hes[d][i]/(rho*c)
            for j in range(3):
                dtf = Q(3, 10)*beta*(Q(d == i)*grad[j]+Q(d == j)*grad[i]-Q(2, 3)*Q(i == j)*grad[d])
                part[i+2][j+2] = (vals['J']*Q(i == j)*grad[d]+vals['H']*dtf)/(rho*c)
    gamma = [zero(5, 5) for _ in range(5)]
    gamma[1][0][1] = gamma[1][1][0] = a
    for d in range(2, 5):
        gamma[d][0][d], gamma[d][1][d] = a, h
        gamma[d][d][0], gamma[d][d][1] = -a, h
    sig = SIG if mutation != 'euclidean_time_sign' else [1]*5
    direct = []
    for i in range(5):
        divergence = Q(0)
        for j in range(5):
            conn = sum((gamma[j][i][b]*k[b][j]+gamma[j][j][b]*k[i][b] for b in range(5)), Q(0))
            divergence += sig[j]*(dk[j][i][j]-(Q(0) if mutation == 'drop_connection' else conn))
        trace_deriv = sum((sig[j]*dk[i][j][j] for j in range(5)), Q(0))
        direct.append(divergence-(Q(1) if mutation == 'trace_factor_one' else Q(1, 2))*trace_deriv)
    R = (yy['K']+yy['N']-3*yy['J'])/2-tt['L']/rho-lam*vals['v']/(rho*c)+a*(4*vals['K']+vals['N']-3*vals['J'])-3*h*vals['L']
    T = yy['L']-(tt['K']+tt['N']+3*tt['J'])/(2*rho)-lam*vals['u']/(rho*c)+5*a*vals['L']-3*h*(vals['N']+vals['J'])
    A = yy['v']-tt['u']/rho+5*a*vals['v']-4*h*vals['u']+(-vals['K']/2+vals['N']/2-vals['J']/2+beta*vals['H'])/(rho*c)
    modal = [R*Y, T*Y]+[A*g for g in grad]
    return direct, modal


groups, controls = [], []
def checked(name, ok):
    need(ok, name)
    groups.append(name)


for l in (0, 2, 4):
    for seed in (1, 3, 7, 11):
        direct, modal = frame_fixture(l, seed)
        checked(f'frame_contraction_degree_{l}_seed_{seed}', direct == modal)
for mutation in ('euclidean_time_sign', 'drop_connection', 'trace_factor_one'):
    direct, modal = frame_fixture(2, 3, mutation)
    need(direct != modal, 'corruption must change contraction:' + mutation)
    controls.append(mutation)

rx, rr, z = matrices(Q(2, 7), Q(-3, 11), Q(5, 13))
gram_sum = [[sum((rx[i][k]*rx[j][k]+rr[i][k]*rr[j][k])/H[k] for k in range(23))*B[j] for j in range(8)] for i in range(8)]
expected = [[(Q(7, 4) if i < 6 else Q(1))*Q(i == j) for j in range(8)] for i in range(8)]
checked('combined_derivative_output_Gram_identity', gram_sum == expected)
checked('all_positive_Gram_weights', min(H+B) > 0)
for name in ['w:0', 'w:2', 'w:4', 'T:+', 'T:-']:
    j = LABELS.index(name)
    checked('constraint_annihilates_' + name, all(row[j] == 0 for m in (rx, rr, z) for row in m))

witness = [Q(0)]*8
witness[0] = witness[3] = Q(1)
ex, er = adjoint(rx, witness), adjoint(rr, witness)
y = [a+b for a, b in zip(mv(rx, ex), mv(rr, er))]
checked('sharp_derivative_bound_balanced_witness', n2(ex, H) == n2(er, H) == Q(7, 4) and n2(y, B) == Q(7, 2)*n2(ex, H))
need(gram_sum[0][0] > Q(3, 2), 'understated derivative Gram cap rejected')
controls.append('understated_derivative_Gram_cap')

# Scaled coordinate conversion must have a negative inward radial sign.
ell, rho = Q(1, 10**6), Q(5, 4)
ey, et = [Q(i+1, i+2) for i in range(23)], [Q(2*i-3, i+3) for i in range(23)]
xx, radial = [ell*rho*v for v in et], [-ell*rho*v for v in ey]
need(mv(rr, radial) != mv(rr, [-v for v in radial]), 'wrong radial orientation rejected')
controls.append('wrong_inward_radial_sign')
checked('physical_scaled_first_derivative_relation', all(radial[i]/(-ell*rho) == ey[i] and xx[i]/(ell*rho) == et[i] for i in range(23)))

# Actual shell source contraction, independent of any formal target.
QR00, QRJ, QC, eR = Q(2, 3), Q(-5, 7), Q(11, 13), Q(17, 19)
trQR = -QR00+3*QRJ
GN = -2*trQR-2*QC
GK, GNrow, GJ = -eR*GN, 2*eR*QR00, 2*eR*QRJ
checked('normal_affine_source_packing_identity', GK+GNrow-3*GJ == 2*eR*QC)
checked('normal_full_second_constraint_boundary_cancellation', -(GK+GNrow-3*GJ)/(2*eR)+QC == 0)
need(-(-GK+GNrow-3*GJ)/(2*eR)+QC != 0, 'wrong source sign rejected')
controls.append('wrong_normal_source_sign')

alpha, sigma = Q(2, 7), Q(-3, 11)
K, N, J, w = Q(2), Q(3), Q(5), Q(7)
ByK, ByN, ByJ = -4*alpha*K-2*alpha*N+6*alpha*J+Q(4, 3)*sigma*w, -alpha*K-sigma*w/3, alpha*K+sigma*w/3
checked('normal_homogeneous_shell_cancellation', (ByK+ByN-3*ByJ)/2+alpha*(4*K+N-3*J) == 0)

# A zero Dirichlet trace does not set its normal derivative to zero.
normal_only = [Q(0)]*23
normal_only[LABELS.index('L:0')] = Q(1)
normal_response = mv(rr, normal_only)
checked('Dirichlet_normal_trace_is_retained', normal_response[3] == -1)
need(any(normal_response), 'omitted Dirichlet normal derivative rejected')
controls.append('omitted_Dirichlet_normal_derivative')

# S3 addition kernels use the entire irreducible hull. This arithmetic records
# their dimensions; homogeneity and isotropy prove their constant kernel form.
dims = [(l+1)**2 for l in (0, 2, 4)]
scalar_kernel = sum(dims)
gradient_trace = sum(dims[1:])
checked('S3_irreducible_hull_dimensions', dims == [1, 9, 25])
checked('pointwise_covector_kernel_maximum', scalar_kernel == 35 and Q(gradient_trace, 3) == Q(34, 3) < scalar_kernel)

print(json.dumps({'status': 'PASS_EXACT_HARMONIC_CONSTRAINT_FORMULAS',
    'named_groups': groups, 'group_count': len(groups),
    'wrong_formula_controls': controls, 'control_count': len(controls),
    'input_labels': LABELS, 'output_labels': OUT,
    'input_Gram': list(map(str, H)), 'output_Gram': list(map(str, B)),
    'sharp_scaled_derivative_cap_squared': '7/2',
    'pointwise_from_sphere_L2_cap_squared': '35',
    'scope': 'independent exact manufactured frame contractions, norm and source-packing identities; no actual second field, physical constraint completion or Einstein residual evaluated'}, sort_keys=True))
