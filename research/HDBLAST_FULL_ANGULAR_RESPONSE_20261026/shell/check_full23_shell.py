#!/usr/bin/env python3
"""Authenticate frozen inputs, then check new exact full23 shell algebra.

This evaluates no old campaign or actual first-field waveform. All scientific
code imports follow complete manifest authentication. Output paths must be new.
"""
import argparse
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction as Q
from pathlib import Path, PurePosixPath

sys.dont_write_bytecode = True
LABELS = {
    'current': 'HDBLAST_COUPLED_MEAN_PILOT_20261025',
    'I': 'HDBLAST_INPUT_PRECISION_20261023',
    'P': 'HDBLAST_COMMON_GAUGE_FORCING_20261016',
    'N': 'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016',
    'G': 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
    'C': 'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023',
    'T': 'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023',
    'latest': 'HDBLAST_SHELL_FEEDBACK_PRECISION_20261025',
}
PINS = {
    'current': 'b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de',
    'I': 'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61',
    'P': '8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6',
    'N': '928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44',
    'G': '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2',
    'C': '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76',
    'T': 'cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72',
    'latest': '1a5ee8aa22953f949da000ba03da3b5a0872e4467a0fbaec9f2d3fe7fba11ec4',
}
HELPER_PIN = 'c01a8f5f95737dcacebdc26379ae0ee5863af9c1977890b3901ad6113103e619'


def need(ok, label):
    if not ok:
        raise ValueError(label)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def authenticate(root, pin):
    need(not root.is_symlink(), 'ROOT_SYMLINK')
    manifest = root / 'MANIFEST.json'
    need(manifest.is_file() and not manifest.is_symlink() and sha(manifest) == pin,
         'MANIFEST_PIN:' + str(root))
    expected = {'MANIFEST.json'}
    for row in json.loads(manifest.read_text())['files']:
        name = row['path']
        rel = PurePosixPath(name)
        need(name == rel.as_posix() and not rel.is_absolute() and '..' not in rel.parts
             and '.' not in rel.parts and '\\' not in name and name not in expected,
             'MEMBER_PATH')
        expected.add(name)
        member = root / name
        need(member.is_file() and not member.is_symlink(), 'MEMBER_TYPE')
        need(member.stat().st_size == row['bytes'] and sha(member) == row['sha256'],
             'MEMBER_BYTES:' + name)
    actual = {'MANIFEST.json'}
    for member in root.rglob('*'):
        need(not member.is_symlink(), 'TREE_SYMLINK')
        if member.is_file():
            actual.add(member.relative_to(root).as_posix())
    need(actual == expected, 'COMPLETE_CLOSURE')
    return {'manifest_sha256': pin, 'files_including_manifest': len(expected)}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rank(rows):
    a = [list(map(Q, row)) for row in rows]
    if not a:
        return 0
    k = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(k, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[k], a[pivot] = a[pivot], a[k]
        z = a[k][j]
        a[k] = [v/z for v in a[k]]
        for i in range(len(a)):
            if i != k and a[i][j]:
                z = a[i][j]
                a[i] = [v-z*w for v, w in zip(a[i], a[k])]
        k += 1
        if k == len(a):
            break
    return k


def encode(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inject-error', choices=('boundary-sign', 'missing-TT',
                                                  'normal-trace', 'wrong-Dirichlet'))
    args = parser.parse_args()
    roots = {key: args.repo_root.resolve()/'research'/name for key, name in LABELS.items()}
    need(not args.output.exists(), 'FRESH_OUTPUT_REQUIRED')
    need(all(root != args.output.resolve() and root not in args.output.resolve().parents
             for root in roots.values()), 'EXTERNAL_OUTPUT_REQUIRED')
    before = {key: authenticate(root, PINS[key]) for key, root in roots.items()}
    helper_path = roots['latest']/'balanced/exact_boundary.py'
    need(sha(helper_path) == HELPER_PIN, 'HELPER_PIN')
    # No project import has occurred before this complete authentication.
    helper = load('shell_authenticated_helper', helper_path)
    inherited = helper.authenticate({key: root for key, root in roots.items() if key != 'latest'})
    need(inherited == {key: row for key, row in before.items() if key != 'latest'},
         'HELPER_AUTHENTICATION_AGREEMENT')
    op = load('shell_authenticated_operator', roots['G']/'projection/inputs/linear_operator.py')
    source = load('shell_authenticated_sources', roots['G']/'projection/inputs/source_adapters.py')
    ring = load('shell_authenticated_ring', roots['N']/'angular/source-audit/source_laurent.py')
    syms = dict(zip(ring.NAMES, ring.variables()))
    a, p, s1, s2, s3, rho = (syms[k] for k in ('a', 'p', 'sigma1', 'sigma2', 'sigma3', 'rho'))
    h, X, V, L, d, v = (syms[k] for k in ('A', 'B', 'C', 'D', 'E', 'v'))
    ell = Q(1, 10**6)
    labels, idx = op.LABELS, op.IDX
    checks = []

    def equal(actual, expected, name):
        need(actual == expected, name)
        checks.append(name)

    shell = source.shell_sources(a=a, p=p, sigmap=s1, sigmapp=s2, sigmappp=s3,
                                 h=h, X=X, V=V, L=L, d=d, v=v)
    G = op.boundary_forcing(ell, rho, shell)
    if args.inject_error == 'boundary-sign':
        G[idx['w:2']] = -G[idx['w:2']]
    if args.inject_error == 'missing-TT':
        G[idx['T:-']] = 0
    if args.inject_error == 'normal-trace':
        G[idx['K:0']] += ell*rho*a*V*V

    # Independent all-point affine boundary-Hessian formulas, no angular samples.
    y2 = (Q(1, 6), Q(1, 2), Q(1))
    grad2 = (Q(4, 3), Q(2), Q(-4))
    hess2 = (Q(8), Q(4), Q(24))
    yhess = (Q(-4, 9), Q(-4, 3), Q(-8, 3), Q(-1, 4), Q(1, 6), Q(2, 3), Q(2, 3))
    tau = 2*X+6*L-16*d
    alpha = a*h*h/4-s1*h*v/6-s2*v*v/6
    beta = a*h+s1*v/3
    pi = p*h*h/4+s2*h*v/2+s3*v*v/2
    sq = (4*X*X+12*L*L-64*L*d, -2*V*V, 4*d*d)
    gn = (4*s2*v*v/3+4*s1*h*v-10*a*h*h+8*a*h*tau, 0, 0)
    gn = [gn[j]-4*a*sq[j] for j in range(3)]
    qc = (4*a*h*h-3*a*h*tau-4*s1*h*v/3+s1*v*tau/3, 0, 0)
    qc = [qc[j]+2*a*sq[j] for j in range(3)]

    def project(q):
        return [q[0]*y2[j]+q[1]*grad2[j]+q[2]*hess2[j] for j in range(3)]

    expected = [ring.P(0) for _ in labels]
    for j, degree in enumerate((0, 2, 4)):
        expected[idx[f'K:{degree}']] = -ell*rho*project(gn)[j]
        expected[idx[f'N:{degree}']] = 2*ell*rho*(-alpha+2*beta*X)*y2[j]
        expected[idx[f'w:{degree}']] = ell*rho*pi*y2[j]
    tensor = [(alpha-2*beta*L)*(y2[j] if j < 3 else 0)-2*beta*d*yhess[j]
              for j in range(7)]
    for name, value in zip(('J:0', 'J:2', 'J:4', 'H:2', 'H:4', 'T:+', 'T:-'), tensor):
        expected[idx[name]] = 2*ell*rho*value
    expected[idx['u:2']] = -ell*rho*beta*V/2
    expected[idx['u:4']] = -ell*rho*beta*V
    for name, actual, target in zip(labels, G, expected):
        equal(actual, target, 'FULL_BOUNDARY_SOURCE:'+name)
    traceQR = [-shell['QK00'][j]+3*shell['QKab'][j] for j in range(3)]
    for j, degree in enumerate((0, 2, 4)):
        equal(shell['GN'][j], -2*traceQR[j]-2*project(qc)[j],
              'NORMAL_GAUGE_TRACE_JOIN:'+str(degree))

    mats = op.matrices()
    PD, PN = mats['PD'], mats['PN']
    if args.inject_error == 'wrong-Dirichlet':
        PD[idx['u:2']][idx['u:2']] = 1
    dirichlet = [name for name in labels if PD[idx[name]][idx[name]]]
    equal(dirichlet, ['L:0', 'L:2', 'L:4', 'v:2', 'v:4'], 'EXACT_DIRICHLET_SLOTS')
    M = op.scale(-ell*rho, op.add(op.scale(a, mats['By_a']), op.scale(p, mats['By_p']),
                                 op.scale(s1, mats['By_sigma1']), op.scale(s2, mats['By_sigma2'])))
    By = op.add(op.scale(a, mats['By_a']), op.scale(p, mats['By_p']),
                op.scale(s1, mats['By_sigma1']), op.scale(s2, mats['By_sigma2']))
    # Orthogonal scalar trace/null decomposition. This never divides by M.
    transform = [[Q(0)]*23 for _ in range(18)]
    inverse = [[Q(0)]*18 for _ in range(23)]
    split_names = []
    for j, degree in enumerate((0, 2, 4)):
        srow, krow, wrow, qrow = 3*j, 3*j+1, 3*j+2, 9+j
        transform[srow][idx[f'N:{degree}']] = -Q(1, 2)
        transform[srow][idx[f'J:{degree}']] = Q(3, 2)
        transform[krow][idx[f'K:{degree}']] = 1
        transform[wrow][idx[f'w:{degree}']] = 1
        transform[qrow][idx[f'N:{degree}']] = transform[qrow][idx[f'J:{degree}']] = 1
        inverse[idx[f'N:{degree}']][srow], inverse[idx[f'N:{degree}']][qrow] = -Q(1, 2), Q(3, 4)
        inverse[idx[f'J:{degree}']][srow], inverse[idx[f'J:{degree}']][qrow] = Q(1, 2), Q(1, 4)
        inverse[idx[f'K:{degree}']][krow] = inverse[idx[f'w:{degree}']][wrow] = 1
        split_names += [f's:{degree}', f'K:{degree}', f'w:{degree}']
    split_names += [f'q:{degree}' for degree in (0, 2, 4)]
    for row, name in enumerate(('u:2', 'u:4', 'H:2', 'H:4', 'T:+', 'T:-'), 12):
        transform[row][idx[name]] = inverse[idx[name]][row] = 1
        split_names.append(name)
    equal(op.mm(transform, inverse), op.ident(18), 'ROBIN_TRACE_SPLIT_RIGHT_INVERSE')
    equal(op.mm(inverse, transform), PN, 'ROBIN_TRACE_SPLIT_LEFT_INVERSE')
    active = [[0, 2*a, 2*s1/3], [4*a, -4*a, 4*s1/3], [0, p/2, -s2/2]]
    equal(sum((entry*entry for row in active for entry in row), ring.P(0)),
          36*a*a+Q(20, 9)*s1*s1+(p*p+s2*s2)/4,
          'THREE_BY_THREE_FROBENIUS_BOUND')
    split_By = [[ring.P(0)]*18 for _ in range(18)]
    for block in range(3):
        for i in range(3):
            for j in range(3):
                split_By[3*block+i][3*block+j] = active[i][j]
    equal(By, op.mm(inverse, op.mm(split_By, transform)), 'ROBIN_THREE_REPEATED_BLOCKS_PLUS_NINE_ZERO')
    split_gram = [[sum((op.gram()[i]*inverse[i][j]*inverse[i][k] for i in range(23)), Q(0))
                   for k in range(18)] for j in range(18)]
    equal(all(split_gram[i][j] == 0 for i in range(18) for j in range(18) if i != j),
          True, 'ROBIN_TRACE_SPLIT_HILBERT_ORTHOGONAL')
    null_vectors = []
    left_vectors = []
    for degree in (0, 2, 4):
        right = [Q(0)]*23
        right[idx[f'N:{degree}']], right[idx[f'J:{degree}']] = 3, 1
        null_vectors.append((f'tangent-tracefree:{degree}', right))
        left = [Q(0)]*23
        left[idx[f'N:{degree}']] = left[idx[f'J:{degree}']] = 1
        left_vectors.append((f'N+J:{degree}', left))
    for name in ('u:2', 'u:4', 'H:2', 'H:4', 'T:+', 'T:-'):
        unit = [Q(i == idx[name]) for i in range(23)]
        null_vectors.append((name, unit))
        left_vectors.append((name, unit))
    for name, vector in null_vectors:
        for i in range(23):
            equal(sum((M[i][j]*vector[j] for j in range(23)), ring.P(0)), 0,
                  'ROBIN_RIGHT_NULL:'+name+':'+str(i))
    for name, vector in left_vectors:
        for j in range(23):
            equal(sum((vector[i]*M[i][j] for i in range(23)), ring.P(0)), 0,
                  'ROBIN_LEFT_NULL:'+name+':'+str(j))
    witness_M = op.scale(-ell, op.add(op.scale(Q(2), mats['By_a']),
                                     op.scale(Q(3), mats['By_p']),
                                     op.scale(Q(-6), mats['By_sigma1']),
                                     op.scale(Q(5), mats['By_sigma2'])))
    equal(rank(witness_M), 9, 'GENERIC_ROBIN_RANK_WITNESS')
    monomials = sorted({term for item in G for term in ring.P(item).terms})
    source_matrix = [[ring.P(item).terms.get(term, Q(0)) for term in monomials] for item in G]
    source_rank = rank(source_matrix)
    equal(source_rank, 8, 'UNIVERSAL_BOUNDARY_SOURCE_ANGULAR_RANK')
    # A reusable eight-profile source factorization, with its exact H Gram.
    basis = [[Q(0)]*8 for _ in labels]
    for j, degree in enumerate((0, 2, 4)):
        for col, angular in enumerate((y2, grad2, hess2)):
            basis[idx[f'K:{degree}']][col] = angular[j]
        basis[idx[f'N:{degree}']][3] = y2[j]
        basis[idx[f'w:{degree}']][4] = y2[j]
        basis[idx[f'J:{degree}']][6] = y2[j]
    basis[idx['u:2']][5], basis[idx['u:4']][5] = Q(1, 4), Q(1, 2)
    for j, name in enumerate(('J:0', 'J:2', 'J:4', 'H:2', 'H:4', 'T:+', 'T:-')):
        basis[idx[name]][7] = yhess[j]
    profiles = [-ell*rho*g for g in gn]
    profiles += [2*ell*rho*(-alpha+2*beta*X), ell*rho*pi,
                 -2*ell*rho*beta*V, 2*ell*rho*(alpha-2*beta*L), -4*ell*rho*beta*d]
    for i, name in enumerate(labels):
        equal(G[i], sum((basis[i][j]*profiles[j] for j in range(8)), ring.P(0)),
              'EIGHT_PROFILE_FACTORIZATION:'+name)
    equal(rank(basis), 8, 'EIGHT_PROFILE_BASIS_RANK')
    source_gram = [[sum((op.gram()[i]*basis[i][j]*basis[i][k] for i in range(23)), Q(0))
                    for k in range(8)] for j in range(8)]
    equal(G[idx['T:+']], G[idx['T:-']], 'TT_BOUNDARY_FORCING_EQUAL')
    equal(G[idx['T:+']], -Q(8, 3)*ell*rho*beta*d, 'TT_BOUNDARY_FORCING_EXPLICIT')
    equal(G[idx['u:4']], 2*G[idx['u:2']], 'VECTOR_BOUNDARY_SOURCE_RATIO')
    equal(G[idx['H:4']], -Q(2, 3)*G[idx['H:2']], 'HESSIAN_BOUNDARY_SOURCE_RATIO')

    # Exact value/first derivative boundary-lift check, with no g derivative.
    lift_checks = []
    for x, r in ((Q(1), Q(0)), (Q(3, 2), Q(1, 2)), (Q(1), Q(2)), (Q(1), Q(1))):
        t = max(x-r, Q(0))
        value, dx, dr = -t**3/3, -t*t, t*t
        need(dx == -dr, 'BOUNDARY_LIFT_CHARACTERISTIC_DERIVATIVES')
        lift_checks.append({'x': x, 'r': r, 'Lg_for_g_s_squared': value, 'dx': dx, 'dr': dr})

    after = {key: authenticate(root, PINS[key]) for key, root in roots.items()}
    need(after == before, 'FROZEN_CLOSURES_CHANGED')
    result = {
        'schema': 'HDBLAST_FULL23_SHELL_COMPATIBILITY_AUDIT_V1', 'status': 'PASS',
        'authentication_before_and_after': before, 'helper_sha256': HELPER_PIN,
        'complete_closures_authenticated_before_all_project_imports': True,
        'exact_identity_checks': len(checks), 'check_names': checks,
        'dirichlet_slots': dirichlet, 'robin_dimension': 18,
        'generic_robin_rank': 9, 'universal_robin_right_nullity_at_least': 9,
        'robin_split_names': split_names,
        'robin_split_definitions': {'s': '(-N+3*J)/2', 'q': 'N+J',
                                    'N': '-s/2+3*q/4', 'J': 's/2+q/4'},
        'robin_scalar_By_block_on_s_K_w': [['0', '2*a', '2*sigma1/3'],
                                          ['4*a', '-4*a', '4*sigma1/3'],
                                          ['0', 'p/2', '-sigma2/2']],
        'robin_split_gram_diagonal': [split_gram[j][j] for j in range(18)],
        'robin_split_identity': 'M is similar on PN to three repeated -ell*rho*By3 blocks and nine exact zeros; ||M||_H is exactly the Euclidean operator norm of -ell*rho*By3.',
        'robin_squared_norm_upper': '(ell*rho)^2*(36*a^2+20*sigma1^2/9+p^2/4+sigma2^2/4); with sigma1=-2*p: (ell*rho)^2*(36*a^2+329*p^2/36+sigma2^2/4).',
        'robin_right_null_vectors': dict(null_vectors), 'robin_left_null_vectors': dict(left_vectors),
        'universal_boundary_source_rank': source_rank,
        'source_monomials': len(monomials),
        'source_basis_B_rows_in_LABELS_order': basis, 'labels': labels,
        'source_profile_names': ['-ell*rho*gn_Y2', '-ell*rho*gn_grad2', '-ell*rho*gn_hess2',
                                 '2*ell*rho*(-alpha+2*beta*X)', 'ell*rho*pi',
                                 '-2*ell*rho*beta*V', '2*ell*rho*(alpha-2*beta*L)',
                                 '-4*ell*rho*beta*d'],
        'source_profile_definitions': {'alpha': 'a*h*h/4-sigma1*h*v/6-sigma2*v*v/6',
                                       'beta': 'a*h+sigma1*v/3',
                                       'pi': 'p*h*h/4+sigma2*h*v/2+sigma3*v*v/2',
                                       'tau': '2*X+6*L-16*d',
                                       'gn_Y2': '4*sigma2*v*v/3+4*sigma1*h*v-10*a*h*h+8*a*h*tau-4*a*(4*X*X+12*L*L-64*L*d)',
                                       'gn_grad2': '8*a*V*V', 'gn_hess2': '-16*a*d*d'},
        'source_gram_Btranspose_H_B': source_gram,
        'source_norm_identity': 'For G=B*c, ||G||_H^2=c^T*(B^T H B)*c, exactly.',
        'source_formula_symbols': {'h': 'normal first metric amplitude', 'X': 'minus half first 00 amplitude',
                                   'V': 'first 0a gradient amplitude', 'L': 'half angular Y metric coefficient',
                                   'd': 'half angular HessY coefficient', 'v': 'first scalar amplitude',
                                   'beta': 'a*h+sigma1*v/3', 'ell': ell},
        'exact_source_relations': {'Tplus=Tminus': '-8*ell*rho*beta*d/3',
                                   'u4': '2*u2', 'H4': '-2*H2/3'},
        'boundary_lift_identity': 'Lg=-integral_0^max(x-r,0) g(s)ds; (Lg)_x=-g(x-r); (Lg)_r=+g(x-r)',
        'boundary_lift_checks': lift_checks,
        'regularity_requirement': 'Continuous G, G(0)=0 and continuous compatible initial data suffice for the C1 lift; no G_x or G_r bound is required.',
        'physical_constraint_scope': 'Closure preserves the angular subspace; physical constraints require the complete natural/reduced Hessian identities, compatible initial data and subsidiary uniqueness.',
        'limits': ['No old campaign was rerun and no frozen file was modified.',
                   'No actual first-field or second-response numerical waveform was evaluated.',
                   'The rank-8 result concerns the shell affine source, not an invariant bulk/response subspace.',
                   'The universal rank-9 witness establishes generic rank only; special background parameters may lower rank.',
                   'M is not inverted; its universal nine-dimensional Robin nullspace is retained.',
                   'No new numerical C2/higher-source-derivative bound or nonlinear-family claim.'],
    }
    with args.output.open('x') as handle:
        json.dump(encode(result), handle, sort_keys=True, indent=2)
        handle.write('\n')
    print(json.dumps({'status': 'PASS', 'checks': len(checks), 'source_rank': source_rank,
                      'robin_generic_rank': 9, 'output_sha256': sha(args.output)}, sort_keys=True))


if __name__ == '__main__':
    main()
