#!/usr/bin/env python3
"""Authenticate and bind the actual first master to natural field residuals.

Exact Laurent coefficients, arbitrary longitudinal jets, independent linear
connection/Ricci contraction and direct scalar variation. Standard library only.
"""
import argparse
import ast
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
need_helper_pin = '99231502fd782ebb5d4c4a5a684370dd61ad6b2b8dc1afa5806dc2bacaa5d829'
if hashlib.sha256((HERE / 'natural_algebra.py').read_bytes()).hexdigest() != need_helper_pin:
    raise ValueError('natural algebra helper authentication')
spec = importlib.util.spec_from_file_location('natural_algebra', HERE / 'natural_algebra.py')
algebra = importlib.util.module_from_spec(spec)
spec.loader.exec_module(algebra)
P, derivation, master_derivative = algebra.P, algebra.derivation, algebra.master_derivative

PINS = {
    'HDBLAST_REGULAR_CONE_EVOLUTION_20261013': 'a873132efe11755217ddb865fd0156041b8c1ad21ae1c45310183da4fb6140ae',
    'HDBLAST_REGULAR_SHELL_SPECTRA_20261013': 'd1d47c66cef1ed79979c22ae7145dc66d35fa35be69dba900f9912ac01f31e2f',
    'HDBLAST_INTERIOR_MASTER_JETS_20261015': 'b043a13a228d11296fbbde292aa15f7a024deab9ff915bdce452c52c73745bbb',
    'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016': '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2',
    'HDBLAST_BULK_SHELL_CONSTRAINT_BINDING_20261030': '7eae74dedfef544199d410665deb5ebf013d4cf484745428aba37dac7ba0671e',
}
CONTROLS = (
    'ricci_upper_connection', 'ricci_time_signature', 'omit_potential_metric',
    'wrong_potential_scalar_factor', 'wrong_chi_sign', 'omit_scalar_gauge_term',
    'omit_hamiltonian', 'master_plus4', 'wrong_angular_eigenvalue',
    'omit_temporal_theta_derivative', 'lose_radial_constant',
    'wrong_potential_p2', 'wrong_Vpp_sign', 'reverse_radial_sign',
    'freeze_radial_rho', 'wrong_scalar_master_derivative',
    'claim_point_scalar_equation', 'wrong_source_scalar_shift',
    'source_physical_yy_without_rho_term', 'source_temporal_theta_wrong_factor',
    'claim_comparator_is_actual', 'scalar_derivative_missing_ell',
)
EXPECTED_FAILURES = {
    'ricci_upper_connection': 'independent_connection_Ricci_(0, 0)',
    'ricci_time_signature': 'independent_connection_Ricci_(0, 0)',
    'omit_potential_metric': 'natural_Einstein_factor_(0, 0)',
    'wrong_potential_scalar_factor': 'natural_Einstein_factor_(0, 0)',
    'wrong_chi_sign': 'mixed_reconstruction',
    'omit_scalar_gauge_term': 'direct_scalar_variation_formula',
    'omit_hamiltonian': 'natural_Einstein_factor_(0, 0)',
    'master_plus4': 'independent_temporal_constant',
    'wrong_angular_eigenvalue': 'independent_temporal_constant',
    'omit_temporal_theta_derivative': 'independent_temporal_constant',
    'lose_radial_constant': 'actual_separated_master_factor',
    'wrong_potential_p2': 'actual_separated_master_factor',
    'wrong_Vpp_sign': 'actual_separated_master_factor',
    'reverse_radial_sign': 'actual_scalar_reconstruction',
    'freeze_radial_rho': 'variable_rho_second_radial_chain',
    'wrong_scalar_master_derivative': 'direct_scalar_master_factor',
    'claim_point_scalar_equation': 'invalid_point_only_scalar_claim',
    'wrong_source_scalar_shift': 'actual_scalar_Lie_pullback_plus',
    'source_physical_yy_without_rho_term': 'actual_physical_jet_adapter_yy',
    'source_temporal_theta_wrong_factor': 'actual_source_vartheta_scaled_log_derivative',
    'claim_comparator_is_actual': 'invalid_comparator_for_actual_equation',
    'scalar_derivative_missing_ell': 'actual_scalar_residual_and_derivative_factor',
}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def authenticate(repo):
    result = {}
    for name, pin in PINS.items():
        root = repo / 'research' / name
        need(not root.is_symlink(), 'symlink capsule root')
        raw = (root / 'MANIFEST.json').read_bytes()
        need(digest(raw) == pin, 'manifest authentication: ' + name)
        manifest = json.loads(raw)
        members = [row['path'] for row in manifest['files']]
        need(len(set(members)) == len(members) and 'MANIFEST.json' not in members,
             'nonunique/self-referential source manifest')
        inventory = set()
        for candidate in root.rglob('*'):
            need(not candidate.is_symlink(), 'symlink in source closure')
            if candidate.is_file():
                inventory.add(str(candidate.relative_to(root)))
        need(inventory == set(members) | {'MANIFEST.json'}, 'source closure inventory: ' + name)
        total = 0
        for row in manifest['files']:
            rel = Path(row['path'])
            need(not rel.is_absolute() and '..' not in rel.parts, 'unsafe manifest member')
            data = (root / rel).read_bytes()
            need(len(data) == row['bytes'] and digest(data) == row['sha256'],
                 'member authentication: ' + name + '/' + row['path'])
            total += len(data)
        result[name] = {'manifest_sha256': pin, 'members': len(manifest['files']), 'bytes': total}
    return result


def vars_named(names):
    return [P.var(s) for s in names.split()]


rho, a, p, Vp, Vpp, Vppp, V, psi = vars_named('rho a p Vp Vpp Vppp V psi')
ay = -rho ** -2 - p ** 2 / 3
py = Vp - 4 * a * p
Vconstraint = 6 / rho ** 2 + p ** 2 / 2 - 6 * a ** 2


def d(poly, index):
    rules = {}
    if index == 0:
        rules = {'rho': a * rho, 'a': ay, 'p': py, 'Vp': Vpp * p,
                 'Vpp': Vppp * p, 'V': Vp * p}
    for name in P.cv(poly).names():
        if name.startswith('psi'):
            rules[name] = P.var('psi' + ''.join(sorted(name[3:] + str(index))))
    return derivation(poly, rules)


def box_slice(poly, signature=(-1, 1, 1, 1)):
    return sum((signature[i - 1] * d(d(poly, i), i) for i in range(1, 5)), P())


class Audit:
    def __init__(self, injection):
        self.injection = injection
        self.rows = []

    def equal(self, name, left, right):
        remainder = P.cv(left) - right
        need(not remainder, 'mathematical mismatch ' + name + ': ' + str(remainder.encode()[:3]))
        self.rows.append({'name': name, 'zero_collected_coefficients': True})

    def nonzero(self, name, expression):
        need(bool(P.cv(expression)), 'expected nonzero witness: ' + name)
        self.rows.append({'name': name, 'nonzero_witness': P.cv(expression).encode()})


def derive_longitudinal(audit):
    """Use slice Riemann normal coordinates at an arbitrary point.

    gamma=(-,+,+,+), d gamma=0 at that point. P is a tensor built from
    gamma, dPsi and warp coefficients, so its first derivative needs no
    derivative of the slice Christoffel. The calculation is tensorial at every
    point and leaves all first/second scalar-function jets independent.
    """
    inj = audit.injection
    eta = [-1, 1, 1, 1]
    gdiag = [P(1)] + [rho ** 2 * s for s in eta]
    inverse = [P(1)] + [s / rho ** 2 for s in eta]
    hdiag = [-4 * psi] + [2 * rho ** 2 * psi * s for s in eta]

    def h(i, j):
        return hdiag[i] if i == j else P()

    def gamma(up, i, j):
        if up == 0 and i == j and i > 0:
            return -a * gdiag[i]
        if up > 0 and ((i == 0 and j == up) or (j == 0 and i == up)):
            return a
        return P()

    def nabla_h(k, i, j):
        return d(h(i, j), k) - sum((gamma(v, k, i) * h(v, j) +
                gamma(v, k, j) * h(i, v) for v in range(5)), P())

    conn = {(u, i, j): inverse[u] * (nabla_h(i, j, u) + nabla_h(j, i, u) -
             nabla_h(u, i, j)) / 2
            for u in range(5) for i in range(5) for j in range(5)}
    # Reduce the independently generated P tensor to a universal expression
    # in gamma, grad(Psi), and warp coefficients BEFORE differentiating it.
    # At an arbitrary slice RNC point d gamma=0; the following expression
    # contains no slice Christoffel to differentiate. Its covariant derivative
    # below retains each of P's three tensor-slot connection terms.
    universal = {}
    for u in range(5):
        for i in range(5):
            for j in range(5):
                if u == 0:
                    if i == j == 0:
                        value = -2 * d(psi, 0)
                    elif i == 0 or j == 0:
                        value = -2 * d(psi, j if i == 0 else i)
                    else:
                        value = -rho ** 2 * (d(psi, 0) + 6 * a * psi) * eta[i - 1] if i == j else P()
                elif i == j == 0:
                    value = 2 * d(psi, u) * eta[u - 1] / rho ** 2
                elif i == 0 or j == 0:
                    value = d(psi, 0) if u == (j if i == 0 else i) else P()
                else:
                    value = (d(psi, i) if u == j else P()) + (d(psi, j) if u == i else P())
                    if i == j:
                        value -= eta[i - 1] * eta[u - 1] * d(psi, u)
                universal[u, i, j] = value
                need(conn[u, i, j] == value, 'mathematical mismatch universal_first_connection')
    audit.rows.append({'name': 'universal_first_connection_tensor', 'zero_collected_coefficients': True,
                       'ordered_tensor_entries': 125, 'independent_lower_slots': True})

    def nabla_conn(k, u, i, j):
        out = d(conn[u, i, j], k)
        for v in range(5):
            if inj != 'ricci_upper_connection':
                out += gamma(u, k, v) * conn[v, i, j]
            out -= gamma(v, k, i) * conn[u, v, j]
            out -= gamma(v, k, j) * conn[u, i, v]
        return out

    ric = {(i, j): sum((nabla_conn(u, u, i, j) - nabla_conn(j, u, i, u)
                         for u in range(5)), P())
           for i in range(5) for j in range(i, 5)}
    bs = box_slice(psi, (1, 1, 1, 1) if inj == 'ricci_time_signature' else (-1, 1, 1, 1))
    expected = {(0, 0): -4 * d(d(psi, 0), 0) - 16 * a * d(psi, 0) + 2 * bs / rho ** 2}
    for i in range(1, 5):
        expected[0, i] = -3 * d(d(psi, 0) + 2 * a * psi, i)
        for j in range(i, 5):
            expected[i, j] = eta[i - 1] * (-bs - rho ** 2 *
                 (d(d(psi, 0), 0) + 10 * a * d(psi, 0) + 6 * (ay + 4 * a ** 2) * psi)) if i == j else P()
    for key in sorted(ric):
        audit.equal('independent_connection_Ricci_' + str(key), ric[key], expected[key])
    # C^a = g^{bc} P^a_bc, obtained from the same independent connection.
    C = [sum((inverse[i] * conn[u, i, i] for i in range(5)), P()) for u in range(5)]
    audit.equal('nonharmonic_C_y', C[0], -6 * d(psi, 0) - 24 * a * psi)
    for i in range(1, 5):
        audit.equal('nonharmonic_C_slice_' + str(i), C[i], P())
    chi = (-1 if inj == 'wrong_chi_sign' else 1) * -3 * (d(psi, 0) + 2 * a * psi) / p
    Pmaster = rho ** 2 * psi
    M = d(d(Pmaster, 0), 0) - (2 * a + 2 * py / p) * d(Pmaster, 0) + \
        (box_slice(Pmaster) + 4 * Pmaster) / rho ** 2 - F(2, 3) * p ** 2 * Pmaster
    audit.equal('mixed_reconstruction', p * chi, -3 * (d(psi, 0) + 2 * a * psi))
    for (i, j), Rij in ric.items():
        metric_term = 0 if inj == 'omit_potential_metric' else F(2, 3) * V * h(i, j)
        scalarfactor = F(1, 3) if inj == 'wrong_potential_scalar_factor' else F(2, 3)
        Eij = Rij - (p * d(chi, j) if i == 0 else 0) - (p * d(chi, i) if j == 0 else 0)
        if i == j:
            Eij -= scalarfactor * Vp * chi * gdiag[i] + metric_term
        Eij = Eij.substitute({} if inj == 'omit_hamiltonian' else {'V': Vconstraint})
        goal = 2 * M / rho ** 2 if (i, j) == (0, 0) else (-eta[i - 1] * M if i == j else P())
        audit.equal('natural_Einstein_factor_' + str((i, j)), Eij, goal)
    # Scalar variation independently evaluated. Hess(phi0) has yy=py,
    # mu,nu=a*p*g_mu,nu; no scalar Bianchi implication is used.
    hHess = hdiag[0] * py + sum((inverse[i] ** 2 * hdiag[i] * a * p * gdiag[i]
                                  for i in range(1, 5)), P())
    scalar_gauge = P() if inj == 'omit_scalar_gauge_term' else C[0] * p
    scalar = d(d(chi, 0), 0) + 4 * a * d(chi, 0) + box_slice(chi) / rho ** 2 - \
             hHess - scalar_gauge - Vpp * chi
    audit.equal('direct_scalar_variation_formula', scalar,
                d(d(chi, 0), 0) + 4 * a * d(chi, 0) + box_slice(chi) / rho ** 2 -
                Vpp * chi + 6 * p * d(psi, 0) + 4 * Vp * psi)
    factor = -3 / (rho ** 2 * p) * (d(M, 0) + 2 * a * M)
    if inj == 'wrong_scalar_master_derivative':
        factor = -3 / (rho ** 2 * p) * d(M, 0)
    audit.equal('direct_scalar_master_factor', scalar, factor)
    # Values of the master residual at a point do not give its derivative.
    R, Ry = vars_named('residual residual_y')
    point_claim = -3 / (rho ** 2 * p) * Ry
    audit.nonzero('point_value_zero_still_leaves_scalar_derivative', point_claim)
    if inj == 'claim_point_scalar_equation':
        audit.equal('invalid_point_only_scalar_claim', point_claim, 0)
    # Constraint consistency follows with the same full background relations.
    audit.equal('background_V_constraint_derivative', d(Vconstraint, 0), Vp * p)
    return {'Ricci': 15, 'Einstein': 15, 'scalar_direct': True,
            'arbitrary_slice_hessian': True, 'full_parameter_functions': True}


def safe_expression(node, env):
    allowed = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Name, ast.Load, ast.Constant,
               ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub, ast.UAdd,
               ast.Call, ast.Tuple, ast.List, ast.Dict, ast.Subscript)
    tree = ast.Expression(body=node)
    need(all(isinstance(n, allowed) for n in ast.walk(tree)), 'unapproved source expression')
    return eval(compile(ast.fix_missing_locations(tree), '<authenticated-source-expression>', 'eval'),
                {'__builtins__': {}}, env)


def bind_sources(repo, audit):
    geom = repo / 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric'
    tree = ast.parse((geom / 'actual_first_jets.py').read_text())
    function, = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'frame_amplitudes']
    namespace = {'F': F, 'ELL': P.var('ell')}
    module = ast.Module(body=[function], type_ignores=[])
    exec(compile(ast.fix_missing_locations(module), '<authenticated-frame-function>', 'exec'),
         {'__builtins__': {}, **namespace}, namespace)
    alpha, th, ell, W, Wr, H, c = vars_named('alpha vartheta ell W Wr H c')
    temporal, = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'temporal']
    temporal_rhs, = [n.value for n in temporal.body if isinstance(n, ast.Assign) and
                    any(isinstance(t, ast.Name) and t.id == 'rhs' for t in n.targets)]
    th_rhs = safe_expression(temporal_rhs, {'F': F, 'ELL': ell, 'h': H})
    target_temporal_rhs = -F(3, 2) * ell * H
    if audit.injection == 'source_temporal_theta_wrong_factor':
        target_temporal_rhs = -3 * ell * H
    audit.equal('actual_source_vartheta_scaled_log_derivative', th_rhs, target_temporal_rhs)
    B = 3 * alpha / (rho ** 3 * p)
    R = alpha / rho ** 2
    k = rho * (F(5, 2) * a - Vp / p)
    Pmaster = alpha * th * W
    out, scalar = namespace['frame_amplitudes'](W, [P(), P(), P()], rho, p, a, R, B, k, c, H, th)
    audit.equal('actual_source_longitudinal_A', out['A'], -4 * R * th * W)
    audit.equal('actual_source_longitudinal_D', out['D'], -2 * R * th * W)
    audit.equal('actual_source_longitudinal_L', out['L'], 2 * R * th * W)
    for name in ('B', 'C', 'E', 'M'):
        audit.equal('actual_source_longitudinal_zero_' + name, out[name], 0)
    # Coordinate alpha_y follows from the actual z-recurrence read below.
    coretree = ast.parse((geom / 'interval_core.py').read_text())
    background, = [n for n in coretree.body if isinstance(n, ast.FunctionDef) and n.name == 'background']
    rhs_nodes = [n.value for n in ast.walk(background) if isinstance(n, ast.Assign) and
                 isinstance(n.value, ast.Dict) and
                 any(isinstance(t, ast.Name) and t.id == 'rhs' for t in n.targets)]
    need(len(rhs_nodes) == 1 and isinstance(rhs_nodes[0], ast.Dict), 'unique background RHS')
    rhs_nodes = {k0.value: v for k0, v in zip(rhs_nodes[0].keys, rhs_nodes[0].values)}
    env = {'rho': rho, 'a': a, 'p': p, 'alpha': alpha, 'h': Vp / p, 'F': F,
           'poly': lambda coefficients, x: Vp, 'UP': [None, None], 'e': P.var('phi')}
    h_nodes = [n.value for n in ast.walk(background) if isinstance(n, ast.Assign) and
               any(isinstance(t, ast.Name) and t.id == 'h' for t in n.targets)]
    need(len(h_nodes) == 2, 'background h-ratio source shape')
    for i, node in enumerate(h_nodes):
        audit.equal('actual_source_background_potential_ratio_' + str(i), safe_expression(node, env), Vp / p)
    alphaz = safe_expression(rhs_nodes['alpha'], env)
    alphay = -alphaz / rho
    audit.equal('actual_source_alpha_y', alphay, alpha * (Vp / p - F(5, 2) * a))
    audit.equal('actual_source_rho_z', safe_expression(rhs_nodes['rho'], env), -rho ** 2 * a)
    audit.equal('actual_source_a_z', safe_expression(rhs_nodes['a'], env), -rho * ay)
    audit.equal('actual_source_p_z', safe_expression(rhs_nodes['p'], env), -rho * py)
    py_master = alphay * th * W - alpha * th * Wr / (ell * rho)
    if audit.injection == 'reverse_radial_sign':
        py_master = alphay * th * W + alpha * th * Wr / (ell * rho)
    audit.equal('actual_scalar_reconstruction', scalar, -3 * py_master / (rho ** 2 * p))
    # The source's scalar plus-sign Lie shift is checked with free f.
    f = P.var('f')
    _, pulled = namespace['frame_amplitudes'](W, [f, P(), P()], rho, p, a, R, B, k, c, H, th)
    shift = -p * f if audit.injection == 'wrong_source_scalar_shift' else p * f
    audit.equal('actual_scalar_Lie_pullback_plus', pulled - scalar, shift)
    # Execute just the authenticated physical-derivative adapter with an exact
    # ordinary-jet interface. No numerical Taylor remainder is differentiated.
    physical, = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'physical_derivatives']
    physical_namespace = {}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[physical], type_ignores=[])),
                 '<authenticated-physical-derivatives>', 'exec'),
         {'__builtins__': {}, 'ELL': ell}, physical_namespace)

    class Value:
        def __init__(self, value):
            self.value = value

        def val(self):
            return self.value

    class Jet(Value):
        def derivative(self, i, j):
            return P.var('Z' + 'x' * i + 'r' * j)

    actual_derivatives = physical_namespace['physical_derivatives'](Jet(P.var('Z')), Value(rho), Value(a))
    z, zx, zr, zxx, zxr, zrr = vars_named('Z Zx Zr Zxx Zxr Zrr')
    expected_derivatives = {'': z, 'y': -zr / (ell * rho), 't': zx / ell,
                            'yy': zrr / (ell * rho) ** 2 + a * zr / (ell * rho),
                            'yt': -zxr / (ell ** 2 * rho), 'tt': zxx / ell ** 2}
    if audit.injection == 'source_physical_yy_without_rho_term':
        expected_derivatives['yy'] = zrr / (ell * rho) ** 2
    for suffix, value in expected_derivatives.items():
        audit.equal('actual_physical_jet_adapter_' + (suffix or 'value'), actual_derivatives[suffix], value)
    # Bind potential Q directly to its authenticated retained expression.
    potentialpath = repo / 'research/HDBLAST_INTERIOR_MASTER_JETS_20261015/potential/check_actual_potential.py'
    pottree = ast.parse(potentialpath.read_text())
    formal, = [n for n in pottree.body if isinstance(n, ast.FunctionDef) and n.name == 'formal_identities']
    Qnode, = [n.value for n in formal.body if isinstance(n, ast.Assign) and
              any(isinstance(t, ast.Name) and t.id == 'Q' for t in n.targets)]
    Qsource = safe_expression(Qnode, {'h': Vp / p, 'a': a, 'u2': Vpp, 'p': p, 'F': F})
    Q = 2 * (Vp / p) ** 2 - 10 * a * Vp / p + F(35, 4) * a ** 2 - Vpp - p ** 2 / 6
    audit.equal('actual_registered_potential_Q', Qsource, Q)
    return {'alpha': alpha, 'vartheta': th, 'ell': ell, 'W': W, 'Wr': Wr,
            'H': H, 'c': c, 'Q': Q, 'alphay': alphay}


def separated_factorization(values, audit):
    alpha, th, ell, W, Wr, H, c, Q, alphay = [values[n] for n in
        ('alpha', 'vartheta', 'ell', 'W', 'Wr', 'H', 'c', 'Q', 'alphay')]
    inj = audit.injection
    logalpha = alphay / alpha
    radial_constant = d(logalpha, 0) + logalpha ** 2 - (2 * a + 2 * py / p) * logalpha - F(2, 3) * p ** 2
    audit.equal('independent_alpha_radial_constant', radial_constant, -Q + F(5, 2) / rho ** 2)
    drho = 0 if inj == 'freeze_radial_rho' else a * rho
    # Derive W_yy from W_y=-Wr/(ell*rho), with variable rho retained.
    Wy = -Wr / (ell * rho)
    Wrr, Wxx, Wx = vars_named('Wrr Wxx Wx')
    Wyy = Wrr / (ell * rho) ** 2 + drho * Wr / (ell * rho ** 2)
    audit.equal('variable_rho_second_radial_chain', Wyy,
                Wrr / (ell * rho) ** 2 + a * Wr / (ell * rho))
    radial = alpha * th * (Wyy + (2 * logalpha - 2 * a - 2 * py / p) * Wy + radial_constant * W)
    thlog = P() if inj == 'omit_temporal_theta_derivative' else -F(3, 2) * H
    thsecond = (-F(3, 2) * (1 - H ** 2) + thlog ** 2)
    eigenvalue = 6 if inj == 'wrong_angular_eigenvalue' else 8
    plus4 = 0 if inj == 'master_plus4' else 4
    boxterm = alpha * th / rho ** 2 * (-Wxx / ell ** 2 -
        (2 * thlog + 3 * H) * Wx / ell + (-thsecond - 3 * H * thlog - eigenvalue / c ** 2 + plus4) * W)
    # The dS relation c^-2=1-H^2 is exact; do not select interval endpoints.
    full = (radial + boxterm).substitute({'c': P.var('c')})
    def c_quotient(poly):
        out = P()
        for m, coeff in poly.terms.items():
            powers = dict(m)
            cpow = powers.pop('c', 0)
            need(cpow in (0, -2), 'unexpected cosh exponent')
            term = P({tuple(sorted(powers.items())): coeff})
            out += term * ((1 - H ** 2) if cpow == -2 else 1)
        return out
    full = c_quotient(full)
    temporal_constant = c_quotient(-thsecond - 3 * H * thlog - eigenvalue / c ** 2 + plus4)
    audit.equal('independent_temporal_constant', temporal_constant, -F(5, 2) + F(35, 4) * H ** 2)
    if inj == 'lose_radial_constant':
        full -= alpha * th * F(5, 2) * W / rho ** 2
    Qtarget = Q
    if inj == 'wrong_potential_p2':
        Qtarget += p ** 2 / 3
    if inj == 'wrong_Vpp_sign':
        Qtarget += 2 * Vpp
    residual = Wxx - Wrr + ell ** 2 * (rho ** 2 * Qtarget - F(35, 4) * H ** 2) * W
    audit.equal('actual_separated_master_factor', full, -alpha * th * residual / (ell ** 2 * rho ** 2))
    comparator_defect = ell ** 2 * (rho ** 2 * Q - F(35, 4) * H ** 2) * W
    audit.nonzero('free_comparator_retains_actual_potential_defect', comparator_defect)
    if inj == 'claim_comparator_is_actual':
        audit.equal('invalid_comparator_for_actual_equation', comparator_defect, 0)
    # Differentiate the complete residual identity, treating residual and its
    # radial derivative as free independent symbols (master order three).
    Rm, Rmr = vars_named('Rmaster Rmaster_r')
    multiplier = -alpha * th / (ell ** 2 * rho ** 2)
    DyM = multiplier * ((logalpha - 2 * a) * Rm - Rmr / (ell * rho))
    scalar = -3 / (rho ** 2 * p) * (DyM + 2 * a * multiplier * Rm)
    scalar_radial_multiplier = rho if inj == 'scalar_derivative_missing_ell' else ell * rho
    audit.equal('actual_scalar_residual_and_derivative_factor', scalar,
                3 * alpha * th / (ell ** 2 * rho ** 4 * p) * (logalpha * Rm - Rmr / scalar_radial_multiplier))
    audit.nonzero('master_residual_derivative_cannot_be_dropped',
                  -3 * alpha * th * Rmr / (ell ** 3 * rho ** 5 * p))
    audit.equal('opposite_five_halves_cancel', F(5, 2) + -F(5, 2), 0)


def run(repo, injection):
    pins_before = authenticate(repo)
    audit = Audit(injection)
    tensor = derive_longitudinal(audit)
    values = bind_sources(repo, audit)
    separated_factorization(values, audit)
    pins_after = authenticate(repo)
    need(pins_before == pins_after, 'source closure changed')
    return {
        'status': 'PASS_ACTUAL_FIRST_NATURAL_MASTER_BINDING',
        'arithmetic': 'exact rational Laurent polynomial coefficients; no floats or sampling',
        'input_closures': pins_before,
        'conventions': {'E': 'Ric-dphi*dphi-(2/3)Vg', 'S': 'Box(phi)-Vprime',
            'slice_signature': '-+++', 'outward_y': 'dy=-ell*rho*dr',
            'master': 'W_xx-W_rr+ell^2*(rho^2 Q-(35/4)tanh(t)^2)*W',
            'metric': 'h_yy=-4Psi; h_mu,nu=2rho^2 Psi gamma_mu,nu',
            'scalar': '-3*(rho^2 Psi)_y/(rho^2 p)',
            'harmonic_transfer': 'prior exact implemented plus-sign Lie pullback and natural covariance'},
        'retained_analytic_master_PDE_source': {
            'path': 'research/HDBLAST_INTERIOR_MASTER_JETS_20261015/jets/HIGHER_JETS.md',
            'meaning': 'actual smooth master solves the displayed PDE as a neighbourhood identity',
            'treatment': 'authenticated analytic input, not proof of solution existence from coefficient algebra'},
        'independent_tensor_route': tensor,
        'checks': audit.rows,
        'scope': {'background_field_equations_retained': True, 'p_and_rho_nonzero_collar': True,
                  'actual_master_PDE_and_regular_solution_are_analytic_inputs': True,
                  'first_scalar_equation_requires_master_neighbourhood_identity': True,
                  'existing_analytic_reconstruction_proof_already_present': True,
                  'new_dynamic_shell_trace_join': False,
                  'quadratic_bulk_Bianchi_map_bound': False,
                  'finite_amplitude_nonlinear_family_proved': False,
                  'cosmological_origin_evidence': False},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inject-error', choices=('none',) + CONTROLS, default='none')
    parser.add_argument('--run-controls', action='store_true')
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse overwrite')
    report = run(args.repo.resolve(), args.inject_error)
    if args.run_controls:
        records = []
        with tempfile.TemporaryDirectory(prefix='hdblast-natural-controls-') as temp:
            for control in CONTROLS:
                out = Path(temp) / (control + '.json')
                command = [sys.executable, '-I', '-S', '-B']
                if sys.flags.optimize:
                    command.append('-O')
                command += [str(Path(__file__).resolve()), '--repo', str(args.repo.resolve()),
                            '--output', str(out), '--inject-error', control]
                completed = subprocess.run(command, capture_output=True, text=True, timeout=180)
                expected = 'mathematical mismatch ' + EXPECTED_FAILURES[control] + ':'
                need(completed.returncode != 0 and not out.exists() and
                     expected in completed.stderr,
                     'corruption was not mathematically rejected: ' + control)
                records.append({'control': control, 'fresh_process_rejected': True,
                                'expected_failed_identity': EXPECTED_FAILURES[control], 'output_absent': True})
        report['wrong_formula_controls'] = records
    args.output.write_text(json.dumps(report, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'checks': len(report['checks']),
                      'controls': len(report.get('wrong_formula_controls', [])),
                      'report_sha256': digest(args.output.read_bytes())}, sort_keys=True))


if __name__ == '__main__':
    main()
