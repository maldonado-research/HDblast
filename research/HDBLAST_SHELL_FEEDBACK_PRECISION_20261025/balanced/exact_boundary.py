#!/usr/bin/env python3
"""Fixed-budget, outward-rational evaluation of an exact-u Robin repair.

Only this new candidate is evaluated. No claim about the exact response sign
is made, and no frozen numerical producer is executed.
"""
import argparse
import ast
import copy
import hashlib
import importlib.util
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
LABELS = {
    'current': 'HDBLAST_COUPLED_MEAN_PILOT_20261025',
    'I': 'HDBLAST_INPUT_PRECISION_20261023',
    'P': 'HDBLAST_COMMON_GAUGE_FORCING_20261016',
    'N': 'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016',
    'G': 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
    'C': 'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023',
    'T': 'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023',
}
PINS = {
    'current': 'b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de',
    'I': 'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61',
    'P': '8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6',
    'N': '928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44',
    'G': '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2',
    'C': '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76',
    'T': 'cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72',
}
SEGMENTS = ((F(65, 64), F(66, 64), 16),
            (F(66, 64), F(9, 8), 112),
            (F(9, 8), F(5, 4), 128))


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def authenticate(roots):
    result = {}
    for key, root in roots.items():
        need(not root.is_symlink(), 'ROOT_SYMLINK:' + key)
        manifest = root / 'MANIFEST.json'
        need(sha(manifest) == PINS[key], 'MANIFEST:' + key)
        seen = {'MANIFEST.json'}
        for entry in json.loads(manifest.read_text())['files']:
            rel = Path(entry['path'])
            need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,
                 'MEMBER_PATH:' + key)
            seen.add(str(rel))
            path = root / rel
            need(path.is_file() and not path.is_symlink(), 'MEMBER_TYPE:' + str(path))
            need(path.stat().st_size == entry['bytes'] and sha(path) == entry['sha256'],
                 'MEMBER_BYTES:' + str(path))
        for path in root.rglob('*'):
            need(not path.is_symlink(), 'TREE_SYMLINK:' + str(path))
            if path.is_file():
                need(str(path.relative_to(root)) in seen, 'UNLISTED_FILE:' + str(path))
        result[key] = {'manifest_sha256': PINS[key], 'files_including_manifest': len(seen)}
    return result


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    roots = {key: args.repo_root.resolve()/'research'/name for key, name in LABELS.items()}
    need(not args.output.exists(), 'FRESH_OUTPUT_REQUIRED')
    need(all(root not in args.output.resolve().parents for root in roots.values()),
         'OUTPUT_INSIDE_FROZEN_CAPSULE')
    before = authenticate(roots)  # Before every project import.
    native = load('balanced_frozen_native', roots['current']/'native_pilot.py')
    init = load('balanced_frozen_intervals', roots['N']/'numeric/evaluate_initial_profiles.py')
    kernel = load('balanced_frozen_kernel', roots['current']/'reflected_ridge.py')
    direct = load('balanced_frozen_direct', roots['I']/'initial/direct_difference_pilot.py')
    depmap = {key: str(roots[key]) for key in ('P', 'N', 'G', 'C', 'T')}
    depmap['B'] = str(roots['I']/'branch/BRANCH_COMPARISON_REPORT.json')
    dependency_before = direct.authenticate(depmap)
    need(native.authenticate_precision(roots['I'])['manifest_sha256'] == PINS['I'],
         'PRECISION_NATIVE_AUTHENTICATION')
    I, J, ell = init.I, init.J, init.ELL
    need(ell == F(1, 10**6) and init.GRID == 1 << 112, 'ARITHMETIC_CONVENTION')
    frozen = json.loads((roots['I']/'initial/RECENTERED_CANDIDATE_NORMAL.json').read_text())
    coefficients = frozen['representation_background_coefficients_ascending']
    digest = hashlib.sha256(json.dumps(coefficients, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    need(digest == native.COEFFICIENT_PIN == frozen['representation_background_coefficients_sha256'],
         'FIXED_REGISTERED_COEFFICIENTS')
    cs = {key: list(map(F, row)) for key, row in coefficients.items()}
    data = json.loads((roots['current']/'results/PROFILES.json').read_text())
    old = json.loads((roots['current']/'results/EVALUATED_CANDIDATE.json').read_text())
    source = kernel.from_data(data)
    M = [list(map(F, row)) for row in data['matrices']['M']]
    m = M[3][4]
    need(all(M[3][j] == 0 for j in (1, 2, 3)), 'J_ROW_STRUCTURE')
    need(m < 0, 'J_W_ROBIN_SIGN')
    need(all(F(row[0]) == 0 for row in data['profiles']['u']), 'P1_K_DISPLACEMENT_ZERO')

    # Extract unchanged registered f expressions. Omit later unused velocity
    # calculations and replace the zero-collar return by its f component.
    tree = ast.parse((roots['current']/'native_pilot.py').read_text())
    main_node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
    wanted = ('constant', 'polynomial', 'forty_squares', 'nominal_profiles')
    nodes = [copy.deepcopy(n) for n in main_node.body if isinstance(n, ast.FunctionDef) and n.name in wanted]
    need(tuple(n.name for n in nodes) == wanted, 'AST_REGISTERED_FUNCTIONS')
    nominal = nodes[-1]
    stop = next(i for i, n in enumerate(nominal.body)
                if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'f' for t in n.targets))
    nominal.body = nominal.body[:stop+1] + [ast.Return(value=ast.Name(id='f', ctx=ast.Load()))]
    nominal.body[0].body[0] = ast.Return(value=ast.Call(func=ast.Name(id='constant', ctx=ast.Load()), args=[ast.Constant(0)], keywords=[]))
    extracted = ast.fix_missing_locations(ast.Module(body=nodes, type_ignores=[]))
    ns = {'I': I, 'J': J, 'F': F, 'init': init, 'ell': ell, 'cs': cs,
          'order': 6, 'power': 1 << 40, 'need': need}
    exec(compile(extracted, str(roots['current']/'native_pilot.py')+':extracted_f', 'exec'), ns)
    f_source_hash = hashlib.sha256(ast.dump(extracted).encode()).hexdigest()

    def f_jet(r, order):
        ns['order'] = order
        return ns['nominal_profiles'](r)

    def C_jet(r, order):
        rr = J([r, 1/ell]+[0]*(order-1))
        z = ell*rr
        rho = ns['polynomial'](cs['rho'], z)
        pp = ns['polynomial'](cs['p'], z)
        return 3/(2*ell*rho*pp)

    # m*C' f is differentiated in z by the extracted jet implementation;
    # c4 * 4! * ell^4 converts its fourth derivative back to r.
    def integrand(r, order):
        ff = f_jet(r, order)
        cc = C_jet(r, order)
        return m*ell*cc.dz()*ff

    point_cache = {}
    def point(r):
        if r not in point_cache:
            point_cache[r] = integrand(I(r), 2).val()
        return point_cache[r]

    # Exact polynomial numerator avoids artificial cancellation in 1-m*C.
    rp = [sum((cs['rho'][i]*cs['p'][k-i] for i in range(len(cs['rho']))
               if 0 <= k-i < len(cs['p'])), F(0))
          for k in range(len(cs['rho'])+len(cs['p'])-1)]
    numerator = [2*ell*c for c in rp]
    numerator[0] -= 3*m
    def balanced_endpoint(x):
        z = I(ell*x)
        rho_p = init.poly(rp, z)
        factor = init.poly(numerator, z)/(2*ell*rho_p)
        return factor*f_jet(I(x), 2).val(), factor

    need(sum(n for a, b, n in SEGMENTS) == 256 and all(n % 2 == 0 for a, b, n in SEGMENTS),
         'FIXED_QUADRATURE_BUDGET')
    integral_sum = I(0)
    error_sum = F(0)
    certificates = []
    samples = []
    for segment_index, (left, right, count) in enumerate(SEGMENTS):
        h = (right-left)/count
        for k in range(0, count, 2):
            a, mid, b = left+k*h, left+(k+1)*h, left+(k+2)*h
            simpson = (b-a)/6*(point(a)+4*point(mid)+point(b))
            fourth = math.factorial(4)*ell**4*integrand(I(a, b), 6).c[4]
            radius = (b-a)**5*fourth.mag()/2880
            integral_sum = integral_sum+simpson
            error_sum += radius
            certificates.append({'r': [a, b], 'midpoint': mid,
                                 'simpson_interval': simpson.row(),
                                 'fourth_derivative_interval': fourth.row(),
                                 'remainder_upper': radius})
            if len(certificates) % 16 == 0:
                print(json.dumps({'completed_simpson_panels': len(certificates),
                                  'fixed_total_simpson_panels': 128}), flush=True)
        if right not in (F(9, 8), F(5, 4)):
            continue
        x = right
        integrated = integral_sum + I(-error_sum, error_sum)
        endpoint, factor = balanced_endpoint(x)
        bulk = kernel.reflected_kernel(source.F, source.G, x, F(0))['value'][3]
        pv = source.v.primitive(x, 1)
        qv = source.v.primitive(x, 2)
        velocity_and_bulk = pv[3]+bulk-M[3][0]*qv[0]-m*qv[4]
        candidate = endpoint+integrated+velocity_and_bulk
        old_sample = next(s for s in old['samples'] if F(s['x']) == x and F(s['r']) == 0)
        old_candidate = I(*old_sample['complete_candidate_interval_evaluation']['value'][3])
        samples.append({'x': x, 'r': F(0), 'new_J_candidate_interval': candidate.row(),
                        'balanced_endpoint_interval': endpoint.row(),
                        'one_minus_mC_interval': factor.row(),
                        'm_Cprime_f_integral_interval': integrated.row(),
                        'simpson_sum_interval': integral_sum.row(),
                        'simpson_remainder_upper': error_sum,
                        'old_velocity_and_bulk_exact': velocity_and_bulk,
                        'old_bulk_J_exact': bulk,
                        'old_velocity_primitive_J_exact': pv[3],
                        'old_boundary_velocity_J_exact': -M[3][0]*qv[0]-m*qv[4],
                        'old_J_candidate_interval': old_candidate.row(),
                        'candidate_change_interval': (candidate-old_candidate).row(),
                        'display_only': {'new_J_midpoint': float((candidate.lo+candidate.hi)/2),
                                         'new_J_radius': float((candidate.hi-candidate.lo)/2),
                                         'quadrature_remainder': float(error_sum),
                                         'balanced_endpoint': float((endpoint.lo+endpoint.hi)/2),
                                         'balanced_integral': float((integrated.lo+integrated.hi)/2),
                                         'old_velocity_and_bulk': float(velocity_and_bulk)}})
        print(json.dumps({'x': str(x), **samples[-1]['display_only']}), flush=True)

    need(len(point_cache) == 257 and len(certificates) == 128, 'EXACT_FIXED_EVALUATION_COUNTS')
    need(authenticate(roots) == before, 'FROZEN_CLOSURE_CHANGED')
    need(direct.authenticate(depmap) == dependency_before, 'UPSTREAM_CLOSURE_CHANGED')
    result = {'schema': 'HDBLAST_EXACT_BOUNDARY_CANDIDATE_QUADRATURE_V1',
              'status': 'EVALUATED_NEW_EXACT_BOUNDARY_CANDIDATE_ONLY',
              'authentication': before,
              'complete_closures_authenticated_before_and_after': True,
              'candidate': 'B0=free(exact registered u0,old vP1); B1=old bulk feedback + L(Mstar PN(exact u0+primitive(old vP1)))',
              'J_shell_formula': '(1-m C(x))*f(x) + integral_0^x m*Cprime*f + Pv_J + old_bulk_J - M_JK*Qv_K - m*Qv_w',
              'definitions': {'f': 'registered u0_J', 'C': '3/(2*ell*rho_nom*p_nom)',
                              'm': m, 'M_JK': M[3][0], 'Cprime': 'ordinary r derivative of nominal degree-four background polynomials'},
              'arithmetic': {'outward_grid_denominator': str(init.GRID), 'decisions': 'exact Fractions only',
                             'AST_registered_f_sha256': f_source_hash,
                             'coefficient_sha256': digest,
                             'point_jet_order': 2, 'whole_panel_jet_order': 6,
                             'quadrature': 'composite Simpson; per-panel error <= width^5 sup|h4|/2880',
                             'fixed_subintervals': 256, 'simpson_panels': 128,
                             'distinct_point_evaluations': 257, 'whole_panel_evaluations': 128},
              'fixed_segments': SEGMENTS, 'samples': samples, 'panel_certificates': certificates,
              'limits': ['Candidate evaluation only; total-response and observable errors must be supplied separately.',
                         'No actual-response, observable-sign, nonlinear-family or physical-origin claim.',
                         'No frozen mesh or old producer was rerun or modified.']}
    with args.output.open('x') as handle:
        json.dump(encode(result), handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({'output': str(args.output), 'sha256': sha(args.output)}), flush=True)


if __name__ == '__main__':
    main()
