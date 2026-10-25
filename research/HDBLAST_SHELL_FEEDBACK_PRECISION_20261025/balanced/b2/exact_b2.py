#!/usr/bin/env python3
"""Bounded exact-u constant-operator two-feedback shell evaluation."""
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


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def need(ok, msg):
    if not ok:
        raise ValueError(msg)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encoded(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {k: encoded(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encoded(v) for v in x]
    if hasattr(x, 'row'):
        return x.row()
    return x


def plus(a, b):
    return [x+y for x, y in zip(a, b)]


def minus(a, b):
    return [x-y for x, y in zip(a, b)]


def scale(c, a):
    return [c*x for x in a]


def mv(m, a):
    return [sum((x*y for x, y in zip(row, a)), 0) for row in m]


def add(*vectors):
    out = [0]*5
    for row in vectors:
        out = plus(out, row)
    return out


E = (1, -1, 1, 1, 1)
def reflect(row):
    return [e*x for e, x in zip(E, row)]


def project(row):
    return [x if e == 1 else 0*x for e, x in zip(E, row)]


def shell_algebra(T, U, V, matrices):
    """Linear in exact moments U0..U4 and V1..V5; scalar type agnostic."""
    A, B, D, M = [matrices[k] for k in ('A', 'B', 'D', 'M')]
    ApB = [[a+b for a, b in zip(ar, br)] for ar, br in zip(A, B)]
    BmA = [[b-a for a, b in zip(ar, br)] for ar, br in zip(A, B)]
    def moments(k, negative=False):
        uu, up, vv, vp = U[k], U[k-1], V[k], V[k+1]
        if negative:
            uu = scale((-1)**k, reflect(uu))
            up = scale((-1)**(k-1), reflect(up))
            vv = scale((-1)**k, reflect(vv))
            vp = scale((-1)**(k+1), reflect(vp))
        ff = scale(F(1, 2), plus(mv(ApB, plus(up, vv)), mv(D, plus(uu, vp))))
        gg = scale(F(1, 2), plus(mv(BmA, minus(up, vv)), mv(D, minus(uu, vp))))
        return ff, gg
    fg = {k: moments(k) for k in range(1, 5)}
    gm = {k: moments(k, True)[1] for k in range(1, 5)}
    ff = {k: fg[k][0] for k in fg}
    gg = {k: fg[k][1] for k in fg}
    g = {k: mv(M, project(plus(U[k], V[k+1]))) for k in range(1, 4)}
    b0 = project(plus(U[0], V[1]))
    bulk1 = project(add(scale(T, ff[1]), scale(-1, ff[2]),
                        scale(F(1, 2), plus(gg[2], gm[2]))))
    b1 = minus(bulk1, g[1])
    H = add(scale(T, ff[2]), scale(-2, ff[3]), scale(F(1, 2), minus(gg[3], gm[3])))
    bh = minus(project(H), g[2])
    ch = minus(scale(F(1, 4), add(scale(T*T, ff[1]), scale(-T, ff[2]),
                scale(T, gg[2]), scale(-1, plus(gg[3], gm[3])), reflect(H))),
                scale(F(1, 2), g[2]))
    vh = minus(scale(F(1, 4), add(scale(T*T, ff[2]), scale(-3*T, ff[3]),
                scale(3, ff[4]), scale(T, gg[3]), scale(-2, gg[4]), gm[4],
                reflect(add(scale(T, ff[3]), scale(-3, ff[4]),
                            scale(F(1, 2), plus(gg[4], gm[4])))))), scale(F(1, 2), g[3]))
    b2 = minus(project(add(mv(ApB, ch), scale(-1, mv(B, bh)), mv(D, vh))),
               mv(M, project(bh)))
    return {'B0': b0, 'B1star_bulk': bulk1, 'B1star_boundary': scale(-1, g[1]),
            'B1star': b1, 'B2star': b2, 'two_term': plus(b0, b1),
            'three_term': add(b0, b1, b2),
            'B1star_shell_history': bh, 'B1star_characteristic_history': ch,
            'B1star_triangle_history': vh}


def velocity_moments(T, data):
    """Exact direct integration of frozen P1 segments, through order five."""
    nodes = list(map(F, data['nodes']))
    vals = [list(map(F, row)) for row in data['profiles']['v']]
    out = {k: [F(0)]*5 for k in range(1, 6)}
    for a, b, va, vb in zip(nodes, nodes[1:], vals, vals[1:]):
        if T <= a:
            break
        end = min(T, b)
        d, e = T-a, T-end
        slopes = [(y-x)/(b-a) for x, y in zip(va, vb)]
        for k in range(1, 6):
            for j in range(5):
                out[k][j] += ((va[j]+slopes[j]*d)*(d**k-e**k)/math.factorial(k)
                              -k*slopes[j]*(d**(k+1)-e**(k+1))/math.factorial(k+1))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    need(not args.output.exists(), 'FRESH_OUTPUT')
    auth_path = Path(__file__).resolve().parent.parent/'exact_boundary.py'
    need(sha(auth_path) == 'c01a8f5f95737dcacebdc26379ae0ee5863af9c1977890b3901ad6113103e619', 'AUTHENTICATOR_SOURCE')
    auth = load('b2_boundary_auth', auth_path)
    roots = {k: args.repo_root.resolve()/'research'/name for k, name in auth.LABELS.items()}
    need(all(root not in args.output.resolve().parents for root in roots.values()), 'OUTPUT_INSIDE_FROZEN')
    before = auth.authenticate(roots)
    native = load('b2_frozen_native', roots['current']/'native_pilot.py')
    init = load('b2_frozen_interval', roots['N']/'numeric/evaluate_initial_profiles.py')
    direct = load('b2_frozen_direct', roots['I']/'initial/direct_difference_pilot.py')
    depmap = {key: str(roots[key]) for key in ('P', 'N', 'G', 'C', 'T')}
    depmap['B'] = str(roots['I']/'branch/BRANCH_COMPARISON_REPORT.json')
    dependency_before = direct.authenticate(depmap)
    I, J, ell = init.I, init.J, init.ELL
    need(ell == F(1, 10**6) and init.GRID == 1 << 112, 'ARITHMETIC_CONVENTION')
    frozen = json.loads((roots['I']/'initial/RECENTERED_CANDIDATE_NORMAL.json').read_text())
    coefficients = frozen['representation_background_coefficients_ascending']
    digest = hashlib.sha256(json.dumps(coefficients, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    need(digest == native.COEFFICIENT_PIN == frozen['representation_background_coefficients_sha256'], 'NOMINAL_COEFFICIENTS')
    cs = {k: list(map(F, row)) for k, row in coefficients.items()}
    data = json.loads((roots['current']/'results/PROFILES.json').read_text())
    matrices = {k: [list(map(F, row)) for row in m] for k, m in data['matrices'].items()}
    need(all(x == 0 for x in matrices['M'][1]), 'BOUNDARY_NEUMANN_OUTPUT')
    tree = ast.parse((roots['current']/'native_pilot.py').read_text())
    mn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
    names = ('constant', 'polynomial', 'forty_squares', 'nominal_profiles')
    nsnodes = [copy.deepcopy(n) for n in mn.body if isinstance(n, ast.FunctionDef) and n.name in names]
    need(tuple(n.name for n in nsnodes) == names, 'AST_FUNCTION_ORDER')
    fn = nsnodes[-1]
    stop = next(i for i, n in enumerate(fn.body) if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'f' for t in n.targets))
    fn.body = fn.body[:stop+1]+ast.parse('return f, -u*u/3').body
    fn.body[0].body = ast.parse('return constant(0), constant(0)').body
    module = ast.fix_missing_locations(ast.Module(body=nsnodes, type_ignores=[]))
    ns = {'I': I, 'J': J, 'F': F, 'init': init, 'ell': ell, 'cs': cs,
          'order': 6, 'power': 1 << 40, 'need': need}
    exec(compile(module, str(roots['current']/'native_pilot.py')+':b2_extracted', 'exec'), ns)
    asthash = hashlib.sha256(ast.dump(module).encode()).hexdigest()
    targets = (F(9, 8), F(5, 4))
    def base(r, order):
        ns['order'] = order
        ff, nn = ns['nominal_profiles'](r)
        rr = J([r, 1/ell]+[0]*(order-1))
        rho = ns['polynomial'](cs['rho'], ell*rr)
        pp = ns['polynomial'](cs['p'], ell*rr)
        cc = 3/(2*ell*rho*pp)
        return ff, nn, cc, ell*cc.dz(), rr
    cache = {}
    def node(r):
        if r not in cache:
            cache[r] = base(I(r), 3 if r in targets else 2)
        return cache[r]
    def integrands(T, values):
        ff, nn, cc, cp, rr = values
        weights = {k: (T-rr)**(k-1)/math.factorial(k-1) for k in range(1, 5)}
        rows = {}
        for k in range(1, 5):
            ww = -weights[k]*cp*ff if k == 1 else (weights[k-1]*cc-weights[k]*cp)*ff
            rows[k] = [weights[k]*nn, weights[k]*ff, ww]
        return rows
    sums = {T: {k: [I(0)]*3 for k in range(1, 5)} for T in targets}
    errors = {T: {k: [F(0)]*3 for k in range(1, 5)} for T in targets}
    panel_certificates = []
    need(sum(n for a, b, n in auth.SEGMENTS) == 256, 'FIXED_MESH_CAP')
    for left, right, count in auth.SEGMENTS:
        step = (right-left)/count
        for index in range(0, count, 2):
            a, mid, b = left+index*step, left+(index+1)*step, left+(index+2)*step
            av, midv, bv = node(a), node(mid), node(b)
            whole = base(I(a, b), 6)
            row = {'r': [a, b], 'midpoint': mid, 'targets': {}}
            for T in targets:
                if a >= T:
                    continue
                aa, mm, bb = [integrands(T, v) for v in (av, midv, bv)]
                jj = integrands(T, whole)
                moments = {}
                for k in range(1, 5):
                    entries = []
                    for j in range(3):
                        simpson = (b-a)/6*(aa[k][j].val()+4*mm[k][j].val()+bb[k][j].val())
                        fourth = math.factorial(4)*ell**4*jj[k][j].c[4]
                        rem = (b-a)**5*fourth.mag()/2880
                        sums[T][k][j] = sums[T][k][j]+simpson
                        errors[T][k][j] += rem
                        entries.append({'simpson_interval': simpson, 'fourth_derivative_interval': fourth,
                                        'remainder_upper': rem})
                    moments[k] = entries
                row['targets'][str(T)] = moments
            panel_certificates.append(row)
            if len(panel_certificates) % 16 == 0:
                print(json.dumps({'completed_panels': len(panel_certificates), 'fixed_total': 128}), flush=True)
    need(len(cache) == 257 and len(panel_certificates) == 128, 'FIXED_COUNTS')
    samples = []
    for T in targets:
        ff, nn, cc, cp, rr = node(T)
        endpoint = [I(0), I(0), nn.val(), ff.val(), (ell*cc*ff.dz()).val()]
        U = {0: endpoint}
        Ucentral = {0: endpoint}
        Uerror = {0: [I(0)]*5}
        for k in range(1, 5):
            U[k], Ucentral[k], Uerror[k] = [I(0)]*5, [I(0)]*5, [I(0)]*5
            for j, component in enumerate((2, 3, 4)):
                center = sums[T][k][j]
                if k == 1 and component == 4:
                    center = center+(cc*ff).val()
                radius = errors[T][k][j]
                Ucentral[k][component] = center
                Uerror[k][component] = I(-radius, radius)
                U[k][component] = center+Uerror[k][component]
        V = velocity_moments(T, data)
        result = shell_algebra(T, U, V, matrices)
        central_result = shell_algebra(T, Ucentral, V, matrices)
        quadrature = shell_algebra(T, Uerror, {k: [F(0)]*5 for k in V}, matrices)
        quadratic_caps = {name: [x.mag() if hasattr(x, 'mag') else abs(x) for x in vec]
                          for name, vec in quadrature.items()}
        need(result['three_term'][1].lo == result['three_term'][1].hi == 0, 'EXACT_DIRICHLET_SHELL')
        sample = {'x': T, 'r': F(0), 'exact_u_moment_enclosures': U,
                  'moment_simpson_sum_enclosures_before_scalar_endpoint': sums[T],
                  'moment_remainder_uppers_N_J_w': errors[T],
                  'velocity_moments_exact': V,
                  'evaluation': result,
                  'evaluation_omitting_quadrature_remainders': central_result,
                  'quadrature_error_component_uppers': quadratic_caps,
                  'display_only': {key: {'midpoints': [float((v.lo+v.hi)/2) if hasattr(v, 'lo') else float(v) for v in row],
                                         'radii': [float((v.hi-v.lo)/2) if hasattr(v, 'lo') else 0.0 for v in row]}
                                   for key, row in result.items()}}
        samples.append(sample)
        print(json.dumps({'x': str(T), 'B1star': sample['display_only']['B1star'],
                          'B2star': sample['display_only']['B2star'],
                          'three_term': sample['display_only']['three_term']}), flush=True)
    need(auth.authenticate(roots) == before, 'FROZEN_CLOSURES_CHANGED')
    need(direct.authenticate(depmap) == dependency_before, 'DEPENDENCY_CLOSURES_CHANGED')
    out = {'schema': 'HDBLAST_EXACT_U_B2_SHELL_CANDIDATE_V1',
           'status': 'EVALUATED_BOUNDED_CONSTANT_OPERATOR_THREE_TERM_CANDIDATE',
           'mean_order': ['K0', 'L0', 'N0', 'J0', 'w0'],
           'candidate': 'B0=free(exact_registered_u0,old_vP1); S=constant frozen bulk and Robin operator; V=B0+S B0+S^2 B0',
           'authentication': before, 'complete_closures_authenticated_before_and_after': True,
           'source_sha256': sha(Path(__file__)), 'authenticator_sha256': sha(auth_path),
           'registered_f_N_AST_sha256': asthash, 'candidate_coefficients_sha256': digest,
           'frozen_source_pins': {str(path.relative_to(args.repo_root.resolve())): sha(path)
                                  for path in (roots['current']/'native_pilot.py',
                                               roots['N']/'numeric/evaluate_initial_profiles.py',
                                               roots['current']/'results/PROFILES.json',
                                               roots['I']/'initial/RECENTERED_CANDIDATE_NORMAL.json')},
           'matrices': matrices,
           'arithmetic': {'outward_grid_denominator': str(init.GRID), 'only_exact_rational_decisions': True,
                          'fixed_subintervals': 256, 'distinct_point_nodes': 257, 'whole_panel_evaluations': 128,
                          'point_jet_order': 2, 'endpoint_jet_order': 3, 'panel_jet_order': 6,
                          'moment_component_order': ['N0', 'J0', 'w0'],
                          'Simpson_bound': 'width^5 * sup|h4| / 2880; h4=24*ell^4*jet.c[4]',
                          'w_moments': 'W1=C(T)f(T)-int Cprime f; Wk=int[(p{k-1} C-pk Cprime)f], k>=2'},
           'segments': auth.SEGMENTS, 'samples': samples, 'panel_certificates': panel_certificates,
           'limits': ['Candidate evaluation only. Full actual-coefficient defects and T^3 total-response tail require a separate ledger.',
                      'No exact-response sign, observable sign, nonlinear family, or physical-origin conclusion.',
                      'Frozen P1 displacement and derivative are absent from both feedback terms; frozen velocity is retained.',
                      'No old numerical stage was rerun, and no frozen capsule was modified.']}
    with args.output.open('x') as handle:
        json.dump(encoded(out), handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({'output': str(args.output), 'sha256': sha(args.output)}), flush=True)


if __name__ == '__main__':
    main()
