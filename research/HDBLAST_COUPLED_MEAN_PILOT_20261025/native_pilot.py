#!/usr/bin/env python3
"""Bounded coupled five-mean feedback pilot; all scientific decisions use Fractions.

The exact registered rational displacement remains in B0. The three P1 tables
approximate u and du only inside feedback; v also defines B0's velocity.
A failed accuracy gate is a finite scientific result, not a process failure.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True

PRECISION_PIN = 'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61'
COEFFICIENT_PIN = '739240283c2bee53afa559c4ce71c1e3356cf2d0ad4e4b0f18b9621c077e20d3'
OLD_COMMIT = '2fb40bdc5a1785ef6a55bc50e6c9d81a3e3e9978'
WEIGHTS = (1, 2, 1, 3, 1)
LEVELS = (64, 128, 256)
Q = F(3409, 8000000)
UCAP = F(103 * 10**31)
SOURCE_CAP = F(52162547172519, 2)
TARGET = F(10**29)


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def authenticate_precision(root):
    need(not root.is_symlink(), 'PRECISION_ROOT_SYMLINK')
    need(sha(root / 'MANIFEST.json') == PRECISION_PIN, 'PRECISION_MANIFEST')
    entries = json.loads((root / 'MANIFEST.json').read_text())['files']
    seen = set()
    for entry in entries:
        rel = Path(entry['path'])
        need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,
             'PRECISION_MANIFEST_PATH')
        seen.add(str(rel))
        path = root / rel
        need(not path.is_symlink() and path.is_file(), 'PRECISION_MEMBER_TYPE')
        need(path.stat().st_size == entry['bytes'] and sha(path) == entry['sha256'],
             'PRECISION_MEMBER_HASH:' + str(rel))
    for path in root.rglob('*'):
        need(not path.is_symlink(), 'PRECISION_SYMLINK')
        if path.is_file():
            need(str(path.relative_to(root)) in seen | {'MANIFEST.json'},
                 'PRECISION_UNLISTED_FILE')
    return {'manifest_sha256': PRECISION_PIN, 'files': len(entries)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    output = args.output_dir.resolve()
    research = repo / 'research'
    labels = {
        'P': 'HDBLAST_COMMON_GAUGE_FORCING_20261016',
        'N': 'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016',
        'G': 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
        'C': 'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023',
        'T': 'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023',
    }
    roots = {key: research / name for key, name in labels.items()}
    precision = research / 'HDBLAST_INPUT_PRECISION_20261023'
    need(not output.exists(), 'OUTPUT_DIRECTORY_MUST_BE_NEW')
    need(precision not in output.parents and all(root not in output.parents for root in roots.values()),
         'OUTPUT_INSIDE_FROZEN_CAPSULE')
    current_before = authenticate_precision(precision)
    direct = load_module('frozen_direct_difference', precision / 'initial/direct_difference_pilot.py')
    branch = precision / 'branch/BRANCH_COMPARISON_REPORT.json'
    depmap = {key: str(value) for key, value in roots.items()}
    depmap['B'] = str(branch)
    before = direct.authenticate(depmap)
    need(sha(branch) == direct.BRANCH_PIN, 'BRANCH_REPORT_PIN')
    branch_report = json.loads(branch.read_text())
    need(branch_report['status'] == 'PASS_RETAINED_BRANCH_ANALYTIC_COMPARISON_ENCLOSURE',
         'BRANCH_STATUS')
    frozen_path = precision / 'initial/RECENTERED_CANDIDATE_NORMAL.json'
    frozen = json.loads(frozen_path.read_text())
    coefficients = frozen['representation_background_coefficients_ascending']
    digest = hashlib.sha256(json.dumps(coefficients, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    need(digest == COEFFICIENT_PIN == frozen['representation_background_coefficients_sha256'],
         'FIXED_CANDIDATE_COEFFICIENTS')
    need(frozen['passes_1e28'] and frozen['full_footprint_covered'], 'RETAINED_INITIAL_GATE')
    efree = F(frozen['free_wave_X_cap'])
    need(efree < 10**28, 'RETAINED_INITIAL_CAP')
    init = load_module('frozen_interval_initial', roots['N'] / 'numeric/evaluate_initial_profiles.py')
    op = load_module('frozen_linear_operator', roots['G'] / 'projection/inputs/linear_operator.py')
    I, J, ell = init.I, init.J, init.ELL
    need(ell == F(1, 10**6) and init.GRID == 1 << 112, 'ARITHMETIC_CONVENTION')
    cs = {key: list(map(F, value)) for key, value in coefficients.items()}
    power = 1 << 40
    order = 6

    def constant(value, n=order):
        return J([value] + [0] * n)

    def polynomial(values, z):
        result = constant(0, z.n)
        for value in reversed(values):
            result = result * z + value
        return result

    def forty_squares(value):
        for _ in range(40):
            value = value * value
        return value

    def nominal_profiles(r):
        if r.hi <= F(65, 64):
            return [constant(0) for _ in range(5)], [constant(0) for _ in range(5)]
        need(r.lo >= F(65, 64) and
             not (r.lo < F(66, 64) < r.hi), 'PROFILE_CELL_CROSSES_CUTOFF')
        rr = J([r, 1 / ell] + [0] * (order - 1))
        z = ell * rr
        bg = {key: polynomial(value, z) for key, value in cs.items()}
        rho, eta, p, a, alpha = [bg[key] for key in ('rho', 'eta', 'p', 'a', 'alpha')]
        need(rho.val().lo > 0 and p.val().hi < 0, 'NOMINAL_BACKGROUND_SIGNS')
        h = init.poly(init.UP[1], eta) / p
        k = rho * (F(5, 2) * a - h)
        C, B = alpha / rho**2, 3 * alpha / (rho**3 * p)
        dr = rr - F(3, 2)
        denominator = F(1, 4) - dr * dr
        denominator.c[0] = F(1, 4) - (r - F(3, 2))**2
        t = 1 / denominator - 4
        packet = forty_squares(1 / (1 + t / power))
        if r.lo < F(66, 64):
            vcut = 64 * rr - 65
            packet = packet * (126*vcut**5 - 420*vcut**6 + 540*vcut**7
                               - 315*vcut**8 + 70*vcut**9)
        E = forty_squares(1 / (1 + z / power))
        ir, ip = 1 / rho, 1 / p
        dy = lambda value: -ir * value.dz()
        u = C * packet
        v = B * (packet.dz() + k * packet)
        vy, ut = dy(v), E * u
        d = 3 * ir * ip
        vt = E * (v - d * u)
        U2 = init.poly(init.UP[2], eta)
        aa = (-2*vy**2 - 8*p*u*vy + 8*p*(h+2*a)*u*v
              + (16*p*p-192*ir**2)*u*u - (2*U2+F(8, 3)*p*p)*v*v
              - ir**2*(12*ut*ut+2*vt*vt))
        dd = ir**2 * (36*u*u - 2*v*v)
        f = rho**2 * (aa / 6 + F(4, 3)*dd) / 6
        fy = dy(f)
        K = E*ir/3 * (18*u*u-v*v+d*u*v)
        Pi = E*ir*ip * (d*u*vy+(d*(h-ir)+10*p)*u*v
                        +(30*a-6*ir)*u*u+(2*ir-a)*v*v)
        zero = constant(0)
        displacement = [zero, zero, -u*u/3, f, -F(3, 2)*ip*fy]
        velocity = [
            ell*(2*rho*K-8*u*ut)/6,
            ell*rho*((-26*u*dy(u)-2*a*u*u)/6-F(3, 2)*fy-3*a*f),
            ell*(-8*rho*K+44*u*ut)/6,
            ell*(2*rho*K+4*u*ut)/6,
            ell*(rho*Pi+2*u*vt)/6,
        ]
        return displacement, velocity

    def profiles(r):
        u, v = nominal_profiles(r)
        values = [[item.val() for item in u],
                  [(ell*item.dz()).val() for item in u],
                  [item.val() for item in v]]
        seconds = [[2*ell**2*item.c[2] for item in u],
                   [6*ell**3*item.c[3] for item in u],
                   [2*ell**2*item.c[2] for item in v]]
        return values, seconds

    def sqrt_upper(value):
        need(value >= 0, 'NEGATIVE_NORM_SQUARE')
        if not value:
            return F(0)
        grid = 1 << 112
        return F(math.isqrt((value*grid*grid).__ceil__()) + 1, grid)

    def hnorm(values):
        return sqrt_upper(sum(F(w)*value*value for w, value in zip(WEIGHTS, values)))

    indices = [op.IDX[name + ':0'] for name in ('K', 'L', 'N', 'J', 'w')]
    need(indices == [0, 3, 6, 9, 12], 'MEAN_INDICES')
    primitive_full = op.matrices()
    for name, matrix in primitive_full.items():
        need(all(matrix[i][j] == matrix[j][i] == 0
                 for i in indices for j in range(op.N) if j not in indices),
             'MEAN_BLOCK_NOT_CLOSED:' + name)
    primitive = {name: [[matrix[i][j] for j in indices] for i in indices]
                 for name, matrix in primitive_full.items()}
    need([op.gram()[i] for i in indices] == list(map(F, WEIGHTS)), 'MEAN_GRAM')
    for name in ('Delta', 'Ca', 'Ch', 'Ephi'):
        need(not any(value for row in primitive[name] for value in row), 'MEAN_ANGULAR_TERM')

    def matrix_sum(*terms):
        return [[sum((coefficient*matrix[i][j] for coefficient, matrix in terms), I(0))
                 for j in range(5)] for i in range(5)]

    def pmat(name):
        return primitive[name]

    def actual_background(r):
        z = ell*r
        result = {}
        for key in ('rho', 'eta', 'p', 'a'):
            coeff = [I(*value) for value in frozen['shell_Taylor_coefficients'][key]]
            value = I(0)
            for entry in reversed(coeff):
                value = value*z + entry
            result[key] = value + I(*frozen['uniform_fifth_normalized_ODE_derivative'][key])*z**5
        return result

    def assemble(x, r):
        state = actual_background(r)
        rho, eta, p, a = [state[key] for key in ('rho', 'eta', 'p', 'a')]
        need(rho.lo > 0 and p.hi < 0, 'ACTUAL_BACKGROUND_SIGNS')
        t = ell*x
        # tanh(t) lies in [t-t^3/3,t] for nonnegative t; use monotonicity
        # of the lower polynomial on the tiny retained time interval.
        ht = I(t.lo-t.lo**3/3, t.hi)
        h = ht/rho
        V, V1, V2 = [init.poly(init.UP[j], eta) for j in range(3)]
        curvature_d = 1/rho**2-a*a
        curvature_b = curvature_d+p*p/3
        A = matrix_sum((-3*ell*ht, pmat('I')), (-2*ell*rho*a, pmat('O')),
                       (ell*rho*p, pmat('Aphi')))
        B = matrix_sum((-3*ell*rho*a, pmat('I')), (ell*rho*p, pmat('Bphi')))
        C = matrix_sum((a*a, pmat('Sa')), (-a*a, pmat('O2')),
                       (a*h, pmat('Sah')), (-3*a*h, pmat('O')),
                       (h*h, pmat('Sh')), (-2*curvature_b, pmat('Rb')),
                       (-2*curvature_d, pmat('Rd')), (F(4, 3)*V, pmat('Pmet')),
                       (F(4, 3)*V1, pmat('Pg')), (-V2, pmat('Pww')),
                       (-V1, pmat('W_K')), (3*p*h, pmat('W_L')))
        D = [[ell**2*rho**2*entry for entry in row] for row in C]
        M = matrix_sum((-ell*rho*a, pmat('By_a')), (-ell*rho*p, pmat('By_p')),
                       (2*ell*rho*p, pmat('By_sigma1')),
                       (-4*ell*rho*(1+eta), pmat('By_sigma2')))
        return {'A': A, 'B': B, 'D': D, 'M': M}

    local = assemble(I(0), I(0))
    global_boxes = assemble(I(0, F(3, 2)), I(0, F(7, 4)))
    global_boxes['M'] = local['M']  # Robin coefficients are shell values.
    matrices = {name: [[(value.lo+value.hi)/2 for value in row] for row in matrix]
                for name, matrix in local.items()}

    def matrix_norm(matrix):
        return sqrt_upper(sum(F(WEIGHTS[i], WEIGHTS[j])*matrix[i][j]**2
                              for i in range(5) for j in range(5)))

    matrix_caps = {name: matrix_norm(value) for name, value in matrices.items()}
    coefficient_errors = {
        name: matrix_norm([[(box[i][j]-matrices[name][i][j]).mag() for j in range(5)]
                           for i in range(5)])
        for name, box in global_boxes.items()
    }
    need(all(matrices['M'][1][j] == matrices['M'][j][1] == 0 for j in range(5)),
         'ROBIN_DIRICHLET_MIXING')
    cache = {}

    def node(r):
        if r not in cache:
            vals, unused = profiles(I(r))
            centers = [[(value.lo+value.hi)/2 for value in row] for row in vals]
            radii = [[(value.hi-value.lo)/2 for value in row] for row in vals]
            cache[r] = (centers, radii)
        return cache[r]

    level_reports = []
    selected_table = None
    for level in LEVELS:
        nodes = [F(0)] + [F(k, level) for k in range(65*level//64, 7*level//4+1)]
        need(F(65, 64) in nodes and F(66, 64) in nodes and nodes[-1] == F(7, 4),
             'MESH_COVERAGE')
        for r in nodes:
            node(r)
        errors = [F(0)]*3
        rows = []
        for left, right in zip(nodes, nodes[1:]):
            if right <= F(65, 64):
                caps = [[F(0)]*5 for _ in range(3)]
            else:
                unused, seconds = profiles(I(left, right))
                caps = [[max(node(left)[1][kind][i], node(right)[1][kind][i])
                         +(right-left)**2*seconds[kind][i].mag()/8 for i in range(5)]
                        for kind in range(3)]
            cell_norms = [hnorm(row) for row in caps]
            errors = [max(old, value) for old, value in zip(errors, cell_norms)]
            rows.append({'r': [str(left), str(right)],
                         'component_error_caps': [[str(v) for v in row] for row in caps],
                         'H_error_caps': list(map(str, cell_norms))})
        eu, edu, ev = errors
        e0 = efree+ev+SOURCE_CAP
        x0 = (1+Q)*UCAP+e0
        ea, eb, ed, em = [coefficient_errors[key] for key in ('A', 'B', 'D', 'M')]
        na, nb, nd, nm = [matrix_caps[key] for key in ('A', 'B', 'D', 'M')]
        df_coeff = (ea+eb+F(3, 4)*ed)*x0
        dg_coeff = F(1, 4)*em*x0
        df_profile = (na+nb)*edu+nd*eu
        dg_profile = nm*eu
        delta_f, delta_g = df_coeff+df_profile, dg_coeff+dg_profile
        e1 = F(3, 2)*delta_f+delta_g
        tail = Q*Q*UCAP
        total = tail+(1+Q)*e0+e1
        ledger = {'epsilon_u': eu, 'epsilon_du': edu, 'epsilon_v': ev,
                  'e_free_retained': efree, 'omitted_forcing': SOURCE_CAP,
                  'e0': e0, 'B0_X_upper': x0, 'delta_f_coefficient': df_coeff,
                  'delta_f_profile': df_profile, 'delta_g_coefficient': dg_coeff,
                  'delta_g_profile': dg_profile, 'delta_f': delta_f, 'delta_g': delta_g,
                  'e1': e1, 'two_term_tail': tail, 'total_X_error_upper': total}
        entry = {'mesh_denominator': level, 'nodes': len(nodes), 'cells': len(rows),
                 'ledger': {key: str(value) for key, value in ledger.items()},
                 'display_only': {key: float(value) for key, value in ledger.items()},
                 'passes_target': total < TARGET, 'cell_certificates': rows}
        level_reports.append(entry)
        selected_table = {'nodes': list(map(str, nodes)),
                          **{name: [[str(value) for value in node(r)[0][kind]] for r in nodes]
                             for kind, name in enumerate(('u', 'du', 'v'))}}
        print(json.dumps({'mesh_denominator': level, 'total_X_error_upper': str(total),
                          'display_only': float(total), 'passes_target': total < TARGET}),
              flush=True)
        if total < TARGET:
            break

    passed = level_reports[-1]['passes_target']

    def velocity_at_positive(r):
        for i, (left, right) in enumerate(zip(nodes, nodes[1:])):
            if left <= r <= right:
                weights = ((right-r)/(right-left), (r-left)/(right-left))
                return [weights[0]*node(left)[0][2][j]+weights[1]*node(right)[0][2][j]
                        for j in range(5)]
        raise ValueError('VELOCITY_EVALUATION_OUTSIDE_FOOTPRINT')

    def velocity_primitive_positive(r):
        result = [F(0)]*5
        for left, right in zip(nodes, nodes[1:]):
            if r <= left:
                break
            end = min(r, right)
            width = end-left
            for j in range(5):
                vl, vr = node(left)[0][2][j], node(right)[0][2][j]
                result[j] += width*vl+width**2*(vr-vl)/(2*(right-left))
        need(0 <= r <= F(7, 4), 'PRIMITIVE_OUTSIDE_FOOTPRINT')
        return result

    parity = (1, -1, 1, 1, 1)

    def reflected_initial(r):
        vals, unused = profiles(I(abs(r)))
        negative = r < 0
        u = [value*(parity[j] if negative else 1) for j, value in enumerate(vals[0])]
        du = [value*(-parity[j] if negative else 1) for j, value in enumerate(vals[1])]
        vv = velocity_at_positive(abs(r))
        pv = velocity_primitive_positive(abs(r))
        if negative:
            vv = [parity[j]*value for j, value in enumerate(vv)]
            pv = [-parity[j]*value for j, value in enumerate(pv)]
        return u, du, vv, pv

    samples = []
    for x, r in ((F(1), F(0)), (F(9, 8), F(0)), (F(5, 4), F(0)),
                 (F(3, 2), F(1, 4)), (F(3, 4), F(3, 4)), (F(1, 2), F(5, 4))):
        ua, da, va, pa = reflected_initial(r+x)
        ub, db, vb, pb = reflected_initial(r-x)
        value = [(ua[j]+ub[j]+pa[j]-pb[j])/2 for j in range(5)]
        dx = [(da[j]-db[j]+va[j]+vb[j])/2 for j in range(5)]
        dr = [(da[j]+db[j]+va[j]-vb[j])/2 for j in range(5)]
        samples.append({'x': str(x), 'r': str(r),
                        'B0': {name: [entry.row() for entry in entries]
                               for name, entries in (('value', value), ('dx', dx), ('dr', dr))}})

    profiles_output = {
        'schema': 'HDBLAST_COUPLED_RIDGE_PILOT_PROFILES_V1',
        'mean_order': ['K0', 'L0', 'N0', 'J0', 'w0'], 'gram': list(WEIGHTS),
        'ell': str(ell), 'candidate_coefficients_sha256': COEFFICIENT_PIN,
        'initial_displacement': 'exact registered rational displacement, not P1 u',
        'nodes': selected_table['nodes'],
        'profiles': {name: selected_table[name] for name in ('u', 'du', 'v')},
        'free_wave_samples': samples,
        'matrices': {name: [[str(value) for value in row] for row in matrix]
                     for name, matrix in matrices.items()},
        'parity': {'u': [1, -1, 1, 1, 1], 'du': [-1, 1, -1, -1, -1],
                     'v': [1, -1, 1, 1, 1]},
        'candidate_definition': 'B0=free(exact_registered_u,vP1); B1=K(ridge(uP1,duP1,vP1))+L(M PN(uP1+primitive(vP1)))',
    }
    report = {
        'status': 'PASS_BOUNDED_COUPLED_MEAN_X_ERROR' if passed else 'BOUNDED_PILOT_INSUFFICIENT_UPPER_BOUND',
        'passes_target': passed, 'target_X_error': str(TARGET),
        'fixed_mesh_cap': list(LEVELS), 'levels': level_reports,
        'coefficient_error_H_operator_upper': {key: str(value) for key, value in coefficient_errors.items()},
        'constant_matrix_H_operator_upper': {key: str(value) for key, value in matrix_caps.items()},
        'matrix_norm_method': 'Gram-weighted Frobenius upper bound, rational outward square root',
        'profile_error_method': 'max endpoint coefficient error plus h^2/8 times whole-cell second derivative bound; separate analytic cutoff branches',
        'input_commit_reference': OLD_COMMIT,
        'input_hashes': {'precision_manifest': PRECISION_PIN, 'candidate_coefficients': COEFFICIENT_PIN,
                        'frozen_initial_report': sha(frozen_path), 'branch_report': sha(branch),
                        **{key+'_manifest': value['manifest_sha256'] for key, value in before['closures'].items()}},
        'complete_closures_authenticated_before_and_after': True,
        'kernel_and_final_candidate_evaluation': 'separate exact reflected-ridge helper; this producer certifies its stored inputs and complete error ledger',
        'limits': ['Conditional on retained exact branch/operator/source/total-response premises.',
                   'Five mean coefficient response only; no full angular waveform.',
                   'No pointwise second-derivative Einstein residual or invariant observable evaluation.',
                   'No nonlinear family, finite-amplitude remainder or physical-origin inference.',
                   'A failed upper bound proves no lower bound on actual approximation error.'],
    }
    after = direct.authenticate(depmap)
    need(before == after, 'UPSTREAM_INPUT_CHANGED')
    need(authenticate_precision(precision) == current_before, 'PRECISION_INPUT_CHANGED')
    need(sha(branch) == direct.BRANCH_PIN, 'BRANCH_INPUT_CHANGED')
    output.mkdir(parents=True)
    for filename, data in (('PROFILES.json', profiles_output), ('REPORT.json', report)):
        with (output/filename).open('x', encoding='utf-8') as handle:
            json.dump(data, handle, indent=2, sort_keys=True)
            handle.write('\n')
    print(json.dumps({'status': report['status'], 'passes_target': passed,
                      'REPORT_sha256': sha(output/'REPORT.json'),
                      'PROFILES_sha256': sha(output/'PROFILES.json')}), flush=True)


if __name__ == '__main__':
    main()
