#!/usr/bin/env python3
"""Read-only component comparison from the authenticated frozen pilot.

This evaluates finite rational ledgers, not the old numerical campaign.
No frozen source, mesh, profile, or candidate is changed.
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
MEANS = ('K0', 'L0', 'N0', 'J0', 'w0')
WEIGHTS = (1, 2, 1, 3, 1)
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
SOURCE = F(52162547172519, 2)
UCAP = F(103 * 10**31)
Q = F(3409, 8000000)


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


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sqrt_upper(value):
    need(value >= 0, 'NEGATIVE_SQUARE')
    if not value:
        return F(0)
    grid = 1 << 112
    return F(math.isqrt((value * grid * grid).__ceil__()) + 1, grid)


def matvec(matrix, vector):
    return [sum((a*b for a, b in zip(row, vector)), F(0)) for row in matrix]


def matmul(left, right):
    return [[sum((left[i][k]*right[k][j] for k in range(len(right))), F(0))
             for j in range(len(right[0]))] for i in range(len(left))]


def add(*vectors):
    return [sum(column, F(0)) for column in zip(*vectors)]


def positive_resolvent(gamma):
    need(all(value >= 0 for row in gamma for value in row), 'NEGATIVE_COMPARISON')
    contraction = max(map(sum, gamma))
    need(contraction < 1, 'COMPARISON_ROW_CONTRACTION')
    size = len(gamma)
    identity = [[F(i == j) for j in range(size)] for i in range(size)]
    system = [[identity[i][j]-gamma[i][j] for j in range(size)] + identity[i][:]
              for i in range(size)]
    for j in range(size):
        pivot = system[j][j]
        need(pivot != 0, 'RESOLVENT_ZERO_PIVOT')
        system[j] = [entry/pivot for entry in system[j]]
        for i in range(size):
            if i != j:
                multiplier = system[i][j]
                system[i] = [a-multiplier*b for a, b in zip(system[i], system[j])]
    inverse = [row[size:] for row in system]
    need(all(entry >= 0 for row in inverse for entry in row), 'NONPOSITIVE_RESOLVENT')
    iminus = [[identity[i][j]-gamma[i][j] for j in range(size)] for i in range(size)]
    need(matmul(iminus, inverse) == identity and matmul(inverse, iminus) == identity,
         'EXACT_RESOLVENT_IDENTITY')
    return inverse, contraction


def outward(value, digits=9):
    """Positive scientific upper bound with digits significant decimal digits."""
    if value == 0:
        return '0'
    need(value > 0, 'OUTWARD_POSITIVE')
    exponent = len(str(value.numerator))-len(str(value.denominator))
    while value < F(10)**exponent:
        exponent -= 1
    while value >= F(10)**(exponent+1):
        exponent += 1
    units = F(10)**(exponent-digits+1)
    rounded = (value/units).__ceil__()
    return str(rounded) + 'e' + str(exponent-digits+1)


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def main(context_handler=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    research = args.repo_root.resolve() / 'research'
    roots = {key: research/name for key, name in LABELS.items()}
    need(all(root not in args.output.resolve().parents for root in roots.values()),
         'OUTPUT_INSIDE_FROZEN_SOURCE')
    before = authenticate(roots)  # Complete independent closure before any project import.
    native = load_module('frozen_native_ledger', roots['current']/'native_pilot.py')
    need(native.authenticate_precision(roots['I'])['manifest_sha256'] == PINS['I'],
         'NATIVE_PRECISION_AUTHENTICATION')
    direct = load_module('frozen_direct', roots['I']/'initial/direct_difference_pilot.py')
    depmap = {key: str(roots[key]) for key in ('P', 'N', 'G', 'C', 'T')}
    dependency_before = direct.authenticate(depmap)
    init = load_module('frozen_interval', roots['N']/'numeric/evaluate_initial_profiles.py')
    op = load_module('frozen_operator', roots['G']/'projection/inputs/linear_operator.py')
    I, ell = init.I, init.ELL
    need(ell == F(1, 10**6), 'ELL')
    frozen = json.loads((roots['I']/'initial/RECENTERED_CANDIDATE_NORMAL.json').read_text())
    profiles = json.loads((roots['current']/'results/PROFILES.json').read_text())
    report = json.loads((roots['current']/'results/REPORT.json').read_text())
    evaluated = json.loads((roots['current']/'results/EVALUATED_CANDIDATE.json').read_text())
    selected = report['levels'][-1]
    need(selected['mesh_denominator'] == 256 and selected['passes_target'], 'FROZEN_SELECTED_LEVEL')
    matrices = {key: [list(map(F, row)) for row in matrix]
                for key, matrix in profiles['matrices'].items()}
    ids = [op.IDX[name+':0'] for name in ('K', 'L', 'N', 'J', 'w')]
    all_primitive = op.matrices()
    primitive = {name: [[matrix[i][j] for j in ids] for i in ids]
                 for name, matrix in all_primitive.items()}
    need(ids == [0, 3, 6, 9, 12], 'MEAN_ORDER')
    need([op.gram()[i] for i in ids] == list(map(F, WEIGHTS)), 'GRAM')
    need(all(matrix[i][j] == matrix[j][i] == 0 for matrix in all_primitive.values()
             for i in ids for j in range(op.N) if j not in ids), 'MEAN_CLOSURE')

    def matrix_sum(*terms):
        return [[sum((coefficient*matrix[i][j] for coefficient, matrix in terms), I(0))
                 for j in range(5)] for i in range(5)]

    def actual_background(r):
        z = ell*r
        result = {}
        for key in ('rho', 'eta', 'p', 'a'):
            value = I(0)
            for entry in reversed(frozen['shell_Taylor_coefficients'][key]):
                value = value*z + I(*entry)
            result[key] = value+I(*frozen['uniform_fifth_normalized_ODE_derivative'][key])*z**5
        return result

    def assemble(x, r):
        state = actual_background(r)
        rho, eta, p, a = [state[key] for key in ('rho', 'eta', 'p', 'a')]
        need(rho.lo > 0 and p.hi < 0, 'BACKGROUND_SIGNS')
        t = ell*x
        ht = I(t.lo-t.lo**3/3, t.hi)
        h = ht/rho
        V, V1, V2 = [init.poly(init.UP[j], eta) for j in range(3)]
        cd = 1/rho**2-a*a
        cb = cd+p*p/3
        A = matrix_sum((-3*ell*ht, primitive['I']), (-2*ell*rho*a, primitive['O']),
                       (ell*rho*p, primitive['Aphi']))
        B = matrix_sum((-3*ell*rho*a, primitive['I']), (ell*rho*p, primitive['Bphi']))
        C = matrix_sum((a*a, primitive['Sa']), (-a*a, primitive['O2']),
                       (a*h, primitive['Sah']), (-3*a*h, primitive['O']),
                       (h*h, primitive['Sh']), (-2*cb, primitive['Rb']),
                       (-2*cd, primitive['Rd']), (F(4, 3)*V, primitive['Pmet']),
                       (F(4, 3)*V1, primitive['Pg']), (-V2, primitive['Pww']),
                       (-V1, primitive['W_K']), (3*p*h, primitive['W_L']))
        D = [[ell**2*rho**2*entry for entry in row] for row in C]
        M = matrix_sum((-ell*rho*a, primitive['By_a']), (-ell*rho*p, primitive['By_p']),
                       (2*ell*rho*p, primitive['By_sigma1']),
                       (-4*ell*rho*(1+eta), primitive['By_sigma2']))
        return {'A': A, 'B': B, 'D': D, 'M': M}

    local = assemble(I(0), I(0))
    need({name: [[(value.lo+value.hi)/2 for value in row] for row in matrix]
          for name, matrix in local.items()} == matrices, 'FROZEN_NOMINAL_MATRIX_MATCH')

    def norms(matrix):
        return sqrt_upper(sum(F(WEIGHTS[i], WEIGHTS[j])*matrix[i][j]**2
                              for i in range(5) for j in range(5)))

    def scope(label, time, footprint, zero_width):
        boxes = assemble(I(0, time), I(0, footprint))
        boxes['M'] = local['M']
        absolute = {name: [[entry.mag() for entry in row] for row in matrix]
                    for name, matrix in boxes.items()}
        errors = {name: [[(boxes[name][i][j]-matrices[name][i][j]).mag() for j in range(5)]
                         for i in range(5)] for name in matrices}
        if label == 'whole_cone':
            need({key: str(norms(value)) for key, value in errors.items()}
                 == report['coefficient_error_H_operator_upper'], 'FROZEN_COEFFICIENT_NORM_MATCH')
            need({key: str(norms(value)) for key, value in matrices.items()}
                 == report['constant_matrix_H_operator_upper'], 'FROZEN_MATRIX_NORM_MATCH')
        initial_rows = [row for row in frozen['rows'] if F(row['r'][0]) < footprint]
        initial = [[max(sum((F(row['channel_component_error_caps'][channel][kind][i])
                              for channel in range(3)), F(0)) for row in initial_rows)
                    for i in range(5)] for kind in range(3)]
        cells = [cell for cell in selected['cell_certificates'] if F(cell['r'][0]) < footprint]
        eps = [[max(F(cell['component_error_caps'][kind][i]) for cell in cells)
                for i in range(5)] for kind in range(3)]
        nodes = [k for k, r in enumerate(profiles['nodes']) if F(r) <= footprint]
        need(F(profiles['nodes'][nodes[-1]]) == footprint, 'FOOTPRINT_IS_SELECTED_KNOT')
        sup = [[max(abs(F(profiles['profiles'][name][k][i])) for k in nodes)
                for i in range(5)] for name in ('u', 'du', 'v')]
        beta = [max((sup[0][i]+eps[0][i]+time*sup[2][i])/3,
                    sup[1][i]+eps[1][i]+sup[2][i]) for i in range(5)]
        efree = [max((initial[0][i]+time*initial[2][i])/3,
                     initial[1][i]+initial[2][i]) for i in range(5)]
        source = [sqrt_upper(SOURCE**2/F(weight)) for weight in WEIGHTS]
        e0 = [max((initial[0][i]+time*(initial[2][i]+eps[2][i]))/3,
                  initial[1][i]+initial[2][i]+eps[2][i])+source[i] for i in range(5)]
        gamma = [[time*(absolute['A'][i][j]+absolute['B'][i][j]
                        +zero_width*absolute['D'][i][j])
                  +zero_width*absolute['M'][i][j] for j in range(5)] for i in range(5)]
        inverse, contraction = positive_resolvent(gamma)
        df_coeff_matrix = [[errors['A'][i][j]+errors['B'][i][j]+zero_width*errors['D'][i][j]
                            for j in range(5)] for i in range(5)]
        df_coefficient = matvec(df_coeff_matrix, beta)
        df_profile = add(matvec([[abs(matrices['A'][i][j])+abs(matrices['B'][i][j])
                                 for j in range(5)] for i in range(5)], eps[1]),
                         matvec([[abs(a) for a in row] for row in matrices['D']], eps[0]))
        dg_coefficient = [zero_width*item for item in matvec(errors['M'], beta)]
        dg_profile = matvec([[abs(a) for a in row] for row in matrices['M']], eps[0])
        e1 = [time*(a+b)+c+d for a, b, c, d in
              zip(df_coefficient, df_profile, dg_coefficient, dg_profile)]
        full = matvec(inverse, add(beta, e0))
        propagated_input = matvec(inverse, e0)
        two_feedback_tail = matvec(matmul(gamma, gamma), matvec(inverse, beta))
        E = add(propagated_input, two_feedback_tail, e1)
        identity_E = add(matvec(matmul(gamma, gamma), full), e0, matvec(gamma, e0), e1)
        need(E == identity_E, 'TWO_EXACT_LEDGER_FORMS')
        old_cap = [sqrt_upper(UCAP**2/F(weight)) for weight in WEIGHTS]
        old_route = add(matvec(matmul(gamma, gamma), old_cap), e0, matvec(gamma, e0), e1)
        mixed_full_cap = [min(a, b) for a, b in zip(full, old_cap)]
        mixed_route = add(matvec(matmul(gamma, gamma), mixed_full_cap), e0, matvec(gamma, e0), e1)
        chosen = [min(a, b, c) for a, b, c in zip(E, old_route, mixed_route)]
        sources = {'A_abs': absolute['A'], 'B_abs': absolute['B'], 'D_abs': absolute['D'],
                   'M_abs': absolute['M'], 'coefficient_entry_errors': errors}
        return {'label': label, 'time_limit': time, 'footprint': footprint,
                'zero_wedge_value_factor': zero_width, 'comparison_matrix': gamma,
                'comparison_max_row_sum': contraction, 'positive_resolvent': inverse,
                'matrix_entries': sources, 'initial_component_errors_u_du_v': initial,
                'p1_component_errors_u_du_v': eps, 'initial_sup_u_du_v': sup,
                'B0_component_X_upper': beta, 'e_free_component': efree,
                'omitted_forcing_component_upper': source, 'e0': e0,
                'delta_f_coefficient': df_coefficient, 'delta_f_profile': df_profile,
                'delta_g_coefficient': dg_coefficient, 'delta_g_profile': dg_profile,
                'e1': e1, 'full_response_component_X_upper': full,
                'resolvent_propagated_input': propagated_input,
                'resolvent_two_feedback_tail': two_feedback_tail,
                'resolvent_total_component_X_error_upper': E,
                'retained_global_U_route_error_upper': old_route,
                'mixed_full_response_component_X_upper': mixed_full_cap,
                'mixed_U_route_error_upper': mixed_route,
                'chosen_component_X_error_upper': chosen,
                'chosen_component_X_error_outward': list(map(outward, chosen))}

    ledgers = [scope('whole_cone', F(3, 2), F(7, 4), F(3, 4)),
               scope('shell_history_9_8', F(9, 8), F(9, 8), F(1, 8)),
               scope('shell_history_5_4', F(5, 4), F(5, 4), F(1, 4))]

    def rounded_upper(values):
        grid = 1 << 112
        return [F((item*grid).__ceil__(), grid) for item in values]

    def ten_from(ledger):
        tau, c = ledger['time_limit'], ledger['zero_wedge_value_factor']
        kappa = (tau*tau-1)/2
        entry = ledger['matrix_entries']
        A, B, D, M = [entry[key+'_abs'] for key in ('A', 'B', 'D', 'M')]
        C = [[A[i][j]+B[i][j] for j in range(5)] for i in range(5)]
        H = [[kappa*D[i][j]+c*M[i][j] for j in range(5)]
             +[kappa*C[i][j] for j in range(5)] for i in range(5)]
        H += [[tau*D[i][j]+M[i][j] for j in range(5)]
              +[tau*C[i][j] for j in range(5)] for i in range(5)]
        R10, contraction = positive_resolvent(H)
        eps, sup = ledger['p1_component_errors_u_du_v'], ledger['initial_sup_u_du_v']
        beta_a = [sup[0][i]+eps[0][i]+tau*sup[2][i] for i in range(5)]
        beta_d = [sup[1][i]+eps[1][i]+sup[2][i] for i in range(5)]
        beta = beta_a+beta_d
        initial = ledger['initial_component_errors_u_du_v']
        fc = [sqrt_upper(F(17387513156573)**2/F(weight)) for weight in WEIGHTS]
        gc = [sqrt_upper(F(3851400)**2/F(weight)) for weight in WEIGHTS]
        need(F(3, 2)*F(17387513156573)+F(3851400) == SOURCE, 'FORCING_SPLIT')
        e0_a = [initial[0][i]+tau*(initial[2][i]+eps[2][i])+kappa*fc[i]+c*gc[i]
                for i in range(5)]
        e0_d = [initial[1][i]+initial[2][i]+eps[2][i]+tau*fc[i]+gc[i] for i in range(5)]
        e0 = e0_a+e0_d
        delta = entry['coefficient_entry_errors']
        delta_C = [[delta['A'][i][j]+delta['B'][i][j] for j in range(5)] for i in range(5)]
        df = add(matvec(delta_C, beta_d), matvec(delta['D'], beta_a), ledger['delta_f_profile'])
        dg_coefficient = matvec(delta['M'], beta_a)
        dg_current = add(dg_coefficient, ledger['delta_g_profile'])
        def e1_for(dg):
            return [kappa*df[i]+c*dg[i] for i in range(5)]+[tau*df[i]+dg[i] for i in range(5)]
        inherited_ray = [sqrt_upper(UCAP**2/F(weight)) for weight in WEIGHTS]
        full = rounded_upper(matvec(R10, add(beta, e0)))
        full_cap = [min(a, b) for a, b in zip(full, [c*x for x in inherited_ray]+inherited_ray)]
        tail = rounded_upper(matvec(matmul(H, H), full_cap))
        input_error = rounded_upper(add(e0, matvec(H, e0)))
        current_error = rounded_upper(add(tail, input_error, e1_for(dg_current)))
        exact_g_error = rounded_upper(add(tail, input_error, e1_for(dg_coefficient)))
        return {'time': tau, 'kappa': kappa, 'value_then_derivative_order': MEANS+MEANS,
                          'comparison_matrix': H, 'comparison_max_row_sum': contraction,
                          'exact_positive_inverse_and_identity_verified': True,
                          'beta_value_derivative': beta, 'e0_value_derivative': e0,
                          'full_response_value_derivative_upper': full_cap,
                          'tail_value_derivative': tail, 'input_value_derivative': input_error,
                          'delta_f': df, 'delta_g_coefficient': dg_coefficient,
                          'current_error_value_derivative': current_error,
                          'hypothetical_exact_boundary_error_value_derivative': exact_g_error,
                          'current_error_outward': list(map(outward, current_error)),
                          'hypothetical_exact_boundary_error_outward': list(map(outward, exact_g_error))}

    if context_handler is not None:
        # Reuse the authenticated finite-ledger context in a later bounded
        # operator experiment. No old numerical producer is invoked here.
        context_handler(locals())
        return

    ten_state = [ten_from(item) for item in ledgers[1:]]

    def integral_ten_tail(absolute, left, right, time):
        C = [[absolute['A_abs'][i][j]+absolute['B_abs'][i][j] for j in range(5)] for i in range(5)]
        D, M = absolute['D_abs'], absolute['M_abs']
        zeros = [[F(0)]*5 for _ in range(5)]
        H0 = [[-D[i][j]/2-M[i][j] for j in range(5)]+[-C[i][j]/2 for j in range(5)] for i in range(5)]
        H0 += [M[i][:]+zeros[i][:] for i in range(5)]
        H1 = [M[i][:]+zeros[i][:] for i in range(5)]+[D[i][:]+C[i][:] for i in range(5)]
        H2 = [[D[i][j]/2 for j in range(5)]+[C[i][j]/2 for j in range(5)] for i in range(5)]
        H2 += [zeros[i][:]+zeros[i][:] for i in range(5)]
        outer = [D[i][:]+C[i][:] for i in range(5)]
        L0 = [M[i][:]+zeros[i][:] for i in range(5)]
        L1 = [D[i][:]+C[i][:] for i in range(5)]
        ML0, ML1 = matmul(M, L0), matmul(M, L1)
        bulk = [matmul(outer, H) for H in (H0, H1, H2)]
        coefficients = [[[time*ML0[i][j] for j in range(10)] for i in range(5)],
                        [[bulk[0][i][j]+time*ML1[i][j]-ML0[i][j] for j in range(10)] for i in range(5)],
                        [[bulk[1][i][j]-ML1[i][j] for j in range(10)] for i in range(5)],
                        bulk[2]]
        result = [[sum((matrix[i][j]*(right**(k+1)-left**(k+1))/(k+1)
                        for k, matrix in enumerate(coefficients)), F(0))
                   for j in range(10)] for i in range(5)]
        need(all(item >= 0 for row in result for item in row), 'TEN_SHELL_TAIL_POSITIVE')
        return result

    def integral_comparison(absolute, left, right, square=False, shell_tail_time=None):
        # Gamma(a)=a(A+B)+a(a-1)D+(a-1)M; polynomial coefficients
        # can be signed, but exact integrated values are nonnegative.
        C = [[absolute['A_abs'][i][j]+absolute['B_abs'][i][j]
              for j in range(5)] for i in range(5)]
        D, M = absolute['D_abs'], absolute['M_abs']
        polynomials = [[[-M[i][j] for j in range(5)] for i in range(5)],
                       [[C[i][j]-D[i][j]+M[i][j] for j in range(5)] for i in range(5)],
                       D]
        if square or shell_tail_time is not None:
            left_polynomials = polynomials
            if shell_tail_time is not None:
                # Bulk part Gamma_bulk(a) Gamma(a), plus exact nested-shell
                # lift integral (tau-a) M Gamma(a).
                left_polynomials = [[[shell_tail_time*M[i][j] for j in range(5)] for i in range(5)],
                                    [[C[i][j]-D[i][j]-M[i][j] for j in range(5)] for i in range(5)], D]
            polynomials = [[[sum((matmul(left_polynomials[a], polynomials[b])[i][j]
                                  for a in range(3) for b in range(3) if a+b == k), F(0))
                             for j in range(5)] for i in range(5)] for k in range(5)]
        result = [[sum((matrix[i][j]*(right**(k+1)-left**(k+1))/(k+1)
                        for k, matrix in enumerate(polynomials)), F(0))
                   for j in range(5)] for i in range(5)]
        need(all(entry >= 0 for row in result for entry in row), 'INTEGRATED_COMPARISON_POSITIVE')
        return result

    integrated = []
    local_cache = {item['time_limit']: item for item in ledgers[1:]}
    ten_cache = {item['time']: item for item in ten_state}
    def upward_vector(values):
        # Charge every accumulation rounding explicitly in the upper bound.
        grid = 1 << 112
        return [F((value*grid).__ceil__(), grid) for value in values]

    def velocity_error_integral(time):
        result = [F(0)]*5
        for rows, channel_sum in ((frozen['rows'], True), (selected['cell_certificates'], False)):
            for row in rows:
                left, right = map(F, row['r'])
                width = max(F(0), min(time, right)-left)
                for i in range(5):
                    cap = (sum((F(row['channel_component_error_caps'][channel][2][i])
                                for channel in range(3)), F(0)) if channel_sum else
                           F(row['component_error_caps'][2][i]))
                    result[i] += width*cap
        return result

    for time in (F(9, 8), F(5, 4)):
        boundaries = [F(1)] + [F(r) for r in profiles['nodes'] if F(1) < F(r) <= time]
        tail, input_feedback, coeff, bulk_profile, boundary_profile = ([F(0)]*5 for _ in range(5))
        sharper_tail, sharper_input_feedback = ([F(0)]*5 for _ in range(2))
        ten_tail, ten_input_feedback, ten_coefficient = ([F(0)]*5 for _ in range(3))
        for left, right in zip(boundaries, boundaries[1:]):
            if right not in local_cache:
                local_cache[right] = scope('integration_cell_upper', right, right, right-1)
            upper = local_cache[right]
            if right not in ten_cache:
                ten_cache[right] = ten_from(upper)
            ten_upper = ten_cache[right]
            entries = upper['matrix_entries']
            int_gamma = integral_comparison(entries, left, right)
            int_gamma_squared = integral_comparison(entries, left, right, square=True)
            tail = upward_vector(add(tail, matvec(int_gamma_squared, upper['mixed_full_response_component_X_upper'])))
            int_shell_tail = integral_comparison(entries, left, right, shell_tail_time=time)
            sharper_tail = upward_vector(add(sharper_tail, matvec(int_shell_tail, upper['mixed_full_response_component_X_upper'])))
            input_feedback = upward_vector(add(input_feedback, matvec(int_gamma, upper['e0'])))
            bulk_entries = dict(entries)
            bulk_entries['M_abs'] = [[F(0)]*5 for _ in range(5)]
            int_bulk_gamma = integral_comparison(bulk_entries, left, right)
            initial_cell_u = [max(sum((F(row['channel_component_error_caps'][channel][0][i])
                                      for channel in range(3)), F(0))
                                  for row in frozen['rows'] if F(row['r'][0]) < right and F(row['r'][1]) > left)
                              for i in range(5)]
            int_shell_input = add([(right-left)*(a+b) for a, b in
                                   zip(initial_cell_u, velocity_error_integral(right))],
                                  [(right-left)*(right+left-2)/2*sqrt_upper(SOURCE**2/F(weight))
                                   for weight in WEIGHTS])
            sharper_input_feedback = upward_vector(add(sharper_input_feedback,
                 matvec(int_bulk_gamma, upper['e0']), matvec(entries['M_abs'], int_shell_input)))
            full10 = ten_upper['full_response_value_derivative_upper']
            full5 = upper['mixed_full_response_component_X_upper']
            minfull10 = [min(a, b) for a, b in zip(full10, [(right-1)*x for x in full5]+full5)]
            ten_tail = upward_vector(add(ten_tail,
                matvec(integral_ten_tail(entries, left, right, time), minfull10)))
            C = [[entries['A_abs'][i][j]+entries['B_abs'][i][j] for j in range(5)] for i in range(5)]
            outer = [entries['D_abs'][i][:]+C[i][:] for i in range(5)]
            ten_input_feedback = upward_vector(add(ten_input_feedback,
                [(right*right-left*left)/2*x for x in matvec(outer, ten_upper['e0_value_derivative'])],
                matvec(entries['M_abs'], int_shell_input)))
            errors = entries['coefficient_entry_errors']
            error_outer = [errors['D'][i][:]+[errors['A'][i][j]+errors['B'][i][j] for j in range(5)] for i in range(5)]
            ten_coefficient = upward_vector(add(ten_coefficient,
                [(right*right-left*left)/2*x for x in matvec(error_outer, ten_upper['beta_value_derivative'])],
                [(right-left)*x for x in matvec(errors['M'], ten_upper['beta_value_derivative'][:5])]))
            coefficient_entries = {name+'_abs': matrix for name, matrix in entries['coefficient_entry_errors'].items()}
            int_error = integral_comparison(coefficient_entries, left, right)
            coeff = upward_vector(add(coeff, matvec(int_error, upper['B0_component_X_upper'])))
            bulk_profile = add(bulk_profile, [(right*right-left*left)/2*item
                                             for item in upper['delta_f_profile']])
            boundary_profile = add(boundary_profile, [(right-left)*item
                                                      for item in upper['delta_g_profile']])
        initial_at_time = [min(sum((F(row['channel_component_error_caps'][channel][0][i])
                                   for channel in range(3)), F(0))
                               for row in frozen['rows'] if F(row['r'][0]) <= time <= F(row['r'][1]))
                           for i in range(5)]
        velocity_integral = velocity_error_integral(time)
        direct_input = add(initial_at_time, velocity_integral,
                           [(time-1)*sqrt_upper(SOURCE**2/F(weight)) for weight in WEIGHTS])
        error_current = add(direct_input, input_feedback, tail, coeff, bulk_profile, boundary_profile)
        error_exact_g = add(direct_input, input_feedback, tail, coeff, bulk_profile)
        sharper_current = add(direct_input, sharper_input_feedback, sharper_tail, coeff, bulk_profile, boundary_profile)
        sharper_exact_g = add(direct_input, sharper_input_feedback, sharper_tail, coeff, bulk_profile)
        ten_exact_g = add(direct_input, ten_input_feedback, ten_tail, ten_coefficient, bulk_profile)
        integrated.append({'time': time, 'cells': len(boundaries)-1,
                           'direct_input': direct_input, 'input_feedback': input_feedback,
                           'two_feedback_tail': tail, 'feedback_coefficient_error': coeff,
                           'nested_shell_two_feedback_tail': sharper_tail,
                           'direct_shell_input_feedback': sharper_input_feedback,
                           'bulk_profile_error': bulk_profile, 'boundary_profile_error': boundary_profile,
                           'current_candidate_shell_value_error_upper': error_current,
                           'current_candidate_shell_value_error_outward': list(map(outward, error_current)),
                           'hypothetical_exact_displacement_boundary_shell_error_upper': error_exact_g,
                           'hypothetical_exact_displacement_boundary_shell_error_outward': list(map(outward, error_exact_g)),
                           'sharper_current_candidate_shell_error_upper': sharper_current,
                           'sharper_current_candidate_shell_error_outward': list(map(outward, sharper_current)),
                           'sharper_hypothetical_exact_boundary_shell_error_upper': sharper_exact_g,
                           'sharper_hypothetical_exact_boundary_shell_error_outward': list(map(outward, sharper_exact_g)),
                           'ten_state_integrated_tail': ten_tail,
                           'ten_state_integrated_input_feedback': ten_input_feedback,
                           'ten_state_integrated_coefficient_error': ten_coefficient,
                           'ten_state_hypothetical_exact_boundary_shell_error_upper': ten_exact_g,
                           'ten_state_hypothetical_exact_boundary_shell_error_outward': list(map(outward, ten_exact_g)),
                           'hypothetical_scope': 'Requires new candidate shell lift using exact registered u0; arithmetic/quadrature error not included here.'})
    shell = []
    for ledger in ledgers[1:]:
        time = ledger['time_limit']
        sample = next(row for row in evaluated['samples'] if F(row['x']) == time and F(row['r']) == 0)
        errors = ledger['chosen_component_X_error_upper']
        values = sample['complete_candidate_interval_evaluation']['value']
        traces = {}
        for i in (2, 3):
            lo, hi = map(F, values[i])
            error = (time-1)*errors[i]
            traces[MEANS[i]] = {'candidate_value': [lo, hi], 'shell_value_error_upper': error,
                               'shell_value_error_outward': outward(error),
                               'conditional_response_value_interval': [lo-error, hi+error],
                               'sign_certified': lo-error > 0 or hi+error < 0}
        lift_integral = [F(0)]*5
        for cell in selected['cell_certificates']:
            left, right = map(F, cell['r'])
            width = max(F(0), min(time, right)-left)
            for i in range(5):
                lift_integral[i] += width*F(cell['component_error_caps'][0][i])
        traces['isolated_J_scalar_P1_lift_error_charge'] = abs(matrices['M'][3][4])*lift_integral[4]
        shell.append({'x': time, 'r': F(0), 'traces': traces})
    need(authenticate(roots) == before, 'CLOSURE_CHANGED_AFTER')
    need(direct.authenticate(depmap) == dependency_before, 'UPSTREAM_CHANGED_AFTER')
    result = {'schema': 'HDBLAST_COMPONENT_COMPARISON_EXPLORATION_V1',
              'mean_order': MEANS, 'gram': WEIGHTS, 'authentication': before,
              'complete_closures_authenticated_before_and_after': True,
              'normalization': 'Xi=max(sup|zi|/3,sup|zi_x|,sup|zi_r|)',
              'ledger_formula': 'E <= R e0 + Gamma^2 R beta + e1; R=(I-Gamma)^(-1)>=0',
              'source_policy': 'Full retained omitted forcing charged once in e0, componentwise via Gram.',
              'scope_premises': ['Same retained exact branch/operator/source and full-response premises.',
                                 'Zero wedge x+r<=1 for exact response and continuous candidate.',
                                 'History comparison restricted to backward causal triangle x+r<=tau.',
                                 'Complete five-by-five absolute comparison retains all couplings.'],
              'ledgers': ledgers, 'shell_sample_intervals': shell,
              'integrated_shell_error_comparison': integrated,
              'value_derivative_ten_state_comparison': ten_state,
              'limits': ['Exploratory component certificate; no full intrinsic-observable sign claim.',
                         'A failed interval sign test does not show an actual sign failure.',
                         'Frozen P1 mesh and candidate unchanged; no old numerical stage rerun.']}
    with args.output.open('x') as handle:
        json.dump(encode(result), handle, sort_keys=True, indent=2)
        handle.write('\n')
    print(json.dumps({'output': str(args.output), 'sha256': sha(args.output),
                      'components_outward': {row['label']: row['chosen_component_X_error_outward']
                                             for row in ledgers},
                      'shell_signs': [{item['x'].__str__(): {key: value['sign_certified']
                                      for key, value in item['traces'].items() if isinstance(value, dict)}}
                                      for item in shell]}, sort_keys=True))


if __name__ == '__main__':
    main()
