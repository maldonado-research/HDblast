#!/usr/bin/env python3
"""Independent exact ledger/matrix/sample review of the full23 response."""
import argparse
import hashlib
import importlib.util
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
COMPUTE = Path(__file__).resolve().parents[1]
ROOT = COMPUTE.parent
DEFAULT_REPO = Path('/workspace/hdblast-research-work/shell-observable-20261025/repo')


def require(ok, message):
    if not ok: raise ValueError(message)


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sqrtup(square):
    if square == 0: return F(0)
    require(square > 0, 'NEGATIVE_SQUARE')
    grid = 2**112
    scaled = square*grid*grid
    upper_integer = -((-scaled.numerator)//scaled.denominator)
    result = F(math.isqrt(upper_integer)+1, grid)
    require(result*result >= square, 'OUTWARD_SQUARE_ROOT')
    return result


def p1_at(nodes, values, parity, argument, primitive=False):
    negative = argument < 0
    endpoint = abs(argument)
    require(0 <= endpoint <= nodes[-1], 'P1_ARGUMENT_DOMAIN')
    out = [F(0)]*23
    for a, b, va, vb in zip(nodes, nodes[1:], values, values[1:]):
        if endpoint < a: break
        end = min(b, endpoint)
        width = end-a
        if primitive:
            for i in range(23):
                out[i] += width*va[i]+width*width*(vb[i]-va[i])/(2*(b-a))
        elif a <= endpoint <= b:
            out = [((b-endpoint)*u+(endpoint-a)*v)/(b-a) for u, v in zip(va, vb)]
            break
        if endpoint <= b: break
    if negative:
        out = [(sign*(-1 if primitive else 1))*v for sign, v in zip(parity, out)]
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=DEFAULT_REPO)
    parser.add_argument('--results-dir', type=Path, default=COMPUTE/'response-normal')
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('RESPONSE_REVIEW.json'))
    args = parser.parse_args()
    paths = {'producer': COMPUTE/'full23_response.py', 'kernel': COMPUTE/'angular_ridge.py',
             'initial_producer': ROOT/'initial/full23_initial.py',
             'initial_report': ROOT/'initial/results/INITIAL_REPORT.json',
             'initial_profiles': ROOT/'initial/results/INITIAL_PROFILES.json',
             'response_report': args.results_dir/'RESPONSE_REPORT.json',
             'response_profiles': args.results_dir/'RESPONSE_PROFILES.json',
             'evaluated': args.results_dir/'EVALUATED_CANDIDATE.json'}
    hashes = {name: sha(path) for name, path in paths.items()}
    initial = load('initial_response_audit_context', paths['initial_producer'])
    context = initial.Context(args.repo_root)
    report = json.loads(paths['response_report'].read_text())
    profiles = json.loads(paths['response_profiles'].read_text())
    evaluated = json.loads(paths['evaluated'].read_text())
    initial_report = json.loads(paths['initial_report'].read_text())
    initial_profiles = json.loads(paths['initial_profiles'].read_text())
    gram = context.weights
    require(report['labels'] == profiles['labels'] == evaluated['labels'] == context.op.LABELS, 'LABELS')
    require(list(map(F, report['gram'])) == list(map(F, profiles['gram'])) == gram, 'GRAM')
    require(F(report['ell']) == F(profiles['ell']) == context.ell, 'ELL')
    for name, item in report['local_dependency_hashes'].items():
        require(sha(ROOT/item['path']) == item['sha256'] and
                (ROOT/item['path']).stat().st_size == item['bytes'], 'PINNED_LOCAL_INPUT:'+name)
    require(profiles['nodes'] == initial_profiles['nodes'] and profiles['profiles'] == initial_profiles['profiles'], 'INITIAL_TABLES_UNCHANGED')
    # Verify that the generic full23 kernel is precisely the declared
    # dimensional extension of the authenticated five-mean retarded kernel.
    transform = json.loads((COMPUTE/'KERNEL_TRANSFORMATION.json').read_text())
    old_kernel_path = context.roots['current']/'reflected_ridge.py'
    require(sha(old_kernel_path) == transform['source_sha256'], 'ORIGINAL_KERNEL_HASH')
    require(hashes['kernel'] == transform['new_source_sha256'], 'NEW_KERNEL_HASH')
    transformed = old_kernel_path.read_text()
    for change in transform['changes']:
        require(transformed.count(change['old']) == change['occurrences'], 'KERNEL_REPLACEMENT_COUNT')
        transformed = transformed.replace(change['old'], change['new'])
    require(transformed == paths['kernel'].read_text(), 'KERNEL_DECLARED_TRANSFORMATION_EXACT')
    kernel = load('response_review_full23_kernel', paths['kernel'])
    require(kernel.GRAM == tuple(gram), 'KERNEL_GRAM')
    blocks = [context.mean_ids,
              [context.op.IDX[name+':2'] for name in ('K', 'L', 'N', 'J', 'w', 'v', 'u', 'H')],
              [context.op.IDX[name+':4'] for name in ('K', 'L', 'N', 'J', 'w', 'v', 'u', 'H')],
              [context.op.IDX['T:+']], [context.op.IDX['T:-']]]
    block_names = ('degree0', 'degree2', 'degree4', 'TTplus', 'TTminus')
    require(sorted(i for block in blocks for i in block) == list(range(23)), 'ORTHOGONAL_BLOCK_PARTITION')
    require(report['orthogonal_blocks'] == {name: [context.op.LABELS[i] for i in block]
                                           for name, block in zip(block_names, blocks)}, 'REPORTED_BLOCK_LABELS')
    block_of = {i: n for n, block in enumerate(blocks) for i in block}
    for name, matrix in context.op.matrices().items():
        require(all(matrix[i][j] == 0 for i in range(23) for j in range(23)
                    if block_of[i] != block_of[j]), 'PRIMITIVE_BLOCK_CLOSURE:'+name)
    matrices = {name: [list(map(F, row)) for row in matrix] for name, matrix in report['constant_matrices'].items()}
    errors = {name: [list(map(F, row)) for row in matrix] for name, matrix in report['coefficient_entry_error_caps'].items()}
    matrix_caps, error_caps = {}, {}
    entry_checks = mean_entries = 0
    for name in ('A', 'B', 'D', 'M'):
        require(profiles['matrices'][name] == report['constant_matrices'][name], 'PROFILE_MATRIX_EQUALITY')
        for i in range(23):
            for j in range(23):
                lo, hi = map(F, report['actual_coefficient_boxes'][name][i][j])
                require(lo <= hi and errors[name][i][j] >= 0, 'MATRIX_BOX_AND_ERROR_SIGN')
                require(errors[name][i][j] >= max(abs(lo-matrices[name][i][j]), abs(hi-matrices[name][i][j])), 'ENTRY_DIFFERENCE_COVER')
                if block_of[i] != block_of[j]:
                    require(lo == hi == matrices[name][i][j] == errors[name][i][j] == 0, 'ACTUAL_CANDIDATE_ERROR_BLOCK_CLOSURE')
                entry_checks += 1
        for row_kind, array, block_key, total_key, store in (
                ('nominal', matrices[name], 'constant_matrix_block_H_operator_upper', 'constant_matrix_H_operator_upper', matrix_caps),
                ('error', errors[name], 'coefficient_error_block_H_operator_upper', 'coefficient_error_H_operator_upper', error_caps)):
            norms = [sqrtup(sum((gram[i]/gram[j]*array[i][j]**2 for i in block for j in block), F(0)))
                     for block in blocks]
            require(norms == [F(report[block_key][name][block_name]) for block_name in block_names], 'BLOCK_GRAM_FROBENIUS:'+row_kind+name)
            require(max(norms) == F(report[total_key][name]), 'MAXIMUM_BLOCK_NORM:'+row_kind+name)
            store[name] = max(norms)
        for i0, i in enumerate(context.mean_ids):
            for j0, j in enumerate(context.mean_ids):
                require(matrices[name][i][j] == F(context.frozen_profiles['matrices'][name][i0][j0]), 'FROZEN_MEAN_MATRIX_ENTRY')
                mean_entries += 1
    require(all(matrices['M'][i][j] == errors['M'][i][j] == 0
                for i in range(23) for j in range(23) if kernel.PARITY[i] == -1 or kernel.PARITY[j] == -1), 'ROBIN_DIRICHLET_ROWS_AND_COLUMNS')
    # Check the saved coefficient boxes against a second exact implementation
    # of assembly at a rational interior point of each retained background box,
    # at time0 and at a positive-time rational point in the valid rectangular
    # tanh/sech enclosures. The latter exercises all time/angular cross terms.
    # These are consistency controls, not a new background solve.
    coefficient_point_checks = 0
    for x, r in ((F(0), F(0)), (F(0), F(7, 4)),
                 (F(3, 2), F(0)), (F(3, 2), F(7, 4))):
        state = context.actual_background(context.I(r))
        rho, eta, p, a = [(state[name].lo+state[name].hi)/2 for name in ('rho', 'eta', 'p', 'a')]
        potential = [context.init.poly(context.init.UP[j], eta) for j in range(3)]
        V, V1, V2 = potential
        d = 1/rho**2-a*a
        t = context.ell*x
        ht, invc = t-t**3/6, 1-t**2/4
        assembled = context.op.assemble(context.ell, rho, a, ht/rho, 1/invc, ht, p,
                                         V1-4*a*p, d+p*p/3, d, V, V1, V2,
                                         -2*p, 4*(1+eta))
        for name in ('A', 'B', 'D', 'M'):
            if name == 'M' and r != 0: continue
            for i in range(23):
                for j in range(23):
                    lo, hi = map(F, report['actual_coefficient_boxes'][name][i][j])
                    require(lo <= assembled[name][i][j] <= hi, 'INDEPENDENT_ASSEMBLY_POINT_CONTAINMENT')
                    coefficient_point_checks += 1
    # Recompute the complete native ledger independently.
    nu, ndu, nv = map(F, initial_report['nominal_H_caps'])
    du, ddu, dv = map(F, initial_report['total_H_caps'])
    eu, edu, ev = map(F, initial_report['P1_H_error_caps'])
    q, source = F(3409, 8000000), F(52162547172519, 2)
    beta = max((nu+F(3, 2)*nv)/3, ndu+nv)
    x0 = max((nu+F(3, 2)*(nv+ev))/3, ndu+nv+ev)
    efree = max((du+F(3, 2)*dv)/3, ddu+dv)
    e0 = efree+ev+source
    ucap = (x0+e0)/(1-q)
    da, db, dd, dm = [error_caps[name] for name in ('A', 'B', 'D', 'M')]
    na, nb, nd, nm = [matrix_caps[name] for name in ('A', 'B', 'D', 'M')]
    dfc, dfp = (da+db+F(3, 4)*dd)*x0, (na+nb)*edu+nd*eu
    dgc, dgp = F(1, 4)*dm*x0, nm*eu
    df, dg = dfc+dfp, dgc+dgp
    e1, tail = F(3, 2)*df+dg, q*q*ucap
    total = tail+(1+q)*e0+e1
    point_factor = sqrtup(F(321, 19))
    expected = {
        'nominal_u_H': nu, 'nominal_du_H': ndu, 'nominal_v_H': nv,
        'delta_u_H': du, 'delta_du_H': ddu, 'delta_v_H': dv,
        'epsilon_u_H': eu, 'epsilon_du_H': edu, 'epsilon_v_H': ev,
        'beta_nominal_free_X_upper': beta, 'B0_X_upper': x0,
        'retained_full_H_contraction': q, 'free_input_error_X_upper': efree,
        'omitted_full_H_forcing_X_upper': source, 'e0': e0,
        'actual_response_X_upper': ucap, 'delta_f_coefficient': dfc,
        'delta_f_profile': dfp, 'delta_g_coefficient': dgc,
        'delta_g_profile': dgp, 'delta_f': df, 'delta_g': dg, 'e1': e1,
        'two_term_tail': tail, 'total_X_error_upper': total,
        'selected_evaluation_factor_upper': point_factor,
        'zero_wedge_pointwise_value_error_upper': F(3, 4)*point_factor*total,
        'pointwise_scaled_derivative_error_upper': point_factor*total,
    }
    require(expected == {key: F(value) for key, value in report['ledger'].items()}, 'COMPLETE_LEDGER_EXACT_EQUALITY')
    require(report['passes_target'] == evaluated['accuracy_target_passed'] == (total < F(10**29)), 'TARGET_GATE')
    require(F(evaluated['total_X_error_upper']) == total and F(report['target_X_error']) == 10**29, 'TARGET_VALUE')
    actual_caps = {name: matrix_caps[name]+error_caps[name] for name in matrix_caps}
    q_concrete = F(3, 2)*(actual_caps['A']+actual_caps['B']+F(3, 4)*actual_caps['D'])+actual_caps['M']/4
    require(0 <= q_concrete <= q < 1, 'CONCRETE_FULL23_CONTRACTION_JOIN')
    candidate = kernel.from_data(profiles)
    nodes = list(map(F, profiles['nodes']))
    velocities = [list(map(F, row)) for row in profiles['profiles']['v']]
    old_samples = json.loads((context.roots['current']/'results/EVALUATED_CANDIDATE.json').read_text())['samples']
    sample_feedback_checks = mean_feedback_checks = free_checks = component_checks = tt_checks = 0
    point_cache = {}
    for k, sample in enumerate(evaluated['samples']):
        x, r = F(sample['x']), F(sample['r'])
        require((x, r) == (F(old_samples[k]['x']), F(old_samples[k]['r'])), 'FROZEN_SAMPLE_LOCATIONS')
        recomputed_feedback = candidate.feedback(x, r, matrices['M'])
        ua, da, ub, dbb = [], [], [], []
        for argument, uout, dout in ((r+x, ua, da), (r-x, ub, dbb)):
            if abs(argument) not in point_cache:
                point_cache[abs(argument)] = context.profiles(context.I(abs(argument)))[0]
            values = point_cache[abs(argument)]
            for i, parity in enumerate(kernel.PARITY):
                uout.append(values[0][i]*(parity if argument < 0 else 1))
                dout.append(values[1][i]*(-parity if argument < 0 else 1))
        va = p1_at(nodes, velocities, kernel.PARITY, r+x)
        vb = p1_at(nodes, velocities, kernel.PARITY, r-x)
        pa = p1_at(nodes, velocities, kernel.PARITY, r+x, True)
        pb = p1_at(nodes, velocities, kernel.PARITY, r-x, True)
        independent_free = {
            'value': [(ua[i]+ub[i]+pa[i]-pb[i])/2 for i in range(23)],
            'dx': [(da[i]-dbb[i]+va[i]+vb[i])/2 for i in range(23)],
            'dr': [(da[i]+dbb[i]+va[i]-vb[i])/2 for i in range(23)]}
        for name in ('value', 'dx', 'dr'):
            for i in range(23):
                feedback = F(sample['feedback_exact'][name][i])
                require(feedback == recomputed_feedback[name][i], 'FULL23_FEEDBACK_EXACT_REPLAY')
                require(independent_free[name][i].row() == sample['free_wave_interval'][name][i], 'INDEPENDENT_P1_FREE_WAVE')
                complete = independent_free[name][i]+feedback
                require(complete.row() == sample['complete_candidate_interval_evaluation'][name][i], 'CANDIDATE_INTERVAL_ADDITION')
                cap = (3 if name == 'value' else 1)*total*sqrtup(1/gram[i])
                require((complete+context.I(-cap, cap)).row() == sample['conditional_exact_response_component_enclosure'][name][i], 'CONDITIONAL_COMPONENT_ERROR_ENCLOSURE')
                sample_feedback_checks += 1; free_checks += 1; component_checks += 1
                if i in (21, 22):
                    require(feedback == complete.lo == complete.hi == 0, 'TT_CANDIDATE_EXACT_SAMPLE_ZERO')
                    tt_checks += 1
            for j, i in enumerate(context.mean_ids):
                require(F(sample['feedback_exact'][name][i]) == F(old_samples[k]['feedback_exact'][name][j]), 'EXACT_FROZEN_MEAN_FEEDBACK')
                mean_feedback_checks += 1
        for dx_box, dt_row in zip(sample['complete_candidate_interval_evaluation']['dx'], sample['complete_candidate_interval_evaluation']['ordinary_t_velocity']):
            require((context.I(*dx_box)/context.ell).row() == dt_row, 'ORDINARY_TIME_SCALING')
    require(len(evaluated['samples']) == evaluated['sample_count'] == 6, 'SIX_SAMPLES')
    require(report['regression']['frozen_mean_constant_entries_reused'] == mean_entries == 100, 'MEAN_MATRIX_COUNT')
    require(report['regression']['exact_mean_feedback_coefficients_equal'] == mean_feedback_checks == 90, 'MEAN_FEEDBACK_COUNT')
    # The actual TT solution need not vanish: full-H source omission is charged.
    for name in ('u', 'du', 'v'):
        require(all(F(row[i]) == 0 for row in profiles['profiles'][name] for i in (21, 22)), 'TT_INITIAL_EXACT_ZERO')
    require(all(matrices[name][i][j] == 0 for name in matrices for i in (21, 22) for j in range(21)), 'TT_AUTONOMOUS_CANDIDATE_BLOCKS')
    context.reauthenticate()
    require(hashes == {name: sha(path) for name, path in paths.items()}, 'REVIEWED_INPUT_CHANGED')
    result = {
        'status': 'PASS_INDEPENDENT_FULL23_RESPONSE_MATRIX_LEDGER_AND_SAMPLE_REVIEW',
        'reviewed_sha256': hashes,
        'coefficient_entry_error_and_block_closure_checks': entry_checks,
        'matrix_norms_recomputed': '40 weighted block Frobenius norms,8 maximum block operator bounds',
        'independent_exact_operator_assembly_point_containments': coefficient_point_checks,
        'frozen_mean_matrix_entries_equal': mean_entries,
        'complete_ledger_exact_entries_recomputed': len(expected),
        'concrete_combined_actual_operator_q_upper': str(q_concrete),
        'concrete_combined_actual_operator_q_display': float(q_concrete),
        'retained_q_upper': str(q),
        'concrete_q_is_at_most_retained_q': True,
        'constant_plus_error_operator_caps': {name: str(value) for name, value in actual_caps.items()},
        'separate_B_upper_may_exceed_retained_individual_cap': actual_caps['B'] > F(1, 10000),
        'contraction_note': 'The combined coefficient bound verifies q directly; an overestimated individual B Frobenius cap is not a contradiction.',
        'kernel_exact_declared_dimensional_transformation_verified': True,
        'full23_exact_feedback_coefficients_replayed': sample_feedback_checks,
        'independently_integrated_P1_free_wave_components_checked': free_checks,
        'conditional_response_component_enclosures_checked': component_checks,
        'frozen_mean_feedback_values_equal': mean_feedback_checks,
        'TT_candidate_sample_components_exact_zero': tt_checks,
        'TT_candidate_global_zero_follows_from_zero_inputs_and_autonomous_blocks': True,
        'actual_TT_zero_not_asserted': True,
        'actual_response_X_upper_display': float(ucap),
        'total_X_error_upper_display': float(total),
        'target_X_error': '100000000000000000000000000000',
        'passes_target': total < 10**29,
        'pointwise_derivative_scope': 'C*E applies to the reconstructed coordinate-scaled derivatives partial_x and partial_r. Angular covariant derivatives require their finite-basis constants and frame factors.',
        'complete_seven_roots_authenticated_before_and_after': True,
        'proof_review': 'ERROR_IDENTITY_REVIEW.txt',
        'limits': ['Conditional on retained exact branch/operator/source/formal physical-completion premises.',
                   'No full source campaign or initial192-cell campaign rerun.',
                   'Exact ridge feedback replay uses the validated generic kernel; independent transformation checks show all integral arithmetic matches the authenticated five-mean helper.',
                   'No C2 pointwise Einstein residual or nonlinear family claim.'],
    }
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
