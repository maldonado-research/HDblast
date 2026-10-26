#!/usr/bin/env python3
"""Independent exact-rational checks of the saved full23 initial ledger."""
import argparse
import hashlib
import importlib.util
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPO = Path('/workspace/hdblast-research-work/shell-observable-20261025/repo')


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok, message):
    if not ok: raise ValueError(message)


def grid_norm(caps, gram):
    square = sum((w*c*c for w, c in zip(gram, caps)), F(0))
    if square == 0: return F(0)
    denominator = 2**112
    scaled = square*denominator**2
    ceiling = -((-scaled.numerator)//scaled.denominator)
    root = math.isqrt(ceiling)+1
    answer = F(root, denominator)
    require(answer**2 >= square, 'OUTWARD_GRAM_BOUND')
    return answer


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=DEFAULT_REPO)
    parser.add_argument('--results-dir', type=Path, default=ROOT/'results')
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('NUMERIC_REPORT_REVIEW.json'))
    args = parser.parse_args()
    paths = {'producer': ROOT/'full23_initial.py',
             'report': args.results_dir/'INITIAL_REPORT.json',
             'profiles': args.results_dir/'INITIAL_PROFILES.json'}
    before_hashes = {name: digest(path) for name, path in paths.items()}
    producer = load_module('initial_producer_independent_ledger_review', paths['producer'])
    context = producer.Context(args.repo_root)
    report = json.loads(paths['report'].read_text())
    profiles = json.loads(paths['profiles'].read_text())
    gram = [F(value) for value in report['gram']]
    require(gram == context.op.gram() == list(map(F, profiles['gram'])), 'EXACT_GRAM')
    require(report['labels'] == profiles['labels'] == context.op.LABELS, 'EXACT_LABEL_ORDER')
    require(F(report['ell']) == F(profiles['ell']) == F(1, 10**6), 'ELL')
    require(report['candidate_background_coefficients_sha256'] ==
            profiles['candidate_coefficients_sha256'] == producer.COEFFICIENT_PIN, 'BACKGROUND_PIN')
    actual_rows = report['rows']
    require(len(actual_rows) == report['continuous_error_cells'] == 192, 'ACTUAL_CELL_COUNT')
    expected_cover = [(F(k, 256), F(k+1, 256)) for k in range(256, 448)]
    require([tuple(map(F, row['r'])) for row in actual_rows] == expected_cover, 'ACTUAL_GAPFREE_COVER')
    require(sum(row['analytic_flat_endpoint_tail'] for row in actual_rows) == 4, 'EXACT_FOUR_TAIL_CELLS')
    require(all(row['analytic_flat_endpoint_tail'] == (F(row['r'][1]) <= F(65, 64))
                for row in actual_rows), 'TAIL_CLASSIFICATION')
    nodes = list(map(F, profiles['nodes']))
    expected_nodes = [F(0)]+[F(k, 256) for k in range(260, 449)]
    require(nodes == expected_nodes and len(nodes) == report['candidate_nodes'] == 190, 'P1_NODES')
    interpolation = report['interpolation_cells']
    require(len(interpolation) == 189, 'P1_CELL_COUNT')
    require([tuple(map(F, row['r'])) for row in interpolation] == list(zip(nodes, nodes[1:])), 'P1_GAPFREE_COVER')
    require(F(65, 64) in nodes and F(66, 64) in nodes, 'ALIGNED_CUTOFF_NODES')
    require(not any(F(row['r'][0]) < edge < F(row['r'][1])
                    for row in interpolation for edge in (F(65, 64), F(66, 64))), 'NO_CUTOFF_CROSSING')
    channel_components = [[[F(0)]*23 for _ in range(3)] for _ in range(3)]
    mean_error_checks = 0
    old_rows = context.frozen['rows']
    for row in actual_rows:
        left, right = map(F, row['r'])
        wide = [old for old in old_rows if F(old['r'][0]) <= left and right <= F(old['r'][1])]
        require(len(wide) == 1, 'UNIQUE_FROZEN_WIDER_CELL')
        errors = row['channel_component_error_caps']
        require(len(errors) == 3 and all(len(channel) == 3 for channel in errors), 'ERROR_CHANNEL_KIND_SHAPE')
        for channel in range(3):
            for kind in range(3):
                caps = list(map(F, errors[channel][kind]))
                require(len(caps) == 23 and min(caps) >= 0, 'ERROR_COMPONENT_SHAPE_SIGN')
                channel_components[channel][kind] = [max(a, b) for a, b in zip(channel_components[channel][kind], caps)]
                for j, i in enumerate(context.mean_ids):
                    require(caps[i] <= F(wide[0]['channel_component_error_caps'][channel][kind][j]), 'FROZEN_MEAN_DOMINATION')
                    mean_error_checks += 1
    require(channel_components == [[[F(x) for x in row] for row in channel]
                                   for channel in report['channel_component_error_caps']], 'CHANNEL_GLOBAL_MAXIMA')
    channel_norms = [[grid_norm(row, gram) for row in channel] for channel in channel_components]
    require(channel_norms == [list(map(F, row)) for row in report['channel_H_caps']], 'CHANNEL_GRAM_NORMS')
    total = [sum((row[kind] for row in channel_norms), F(0)) for kind in range(3)]
    require(total == list(map(F, report['total_H_caps'])), 'THREE_CHANNEL_TRIANGLE_SUM')
    mean_node_checks = 0
    for kind, name in enumerate(('u', 'du', 'v')):
        table = profiles['profiles'][name]
        old_table = context.frozen_profiles['profiles'][name]
        require(len(table) == 190, 'P1_TABLE_NODE_COUNT')
        for row, old in zip(table, old_table):
            require(len(row) == 23, 'P1_TABLE_COMPONENT_COUNT')
            for j, i in enumerate(context.mean_ids):
                require(F(row[i]) == F(old[j]), 'FROZEN_MEAN_NODE_EQUALITY')
                mean_node_checks += 1
    p1_components = [[F(0)]*23 for _ in range(3)]
    second_components = [[F(0)]*23 for _ in range(3)]
    endpoint_residuals = {}
    cell_norm_checks = 0
    for row in interpolation:
        left, right = map(F, row['r'])
        h = right-left
        for kind in range(3):
            caps = list(map(F, row['component_error_caps'][kind]))
            seconds = list(map(F, row['second_derivative_component_caps'][kind]))
            require(len(caps) == len(seconds) == 23 and min(caps+seconds) >= 0, 'P1_CAP_SHAPE_SIGN')
            require(grid_norm(caps, gram) == F(row['H_error_caps'][kind]), 'P1_CELL_GRAM_NORM')
            residual = [cap-h*h*derivative/8 for cap, derivative in zip(caps, seconds)]
            require(min(residual) >= 0, 'P1_CURVATURE_TERM_RETAINED')
            if right <= F(65, 64): require(not any(caps+seconds), 'NOMINAL_ZERO_REGION_EXACT')
            endpoint_residuals[(left, right, kind)] = residual
            p1_components[kind] = [max(a, b) for a, b in zip(p1_components[kind], caps)]
            second_components[kind] = [max(a, b) for a, b in zip(second_components[kind], seconds)]
            cell_norm_checks += 1
    require(p1_components == [list(map(F, row)) for row in report['P1_component_error_caps']], 'P1_GLOBAL_MAXIMA')
    p1_norms = [grid_norm(row, gram) for row in p1_components]
    require(p1_norms == list(map(F, report['P1_H_error_caps'])), 'P1_GLOBAL_GRAM_NORMS')
    require(second_components == [list(map(F, row)) for row in report['nominal_second_derivative_component_caps']], 'SECOND_DERIVATIVE_GLOBAL_MAXIMA')
    require([grid_norm(row, gram) for row in second_components] ==
            list(map(F, report['nominal_second_derivative_H_caps'])), 'SECOND_DERIVATIVE_GRAM_NORMS')
    require([grid_norm(list(map(F, row)), gram) for row in report['nominal_component_caps']] ==
            list(map(F, report['nominal_H_caps'])), 'NOMINAL_GRAM_NORMS')
    # The endpoint radii are not stored.  Review the producer loop and verify
    # selected adjacent cells by reevaluating only their endpoint enclosures.
    # This is separate from the 192-cell paired error calculation.
    source = paths['producer'].read_text()
    required_formula = 'max(node_cache[left][1][kind][i], node_cache[right][1][kind][i])\n                     +(right-left)**2*seconds[kind][i].mag()/8'
    require(required_formula in source, 'PRODUCER_P1_FORMULA')
    require('max(abs(box.lo-center), abs(box.hi-center))' in source, 'PRODUCER_ENDPOINT_RADIUS_FORMULA')
    selected_cells = [(F(65, 64), F(261, 256)), (F(263, 256), F(66, 64)),
                      (F(66, 64), F(265, 256)), (F(383, 256), F(3, 2)),
                      (F(447, 256), F(7, 4))]
    evaluated_nodes = {}
    endpoint_formula_checks = 0
    for left, right in selected_cells:
        for r in (left, right):
            if r not in evaluated_nodes:
                values, unused = context.profiles(context.I(r))
                centers = [profiles['profiles'][name][nodes.index(r)] for name in ('u', 'du', 'v')]
                evaluated_nodes[r] = [[max(abs(box.lo-F(center)), abs(box.hi-F(center)))
                                       for box, center in zip(row, centers[kind])]
                                      for kind, row in enumerate(values)]
        for kind in range(3):
            expected = [max(a, b) for a, b in zip(evaluated_nodes[left][kind], evaluated_nodes[right][kind])]
            require(expected == endpoint_residuals[(left, right, kind)], 'SELECTED_EXACT_ENDPOINT_RADIUS_CHARGE')
            endpoint_formula_checks += 23
    require(mean_error_checks == 8640 and mean_node_checks == 2850, 'MEAN_CHECK_COUNTS')
    require(report['mean_consistency']['entries_compared'] == 8640 and
            report['mean_consistency']['new_error_entries_dominated_by_frozen_wider_cells'] == 8640 and
            report['mean_consistency']['retained_mean_node_centers_reused'] == 2850, 'REPORTED_MEAN_COUNTS')
    require(total[2] > 0 and p1_norms[2] > 0, 'SEPARATE_NONZERO_VELOCITY_ERRORS_RETAINED')
    require(report['derivative_convention']['v'] == 'partial_x U; x=t/ell', 'VELOCITY_CONVENTION')
    context.reauthenticate()
    require(before_hashes == {name: digest(path) for name, path in paths.items()}, 'REVIEW_INPUT_CHANGED')
    result = {
        'status': 'PASS_INDEPENDENT_FULL23_INITIAL_NUMERIC_LEDGER_REVIEW',
        'reviewed_sha256': before_hashes,
        'actual_error_cells_gapfree': 192, 'actual_cover': ['1', '7/4'], 'actual_exact_zero_cover': ['0', '1'],
        'nominal_P1_cells_gapfree': 189, 'nominal_P1_nodes': 190, 'nominal_cover': ['0', '7/4'],
        'cutoff_nodes_aligned': ['65/64', '66/64'], 'analytic_actual_tail_cells': 4,
        'frozen_mean_error_domination_checks': mean_error_checks,
        'frozen_mean_node_equality_checks': mean_node_checks,
        'P1_cell_Gram_norms_recomputed': cell_norm_checks,
        'channel_global_maxima_and_Gram_norms_recomputed': True,
        'total_initial_error_triangle_sums_recomputed': True,
        'P1_global_maxima_and_Gram_norms_recomputed': True,
        'nominal_and_second_derivative_Gram_norms_recomputed': True,
        'P1_curvature_term_h_squared_over_eight_verified_all_cells': True,
        'P1_endpoint_radius_formula_verified_by_source_and_selected_exact_nodes': endpoint_formula_checks,
        'P1_selected_endpoint_cells': [[str(a), str(b)] for a, b in selected_cells],
        'velocity_input_error_and_P1_error_separately_retained': True,
        'total_initial_H_caps_display': list(map(float, total)),
        'P1_H_error_caps_display': list(map(float, p1_norms)),
        'complete_roots_authenticated_before_and_after': True,
        'limits': ['No complete numerical campaign rerun; saved exact ledgers recomputed and selected endpoint evaluations checked.',
                   'Endpoint radii are not all stored; all-cell formula reviewed in producer source, with345 exact selected component checks.',
                   'Parity review handled separately.'],
    }
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
