#!/usr/bin/env python3
"""Review support, frozen-mean reuse, exact zero channels, and cover ledgers."""
import argparse
import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = Path('/workspace/hdblast-research-work/shell-observable-20261025/repo')


def need(value, message):
    if not value:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=REPO)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('initial23_output_producer', ROOT/'full23_initial.py')
    producer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(producer)
    context = producer.Context(args.repo_root)
    report_path = ROOT/'results/INITIAL_REPORT.json'
    profiles_path = ROOT/'results/INITIAL_PROFILES.json'
    report = json.loads(report_path.read_text())
    profiles = json.loads(profiles_path.read_text())
    labels = context.op.LABELS
    expected_u = [name for name in labels if name.split(':')[0] in ('N', 'J', 'w')]
    expected_v = [name for name in labels if name.split(':')[0] in ('K', 'L', 'N', 'J', 'w', 'u')]
    need([name for name, cap in zip(labels, report['nominal_component_caps'][0]) if F(cap)] == expected_u,
         'NOMINAL_DISPLACEMENT_SUPPORT')
    need([name for name, cap in zip(labels, report['nominal_component_caps'][2]) if F(cap)] == expected_v,
         'NOMINAL_VELOCITY_SUPPORT')
    exact_zero_checks = 0
    for row in report['rows']:
        for c in range(3):
            for kind in range(3):
                allowed = expected_v if kind == 2 else expected_u
                for name, cap in zip(labels, row['channel_component_error_caps'][c][kind]):
                    if name not in allowed:
                        need(F(cap) == 0, 'ERROR_ZERO_CHANNEL:'+name)
                        exact_zero_checks += 1
    rows = report['rows']
    need(len(rows) == 192 and F(rows[0]['r'][0]) == 1 and F(rows[-1]['r'][1]) == F(7, 4), 'ACTUAL_COVER')
    need(all(a['r'][1] == b['r'][0] for a, b in zip(rows, rows[1:])), 'ACTUAL_COVER_GAP')
    cells = report['interpolation_cells']
    need(len(cells) == 189 and cells[0]['r'][0] == '0' and cells[-1]['r'][1] == '7/4', 'P1_COVER')
    need(all(a['r'][1] == b['r'][0] for a, b in zip(cells, cells[1:])), 'P1_COVER_GAP')
    for row in rows+cells:
        left, right = map(F, row['r'])
        need(not any(left < junction < right for junction in (F(65, 64), F(66, 64))), 'CUTOFF_CROSSING')
    mean_checks = 0
    for r, old_r in zip(profiles['nodes'], context.frozen_profiles['nodes']):
        need(r == old_r, 'FROZEN_NODE_LOCATION')
    for name in producer.KINDS:
        for new, old in zip(profiles['profiles'][name], context.frozen_profiles['profiles'][name]):
            for i, retained in zip(context.mean_ids, old):
                need(new[i] == retained, 'FROZEN_MEAN_NODE_CENTER')
                mean_checks += 1
    dominated = 0
    for row in rows:
        left, right = map(F, row['r'])
        old = next(old for old in context.frozen['rows'] if F(old['r'][0]) <= left and right <= F(old['r'][1]))
        for c in range(3):
            for kind in range(3):
                for j, i in enumerate(context.mean_ids):
                    need(F(row['channel_component_error_caps'][c][kind][i]) <= F(old['channel_component_error_caps'][c][kind][j]),
                         'FROZEN_WIDER_CELL_DOMINANCE')
                    dominated += 1
    # Re-evaluate a small set of nodes; compare nonmean P1 coefficients with
    # the expression midpoint and independently retained means with old data.
    selected = [F(0), F(65, 64), F(66, 64), F(9, 8), F(5, 4), F(7, 4)]
    node_checks = 0
    for r in selected:
        index = list(map(F, profiles['nodes'])).index(r)
        values, _ = context.profiles(context.I(r))
        for kind, name in enumerate(producer.KINDS):
            for i, box in enumerate(values[kind]):
                center = F(profiles['profiles'][name][index][i])
                if i not in context.mean_ids:
                    need(center == (box.lo+box.hi)/2, 'NOMINAL_NODE_EXPRESSION')
                node_checks += 1
    context.reauthenticate()
    result = {
        'status': 'PASS_FULL23_INITIAL_OUTPUT_SUPPORT_AND_FROZEN_MEAN_REVIEW',
        'producer_sha256': producer.sha(ROOT/'full23_initial.py'),
        'report_sha256': producer.sha(report_path), 'profiles_sha256': producer.sha(profiles_path),
        'exact_zero_channel_error_checks': exact_zero_checks,
        'frozen_mean_nodal_center_checks': mean_checks,
        'frozen_wider_cell_error_dominance_checks': dominated,
        'sampled_node_expression_slots': node_checks,
        'gap_free_closed_actual_cells': len(rows), 'gap_free_P1_cells': len(cells),
        'complete_roots_authenticated_before_and_after': True,
        'regularity': report['regularity'],
        'metadata_correction': 'The running producer had loaded the prior V>=C2 descriptive sentence; only that report sentence was corrected after the run. All numerical expressions, nodes, caps and decisions are unchanged. Current source emits the corrected C1 sentence.',
    }
    path = HERE/'OUTPUT_REVIEW.json'
    path.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
