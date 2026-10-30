#!/usr/bin/env python3
"""Reject deliberately corrupted refinement certificates with the pinned auditor.

Only the named coverage/aggregation audit functions are loaded from that file.
No project scientific code or interval-leaf producer is imported here.
"""
import argparse
import ast
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def read(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return raw


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', required=True, type=Path)
ap.add_argument('--audit-checker', type=Path)
ap.add_argument('--refined-report', type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
auditor = args.audit_checker or Path(__file__).with_name('check_refined_response.py')
report = args.refined_report or Path(__file__).resolve().parents[1]/'reference/REFINED.json'
source = read(auditor,
    '0e182ce5421b2cb42726ea452ef7bd9b26e643dba445916dbc2aa118567e3240')
refined = json.loads(read(report,
    'd97f7393cb9467e78387368754981e81cf763da043cab0a4df65978f52b72122'))
initial = json.loads(read(root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results/INITIAL_REPORT.json',
    'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'))
names = {'need', 'interval', 'caps_array', 'aggregate', 'verify_refinement'}
parsed = ast.parse(source, filename=str(auditor))
nodes = [node for node in parsed.body
         if isinstance(node, ast.FunctionDef) and node.name in names]
need({node.name for node in nodes} == names and len(nodes) == len(names),
     'exact named coverage-auditor function set')
namespace = {'Q': Q}
exec(compile(ast.Module(body=nodes, type_ignores=[]), str(auditor), 'exec'), namespace)
audit = namespace['verify_refinement']
labels = initial['labels']
weights = list(map(Q, initial['gram']))
totals, stats = audit(refined, initial, labels, weights)
need(stats == {'inherited_cells': 192, 'refined_parent_cells': 1,
               'verified_refined_leaves': 10, 'component_caps_aggregated_per_leaf': 207},
     'successful certificate has expected independently counted cover')
need(totals == list(map(Q, refined['total_H_caps'])), 'baseline totals')


def parent(d):
    return next(row for row in d['coverage_cells'] if row['refined'])


def leaf(d):
    return parent(d)['leaves'][0]


def missing_parent(d):
    d['coverage_cells'].pop(0)


def wrong_parent_identity(d):
    d['coverage_cells'][0]['old_row_index'] = 1


def missing_leaf(d):
    parent(d)['leaves'].pop(1)


def duplicated_leaf(d):
    parent(d)['leaves'].append(copy.deepcopy(leaf(d)))


def invalid_binary_path(d):
    leaf(d)['path'] += 'x'


def wrong_leaf_endpoint(d):
    row = leaf(d)
    row['r'][1] = str((Q(row['r'][0])+Q(row['r'][1]))/2)


def negative_raw_cap(d):
    leaf(d)['raw_channel_component_error_caps'][1][1][12] = '-1'


def missing_channel(d):
    leaf(d)['raw_channel_component_error_caps'].pop()


def incorrect_cap_intersection(d):
    leaf(d)['channel_component_error_caps'][1][1][12] = '0'


def incorrect_global_component_maximum(d):
    d['channel_component_error_caps'][1][1][12] = '0'


def understated_Gram_norm(d):
    d['channel_H_caps'][1][1] = '0'


def understated_channel_sum(d):
    d['total_H_caps'][1] = '0'


def wrong_Gram_weight(d):
    d['gram'][12] = '2'


def changed_candidate(d):
    d['same_candidate'] = False


def wrong_derivative_convention(d):
    d['derivative_convention']['du'] = 'partial_x U'


def changed_producer_pin(d):
    d['input_sha256']['initial_producer'] = '0'*64


controls = [missing_parent, wrong_parent_identity, missing_leaf, duplicated_leaf,
            invalid_binary_path, wrong_leaf_endpoint, negative_raw_cap,
            missing_channel, incorrect_cap_intersection,
            incorrect_global_component_maximum, understated_Gram_norm,
            understated_channel_sum, wrong_Gram_weight, changed_candidate,
            wrong_derivative_convention, changed_producer_pin]
rejected = []
for mutation in controls:
    damaged = copy.deepcopy(refined)
    mutation(damaged)
    try:
        audit(damaged, initial, labels, weights)
    except RuntimeError:
        rejected.append(mutation.__name__)
    else:
        raise RuntimeError('CORRUPTED_CERTIFICATE_ACCEPTED:' + mutation.__name__)
need(len(rejected) == len(controls), 'all deliberate corruptions rejected')
print(json.dumps({'status': 'PASS', 'baseline_coverage_audit': stats,
    'deliberately_corrupted_certificates_rejected': len(rejected),
    'controls': rejected,
    'scope': 'finite guard tests of exact coverage, positivity, max/Gram/triangle aggregation and same-candidate metadata; does not regenerate or independently prove the scientific leaf intervals'}, sort_keys=True))
