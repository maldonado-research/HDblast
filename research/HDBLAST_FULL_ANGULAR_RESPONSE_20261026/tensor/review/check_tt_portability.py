#!/usr/bin/env python3
"""Check an I/O-only TT portability change against retained execution records.

Only Python standard-library modules are used.  No source oracle is executed.
The reference arguments select preserved records outside the publication tree.
"""
import argparse
import ast
import copy
import hashlib
import json
from pathlib import Path

CAPSULES = (
    'HDBLAST_COUPLED_MEAN_PILOT_20261025',
    'HDBLAST_INPUT_PRECISION_20261023',
    'HDBLAST_COMMON_GAUGE_FORCING_20261016',
    'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016',
    'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
    'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023',
    'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023',
    'HDBLAST_SHELL_FEEDBACK_PRECISION_20261025',
)
CORE_FUNCTIONS = ('forcing', 'operator', 'finish', 'probes', 'partition',
                  'authenticate', 'encode')


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def function_ast(tree, name):
    if name in ('forcing', 'operator', 'finish'):
        tree = next(node for node in tree.body
                    if isinstance(node, ast.ClassDef) and node.name == 'Context')
    node = next(node for node in tree.body
                if isinstance(node, ast.FunctionDef) and node.name == name)
    return ast.dump(node, include_attributes=False)


def payload(report):
    value = copy.deepcopy(report)
    value.pop('wall_seconds', None)
    value.pop('evaluator_sha256')
    value['source'].pop('path')
    value['latest_authentication_helper'].pop('path')
    return value


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for flag in ('repo-root', 'reference-report', 'reference-evaluator',
                 'report', 'optimized-report', 'evaluator', 'output'):
        p.add_argument('--' + flag, type=Path, required=True)
    args = p.parse_args()
    target = args.output.resolve()
    protected = [(args.repo_root.resolve()/'research'/name).resolve() for name in CAPSULES]
    need(not args.output.exists() and not args.output.is_symlink(), 'FRESH_OUTPUT_REQUIRED')
    need(all(target != root and root not in target.parents for root in protected),
         'OUTPUT_INSIDE_FROZEN_CAPSULE')
    source_before = ast.parse(args.reference_evaluator.read_text())
    source_after = ast.parse(args.evaluator.read_text())
    function_pins = {}
    for name in CORE_FUNCTIONS:
        before = function_ast(source_before, name)
        after = function_ast(source_after, name)
        need(before == after, 'CORE_FUNCTION_AST_CHANGED:' + name)
        function_pins[name] = hashlib.sha256(after.encode()).hexdigest()
    old = json.loads(args.reference_report.read_text())
    new = json.loads(args.report.read_text())
    optimized = json.loads(args.optimized_report.read_text())
    need(old['evaluator_sha256'] == sha(args.reference_evaluator), 'REFERENCE_SOURCE_BINDING')
    need(new['evaluator_sha256'] == sha(args.evaluator), 'FINAL_SOURCE_BINDING')
    need('wall_seconds' in old and 'wall_seconds' not in new, 'TIMING_METADATA_CHANGE')
    need(args.report.read_bytes() == args.optimized_report.read_bytes(), 'NORMAL_OPTIMIZED_BYTES')
    for section in ('source', 'latest_authentication_helper'):
        path = new[section]['path']
        need(not Path(path).is_absolute() and path.startswith('research/')
             and '..' not in Path(path).parts, 'LOGICAL_REPO_PATH:' + section)
        need(old[section]['path'].endswith('/' + path), 'SOURCE_PATH_IDENTITY:' + section)
    before, after = payload(old), payload(new)
    need(before == after == payload(optimized), 'EXACT_NUMERICAL_PAYLOAD_CHANGED')
    interval = new['partition_runs'][-1]['total_Tplus_response_interval']
    need(interval == old['partition_runs'][-1]['total_Tplus_response_interval'],
         'EXACT_INTERVAL_CHANGED')
    out = {
        'schema': 'HDBLAST_TT_IO_ONLY_PORTABILITY_REVIEW_V1',
        'status': 'PASS_EXACT_NUMERICAL_PAYLOAD_UNCHANGED',
        'all_numerical_payload_fields_exactly_equal': True,
        'exact_total_interval_unchanged': interval,
        'normal_optimized_entire_bytes_identical': True,
        'permitted_metadata_changes': ['evaluator_sha256', 'remove wall_seconds',
                                       'source/path and latest_authentication_helper/path become repo-relative'],
        'normalized_before_after_payload_sha256': canonical_hash(after),
        'core_function_ast_sha256_identical_before_after': function_pins,
        'files_sha256': {
            'reference/tensor/tt_response.py': sha(args.reference_evaluator),
            'reference/tensor/TT_TOTAL_RESPONSE_NORMAL.json': sha(args.reference_report),
            'tensor/tt_response.py': sha(args.evaluator),
            'tensor/TT_TOTAL_RESPONSE_NORMAL.json': sha(args.report),
            'tensor/TT_TOTAL_RESPONSE_OPTIMIZED.json': sha(args.optimized_report),
            'tensor/review/check_tt_portability.py': sha(Path(__file__)),
        },
        'limits': ['Reference artifacts remain unchanged outside the publication tree.',
                   'This review executes no source oracle and grants no new mathematical premise.'],
    }
    with args.output.open('x') as handle:
        json.dump(out, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({'status': out['status'],
                      'payload_sha256': out['normalized_before_after_payload_sha256']}))


if __name__ == '__main__':
    main()
