#!/usr/bin/env python3
"""Seven targeted false-certificate controls, all in the current process.

These cases guard consequential joins and metadata conventions. They do not
establish rejection of every possible fabricated oracle row or reprove every
inherited analytic assumption.
"""
import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('independent_boundary_review', HERE / 'review_boundary_response.py')
reviewer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reviewer)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def changed(report, name):
    r = copy.deepcopy(report)
    if name == 'halve_reflected_overlap_weight':
        row = next(x for x in r['actual_bulk_source_cell_enclosures'] if x['region'] == 'Dminus')
        row['wave_kernel_weight'] = '1/2'
    elif name == 'omit_active_bulk_cell':
        r['actual_bulk_source_cell_enclosures'].pop()
    elif name == 'claim_exterior_errors_on_global_bbox':
        row = next(x for x in r['actual_bulk_source_cell_enclosures'] if x['actual_error_regime'] == 'global')
        row['actual_error_regime'] = 'exterior'; row['bbox_exterior_error_eligible'] = True
    elif name == 'reverse_retarded_boundary_lift':
        row = r['actual_boundary_source_interval_enclosures'][0]
        row['signed_negative_LG_contribution_interval'] = row['G_integral_contribution_interval'][:]
    elif name == 'omit_full_metric_lowered_QC_y':
        r['actual_boundary_source_interval_enclosures'][0]['varying_metric_lowered_QC_y_0'] = ['0', '0']
    elif name == 'double_one_copy_outward_shell_QK':
        row = r['actual_boundary_source_interval_enclosures'][0]
        for field in ('one_copy_outward_QK00_0', 'one_copy_outward_QKab_0'):
            row[field] = [str(2 * reviewer.F(value)) for value in row[field]]
    elif name == 'omit_Robin_feedback_allowance':
        row = r['full_feedback_ledger']
        row['Robin_trace_value_allowance'] = '0'
        row['total_point_mean_allowance'] = row['bulk_value_allowance']
    else:
        raise ValueError('UNKNOWN_CONTROL')
    return r


def main():
    reviewer.need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode, 'USE_PYTHON_I_S_B')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-root', type=Path, required=True)
    ap.add_argument('--previous-source-root', type=Path, required=True)
    ap.add_argument('--producer-script', type=Path, required=True)
    ap.add_argument('--report', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    reviewer.need(args.output.is_absolute() and not args.output.exists(), 'NEW_CONTROL_OUTPUT')
    raw = args.report.read_bytes(); original = reviewer.decode(raw)
    own_pin = sha(Path(__file__).read_bytes()); review_pin = sha((HERE / 'review_boundary_response.py').read_bytes())
    reviewer.review(original, args.source_root.resolve(), args.previous_source_root.resolve(), args.producer_script.resolve(), False)
    names = ('halve_reflected_overlap_weight', 'omit_active_bulk_cell',
             'claim_exterior_errors_on_global_bbox', 'reverse_retarded_boundary_lift',
             'omit_full_metric_lowered_QC_y', 'double_one_copy_outward_shell_QK',
             'omit_Robin_feedback_allowance')
    rejected = []
    for name in names:
        try:
            reviewer.review(changed(original, name), args.source_root.resolve(), args.previous_source_root.resolve(), args.producer_script.resolve(), False)
        except ValueError as error:
            rejected.append({'case': name, 'rejected': True, 'error': str(error)})
        else:
            raise ValueError('FALSE_CERTIFICATE_ACCEPTED_' + name)
    reviewer.need(len(rejected) == 7 and sha(Path(__file__).read_bytes()) == own_pin
                  and sha((HERE / 'review_boundary_response.py').read_bytes()) == review_pin
                  and args.report.read_bytes() == raw, 'CONTROL_INPUTS_OR_PROGRAM_CHANGED')
    result = {'status': 'PASS_SEVEN_TARGETED_BOUNDARY_NUMERICAL_REJECTIONS',
              'baseline_passed': True, 'count': 7, 'rejections': rejected,
              'producer_report_sha256': sha(raw), 'reviewer_source_sha256': review_pin,
              'control_source_sha256': own_pin,
              'scope': 'Seven targeted in-process numerical certificate corruptions per execution mode; no separate child per case and no exhaustive false-source proof.'}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': result['status'], 'count': result['count'], 'sha256': sha(args.output.read_bytes())}, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError) as error:
        sys.exit('BOUNDARY_NUMERICAL_CONTROLS_FAILED:' + str(error))
