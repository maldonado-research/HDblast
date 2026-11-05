"""Exact postprocessing of retained continuous waveform cells.

This does not execute a legacy producer, solve a PDE or reprove retained
continuous-cell/reflection/same-state analytic premises.
"""
import argparse
import csv
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PINS = {
    'OLD_FUTURE_PROOF.md': '7897073c52ee87514b26d9377ac9b3cd88a4b082c47f35bf36164e6f97ad2458',
    'OLD_FUTURE_REPORT.json': 'd2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
    'OLD_FUTURE_CELLS.csv': '2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7',
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checks = []

    def need(ok, name):
        if not ok:
            raise RuntimeError(name)
        checks.append(name)

    output = args.output
    need(output.is_absolute() and not output.exists() and
         all(not p.is_symlink() for p in [output, *output.parents]) and
         not output.resolve().is_relative_to(HERE), 'fresh absolute external output')
    payload = {}
    for name, pin in PINS.items():
        path = HERE / 'inputs' / name
        need(path.is_file() and not path.is_symlink(), 'regular selected source ' + name)
        raw = path.read_bytes()
        need(hashlib.sha256(raw).hexdigest() == pin, 'source SHA256 ' + name)
        payload[name] = raw
    report = json.loads(payload['OLD_FUTURE_REPORT.json'])
    ell = F(report['state']['ell'])
    B = F(264000000001, 13738809600)
    need(report['schema'] == 'hdblast-inertial-shell-waveform-v1' and
         report['status'] == 'PASS', 'retained actual waveform report schema')
    need(ell == F(1, 10**6) and F(report['state']['T']) == 3 * ell and
         report['state']['time_coordinate'] == 'x=t/ell=tau/(rho_b ell)' and
         report['state']['initial_velocity'] == 'exp(-z) Q(z/ell)',
         'registered ell clock and velocity metadata')
    need(B < 20 and B > 0, 'positive exact companion reflection bound below twenty')
    cell = report['cell_data']
    raw = payload['OLD_FUTURE_CELLS.csv']
    need(cell['bytes'] == len(raw) and cell['sha256'] == PINS['OLD_FUTURE_CELLS.csv'],
         'CSV bound to retained actual report')
    dt, dv = cell['time_grid_denominator'], cell['value_denominator']
    need(dt == 4096 and dv == 2**112 == report['arithmetic']['value_denominator'],
         'time and value denominators')
    reader = csv.DictReader(raw.decode('utf-8').splitlines())
    need(reader.fieldnames == ['left_grid', 'right_grid', 'H0_lo', 'H0_hi',
                              'error_upper', 'H_lo', 'H_hi', 'ell_S_lo', 'ell_S_hi'],
         'exact retained cell columns')
    integer_rows = [{key: int(value) for key, value in row.items()} for row in reader]
    need(len(integer_rows) == cell['rows'] == 4353 and
         cell['support_cells'] == 4096 and cell['tail_cells'] == 256,
         'complete continuous cell roster')
    need(integer_rows[0]['left_grid'] == 0 and
         integer_rows[-1]['right_grid'] == 3 * dt and
         all(row['left_grid'] < row['right_grid'] for row in integer_rows) and
         all(a['right_grid'] == b['left_grid'] for a, b in zip(integer_rows, integer_rows[1:])),
         'closed adjacent cell coverage of zero to three')
    need(all(row['H0_lo'] <= row['H0_hi'] and row['H_lo'] <= row['H_hi'] and
             row['ell_S_lo'] <= row['ell_S_hi'] and row['error_upper'] >= 0
             for row in integer_rows), 'all retained interval rows ordered')

    cells = []
    for index, row in enumerate(integer_rows):
        future = (F(row['ell_S_lo'], dv) / ell, F(row['ell_S_hi'], dv) / ell)
        cells.append({'csv_index_zero_based': index, 'raw': row,
                      'left': F(row['left_grid'], dt), 'right': F(row['right_grid'], dt),
                      'future': future, 'past': (future[0] - B, future[1] + B)})
    negatives = [row for row in cells if row['past'][1] < 0 and row['right'] < F(3, 2)]
    positives = [row for row in cells if row['past'][0] > 0 and row['left'] > F(3, 2)]
    need(bool(negatives) and bool(positives), 'strict transported signs on both central sides')
    negative = max(negatives, key=lambda row: row['right'])
    positive = min(positives, key=lambda row: row['left'])
    need(negative['csv_index_zero_based'] == 2047 and
         (negative['left'], negative['right']) == (F(3071, 2048), F(6143, 4096)),
         'latest strict negative whole cell before center')
    need(positive['csv_index_zero_based'] == 2050 and
         (positive['left'], positive['right']) == (F(6145, 4096), F(3073, 2048)),
         'earliest strict positive whole cell after center')
    need(negative['past'][1] < -53 and positive['past'][0] > 63,
         'strict endpoint scalar margins in model units')
    left, right = -positive['left'], -negative['right']
    need((left, right) == (-F(6145, 4096), -F(6143, 4096)) and
         left < right and right - left == F(1, 2048),
         'open past zero bracket and width')
    need(F(1, 2) / (right - left) == 1024,
         'exact improvement relative to prior coarse half-width bracket')
    # The two inner center cells still contain zero after the same transport.
    need(all(row['past'][0] <= 0 <= row['past'][1] for row in cells[2048:2050]),
         'inner cells do not supply an even narrower sign witness')

    def exported(row):
        return {'csv_index_zero_based': row['csv_index_zero_based'], 'raw_integer_row': row['raw'],
                'future_scaled_time_interval': [str(row['left']), str(row['right'])],
                'actual_future_S_interval': [str(x) for x in row['future']],
                'transported_past_S_interval': [str(x) for x in row['past']]}

    result = {
        'schema': 'HDBLAST_PAST_ZERO_CONTINUOUS_CELL_REFINEMENT_V1',
        'status': 'PASS_CONDITIONAL_PAST_ZERO_CELL_REFINEMENT',
        'finite_guards': len(checks), 'check_names': checks,
        'selected_input_sha256': PINS,
        'postprocessor_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'reflection_premise': 'The exact uniform B bound is proved separately in the companion primary proof; this postprocessor does not reprove it.',
        'ell': str(ell), 'uniform_difference_upper': str(B),
        'csv_rows': len(cells), 'negative_transported_cells': sum(row['past'][1] < 0 for row in cells),
        'positive_transported_cells': sum(row['past'][0] > 0 for row in cells),
        'negative_witness_cell': exported(negative), 'positive_witness_cell': exported(positive),
        'open_past_zero_bracket_in_units_of_ell': [str(left), str(right)],
        'bracket_width_in_units_of_ell': str(right-left),
        'width_reduction_factor_from_half_ell_bracket': 1024,
        'strict_endpoint_claims': ['S(-6145ell/4096)>63', 'S(-6143ell/4096)<-53'],
        'conditional_continuity_conclusion': 'At least one zero of the actual first physical scalar amplitude lies in the open stated bracket.',
        'scope': 'Exact postprocessing. Retained continuous-cell, companion reflection, same-state and regularity premises are not rerun or proved by this arithmetic. No uniqueness, full zero count, finer historical future brackets, tail sign, relative error, physical energy or nonlinear/cosmological conclusion.',
        'new_activity': {'legacy_producer_or_checker_executions': 0, 'numerical_PDE_solves': 0,
                         'background_field_source_bump_oracles': 0},
    }
    with output.open('x') as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write('\n')
    print(json.dumps({'status': result['status'], 'finite_guards': len(checks),
                      'bracket_in_units_of_ell': result['open_past_zero_bracket_in_units_of_ell']}, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, ValueError, OSError, KeyError, TypeError) as error:
        print(json.dumps({'status': 'FAIL', 'witness': str(error)}, sort_keys=True))
        sys.exit(2)
