"""Independent whole-cell export, adaptive-error, physical-gain and norm audit.

Uses the separately reviewed exact free-comparator module, never the parent
generator or its error evaluator. It computes adaptive errors from the pinned
independent referee proof receipt and uses exact rational physical products.
"""
import argparse
import csv
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
COMPARATOR_PIN = '27f6185218bd8d2d038530c9c2aadf5a774c38ddf96aee99ffa5f9e6a8902ed4'
ERROR_PIN = '460d2a09aab3f2450772787911667937684f68861a5509a4f441cfdd16a89a51'
GAIN_PIN = '2a2dbc398ac2799418edcddeb917f4dcecce9b9a8049ac73b54541e617b616a4'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=ROOT/'computed')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checks = []

    def check(label, ok):
        if not ok:
            raise RuntimeError(label)
        checks.append(label)

    def authenticate(path, expected):
        raw = path.read_bytes()
        check('authenticated:' + path.relative_to(ROOT).as_posix(), sha256(raw).hexdigest() == expected)
        return raw

    comp_path = ROOT/'sol/inertial_comparator.py'
    comp_raw = authenticate(comp_path, COMPARATOR_PIN)
    error_raw = authenticate(ROOT/'referee/INDEPENDENT_INERTIAL_BOUND_NORMAL.json', ERROR_PIN)
    gain_raw = authenticate(ROOT/'inputs/previous/archive/PHYSICAL_GAIN_ENCLOSURE.json', GAIN_PIN)
    spec = importlib.util.spec_from_file_location('reviewed_free_comparator', comp_path)
    iv = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(iv)
    proof = json.loads(error_raw)
    gain = json.loads(gain_raw)['correlated_retained_eta_box']
    G = tuple(map(F, gain['G']['outward_rational_at_12_decimal_places']))
    K = tuple(map(F, gain['K']['outward_rational_at_12_decimal_places']))
    A, V = F(proof['majorant']['A']), F(proof['majorant']['V'])
    R0, Ev = F(proof['residual']['R0']), F(proof['residual']['boundary_error_velocity'])
    eta = F(proof['residual']['global_H_error'])
    ell, T, M, Gamma = F(1,10**6), F(3,10**6), F(102000), F(250)
    report_path = args.data_dir/'INERTIAL_WAVEFORM_ENCLOSURE.json'
    csv_path = args.data_dir/'INERTIAL_WAVEFORM_CELLS.csv'
    report_raw, csv_raw = report_path.read_bytes(), csv_path.read_bytes()
    report = json.loads(report_raw)
    source_pins = report['source_pins']
    check('producer_provenance_matches_executed_source_pins', set(source_pins) == {'compute_inertial_waveform.py','inertial_error_bounds.py','sol/inertial_comparator.py'} and all(sha256((ROOT/name).read_bytes()).hexdigest() == digest for name,digest in source_pins.items()))
    input_records = report['mathematical_inputs']
    required_inputs = {'background':'fd6fd277176033703a4e75a2f66ab94894318aac49c83264ef979e9a27f0a02c','gain':GAIN_PIN,'shell':'8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8'}
    check('all_mathematical_input_provenance_pins_match_current_bytes', set(input_records) == set(required_inputs) and all(row['sha256'] == required_inputs[name] == sha256((ROOT/row['path']).read_bytes()).hexdigest() and row['bytes'] == (ROOT/row['path']).stat().st_size for name,row in input_records.items()))
    check('approved_state_and_full_grid', report['status'] == 'PASS' and report['state']['ell'] == str(ell) and report['cell_data']['support_cells'] == 4096 and report['cell_data']['tail_cells'] == 256 and report['cell_data']['rows'] == 4353)
    check('authenticated_complete_csv_bytes', report['cell_data']['bytes'] == len(csv_raw) and report['cell_data']['sha256'] == sha256(csv_raw).hexdigest())
    check('root_error_constants_equal_independent_derivation', all(F(report['error_certificate'][k]) == F(proof['residual'][v]) for k,v in [('R0','R0'),('history_integral','history_integral'),('boundary_error_velocity','boundary_error_velocity'),('global_error','global_H_error'),('tail_error','tail_H_error')]))
    check('physical_parameter_boxes_match_authenticated_gain', tuple(map(F,report['physical_gain_interval'])) == G and tuple(map(F,report['physical_norm_prefactor_interval'])) == K)
    mu, gamma, unused = iv.load_shell(ROOT)
    mass = iv.BumpMass(4096)
    N, S = 4096, iv.S
    check('grid_and_value_denominators_correct', report['cell_data']['time_grid_denominator'] == N and report['cell_data']['value_denominator'] == S)
    ranges = [(0,N)]+[(N+j,N+j+1) for j in range(N)]+[(2*N+16*j,2*N+16*(j+1)) for j in range(256)]
    rows = list(csv.DictReader(csv_raw.decode().splitlines()))
    check('whole_grid_exactly_contiguous_4353_cells', len(rows) == len(ranges) and all((int(r['left_grid']),int(r['right_grid'])) == pair for r,pair in zip(rows,ranges)))
    comp_sums, actual_sums = [0,0], [0,0]
    rational_norm = [F(0),F(0)]
    tail_positive, tail_negative_physical = True, True
    for index,(row,(ln,rn)) in enumerate(zip(rows,ranges)):
        left,right = F(ln,N),F(rn,N)
        comp = iv.comparator_cell(left,right,mu,gamma,mass)
        if comp != (int(row['H0_lo']),int(row['H0_hi'])):
            raise RuntimeError('Continuous comparator mismatch at cell'+str(index))
        q_upper = F(iv.q_cell(left,right)[0][1],S)
        t = ell*right
        delta = F(0) if right <= 1 else min(eta,ell*(M*t*q_upper/2+R0+M*V*t*t/4+M*M*A*t**3/6+(1+Gamma*T)*Ev))
        rounded_delta = -((-delta.numerator*S)//delta.denominator)
        if int(row['error_upper']) != rounded_delta:
            raise RuntimeError('Adaptive exact error mismatch at cell'+str(index))
        actual = int(row['H_lo']),int(row['H_hi'])
        if actual != (comp[0]-rounded_delta,comp[1]+rounded_delta):
            raise RuntimeError('Actual interval mismatch at cell'+str(index))
        physical = F(row['ell_S_lo'])/S,F(row['ell_S_hi'])/S
        products = [g*f*h for g in G for f in (1-t*t,F(1)) for h in (F(actual[0],S),F(actual[1],S))]
        if not physical[0] <= min(products) <= max(products) <= physical[1]:
            raise RuntimeError('Physical signed gain/clock enclosure fails at cell'+str(index))
        if left >= 2:
            tail_positive &= actual[0] > 0
            tail_negative_physical &= physical[1] < 0
        for values,total in [(comp,comp_sums),(actual,actual_sums)]:
            lo2 = 0 if values[0] <= 0 <= values[1] else min(v*v for v in values)
            hi2 = max(v*v for v in values)
            total[0] += (rn-ln)*(lo2//S)
            total[1] += (rn-ln)*(-((-hi2)//S))
        low_square = F(0) if actual[0] <= 0 <= actual[1] else min(F(v,S)**2 for v in actual)
        high_square = max(F(v,S)**2 for v in actual)
        rational_norm[0] += F(rn-ln,N)*low_square
        rational_norm[1] += F(rn-ln,N)*high_square
    check('all_4353_comparison_error_actual_and_physical_cells_reconstructed', True)
    check('all_256_tail_cells_positive_H_negative_physical_coefficient', tail_positive and tail_negative_physical)
    check('prefront_export_exactly_zero', all(int(rows[0][k]) == 0 for k in ['H0_lo','H0_hi','error_upper','H_lo','H_hi','ell_S_lo','ell_S_hi']))
    comp_int = [F(v,N*S) for v in comp_sums]
    actual_int = [F(v,N*S) for v in actual_sums]
    check('direct_whole_cell_norm_sums_match_exact_report', comp_int == list(map(F,report['quadrature']['comparison_H0_integral_interval'])) and actual_int == list(map(F,report['quadrature']['actual_H_integral_interval'])))
    check('independent_unrounded_rational_norm_contained', actual_int[0] <= rational_norm[0] <= rational_norm[1] <= actual_int[1])
    C = K[0]*actual_int[0]/ell,K[1]*actual_int[1]/ell
    check('new_wider_C_recomputed_with_correct_inverse_width_jacobian', C == tuple(map(F,report['quadrature']['C_interval'])))
    check('source_and_data_unchanged_throughout_audit', comp_path.read_bytes() == comp_raw and report_path.read_bytes() == report_raw and csv_path.read_bytes() == csv_raw and all(sha256((ROOT/name).read_bytes()).hexdigest() == digest for name,digest in source_pins.items()))
    out = {'status':'PASS_INDEPENDENT_CONTINUOUS_INERTIAL_CELLS_AND_NORM',
        'python_optimization':sys.flags.optimize,'check_count':len(checks),'checks':checks,
        'rows_reconstructed':len(rows),'actual_integral_interval':[str(v) for v in actual_int],
        'C_interval':[str(v) for v in C],
        'display_only_C':[float(v) for v in C],
        'input_pins':{'comparator':COMPARATOR_PIN,'independent_error':ERROR_PIN,'gain':GAIN_PIN,
            'waveform_report':sha256(report_raw).hexdigest(),'waveform_csv':sha256(csv_raw).hexdigest()},
        'own_source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Independent reconstruction of every exported continuous comparison/error/physical cell and the fresh wider-state norm coefficient. The separately reviewed comparator is used; parent generator and error evaluator are not imported. No nonlinear remainder or physical-unit calibration.'}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'checks':len(checks),'rows':len(rows),'C_display_only':out['display_only_C']}))


if __name__ == '__main__':
    main()
