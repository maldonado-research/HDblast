"""Independent exact-rational validation of the explicit material tension enclosure."""
import argparse
import csv
from fractions import Fraction as R
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = 1 << 112
N = 4096
L = R(1, 1000000)
PINS = {'INERTIAL_WAVEFORM_CELLS.csv': '2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7',
        'INERTIAL_WAVEFORM_ENCLOSURE.json': 'd2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
        'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json': '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
        'PRIOR_PHYSICAL_GAIN_ENCLOSURE.json': '2a2dbc398ac2799418edcddeb917f4dcecce9b9a8049ac73b54541e617b616a4'}


def check(data_dir, input_root):
    passed = []

    def gate(ok, name):
        if not ok:
            raise ValueError(name)
        passed.append(name)

    def pair(values):
        return tuple(map(R, values))

    def pow_interval(a, power):
        candidates = [a[0]**power, a[1]**power]
        return (R(0) if a[0] <= 0 <= a[1] else min(candidates), max(candidates))

    def rounded(a):
        return a[0].numerator*D//a[0].denominator, -((-a[1].numerator*D)//a[1].denominator)

    prior_root = input_root/'prior'
    gate(all(sha256((prior_root/name).read_bytes()).hexdigest() == digest for name, digest in PINS.items()),
         'authenticated complete prior scientific inputs')
    old = json.loads((prior_root/'INERTIAL_WAVEFORM_ENCLOSURE.json').read_bytes())
    shell = json.loads((prior_root/'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_bytes())
    report = json.loads((data_dir/'QUADRATIC_FORCING_ENCLOSURE.json').read_bytes())
    gate(report['schema'] == 'hdblast-quadratic-material-tension-forcing-v1' and report['status'] == 'PASS',
         'quadratic material tension report identity')
    state = report['state']
    gate(state['ell'] == str(L) and state['delta'] == '1/1000' and state['domain_x'] == [0, 3],
         'same registered state and full proper shell window')
    gate(state['derivative_convention'] == 'Ueps=U0+eps U1+eps^2 U2/2',
         'second derivative factorial convention')
    gate(state['normal_convention'] == 'outward normal SUM over both copies'
         and report['definitions']['scalar_material_source'] == 'Qphi=-4 psi1^2=-4 S(t)^2 Y^2'
         and report['definitions']['per_copy_scalar_material_source'] == '-2 psi1^2',
         'correct scalar SUM and per-copy tension sources')
    eta = pair(shell['boxes']['eta_b']['exact'])
    sig = tuple(4+4*v for v in eta)
    pars = report['parameters']
    gate(pair(pars['sigma_second_interval']) == sig and pars['sigma_third'] == '4' and sig[0] > 0,
         'actual registered positive tension Hessian and cubic derivative')
    g, k = pair(old['physical_gain_interval']), pair(old['physical_norm_prefactor_interval'])
    gate(pair(pars['physical_gain_interval']) == g and pair(pars['physical_norm_prefactor_interval']) == k,
         'physical gain and proper measure normalization')
    a = report['angular']
    gate(a['sphere_volume_over_pi_squared'] == '2' and a['mean_Y2'] == '1/6' and a['mean_Y4'] == '3/40'
         and a['mean_Z2_squared'] == '1/12' and a['mean_H4_squared'] == '19/720',
         'exact unit three-sphere moments')
    shares = tuple(map(R, a['orthogonal_squared_power_shares_l0_l2_l4']))
    gate(shares == (R(10, 27), R(5, 18), R(19, 54)) and sum(shares) == 1,
         'orthogonal squared angular power shares')
    defs = report['definitions']
    gate(defs['Y_squared_decomposition'] == ['1/6', '1/2', '1']
         and defs['scalar_coefficient_multipliers_of_S_squared'] == ['-2/3', '-2', '-4']
         and defs['tension_coefficient_multipliers_of_sigma_second_S_squared'] == ['1/6', '1/2', '1'],
         'exact projected scalar and constitutive coefficients')
    gate(report['arithmetic']['value_denominator'] == D and report['arithmetic']['time_grid_denominator'] == N
         and report['arithmetic']['fourth_moment_clock_bound'] == '[(1-(ell*right_x)^2)^2,1]',
         'outward arithmetic and retained positive clock inequality')
    cells_path = data_dir/'QUADRATIC_FORCING_CELLS.csv'
    raw = cells_path.read_bytes()
    gate(len(raw) == report['cell_data']['bytes'] and sha256(raw).hexdigest() == report['cell_data']['sha256'],
         'complete quadratic CSV authenticated')
    with (prior_root/'INERTIAL_WAVEFORM_CELLS.csv').open(newline='') as f:
        old_rows = list(csv.DictReader(f))
    with cells_path.open(newline='') as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [[int(v) for v in row] for row in reader]
    expected_header = ['left_grid', 'right_grid']+[f'ell2_{kind}{degree}_{side}'
        for kind in ['Q', 'T'] for degree in [0, 2, 4] for side in ['lo', 'hi']]
    gate(header == expected_header == report['cell_data']['header']
         and len(rows) == len(old_rows) == report['cell_data']['rows'] == 4353,
         'all continuous coefficient cells and semantic columns present')
    j2 = [0, 0]
    j4 = [R(0), R(0)]
    last = 0
    values_ok = True
    for source, row in zip(old_rows, rows):
        left, right = int(source['left_grid']), int(source['right_grid'])
        values_ok = values_ok and row[:2] == [left, right] and left == last and left < right
        last = right
        h = (R(int(source['H_lo']), D), R(int(source['H_hi']), D))
        e = (R(int(source['ell_S_lo']), D), R(int(source['ell_S_hi']), D))
        hsq = rounded(pow_interval(h, 2))
        for side in range(2):
            j2[side] += (right-left)*hsq[side]
        fourth = pow_interval(h, 4)
        j4[0] += R(right-left, N)*(1-(L*R(right, N))**2)**2*fourth[0]
        j4[1] += R(right-left, N)*fourth[1]
        esq = pow_interval(e, 2)
        exact = []
        for factor in [R(-2, 3), R(-2), R(-4)]:
            exact += list(rounded((factor*esq[1], factor*esq[0])))
        for factor in [R(1, 6), R(1, 2), R(1)]:
            exact += list(rounded((factor*sig[0]*esq[0], factor*sig[1]*esq[1])))
        values_ok = values_ok and row[2:] == exact
    gate(values_ok and last == 3*N, 'every quadratic continuous cell independently reconstructed')
    quad = report['quadrature']
    j2 = tuple(R(v, N*D) for v in j2)
    gate(j2 == pair(quad['H_squared_integral_interval']) == pair(old['quadrature']['actual_H_integral_interval']),
         'entire inherited second moment reconstructed')
    gate(tuple(j4) == pair(quad['cosh_minus3_H_fourth_integral_interval']) and j4[0] > 0,
         'entire weighted fourth moment independently reconstructed')
    c = (k[0]*j2[0]/L, k[1]*j2[1]/L)
    gate(c == pair(quad['inherited_C_interval']) == pair(old['quadrature']['C_interval']),
         'first-order proper-time norm identity')
    gate(pair(quad['scalar_signed_integral_interval']) == (-4*c[1], -4*c[0])
         and pair(quad['scalar_L1_interval']) == (4*c[0], 4*c[1]),
         'signed forcing integral and L1 identity')
    g2 = pow_interval(g, 2)
    qnorm = tuple(R(36, 5)*k[i]*g2[i]*j4[i]/L**3 for i in range(2))
    gate(pair(quad['scalar_L2_squared_interval']) == qnorm,
         'fourth-moment proper-time L2 Jacobian and normalization')
    gate(pair(quad['constitutive_tension_integral_interval']) == tuple(sig[i]*c[i] for i in range(2))
         and pair(quad['constitutive_tension_L2_squared_interval']) == tuple(
             R(9, 20)*sig[i]**2*k[i]*g2[i]*j4[i]/L**3 for i in range(2)),
         'constitutive tension integrals retain actual Hessian')
    ht = pair(old['tail_theorem']['actual_H_interval'])
    ep = (g[0]*ht[1], g[1]*(1-(3*L)**2)*ht[0])
    es = pow_interval(ep, 2)
    tail = (-R(2, 3)*es[1], -R(2, 3)*es[0])
    gate(report['tail_theorem']['domain_x'] == [2, 3]
         and pair(report['tail_theorem']['ell2_scalar_monopole_interval']) == tail
         and R('-2.3e-19') < tail[0] <= tail[1] < R('-6.7e-20'),
         'strict monopole bound on full continuous tail')
    gate(defs['shell_measure'] == 'dmu0=rho_b^4 cosh(t)^3 dt dOmega',
         'background proper shell measure explicit')
    gate(report['scope'] == {'explicit_material_tension_contribution': True,
        'complete_bulk_coordinate_quadratic_source': False, 'full_U2_response_computed': False,
        'total_tension_or_stress_second_variation': False, 'nonlinear_remainder_bound': False,
        'novelty_claimed': False, 'cosmological_evidence': False},
        'explicit contribution scope and unresolved second response retained')
    return {'status': 'PASS_QUADRATIC_MATERIAL_TENSION_ENCLOSURE', 'checks': passed,
            'check_count': len(passed), 'csv_rows': len(rows),
            'report_sha256': sha256((data_dir/'QUADRATIC_FORCING_ENCLOSURE.json').read_bytes()).hexdigest(),
            'csv_sha256': sha256(raw).hexdigest(),
            'scope': 'Exact arithmetic validation of the explicit material tension contribution; full U2 and coordinate source remain separate.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=HERE/'computed')
    parser.add_argument('--input-root', type=Path, default=HERE/'inputs')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = check(args.data_dir, args.input_root)
        args.output.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({'status': result['status'], 'check_count': result['check_count']}))
    except (ValueError, KeyError, OSError) as exc:
        print('FAIL: '+str(exc))
        raise SystemExit(1)
