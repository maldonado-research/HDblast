"""Exact outward enclosures of explicit quadratic material tension contributions."""
import argparse
import csv
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEN = 2**112
GRID = 4096
ELL = F(1, 10**6)
PINS = {
    'prior/INERTIAL_WAVEFORM_CELLS.csv': '2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7',
    'prior/INERTIAL_WAVEFORM_ENCLOSURE.json': 'd2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
    'prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json': '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
    'prior/PRIOR_PHYSICAL_GAIN_ENCLOSURE.json': '2a2dbc398ac2799418edcddeb917f4dcecce9b9a8049ac73b54541e617b616a4',
}
HEADER = ['left_grid', 'right_grid', 'ell2_Q0_lo', 'ell2_Q0_hi',
          'ell2_Q2_lo', 'ell2_Q2_hi', 'ell2_Q4_lo', 'ell2_Q4_hi',
          'ell2_T0_lo', 'ell2_T0_hi', 'ell2_T2_lo', 'ell2_T2_hi',
          'ell2_T4_lo', 'ell2_T4_hi']


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def encode(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    return x


def square(a):
    lo, hi = a
    require(lo <= hi, 'ordered interval required')
    return (F(0) if lo <= 0 <= hi else min(lo*lo, hi*hi), max(lo*lo, hi*hi))


def scale(a, c):
    return tuple(sorted((a[0]*c, a[1]*c)))


def positive_product(a, b):
    require(0 <= a[0] <= a[1] and 0 <= b[0] <= b[1], 'positive interval product required')
    return a[0]*b[0], a[1]*b[1]


def floor(x):
    return x.numerator//x.denominator


def ceil(x):
    return -((-x.numerator)//x.denominator)


def fixed(a):
    return floor(a[0]*DEN), ceil(a[1]*DEN)


def calculate(input_root, output):
    input_root, output = Path(input_root), Path(output)
    require(not output.exists(), 'fresh output directory required')
    sources = {}
    for relative, expected in PINS.items():
        raw = (input_root/relative).read_bytes()
        require(sha256(raw).hexdigest() == expected, 'authenticated prior input differs: '+relative)
        sources[relative] = {'bytes': len(raw), 'sha256': expected}
    prior = json.loads((input_root/'prior/INERTIAL_WAVEFORM_ENCLOSURE.json').read_bytes())
    shell = json.loads((input_root/'prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_bytes())
    require(prior['state']['ell'] == str(ELL) and prior['state']['domain'] == [0, 3]
            and prior['arithmetic']['value_denominator'] == DEN, 'same registered wider state required')
    eta = tuple(map(F, shell['boxes']['eta_b']['exact']))
    sigma2 = tuple(4+4*x for x in eta)
    require(sigma2 == tuple(map(F, shell['boxes']['sigma_second_b']['exact']))
            and sigma2[0] > 0, 'actual registered positive tension Hessian required')
    gain = tuple(map(F, prior['physical_gain_interval']))
    k = tuple(map(F, prior['physical_norm_prefactor_interval']))
    gain2 = square(gain)
    j2 = [0, 0]
    j4 = [F(0), F(0)]
    rows = []
    with (input_root/'prior/INERTIAL_WAVEFORM_CELLS.csv').open(newline='') as stream:
        reader = csv.DictReader(stream)
        old_rows = list(reader)
    require(len(old_rows) == 4353, 'complete inherited continuous grid required')
    last = 0
    tail_means = []
    for row in old_rows:
        ln, rn = int(row['left_grid']), int(row['right_grid'])
        require(ln == last and ln < rn <= 3*GRID, 'contiguous whole-cell partition required')
        last = rn
        h = (F(int(row['H_lo']), DEN), F(int(row['H_hi']), DEN))
        e = (F(int(row['ell_S_lo']), DEN), F(int(row['ell_S_hi']), DEN))
        h2, e2 = square(h), square(e)
        h4 = square(h2)
        # This follows by squaring the retained positive bound
        # cosh(t)^(-3/2) >= 1-t^2, valid here for 0 <= t <= 3e-6.
        clock_lo = (1-(ELL*F(rn, GRID))**2)**2
        width = F(rn-ln, GRID)
        for side, value in enumerate(fixed(h2)):
            j2[side] += (rn-ln)*value
        j4[0] += width*clock_lo*h4[0]
        j4[1] += width*h4[1]
        scalar = [scale(e2, F(-2, 3)), scale(e2, -2), scale(e2, -4)]
        tension_base = positive_product(sigma2, e2)
        tension = [scale(tension_base, F(1, 6)), scale(tension_base, F(1, 2)), tension_base]
        rows.append([ln, rn, *(v for iv in scalar+tension for v in fixed(iv))])
        if ln >= 2*GRID:
            tail_means.append(scalar[0])
    require(last == 3*GRID, 'full contact window required')
    j2 = tuple(F(x, GRID*DEN) for x in j2)
    require(j2 == tuple(map(F, prior['quadrature']['actual_H_integral_interval'])),
            'reconstructed second moment equals inherited certificate')
    c = (k[0]*j2[0]/ELL, k[1]*j2[1]/ELL)
    require(c == tuple(map(F, prior['quadrature']['C_interval'])), 'physical proper-time norm identity required')
    scalar_l2sq = scale(positive_product(positive_product(k, gain2), j4), F(36, 5)/ELL**3)
    tension_integral = positive_product(sigma2, c)
    tension_l2sq = scale(positive_product(positive_product(positive_product(square(sigma2), k), gain2), j4),
                         F(9, 20)/ELL**3)
    # A domain-wide tail enclosure uses the retained continuous tail theorem,
    # rather than just endpoint samples or the coarser per-cell union.
    ht = tuple(map(F, prior['tail_theorem']['actual_H_interval']))
    et = (gain[0]*ht[1], gain[1]*(1-(3*ELL)**2)*ht[0])
    qt = scale(square(et), F(-2, 3))
    require(F('-2.3e-19') < qt[0] <= qt[1] < F('-6.7e-20'), 'strict full-tail monopole bounds')
    output.mkdir(parents=True)
    csv_path = output/'QUADRATIC_FORCING_CELLS.csv'
    with csv_path.open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(HEADER)
        writer.writerows(rows)
    report = {
        'schema': 'hdblast-quadratic-material-tension-forcing-v1', 'status': 'PASS',
        'state': {'delta': '1/1000', 'ell': ELL, 'domain_x': [0, 3],
                  'time_coordinate': 'x=t/ell=tau/(rho_b ell)', 'angular_state': 'Y=n1^2-n2^2 on unit S3',
                  'derivative_convention': 'Ueps=U0+eps U1+eps^2 U2/2',
                  'normal_convention': 'outward normal SUM over both copies'},
        'definitions': {
            'scalar_material_source': 'Qphi=-4 psi1^2=-4 S(t)^2 Y^2',
            'per_copy_scalar_material_source': '-2 psi1^2',
            'constitutive_tension_term': 'Tquad=sigma_second_b psi1^2',
            'shell_measure': 'dmu0=rho_b^4 cosh(t)^3 dt dOmega',
            'coefficient_basis': ['1', 'Z2=n1^2+n2^2-1/2', 'H4=Y^2-(n1^2+n2^2)/2+1/12'],
            'Y_squared_decomposition': ['1/6', '1/2', '1'],
            'scalar_coefficient_multipliers_of_S_squared': ['-2/3', '-2', '-4'],
            'tension_coefficient_multipliers_of_sigma_second_S_squared': ['1/6', '1/2', '1']},
        'angular': {'sphere_volume_over_pi_squared': '2', 'mean_Y2': '1/6', 'mean_Y4': '3/40',
                    'mean_Z2_squared': '1/12', 'mean_H4_squared': '19/720',
                    'orthogonal_squared_power_shares_l0_l2_l4': ['10/27', '5/18', '19/54']},
        'arithmetic': {'bits': 112, 'value_denominator': DEN, 'time_grid_denominator': GRID,
                       'norm_arithmetic': 'Exact rational whole-cell powers and sums',
                       'fourth_moment_clock_bound': '[(1-(ell*right_x)^2)^2,1]'},
        'parameters': {'sigma_second_interval': sigma2, 'sigma_third': '4',
                       'physical_gain_interval': gain, 'physical_norm_prefactor_interval': k},
        'quadrature': {'domain_x': [0, 3], 'H_squared_integral_interval': j2,
                       'cosh_minus3_H_fourth_integral_interval': j4,
                       'inherited_C_interval': c,
                       'scalar_signed_integral_interval': scale(c, -4),
                       'scalar_L1_interval': scale(c, 4),
                       'scalar_L2_squared_interval': scalar_l2sq,
                       'constitutive_tension_integral_interval': tension_integral,
                       'constitutive_tension_L2_squared_interval': tension_l2sq},
        'tail_theorem': {'domain_x': [2, 3], 'ell2_scalar_monopole_interval': qt,
                         'published_strict_interval': ['-2.3e-19', '-6.7e-20'],
                         'interpretation': 'Negative angular mean throughout the continuous tail; field vanishes at Y=0.'},
        'cell_data': {'path': csv_path.name, 'rows': len(rows), 'sha256': sha256(csv_path.read_bytes()).hexdigest(),
                      'bytes': csv_path.stat().st_size, 'header': HEADER,
                      'values': 'Every coefficient endpoint is the stored integer divided by 2^112; coefficients are multiplied by ell^2.'},
        'input_pins': sources,
        'scope': {'explicit_material_tension_contribution': True,
                  'complete_bulk_coordinate_quadratic_source': False, 'full_U2_response_computed': False,
                  'total_tension_or_stress_second_variation': False, 'nonlinear_remainder_bound': False,
                  'novelty_claimed': False, 'cosmological_evidence': False}}
    (output/'QUADRATIC_FORCING_ENCLOSURE.json').write_text(json.dumps(encode(report), indent=2)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, default=HERE/'inputs')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = calculate(args.input_root, args.output_dir)
        print(json.dumps({'status': 'PASS', 'rows': result['cell_data']['rows']}))
    except (ValueError, KeyError, OSError) as exc:
        print('FAIL: '+str(exc))
        raise SystemExit(1)
