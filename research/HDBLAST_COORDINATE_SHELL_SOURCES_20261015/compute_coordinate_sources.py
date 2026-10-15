#!/usr/bin/env python3
"""Enclose the registered coordinate shell sources using exact rational cells.

This computes boundary inputs, not the second-order spacetime response.
"""
import argparse
import csv
from fractions import Fraction as F
import hashlib
import io
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
VALUE_DENOMINATOR = 2**144
FIELDS = ('y', 'yprime', 'u', 'utau', 'uy', 'v', 'vy', 'Sp', 'SK', 'EN',
          'dSp', 'dSK', 'dEN')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def interval(lo, hi=None):
    lo, hi = F(lo), F(lo if hi is None else hi)
    require(lo <= hi, 'interval order')
    return lo, hi


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def mul(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def scale(a, factor):
    return mul(a, interval(factor))


def square(a):
    return (F(0) if a[0] <= 0 <= a[1] else min(x*x for x in a),
            max(x*x for x in a))


def reciprocal(a):
    require(a[0] > 0, 'positive reciprocal interval')
    return 1/a[1], 1/a[0]


def maximum_absolute(a):
    return max(map(abs, a))


def outward(a):
    lo, hi = a[0]*VALUE_DENOMINATOR, a[1]*VALUE_DENOMINATOR
    return lo.numerator//lo.denominator, -((-hi.numerator)//hi.denominator)


def rounded(a):
    return tuple(F(v, VALUE_DENOMINATOR) for v in outward(a))


def polynomial_range(p0, p1, value, left, lo, hi):
    def evaluate(x):
        return p0*x-p1+value*(x-left)**2/2
    candidates = [evaluate(lo), evaluate(hi)]
    if value:
        vertex = left-p0/value
        if lo <= vertex <= hi:
            candidates.append(evaluate(vertex))
    return min(candidates), max(candidates)


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def display(value):
    if isinstance(value, F):
        return float(value)
    if isinstance(value, dict):
        return {key: display(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [display(item) for item in value]
    return value


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def atan_bounds(inverse, terms):
    total = sum((F((-1)**k, (2*k+1)*inverse**(2*k+1)) for k in range(terms)), F(0))
    following = total+F((-1)**terms, (2*terms+1)*inverse**(2*terms+1))
    return min(total, following), max(total, following)


def source_coefficients(u, v, ab, pb, sigma2):
    u2, v2, uv = square(u), square(v), mul(u, v)
    dSp = add(scale(mul(pb, u2), -4), scale(mul(sigma2, uv), 2))
    dSK = scale(mul(ab, u2), -12)
    dEN = add(scale(mul(ab, u2), -16), scale(mul(pb, uv), -F(16, 3)))
    material = {'Sp_material': scale(v2, -2),
                'SK_material': scale(mul(sigma2, v2), F(1, 6)),
                'EN_material': scale(mul(sigma2, v2), F(2, 3))}
    return dict(material, dSp=dSp, dSK=dSK, dEN=dEN,
                Sp=add(dSp, material['Sp_material']),
                SK=add(dSK, material['SK_material']),
                EN=add(dEN, material['EN_material']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'computed')
    args = parser.parse_args()
    pins_raw = (ROOT/'inputs/INPUT_PINS.json').read_bytes()
    pins = json.loads(pins_raw)
    raw = {}
    for entry in pins['files']:
        name = entry['path']
        require(name.startswith('inputs/') and '..' not in Path(name).parts,
                'canonical input path')
        value = (ROOT/name).read_bytes()
        require(len(value) == entry['bytes'] and sha(value) == entry['sha256'],
                'input authentication: '+name)
        raw[name] = value
    inherited = json.loads(raw['inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json'])
    shell = json.loads(raw['inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json'])
    require(inherited['state']['ell'] == '1/1000000', 'registered width')
    ell, L = F(1, 10**6), F(3)
    mu = interval(*inherited['shell_parameter_intervals']['mu'])
    def box(name):
        return interval(*shell['boxes'][name]['exact'])
    rho, r2, ab, pb, sigma2 = map(box, ('rho_b', 'H_squared', 'a_b', 'p_b', 'sigma_second_b'))
    require(mu[0] > 0 and rho[0] > 0 and r2[0] > 0 and ab[0] > 0
            and pb[1] < 0 and sigma2[0] > 0, 'background coefficient signs')
    shell_mu = box('mu')
    require(mu[0] <= shell_mu[0] <= shell_mu[1] <= mu[1], 'retained mu enclosure')
    M = F(5, 2)*ell**2
    q = M*L**2/2
    require((ell*L)**2 < F(1, 100) and 0 < q < 1, 'Volterra contraction and clock domain')
    forcing = interval(ell/mu[1], ell/mu[0])
    grid = inherited['cell_data']['time_grid_denominator']
    old_den = inherited['cell_data']['value_denominator']
    require(grid == 4096 and old_den == 2**112, 'retained exact grid')
    source_rows = list(csv.DictReader(io.StringIO(raw['inputs/prior/INERTIAL_WAVEFORM_CELLS.csv'].decode('ascii'))))
    require(len(source_rows) == 4353, 'complete inherited cell count')
    cells, p0, p1, last, BH = [], [F(0), F(0)], [F(0), F(0)], F(0), F(0)
    brackets = [interval(*row) for row in inherited['zero_existence_brackets']]
    bracket_F = [[] for _ in brackets]
    for row in source_rows:
        left, right = F(int(row['left_grid']), grid), F(int(row['right_grid']), grid)
        require(left == last and left < right, 'continuous exact cell partition')
        last = right
        H = interval(F(int(row['H_lo']), old_den), F(int(row['H_hi']), old_den))
        v = interval(F(int(row['ell_S_lo']), old_den)/ell, F(int(row['ell_S_hi']), old_den)/ell)
        H_integral = (p0[0]+min(F(0), H[0]*(right-left)),
                      p0[1]+max(F(0), H[1]*(right-left)))
        Fcell = mul(forcing, (polynomial_range(p0[0], p1[0], H[0], left, left, right)[0],
                              polynomial_range(p0[1], p1[1], H[1], left, left, right)[1]))
        for i, bracket in enumerate(brackets):
            lo, hi = max(left, bracket[0]), min(right, bracket[1])
            if lo <= hi:
                bracket_F[i].append(mul(forcing,
                    (polynomial_range(p0[0], p1[0], H[0], left, lo, hi)[0],
                     polynomial_range(p0[1], p1[1], H[1], left, lo, hi)[1])))
        cells.append({'left': left, 'right': right, 'F': Fcell,
                      'Fprime': mul(forcing, H_integral), 'v': v})
        BH += maximum_absolute(H)*((L-left)**2-(L-right)**2)/2
        for side in (0, 1):
            p0[side] += H[side]*(right-left)
            p1[side] += H[side]*(right**2-left**2)/2
    require(last == L, 'full registered observation window')
    Fsup = max(maximum_absolute(cell['F']) for cell in cells)
    Ybound = Fsup/(1-q)
    E, Eprime = q*Ybound, M*L*Ybound
    Ycoarse = forcing[1]*BH/(1-q)
    Ucoarse, Usigned = r2[1]*Ycoarse, r2[1]*Ybound
    pi = rounded(add(scale(atan_bounds(5, 64), 16), scale(atan_bounds(239, 16), -4)))
    pi2 = square(pi)
    radial_measure = scale(mul(square(square(rho)), pi2), ell)
    L2_prefactor, L1_prefactor = scale(radial_measure, F(3, 20)), scale(radial_measure, F(1, 3))
    EN_L2_prefactor = scale(radial_measure, F(19, 360))
    names = ('Sp', 'SK', 'EN', 'dSp', 'dSK', 'dEN', 'Sp_material', 'SK_material', 'EN_material')
    sums = {name: {'I2': interval(0), 'I1': interval(0), 'signed': interval(0)} for name in names}
    global_envelopes = {}
    output = io.StringIO(newline='')
    writer = csv.DictWriter(output, fieldnames=['left_grid', 'right_grid']+[f'{name}_{side}' for name in FIELDS for side in ('lo', 'hi')], lineterminator='\n')
    writer.writeheader()
    for cell in cells:
        left, right = cell['left'], cell['right']
        clock = interval(1-(ell*right)**2, 1)
        tanh = interval(0, ell*right)
        y, yp = add(cell['F'], interval(-E, E)), add(cell['Fprime'], interval(-Eprime, Eprime))
        u, v = mul(mul(r2, clock), y), cell['v']
        utau = mul(scale(mul(mul(r2, clock), reciprocal(rho)), 1/ell),
                    add(yp, scale(mul(tanh, y), -F(3, 2)*ell)))
        uy = add(scale(mul(ab, u), -2), scale(mul(pb, v), -F(1, 3)))
        vy = add(scale(mul(pb, u), -2), scale(mul(sigma2, v), -F(1, 2)))
        coeff = source_coefficients(u, v, ab, pb, sigma2)
        values = dict(coeff, y=y, yprime=yp, u=u, utau=utau, uy=uy, v=v, vy=vy)
        record = {'left_grid': int(left*grid), 'right_grid': int(right*grid)}
        for name in FIELDS:
            lo, hi = outward(values[name])
            record[name+'_lo'], record[name+'_hi'] = lo, hi
            enclosed = interval(F(lo, VALUE_DENOMINATOR), F(hi, VALUE_DENOMINATOR))
            old = global_envelopes.get(name, enclosed)
            global_envelopes[name] = min(old[0], enclosed[0]), max(old[1], enclosed[1])
        writer.writerow(record)
        weight = interval(1, 1+2*(ell*right)**2)
        width = right-left
        for name in names:
            value = rounded(coeff[name])
            sq = square(value)
            absolute = interval(0 if value[0] <= 0 <= value[1] else min(map(abs, value)), maximum_absolute(value))
            sums[name]['I2'] = add(sums[name]['I2'], scale(mul(weight, sq), width))
            sums[name]['I1'] = add(sums[name]['I1'], scale(mul(weight, absolute), width))
            sums[name]['signed'] = add(sums[name]['signed'], scale(mul(weight, value), width))
    norms = {}
    for name in names:
        degree_four = name in ('EN', 'dEN', 'EN_material')
        norms[name] = {'angular_field': 'H4' if degree_four else 'Y^2',
                       'weighted_x_squared_integral': sums[name]['I2'],
                       'squared_L2': mul(EN_L2_prefactor if degree_four else L2_prefactor, sums[name]['I2'])}
        if not degree_four:
            norms[name].update(L1=mul(L1_prefactor, sums[name]['I1']),
                               signed_integral=mul(L1_prefactor, sums[name]['signed']))
    ratios = {name: norms['d'+name]['squared_L2'][1]/norms[name+'_material']['squared_L2'][0]
              for name in ('Sp', 'SK', 'EN')}
    central = []
    for bracket, Fs in zip(brackets, bracket_F):
        Fc = interval(min(z[0] for z in Fs), max(z[1] for z in Fs))
        yc = add(Fc, interval(-E, E))
        uc = mul(mul(r2, interval(1-(ell*bracket[1])**2, 1)), yc)
        require(uc[0] > 0, 'positive metric trace throughout scalar zero bracket')
        source_zero = source_coefficients(uc, interval(0), ab, pb, sigma2)
        require(source_zero['Sp'][0] > 0 and source_zero['SK'][1] < 0,
                'opposite complete coordinate source signs at actual scalar zeros')
        central.append({'bracket': bracket, 'F': Fc, 'y': yc, 'u': uc,
                        'Sp_at_scalar_zero': source_zero['Sp'], 'SK_at_scalar_zero': source_zero['SK'],
                        'EN_at_scalar_zero': source_zero['EN']})
    gain = interval(*inherited['physical_gain_interval'])
    Htail = interval(*inherited['tail_theorem']['actual_H_interval'])
    vtail = scale(mul(mul(gain, interval(1-(ell*L)**2, 1)), Htail), 1/ell)
    tail = source_coefficients(interval(-Usigned, Usigned), vtail, ab, pb, sigma2)
    require(tail['Sp'][1] < 0 < tail['SK'][0] and tail['EN'][0] > 0, 'uniform tail source signs')
    off = 4*maximum_absolute(pb)
    sym_lower = min(6*ab[0]-off, 4*sigma2[0]-off)
    determinant = add(scale(mul(ab, sigma2), 24), scale(square(pb), -16))
    require(sym_lower > 0 and determinant[0] > 0, 'positive registered symmetrized Robin matrix')
    exact = {'ell': ell, 'L': L, 'mu': mu, 'rho': rho, 'r2': r2, 'a': ab, 'p': pb,
             'sigma2': sigma2, 'M': M, 'q': q, 'Fsup': Fsup, 'Ybound': Ybound,
             'E': E, 'Eprime': Eprime, 'BH': BH, 'Ycoarse': Ycoarse,
             'Ucoarse': Ucoarse, 'Usigned': Usigned, 'pi': pi,
             'L2_prefactor': L2_prefactor, 'L1_prefactor': L1_prefactor,
             'EN_L2_prefactor': EN_L2_prefactor,
             'global_envelopes': global_envelopes, 'central': central,
             'tail': {'domain': interval(2, 3), 'v': vtail, 'coefficients': tail},
             'norms': norms, 'metric_correction_squared_norm_ratio_upper': ratios,
             'symmetrizer': {'weights': [3, 8], 'gershgorin_positive_lower': sym_lower,
                            'determinant_interval': determinant},
             'harmonics': {'basis': ['1', 'Z2', 'H4'], 'source_coefficient_factors': [F(1, 6), F(1, 2), F(1)],
                           'squared_L2_power_shares': [F(10, 27), F(5, 18), F(19, 54)],
                           'EN_scope': 'Retained degree-four scalar Robin row only; not a degree-zero master reduction.'}}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_raw = output.getvalue().encode('ascii')
    csv_name = 'COORDINATE_SHELL_SOURCE_CELLS.csv'
    (args.output_dir/csv_name).write_bytes(csv_raw)
    report = {'schema': 'hdblast-coordinate-shell-sources-enclosure-v1', 'status': 'PASS',
              'runtime': {'python': sys.version.split()[0], 'optimization': sys.flags.optimize},
              'producer_sha256': sha(Path(__file__).read_bytes()),
              'input_pins_sha256': sha(pins_raw),
              'cell_data': {'path': csv_name, 'rows': len(cells), 'time_grid_denominator': grid,
                            'value_denominator': VALUE_DENOMINATOR, 'fields': list(FIELDS),
                            'bytes': len(csv_raw), 'sha256': sha(csv_raw)},
              'exact': encode(exact),
              'display_only': display({'global_envelopes': global_envelopes, 'Usigned': Usigned,
                                       'Ucoarse': Ucoarse, 'norms': {k: v['squared_L2'] for k,v in norms.items()},
                                       'metric_correction_squared_norm_ratio_upper': ratios,
                                       'central': central, 'tail': exact['tail'], 'symmetrizer': exact['symmetrizer']}),
              'scope': {'conditional_on_registered_analytic_inputs': True,
                        'first_representative': 'Retained longitudinal, zero first shell graph and fixed intrinsic identification.',
                        'normal_convention': 'One-copy outward coordinate physical rows; not the outward-normal sum material scalar row.',
                        'taylor_convention': 'Derivative coefficients U=U0+epsilon U1+epsilon^2 U2/2.',
                        'metric_and_scalar_sources': 'Sp and SK multiply Y^2 in the original coordinate scalar and isotropic Israel rows.',
                        'EN_scope': 'Retained degree-four reduced Robin forcing only.',
                        'norms': 'Scalar amplitude norms in background proper shell measure, not an Einstein energy.',
                        'full_U2_evaluated': False, 'nonlinear_remainder_established': False,
                        'coordinate_source_sign_is_gauge_invariant': False, 'cosmological_evidence': False}}
    (args.output_dir/'COORDINATE_SHELL_SOURCE_ENCLOSURE.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'rows': len(cells), 'csv_sha256': sha(csv_raw),
                      'metric_bound': float(Usigned), 'Sp_tail_upper': float(tail['Sp'][1]),
                      'SK_tail_lower': float(tail['SK'][0]), 'correction_squared_norm_ratio_upper': display(ratios)}))


if __name__ == '__main__':
    main()
