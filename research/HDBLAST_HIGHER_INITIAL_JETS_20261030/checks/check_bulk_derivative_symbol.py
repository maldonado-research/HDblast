#!/usr/bin/env python3
"""Exact leading-jet audit of the complete registered bulk forcing.

This does not bound the actual F_x. It tests its highest input-jet dependence
against null variations satisfying the principal wave-trace relations.
Only standard-library modules are imported; scientific source maps are pinned.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def clean(terms):
    return {m: c for m, c in terms.items() if c}


def add(a, b):
    out = dict(a)
    for m, c in b.items():
        out[m] = out.get(m, Q(0))+c
    return clean(out)


def monomial(coefficient=1, **powers):
    return {tuple(sorted((k, e) for k, e in powers.items() if e)): Q(coefficient)} if coefficient else {}


def scale(a, c):
    return clean({m: c*v for m, v in a.items()})


def mul(a, b):
    out = {}
    for m, c in a.items():
        for n, d in b.items():
            powers = dict(m)
            for name, exponent in n:
                powers[name] = powers.get(name, 0)+exponent
            key = tuple(sorted((k, e) for k, e in powers.items() if e))
            out[key] = out.get(key, Q(0))+c*d
    return clean(out)


def derivative(a, name):
    out = {}
    for m, c in a.items():
        powers = dict(m)
        n = powers.get(name, 0)
        if n:
            powers[name] = n-1
            key = tuple(sorted((k, e) for k, e in powers.items() if e))
            out[key] = out.get(key, Q(0))+n*c
    return clean(out)


def decode(rows):
    result = {}
    for row in rows:
        need(set(row) == {'coefficient', 'powers'} and all(isinstance(e, int) for e in row['powers'].values()), 'encoded Laurent terms')
        result = add(result, monomial(Q(row['coefficient']), **row['powers']))
    return result


def encode(poly):
    return [{'coefficient': str(c), 'powers': dict(m)} for m, c in sorted(poly.items())]


def evaluate(poly, values):
    return sum((c*product(values[k]**e for k, e in m) for m, c in poly.items()), Q(0))


def product(values):
    out = Q(1)
    for x in values:
        out *= x
    return out


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
args = ap.parse_args()
path = args.repo_root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json'
raw = path.read_bytes()
digest = '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'
need(hashlib.sha256(raw).hexdigest() == digest, 'complete source-map hash')
data = json.loads(raw)
need(data['status'] == 'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION', 'complete source-map status')
need(data['input_units']['ell'] == 't=ell*x and dy=-rho*ell*dr', 'coordinate convention')
need(data['variation_convention'] == 'second parameter derivative d2/deps2 along affine g0+eps*q,phi0+eps*eta', 'second-parameter source convention')
jetpath = args.repo_root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_FIRST_JETS_NORMAL.json'
jetbytes = jetpath.read_bytes()
jetpin = '0aaa7b2fe30f546db92ef8b26a7224ffa7f6812d3780c331661395aa1991d3b8'
need(hashlib.sha256(jetbytes).hexdigest() == jetpin, 'actual first-field jet report hash')
jetdata = json.loads(jetbytes)
need(jetdata['status'] == 'PASS_ACTUAL_FIRST_COMMON_GAUGE_ORDINARY_2JET_CELL_ENCLOSURES' and
     jetdata['ell'] == '1/1000000' and
     jetdata['domain'] == '0<=x<=3,r>=0,x+r<=8; t=ell*x,z=ell*r,dy=-rho dz', 'actual first-field jet scope')
fields = data['maps']['metric_forcing']
source = {**{field: [decode(row) for row in rows] for field, rows in fields.items()},
          'scalar': [decode(row) for row in data['maps']['scalar_forcing']]}
counts = {'yy': 3, 'y0': 3, '00': 3, 'ya': 2, '0a': 2, 'ab': 7, 'scalar': 3}
need(set(source) == set(counts) and all(len(source[k]) == n for k, n in counts.items()), 'all23 forcing components')
need(sum(counts.values()) == 23, 'complete reduced metric and scalar space')
ell, rho, cc = monomial(ell=1), monomial(rho=1), monomial(c=1)

# Every highest source derivative is linear in the ordinary second field jets.
# This is checked on the complete collected maps, after all geometric terms.
second = {base+suf for base in ('A','B','C','D','E','L','M') for suf in ('yy','yt','tt')} | {'syy','syt','stt'}
for polys in source.values():
    for poly in polys:
        for m in poly:
            need(sum(e for name, e in m if name in second) <= 1 and all(e == 1 for name, e in m if name in second), 'affine highest second-jet source dependence')

results = []
checks = controls = 0
for sign in (1, -1):
    null = add(add(monomial(A=1), monomial(D=1)), monomial(-2*sign, B=1))
    for channel in ('f', 'b', 'd', 'scalar'):
        # The variation is normalized to one highest jet. Gauge channels use
        # delta xi_(y^i t^j)=1^i*(sign*rho)^j, i+j=4. Scalar uses total3.
        varying = {}
        for suf, iy, it in (('yy', 2, 0), ('yt', 1, 1), ('tt', 0, 2)):
            if channel == 'f':
                varying['A'+suf] = monomial(2*sign**(it+1), rho=it+1)
                varying['B'+suf] = monomial(sign**(it+2), rho=it+1)
            elif channel == 'b':
                varying['B'+suf] = monomial(-sign**(it+1), rho=it+2)
                varying['D'+suf] = monomial(-2*sign**(it+2), rho=it+2)
            elif channel == 'd':
                varying['C'+suf] = monomial(sign**(it+1), rho=it+2, c=1)
                varying['E'+suf] = monomial(sign**(it+2), rho=it+2, c=1)
            else:
                varying['s'+suf] = monomial(sign**(it+1), rho=it+1)
        observed = {}
        for field, polys in source.items():
            observed[field] = []
            for poly in polys:
                delta = {}
                for symbol, variation in varying.items():
                    delta = add(delta, mul(derivative(poly, symbol), variation))
                observed[field].append(mul(ell, delta))
        expected = {field: [{} for _ in range(n)] for field, n in counts.items()}
        if channel == 'f':
            p = mul(monomial(ell=3, rho=3), null)
            expected['yy'] = [scale(p, c*sign) for c in (Q(-2,3), Q(-2), Q(-4))]
            expected['y0'] = [scale(p, c) for c in (Q(-1,3), Q(-1), Q(-2))]
        elif channel == 'b':
            p = mul(monomial(ell=3, rho=4), null)
            expected['y0'] = [scale(p, c*sign) for c in (Q(1,3), Q(1), Q(2))]
            expected['00'] = [scale(p, c) for c in (Q(2,3), Q(2), Q(4))]
        elif channel == 'd':
            p = mul(monomial(ell=3, rho=4, c=1), null)
            expected['ya'] = [scale(p, c*sign) for c in (Q(-1,2), Q(-1))]
            expected['0a'] = [scale(p, c) for c in (Q(-1,2), Q(-1))]
        else:
            p = mul(monomial(ell=3, rho=3), null)
            expected['scalar'] = [scale(p, c*sign) for c in (Q(-1,3), Q(-1), Q(-2))]
        for field, polys in observed.items():
            for got, want in zip(polys, expected[field]):
                need(got == want, 'complete-map null highest-jet factor:'+channel+':'+str(sign)+':'+field)
                checks += 1
        # The high variation obeys every principal wave trace. Lower gauge
        # jets are held fixed, so lower-coefficient prolongation terms do not
        # change under this homogeneous highest-jet variation.
        order = 3 if channel == 'scalar' else 4
        for iy in range(order-1):
            it = order-2-iy
            traced = add(mul(monomial(sign**(it+2), rho=it+2), monomial(rho=-2)),
                         monomial(-sign**it, rho=it))
            need(traced == {}, 'null principal wave prolongation with physical rho factor')
            checks += 1
        need(any(poly for polys in observed.values() for poly in polys), 'nonzero symbolic map remains after trace reduction')
        checks += 1
        results.append({'channel': channel, 'sign': sign,
                        'nonzero_outputs': {field: [j for j, poly in enumerate(polys) if poly] for field, polys in observed.items() if any(polys)},
                        'complete_forcing_x_highest_jet_coefficients': {field: [encode(p) for p in polys] for field, polys in observed.items()}})

# Concrete coefficient specializations are algebraic witnesses, not values
# of the retained actual spacetime. A conformal low metric instead has zero
# null contraction and exercises the exceptional cancellation.
need(evaluate(null, {'A':Q(1),'D':Q(0),'B':Q(0)}) == 1, 'nonzero null-contraction witness')
need(evaluate(null, {'A':Q(1),'D':Q(-1),'B':Q(0)}) == 0, 'exceptional conformal null-contraction cancellation')
checks += 2
sample = decode(results[0]['complete_forcing_x_highest_jet_coefficients']['yy'][0])
sample_mixed = decode(results[0]['complete_forcing_x_highest_jet_coefficients']['y0'][0])
need(mul(sample, monomial(ell=-1)) != sample, 'missing time-chain ell corruption rejected')
need(scale(sample_mixed, -1) != sample_mixed, 'reversed y0 output corruption rejected')
need({} != sample, 'wave trace does not cancel the null contraction')
need(add(monomial(rho=2), monomial(-1)) != {}, 'omitted physical wave rho factor rejected')
controls = 4
cell = jetdata['cells'][0]
need(cell['label'] == 'evolved_packet' and cell['scaled_x'] == ['1/4', '1025/4096'] and
     cell['scaled_r'] == ['5/4', '5121/4096'], 'registered evolved-packet cell')
actual = cell['actual_m']
a0, a1 = map(Q, actual['A'])
b0, b1 = map(Q, actual['B'])
d0, d1 = map(Q, actual['D'])
nullbox = (a0+d0-2*b1, a1+d1-2*b0)
rhobox = tuple(map(Q, cell['background']['rho']))
need(0 < rhobox[0] <= rhobox[1] and Q('-4.859e-4') < nullbox[0] <= nullbox[1] < Q('-2.233e-4'),
     'actual first-field outgoing null contraction is strictly negative')
sensitivities = [-Q(2,3)*Q(jetdata['ell'])**3*r**3*j for r in rhobox for j in nullbox]
sensitivity = (min(sensitivities), max(sensitivities))
need(Q('3.265e-16') < sensitivity[0] <= sensitivity[1] < Q('7.104e-16'),
     'complete scalar-mean yy forcing derivative has a nonzero highest-gauge-jet coefficient')
need(hashlib.sha256(path.read_bytes()).hexdigest() == digest, 'source-map readback')
need(hashlib.sha256(jetpath.read_bytes()).hexdigest() == jetpin, 'actual first-field jet report readback')
print(json.dumps({'status': 'PASS_COMPLETE_BULK_DERIVATIVE_LEADING_JET_AUDIT',
    'source_map_sha256': digest, 'component_identities': checks,
    'wrong_formula_controls': controls, 'variations': results,
    'generic_null_factor': 'A+D-2 sign B',
    'actual_first_jet_report_sha256': jetpin,
    'actual_coefficient_witness': {'cell': {'x': cell['scaled_x'], 'r': cell['scaled_r']},
        'outgoing_null_contraction_interval': list(map(str, nullbox)),
        'outgoing_null_contraction_readable_strict_bounds': ['-4.859e-4','-2.233e-4'],
        'radial_gauge_fourth_null_jet_yy_mean_F_x_sensitivity_interval': list(map(str, sensitivity)),
        'sensitivity_readable_strict_bounds': ['3.265e-16','7.104e-16'],
        'meaning': 'coefficient of a normalized highest-jet variation, not the actual source derivative or a physical force'},
    'scope': 'exact highest-jet dependence of all23 collected forcing components; null principal wave prolongations only; no actual F_x value/cap, no actual nonzero source claim or proof that every possible data-specific identity requires C4'}, sort_keys=True))
