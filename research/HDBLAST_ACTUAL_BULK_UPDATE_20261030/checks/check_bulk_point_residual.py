#!/usr/bin/env python3
"""Enclose one actual reduced bulk residual of the prior reconstruction.

The affine source is bounded by authenticated uniform component caps. The
actual coefficients are assembled at the point; this is not a point evaluation
of the affine source or an Einstein/constraint residual.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return raw


def interval(z):
    need(len(z) == 2, 'interval shape')
    a, b = map(F, z)
    need(a <= b, 'interval order')
    return a, b


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def neg(a):
    return -a[1], -a[0]


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    z = [x*y for x in a for y in b]
    return min(z), max(z)


def scale(a, b):
    return mul(a, (b, b))


def total(rows):
    out = (F(0), F(0))
    for z in rows:
        out = add(out, z)
    return out


def mv(matrix, vector):
    return [total(mul(c, v) for c, v in zip(row, vector)) for row in matrix]


def frozen_mv(matrix, vector):
    return [total(scale(v, c) for c, v in zip(row, vector)) for row in matrix]


def encode(row):
    return list(map(str, row))


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', required=True, type=Path)
ap.add_argument('--robin-capsule', type=Path)
ap.add_argument('--contact-capsule', type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
robin = args.robin_capsule or root/'research/HDBLAST_ROBIN_RECONSTRUCTION_20261030'
contact = args.contact_capsule or root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030'
driver = base/'compute/full23_response.py'
DRIVER_PIN = 'efa6a6c2a9c7db596fe024246684f760b807e727fd2ddac8e8126182ba941a73'
pin(driver, DRIVER_PIN)
oldpoint = json.loads(pin(contact/'reference/SHELL.json',
    '9de23a17dfb96361e4d5ad46f31b60cb0da458492d6e21d26c9dd6890e0b1d64'))
point = json.loads(pin(robin/'reference/RESUMMED.json',
    'bc7a0ef6545ddf09e54f16a71f4ac3dc0cf8f1a0a85d37d78a44bf1845ed1f69'))
profile = json.loads(pin(base/'compute/response-normal/RESPONSE_PROFILES.json',
    '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc'))
source_path = root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_SOURCES_NORMAL.json'
SOURCE_PIN = '45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f'
source = json.loads(pin(source_path, SOURCE_PIN))
map_path = root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json'
MAP_PIN = '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'
map_bytes = pin(map_path, MAP_PIN)
need(source['source_map'] == {'bytes': len(map_bytes), 'sha256': MAP_PIN,
    'status': 'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION'},
    'complete geometric source-map premise')
need(point['status'] == oldpoint['status'] == 'PASS' and
     point['point'] == oldpoint['point'] == ['9/8', '0'], 'prior successful point enclosures')
need(source['status'] == 'PASS_CONTINUOUS_ACTUAL_EVOLVED_GEOMETRIC_SOURCE_CELL_ENCLOSURES',
     'actual affine source scope')
need(source['scope']['actual_evolved_common_gauge_bulk_sources_evaluated'] is True,
     'common-gauge source premise')
need(source['domain'] == {'ell':'1/1000000',
    'causal_triangle':'0<=x<=3,r>=0,x+r<=8',
    'angular':'entire unit S3; exact unnormalized registered scalar/vector/tensor coefficient basis'},
    'uniform source domain')
need(source['normalization'] ==
    'Second parameter derivatives along the affine actual-first-pair path; metric_forcing=-2 ell^2 rho^2 QEred; scalar_forcing=ell^2 rho^2 Qphi. QC is the actual-metric-lowered quadratic gauge Hessian.',
    'second-parameter forcing normalization')

# The pinned driver imports only standard library at module scope. Its Response
# first authenticates local dependencies; Context then authenticates all seven
# complete scientific roots before the scientific imports and primitive assembly.
spec = importlib.util.spec_from_file_location('bulk_point_response_driver', driver)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
response = module.Response(root)
ctx = response.context
labels, gram = list(response.labels), list(response.weights)
need(labels == profile['labels'] == point['labels'] == oldpoint['labels'] and
     gram == list(map(F, profile['gram'])) and len(labels) == 23 and min(gram) > 0,
     'same registered modal space')
x, r = F(9, 8), F(0)
actual = response.assemble(response.I(x), response.I(r))
actual = {key: [[(z.lo, z.hi) for z in row] for row in actual[key]]
          for key in ('A', 'B', 'D')}
frozen = {key: [list(map(F, row)) for row in profile['matrices'][key]]
          for key in ('A', 'B', 'D')}

# Map the source's orthogonal physical field blocks to the registered full23
# amplitudes. Exact diagonal Gram and multiplicity equality prove the H norm
# is the same norm, rather than merely assuming that an old response cap is F.
summary = source['global_cover']['summary']
caps = summary['absolute_coefficient_caps']
Hb = summary['positive_bundle_bulk_forcing_H_bound']
mapping = [('yy', 'metric_forcing/yy', ('K:0', 'K:2', 'K:4'), 1),
           ('y0', 'metric_forcing/y0', ('L:0', 'L:2', 'L:4'), 2),
           ('00', 'metric_forcing/00', ('N:0', 'N:2', 'N:4'), 1),
           ('ya', 'metric_forcing/ya', ('v:2', 'v:4'), 2),
           ('0a', 'metric_forcing/0a', ('u:2', 'u:4'), 2),
           ('ab', 'metric_forcing/ab', ('J:0', 'J:2', 'J:4', 'H:2', 'H:4', 'T:+', 'T:-'), 1),
           ('scalar', 'scalar_forcing', ('w:0', 'w:2', 'w:4'), 1)]
Fcap = [F(0)]*23
square = F(0)
gram_entries = 0
fields = {}
for field, prefix, blocklabels, multiplicity in mapping:
    block = Hb['block_Gram_certificates'][field]
    gg = [list(map(F, row)) for row in block['gram']]
    count = len(blocklabels)
    need(block['multiplicity'] == multiplicity and len(gg) == count and
         all(len(row) == count for row in gg), 'source block shape and multiplicity')
    componentcaps = [F(caps[prefix+'/'+str(j)]) for j in range(count)]
    need(min(componentcaps) >= 0, 'nonnegative source magnitude caps')
    fields[field] = componentcaps
    for i, label in enumerate(blocklabels):
        Fcap[labels.index(label)] = componentcaps[i]
        for j in range(count):
            need(multiplicity*gg[i][j] ==
                 (gram[labels.index(label)] if i == j else 0),
                 'exact source/current H Gram identity')
            gram_entries += 1
    blocksquare = multiplicity*sum(abs(gg[i][j])*componentcaps[i]*componentcaps[j]
                                  for i in range(count) for j in range(count))
    need(blocksquare == F(block['squared_norm_upper']), 'source block norm square')
    square += blocksquare
need(ctx.op.pack_fields(fields) == Fcap, 'source packing agrees with authenticated operator')
need(square == F(Hb['squared_upper']) == sum(w*c*c for w, c in zip(gram, Fcap)),
     'whole positive bundle and registered H norm equality')
FH = F(Hb['F_H_upper'])
need(FH == 17387513156573 and (FH-1)**2 <= square < FH**2,
     'genuine uniform pointwise affine F norm bound')

oldrows, rows = oldpoint['full23_point_table'], point['full23_point_table']
need(len(oldrows) == len(rows) == 23, 'complete point tables')
Bs = {key: [interval(row['B0'][key]) for row in oldrows]
      for key in ('value', 'dx', 'dr')}
W = {key: [interval(row['reconstructed_boundary_candidate'][key]) for row in rows]
     for key in ('value', 'dx', 'dr')}
for i in range(23):
    need(oldrows[i]['label'] == rows[i]['label'] == labels[i], 'point row labels')

f0_terms = [frozen_mv(frozen[key], Bs[jet])
            for key, jet in (('A', 'dx'), ('B', 'dr'), ('D', 'value'))]
actual_terms = [mv(actual[key], W[jet])
                for key, jet in (('A', 'dx'), ('B', 'dr'), ('D', 'value'))]
f0 = [total(term[i] for term in f0_terms) for i in range(23)]
LaW = [total(term[i] for term in actual_terms) for i in range(23)]
direct = [sub(LaW[i], f0[i]) for i in range(23)]
# Recombine as coefficient differences and W-Bs before applying matrices. This
# avoids subtracting the two enormous free-wave source terms independently.
residual_terms = []
for key, jet in (('A', 'dx'), ('B', 'dr'), ('D', 'value')):
    delta = [[sub(a, (b, b)) for a, b in zip(ar, br)]
             for ar, br in zip(actual[key], frozen[key])]
    residual_terms.append(mv(delta, Bs[jet]))
    residual_terms.append(mv(actual[key], [sub(a, b) for a, b in zip(W[jet], Bs[jet])]))
lower = [total(term[i] for term in residual_terms) for i in range(23)]
out = []
for i, label in enumerate(labels):
    lo, hi = max(lower[i][0], direct[i][0]), min(lower[i][1], direct[i][1])
    need(lo <= hi, 'two exact residual identities have overlapping enclosures')
    affine = (-Fcap[i], Fcap[i])
    complete = add((lo, hi), affine)
    sign = 'positive' if complete[0] > 0 else 'negative' if complete[1] < 0 else 'undetermined'
    out.append({'label': label,
        'frozen_free_source_f0': encode(f0[i]),
        'actual_lower_operator_LaW': encode(LaW[i]),
        'lower_residual_cancellation_free': encode(lower[i]),
        'lower_residual_direct': encode(direct[i]),
        'lower_residual_intersection': encode((lo, hi)),
        'affine_bulk_source_uniform_interval': encode(affine),
        'complete_reduced_bulk_residual': encode(complete),
        'complete_reduced_bulk_residual_sign': sign})
mean = out[12]
need(mean['label'] == 'w:0' and mean['complete_reduced_bulk_residual_sign'] == 'negative',
     'strict negative registered scalar-mean forcing defect')
mean_bounds = interval(mean['complete_reduced_bulk_residual'])
need(-F(9417*10**19) < mean_bounds[0] <= mean_bounds[1] < -F(9416*10**19),
     'readable strict scalar-mean residual bounds')
lower_square = upper_square = F(0)
for weight, row in zip(gram, out):
    lo, hi = interval(row['complete_reduced_bulk_residual'])
    minimum = F(0) if lo <= 0 <= hi else min(abs(lo), abs(hi))
    lower_square += weight*minimum**2
    upper_square += weight*max(abs(lo), abs(hi))**2
need(lower_square > F(9*10**22)**2, 'nonzero whole residual H norm lower bound')
ctx.reauthenticate()
module.need(module.authenticate_local() == response.local_before, 'local dependency readback')
pin(driver, DRIVER_PIN)
pin(source_path, SOURCE_PIN)
pin(map_path, MAP_PIN)
print(json.dumps({'status':'PASS','point':[str(x),str(r)],
    'residual_convention':'r0=La W0+F-L0 B_s, the forcing defect used by the first actual-bulk update',
    'actual_coefficients':'point-assembled actual background/tanh/sech interval coefficients, not global coefficient boxes',
    'affine_source_method':'authenticated uniform actual common-gauge component bounds; F is not pointwise evaluated',
    'affine_source_global_domain':source['domain'],
    'affine_F_H_uniform_cap':str(FH),
    'source_current_Gram_entries_compared':gram_entries,
    'source_H_norm_square':str(square),
    'complete_residual_H_norm_squared_lower':str(lower_square),
    'complete_residual_H_norm_squared_upper':str(upper_square),
    'strict_primary_residual_bounds':['-9.417e22','-9.416e22'],
    'nonzero_residual_rows':[row['label'] for row in out if row['complete_reduced_bulk_residual_sign']!='undetermined'],
    'full23_point_table':out,
    'input_sha256':{'coefficient_driver':DRIVER_PIN,'actual_affine_source':SOURCE_PIN,
                    'complete_source_map':MAP_PIN},
    'scope':'one conditional registered reduced-evolution residual enclosure of the prior W0, including a uniform affine F bound; not a full pointwise Einstein/constraint residual, bulk-update evaluation, C2 error certificate, nonlinear family or physical observable'},sort_keys=True))
