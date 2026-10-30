#!/usr/bin/env python3
"""Independent arithmetic/source-map audit of one reduced forcing defect.

The coefficient assembler and inherited scientific inputs are authenticated.
No new residual evaluator is imported. A third exact residual decomposition
is compared with the independently produced interval intersection.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def read(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return json.loads(raw)


def iv(value):
    need(len(value) == 2, 'interval shape')
    out = tuple(map(Q, value))
    need(out[0] <= out[1], 'interval order')
    return out


def plus(a, b):
    return a[0] + b[0], a[1] + b[1]


def difference(a, b):
    return a[0] - b[1], a[1] - b[0]


def product(a, b):
    vals = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    return min(vals), max(vals)


def summed(terms):
    result = Q(0), Q(0)
    for term in terms:
        result = plus(result, term)
    return result


def matvec(matrix, values):
    return [summed(product(a, b) for a, b in zip(row, values)) for row in matrix]


def modal_weight(name):
    if name.startswith('T:'):
        return Q(1, 6)
    field, degree = name.split(':')
    degree = int(degree)
    w = {0: Q(1), 2: Q(1, 12), 4: Q(19, 720)}[degree]
    if field == 'L': w *= 2
    if field == 'J': w *= 3
    if field in ('v', 'u'): w *= 2*degree*(degree+2)
    if field == 'H': w *= Q(2, 3)*degree*(degree+2)*(degree*(degree+2)-3)
    return w


def enc(value):
    return list(map(str, value))


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--robin-capsule', type=Path)
ap.add_argument('--contact-capsule', type=Path)
ap.add_argument('--residual-report', type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
robin = args.robin_capsule or root/'research/HDBLAST_ROBIN_RECONSTRUCTION_20261030'
contact = args.contact_capsule or root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030'
old = read(contact/'reference/SHELL.json',
           '9de23a17dfb96361e4d5ad46f31b60cb0da458492d6e21d26c9dd6890e0b1d64')
wave = read(robin/'reference/RESUMMED.json',
            'bc7a0ef6545ddf09e54f16a71f4ac3dc0cf8f1a0a85d37d78a44bf1845ed1f69')
profiles = read(base/'compute/response-normal/RESPONSE_PROFILES.json',
                '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
sourcepath = root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_SOURCES_NORMAL.json'
sourcepin = '45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f'
source = read(sourcepath, sourcepin)
reportpath = args.residual_report or Path(__file__).resolve().parents[1]/'reference/POINT_RESIDUAL.json'
reportpin = '30dc43f1ee1f1830976fefffc589beed42431d8013a1218094f9ecde0327ec69'
reported = read(reportpath, reportpin)
need(all(x['status'] == 'PASS' and x['point'] == ['9/8', '0']
         for x in (old, wave, reported)), 'same successful point')
need(source['status'] == 'PASS_CONTINUOUS_ACTUAL_EVOLVED_GEOMETRIC_SOURCE_CELL_ENCLOSURES'
     and source['scope']['actual_evolved_common_gauge_bulk_sources_evaluated'],
     'inherited actual complete source scope')
need(source['normalization'] ==
     'Second parameter derivatives along the affine actual-first-pair path; metric_forcing=-2 ell^2 rho^2 QEred; scalar_forcing=ell^2 rho^2 Qphi. QC is the actual-metric-lowered quadratic gauge Hessian.',
     'metric/scalar forcing normalization')
need(source['domain']['ell'] == '1/1000000' and
     source['domain']['causal_triangle'] == '0<=x<=3,r>=0,x+r<=8',
     'uniform pointwise source domain covers sample')
map_path = root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json'
map_digest = '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'
map_raw = map_path.read_bytes()
need(hashlib.sha256(map_raw).hexdigest() == map_digest and
     source['source_map'] == {'bytes': len(map_raw), 'sha256': map_digest,
                             'status': 'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION'},
     'authenticated inherited complete source-map premise')
driver = base/'compute/full23_response.py'
driverpin = 'efa6a6c2a9c7db596fe024246684f760b807e727fd2ddac8e8126182ba941a73'
need(hashlib.sha256(driver.read_bytes()).hexdigest() == driverpin, 'coefficient driver hash')
spec = importlib.util.spec_from_file_location('independent_point_coefficient_driver', driver)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
engine = mod.Response(root)  # verifies all local dependencies and seven scientific roots
labels = list(engine.labels)
weights = list(map(modal_weight, labels))
need(len(labels) == 23 and labels == profiles['labels'] == wave['labels'] == old['labels']
     and weights == list(engine.weights) == list(map(Q, profiles['gram'])),
     'independently reconstructed positive Gram convention')

# Independently index the source component arrays by field and label order.
blocks = {'yy': ('K:0', 'K:2', 'K:4'),
          'y0': ('L:0', 'L:2', 'L:4'),
          '00': ('N:0', 'N:2', 'N:4'),
          'scalar': ('w:0', 'w:2', 'w:4'),
          'ya': ('v:2', 'v:4'), '0a': ('u:2', 'u:4'),
          'ab': ('J:0', 'J:2', 'J:4', 'H:2', 'H:4', 'T:+', 'T:-')}
summary = source['global_cover']['summary']
certs = summary['positive_bundle_bulk_forcing_H_bound']
caps = summary['absolute_coefficient_caps']
fcaps = {}
comparisons = 0
for field, names in blocks.items():
    prefix = 'scalar_forcing' if field == 'scalar' else 'metric_forcing/' + field
    values = [Q(caps[prefix + '/' + str(j)]) for j in range(len(names))]
    need(min(values) >= 0, 'source magnitude positivity')
    record = certs['block_Gram_certificates'][field]
    gram = [list(map(Q, row)) for row in record['gram']]
    mult = 2 if field in ('y0', 'ya', '0a') else 1
    need(record['multiplicity'] == mult and len(gram) == len(names)
         and all(len(row) == len(names) for row in gram), 'source block dimensions')
    for i, name in enumerate(names):
        fcaps[name] = values[i]
        for j in range(len(names)):
            need(mult*gram[i][j] == (modal_weight(name) if i == j else 0),
                 'source/current Gram equality')
            comparisons += 1
    need(sum(modal_weight(name)*value*value for name, value in zip(names, values))
         == Q(record['squared_norm_upper']), 'source block norm independently rebuilt')
need(set(fcaps) == set(labels), 'all23 source coefficients mapped once')
norm2 = sum(modal_weight(name)*value*value for name, value in fcaps.items())
need(norm2 == Q(certs['squared_upper']) == Q(reported['source_H_norm_square']),
     'full source norm independently rebuilt')
need(Q(certs['F_H_upper']) == 17387513156573 and norm2 < Q(certs['F_H_upper'])**2,
     'uniform pointwise source H cap, not response SOURCE')

assembled = engine.assemble(engine.I(Q(9, 8)), engine.I(0))
actual = {name: [[(z.lo, z.hi) for z in row] for row in assembled[name]]
          for name in ('A', 'B', 'D')}
frozen = {name: [[(Q(z), Q(z)) for z in row] for row in profiles['matrices'][name]]
          for name in ('A', 'B', 'D')}
W = {jet: [iv(row['reconstructed_boundary_candidate'][jet]) for row in wave['full23_point_table']]
     for jet in ('value', 'dx', 'dr')}
Bs = {jet: [iv(row['B0'][jet]) for row in old['full23_point_table']]
      for jet in ('value', 'dx', 'dr')}

# Third identity: r0 = (La-L0)W + L0(W-Bs) + F.  The primary program uses
# (La-L0)Bs + La(W-Bs) + F and separately the direct LaW-L0Bs identity.
terms = []
f0 = [(Q(0), Q(0)) for _ in labels]
LaW = [(Q(0), Q(0)) for _ in labels]
for name, jet in (('A', 'dx'), ('B', 'dr'), ('D', 'value')):
    delta = [[difference(a, b) for a, b in zip(ar, br)]
             for ar, br in zip(actual[name], frozen[name])]
    terms.append(matvec(delta, W[jet]))
    terms.append(matvec(frozen[name], [difference(a, b) for a, b in zip(W[jet], Bs[jet])]))
    f0 = [plus(a, b) for a, b in zip(f0, matvec(frozen[name], Bs[jet]))]
    LaW = [plus(a, b) for a, b in zip(LaW, matvec(actual[name], W[jet]))]
rows = []
signs = []
for i, name in enumerate(labels):
    primary = reported['full23_point_table'][i]
    need(primary['label'] == name, 'reported row identity')
    need(iv(primary['frozen_free_source_f0']) == f0[i]
         and iv(primary['actual_lower_operator_LaW']) == LaW[i],
         'independent coefficient and point arithmetic')
    cap = fcaps[name]
    need(iv(primary['affine_bulk_source_uniform_interval']) == (-cap, cap),
         'independent same source component cap')
    own = plus(summed(term[i] for term in terms), (-cap, cap))
    other = iv(primary['complete_reduced_bulk_residual'])
    need(max(own[0], other[0]) <= min(own[1], other[1]), 'third identity overlaps primary enclosure')
    sign = 'positive' if own[0] > 0 else 'negative' if own[1] < 0 else 'undetermined'
    if sign != 'undetermined':
        need(primary['complete_reduced_bulk_residual_sign'] in (sign, 'undetermined'),
             'strict residual signs cannot disagree')
        signs.append(name)
    rows.append({'label': name, 'third_identity_residual': enc(own),
                 'independent_residual_sign': sign,
                 'source_component_abs_cap': str(cap)})
need(rows[12]['label'] == 'w:0' and rows[12]['independent_residual_sign'] == 'negative',
     'independently negative scalar-mean forcing defect')
lower, upper = iv(rows[12]['third_identity_residual'])
need(-Q(942*10**20) < lower <= upper < -Q(941*10**20),
     'independent outward scalar bounds')
engine.context.reauthenticate()
mod.need(mod.authenticate_local() == engine.local_before, 'coefficient dependencies unchanged')
need(hashlib.sha256(driver.read_bytes()).hexdigest() == driverpin, 'driver final readback')
need(hashlib.sha256(sourcepath.read_bytes()).hexdigest() == sourcepin, 'source final readback')
print(json.dumps({'status': 'PASS', 'point': ['9/8', '0'],
    'method': 'independent third residual decomposition and direct coefficient/source arithmetic',
    'source_Gram_entries_rebuilt': comparisons,
    'point_row_comparisons': 23,
    'independently_strict_nonzero_rows': signs,
    'full23_point_table': rows,
    'primary_report_sha256': reportpin,
    'scope': 'one conditional reduced-evolution forcing defect; same inherited coefficient assembler and geometric source premises; arithmetic independence only, not new source evaluation, pointwise Einstein residual, C2 error or physical observable'}, sort_keys=True))
