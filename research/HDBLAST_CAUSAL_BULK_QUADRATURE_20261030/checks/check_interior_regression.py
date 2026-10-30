#!/usr/bin/env python3
"""Independent inherited-point, support and interval-inclusion regressions."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
MODULE_PIN = '65879ec37d509f750d70c82a77c40d73ac9d34e0477caef399b71d99fa466e02'
ROBIN_PIN = 'bc7a0ef6545ddf09e54f16a71f4ac3dc0cf8f1a0a85d37d78a44bf1845ed1f69'
RESIDUAL_PIN = '30dc43f1ee1f1830976fefffc589beed42431d8013a1218094f9ecde0327ec69'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return raw


def overlap(a, b):
    return max(a[0], b[0]) <= min(a[1], b[1])


def contains(a, b):
    return a[0] <= b[0] <= b[1] <= a[1]


def encode(z):
    return list(map(str, z))


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--contact-capsule', type=Path, required=True)
ap.add_argument('--robin-capsule', type=Path, required=True)
ap.add_argument('--bulk-capsule', type=Path, required=True)
ap.add_argument('--primitive-panels', type=int, default=16)
args = ap.parse_args()
path = Path(__file__).with_name('interior_evaluator.py')
pin(path, MODULE_PIN)
spec = importlib.util.spec_from_file_location('regression_interior', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
e = module.InteriorEvaluator(args.repo_root, args.primitive_panels,
    args.contact_capsule, args.bulk_capsule)
robin_path = args.robin_capsule/'reference/RESUMMED.json'
robin = json.loads(pin(robin_path, ROBIN_PIN))
residual_path = args.bulk_capsule/'reference/POINT_RESIDUAL.json'
residual_digest = RESIDUAL_PIN
residual = json.loads(pin(residual_path, RESIDUAL_PIN))
need(residual['status'] == robin['status'] == 'PASS' and
     residual['point'] == robin['point'] == ['9/8', '0'] and
     [row['label'] for row in residual['full23_point_table']] == robin['labels'] == e.labels,
     'inherited registered successful point reports')
checks = {}
point = e.evaluate_xr(F(9, 8), 0)
for i, row in enumerate(robin['full23_point_table']):
    for key in ('value', 'dx', 'dr'):
        z = module.add(point['Bs'][key][i], point['delta_W0'][key][i])
        old = tuple(map(F, row['reconstructed_boundary_candidate'][key]))
        need(overlap(z, old), 'PR95 W0 point overlap:' + row['label'] + ':' + key)
checks['PR95_W0_first_jet_interval_overlaps'] = 69
for i, row in enumerate(residual['full23_point_table']):
    need(row['label'] == e.labels[i] and overlap(point['residual'][i],
        tuple(map(F, row['complete_reduced_bulk_residual']))),
        'PR97 complete residual overlap:' + row['label'])
checks['PR97_actual_reduced_residual_overlaps'] = 23

# Exact support identities exercise flat initial data, all normalized primitives
# and all first jets, separately from the magnitude-only F/G caps.
zeros = 0
for interval in ((F(0), F(1)), (F(-1), F(1)), (-e.collar, e.collar)):
    data = e.initial(interval)
    for family in ('U', 'V'):
        for row in data[family]:
            need(all(z == (0, 0) for z in row), 'exact initial collar')
            zeros += len(row)
    need(all(z == (0, 0) for z in data['p']), 'exact derivative collar')
    zeros += 23
for a, theta in (((0, F(1)), (-1, 1)), ((F(1), F(1)), (-1, 1))):
    result = e.evaluate(a, theta)
    for family in ('r0_lower', 'affine_F', 'residual'):
        need(all(z == (0, 0) for z in result[family]), 'inherited zero wedge')
        zeros += 23
    for family in ('Bs', 'delta_W0'):
        for row in result[family].values():
            need(all(z == (0, 0) for z in row), 'zero wedge first jets')
            zeros += 23
checks['exact_support_component_identities'] = zeros

positive = e.initial((F(17, 16), F(35, 32)))
negative = e.initial((-F(35, 32), -F(17, 16)))
parities = 0
for family in ('U', 'V'):
    for j, (p, n) in enumerate(zip(positive[family], negative[family])):
        for i in range(23):
            need(n[i] == module.scale(p[i], (-1)**j*e.parity[i]),
                 'normalized primitive parity')
            parities += 1
for i in range(23):
    need(negative['p'][i] == module.scale(positive['p'][i], -e.parity[i]),
         'initial derivative parity')
    parities += 1
checks['registered_reflection_parity_identities'] = parities

# These are inclusion diagnostics, not the proof of interval enclosure.  They
# cover the shell lift, both reflected branches and a q=0 branch join while the
# mathematical enclosure uses complete panels and analytic remainder bounds.
contained = 0
tests = [((F(9, 8), F(9, 8)), (F(-1), F(-7, 8)), F(-15, 16)),
         ((F(9, 8), F(9, 8)), (F(-1, 16), F(1, 16)), F(0)),
         ((F(9, 8), F(9, 8)), (F(7, 8), F(1)), F(15, 16)),
         ((F(9, 8), F(9, 8)), (-e.collar/F(9, 8)-F(1, 1024),
              -e.collar/F(9, 8)+F(1, 1024)), -e.collar/F(9, 8)),
         ((F(9, 8), F(9, 8)), (e.collar/F(9, 8)-F(1, 1024),
               e.collar/F(9, 8)+F(1, 1024)), e.collar/F(9, 8))]
for alpha, theta, midpoint in tests:
    broad = e.evaluate(alpha, theta)
    sample = e.evaluate(alpha[0], midpoint)
    for family in ('Bs', 'delta_W0'):
        for key in ('value', 'dx', 'dr'):
            for big, small in zip(broad[family][key], sample[family][key]):
                need(contains(big, small), 'whole-panel contains interior first-jet enclosure')
                contained += 1
    for big, small in zip(broad['r0_lower'], sample['r0_lower']):
        need(contains(big, small), 'whole-panel contains interior residual enclosure')
        contained += 1
checks['whole_panel_point_interval_inclusion_diagnostics'] = contained
e.authenticate()
pin(path, MODULE_PIN)
pin(robin_path, ROBIN_PIN)
need(hashlib.sha256(residual_path.read_bytes()).hexdigest() == residual_digest,
     'PR97 residual readback')
print(json.dumps({'status': 'PASS', 'checks': checks,
    'named_check_groups': len(checks),
    'primitive_panels_per_branch': args.primitive_panels,
    'scalar_mean_point_residual': encode(point['residual'][12]),
    'input_sha256': {'interior_evaluator': MODULE_PIN,
                     'PR95_RESUMMED': ROBIN_PIN,
                     'PR97_POINT_RESIDUAL': residual_digest},
    'scope': 'inherited point overlap, exact support/parity and interval-inclusion diagnostics for the authenticated conditional reduced residual evaluator; not independent validation of physical premises'}, sort_keys=True))
