#!/usr/bin/env python3
"""Transfer a validated causal correction to the same conditional target U.

This changes point information, not the inherited whole-domain error theorem.
All arithmetic used to construct the intervals is exact rational arithmetic.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path

QUADRATURE_PIN = '5bd5ce2d40fdf0805d06351160fe7d2e323324983b69aeee5732c1e28013c6f8'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return json.loads(raw)


def iv(z):
    need(len(z) == 2, 'interval shape')
    a, b = map(Q, z)
    need(a <= b, 'ordered interval')
    return a, b


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def mul(a, b):
    v = [x*y for x in a for y in b]
    return min(v), max(v)


def wide(a, e):
    need(e >= 0, 'nonnegative error radius')
    return a[0]-e, a[1]+e


def meet(a, b):
    out = max(a[0], b[0]), min(a[1], b[1])
    need(out[0] <= out[1], 'independently valid boxes must overlap')
    return out


def sqrt_up(q):
    need(q >= 0, 'nonnegative norm square')
    g = 1 << 64
    n = math.isqrt(q.numerator*g*g//q.denominator)
    if Q(n, g)**2 < q:
        n += 1
    out = Q(n, g)
    need(out*out >= q, 'outward rational square root')
    return out


def mv(a, v):
    out = []
    for row in a:
        z = Q(0), Q(0)
        for c, w in zip(row, v):
            z = add(z, mul(c, w))
        out.append(z)
    return out


def enc(z):
    return list(map(str, z))


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', required=True, type=Path)
ap.add_argument('--quadrature-report', required=True, type=Path)
ap.add_argument('--quadrature-sha256')
ap.add_argument('--robin-capsule', type=Path)
ap.add_argument('--contact-capsule', type=Path)
ap.add_argument('--bulk-capsule', type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
robin = args.robin_capsule or root/'research/HDBLAST_ROBIN_RECONSTRUCTION_20261030'
contact = args.contact_capsule or root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030'
bulk = args.bulk_capsule or root/'research/HDBLAST_ACTUAL_BULK_UPDATE_20261030'
need(QUADRATURE_PIN is not None and
     args.quadrature_sha256 in (None, QUADRATURE_PIN), 'registered quadrature digest')
quad = pin(args.quadrature_report, QUADRATURE_PIN)
point = pin(robin/'reference/RESUMMED.json', 'bc7a0ef6545ddf09e54f16a71f4ac3dc0cf8f1a0a85d37d78a44bf1845ed1f69')
shell = pin(contact/'reference/AFFINE_SOURCE.json', 'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e')
bound = pin(bulk/'reference/BULK_BOUND.json', '321a9fc2917f1d6680dab825a087c03858ebf39de7e3fde22507df1dfc4b957f')
old = pin(bulk/'reference/LOCAL.json', '1832fba309ceaf0a138bbddc0aaf14602f81796bc1a7aa8cd938a74470e1b899')
profile = pin(root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/compute/response-normal/RESPONSE_PROFILES.json', '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
response = pin(root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/compute/response-normal/RESPONSE_REPORT.json', 'ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd')
need(all(d['status'] == 'PASS' for d in (quad, point, shell, bound, old)), 'successful reference reports')
need(quad['point'] == point['point'] == old['point'] == ['9/8', '0'], 'same selected shell point')
labels = profile['labels']
weights = list(map(Q, profile['gram']))
parity = list(map(int, profile['parity']['u']))
need(len(labels) == len(weights) == len(parity) == 23 and min(weights) > 0, 'full23 positive Gram')
need(labels == response['labels'] == bound['labels'], 'same registered space')
M = [[iv(z) for z in row] for row in response['actual_coefficient_boxes']['M']]
G = [iv(row['interval']) for row in shell['rows']]
need(len(M) == 23 and all(len(row) == 23 for row in M), 'actual Robin matrix dimensions')
E = Q(bound['one_update_readable_error_strict_upper'])
need(E == Q(3541*10**22) == Q(old['global_one_update_error_readable_strict_upper']) and
     0 < Q(bound['one_update_retained_target_error_X_cap']) < E, 'same inherited global error cap')
d = Q(1, 8)
qrows = quad['full23_point_table']
prows = point['full23_point_table']
orows = old['full23_point_table']
need(len(qrows) == len(prows) == len(orows) == 23, 'complete point tables')
correction = [{k: iv(row['correction'][k]) for k in ('value', 'dx', 'dr')} for row in qrows]
W0 = [{k: iv(row['reconstructed_boundary_candidate'][k]) for k in ('value', 'dx', 'dr')} for row in prows]
W1 = [{k: add(w[k], c[k]) for k in w} for w, c in zip(W0, correction)]
trace = [w['value'] for w in W1]
radial_law = [add(z, g) for z, g in zip(mv(M, trace), G)]
rows = []
strict = 0
signs = {k: [] for k in ('value', 'dx', 'dr')}
for i, label in enumerate(labels):
    need(qrows[i]['label'] == prows[i]['label'] == orows[i]['label'] == shell['rows'][i]['label'] == label, 'point row identity')
    need(prows[i]['parity'] == orows[i]['parity'] == parity[i], 'same output parity')
    factor = sqrt_up(1/weights[i])
    dual = sqrt_up(sum(max(abs(a), abs(b))**2/weights[j] for j, (a, b) in enumerate(M[i])))
    if parity[i] > 0:
        W1[i]['dr'] = meet(W1[i]['dr'], radial_law[i])
        target = {'value': wide(W1[i]['value'], d*E*factor),
                  'dx': wide(W1[i]['dx'], E*factor),
                  'dr': wide(W1[i]['dr'], d*E*dual)}
    else:
        need(W1[i]['value'] == W1[i]['dx'] == correction[i]['value'] == correction[i]['dx'] == (0, 0), 'Dirichlet trace identities')
        need(G[i] == (0, 0) and all(z == (0, 0) for z in M[i]), 'Dirichlet source and matrix rows')
        target = {'value': (Q(0), Q(0)), 'dx': (Q(0), Q(0)), 'dr': wide(W1[i]['dr'], E*factor)}
    previous = {k: iv(orows[i]['conditional_actual_target_U_best_point_boxes'][k]) for k in target}
    best = {k: meet(target[k], previous[k]) for k in target}
    ratios = {}
    for k in target:
        need(previous[k][0] <= best[k][0] <= best[k][1] <= previous[k][1], 'best point containment')
        oldwidth = previous[k][1]-previous[k][0]
        newwidth = best[k][1]-best[k][0]
        strict += newwidth < oldwidth
        ratios[k] = str(oldwidth/newwidth) if newwidth else None
        if correction[i][k][0] > 0 or correction[i][k][1] < 0:
            signs[k].append(label)
    if dual == 0 and parity[i] > 0:
        need(correction[i]['dr'] == (0, 0) and best['dr'] == previous['dr'], 'unchanged zero-M source-only normal box')
    rows.append({'label': label, 'parity': parity[i],
                 'correction': {k: enc(z) for k, z in correction[i].items()},
                 'enclosed_actual_bulk_candidate_W1': {k: enc(z) for k, z in W1[i].items()},
                 'new_target_U_point_enclosure': {k: enc(z) for k, z in target.items()},
                 'conditional_actual_target_U_best_point_boxes': {k: enc(z) for k, z in best.items()},
                 'previous_to_best_width_ratios': ratios})
need(labels[12] == 'w:0', 'registered scalar mean')
print(json.dumps({'status': 'PASS', 'point': ['9/8', '0'],
    'method': 'validated causal correction plus same global C1 error; intersect with previous best target boxes',
    'unchanged_global_error_readable_strict_upper': str(E),
    'point_box_containments': 69, 'strictly_narrower_best_target_boxes': strict,
    'sign_determined_correction_rows': signs,
    'scalar_mean_previous_to_best_width_ratios': rows[12]['previous_to_best_width_ratios'],
    'full23_point_table': rows,
    'input_sha256': {'quadrature': QUADRATURE_PIN},
    'scope': 'same conditional target U and integral-defined finite-iteration W1; numerical correction sign need not be target error sign; no improved global error theorem, quantitative C2 certificate, full Einstein/constraint residual, nonlinear family or physical observable'}, sort_keys=True))
