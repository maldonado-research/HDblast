#!/usr/bin/env python3
"""Independently audit complete-cover receipts, caps and Hilbert conversion.

The source map and interval extension have a separate independent review.
This checker audits every stored cell and summary, computes sphere Gram
matrices by an independent unreduced monomial integral, and optionally
replays selected cell records against the reviewed numeric producer.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import factorial, isqrt
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
CHECKPOINT = ROOT.parent


def need(ok, text):
    if not ok:
        raise ValueError(text)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def interval(row):
    need(isinstance(row, list) and len(row) == 2, 'interval pair')
    lo, hi = map(F, row)
    need(lo <= hi, 'stored interval ordered')
    return lo, hi


def leaves(tree, prefix=()):
    if isinstance(tree, dict):
        for key, value in tree.items():
            yield from leaves(value, prefix+(str(key),))
    elif isinstance(tree, list) and len(tree) == 2 and all(isinstance(v, str) for v in tree):
        yield '/'.join(prefix), interval(tree)
    elif isinstance(tree, list):
        for index, value in enumerate(tree):
            yield from leaves(value, prefix+(str(index),))
    else:
        raise ValueError('unexpected source leaf')


def summarize(rows):
    caps, signs, schema = {}, {}, None
    for row in rows:
        values = dict(leaves(row['sources']))
        if schema is None:
            schema = set(values)
        need(set(values) == schema, 'all source rows retained in every cell')
        for key, (lo, hi) in values.items():
            caps[key] = max(caps.get(key, F(0)), abs(lo), abs(hi))
            if lo > 0 or hi < 0:
                signs[key] = signs.get(key, 0)+1
    return caps, signs


def moment(exponent):
    # Dirichlet(1/2,1/2,1/2,1/2) sphere moment, via half-integer
    # Gamma ratios.  No producer polynomial multiplication/reduction.
    if any(n % 2 for n in exponent):
        return F(0)
    half = [n//2 for n in exponent]
    out = F(1, factorial(sum(half)+1))
    for n in half:
        out *= F(factorial(2*n), 4**n*factorial(n))
    return out


def inner(a, b):
    return sum((v*w*moment(tuple(x+y for x, y in zip(e, f)))
                for component in a for e, v in a[component].items() for f, w in b[component].items()), F(0))


def audit_norm(summary):
    angular = load_module(CHECKPOINT/'projection/angular_generic.py', 'independent_basis_for_norm')
    caps = {key: F(value) for key, value in summary['absolute_coefficient_caps'].items()}
    observed = summary['positive_bundle_bulk_forcing_H_bound']
    total = F(0)
    weights = (('yy', angular.SCALARS, 1), ('y0', angular.SCALARS, 2), ('00', angular.SCALARS, 1),
               ('ya', angular.GRADS, 2), ('0a', angular.GRADS, 2), ('ab', angular.TENSORS, 1),
               ('scalar', angular.SCALARS, 1))
    for name, basis, multiplicity in weights:
        prefix = 'scalar_forcing' if name == 'scalar' else 'metric_forcing/'+name
        cc = [caps[prefix+'/'+str(i)] for i in range(len(basis))]
        gram = [[inner(a, b) for b in basis] for a in basis]
        need(all(gram[i][i] > 0 for i in range(len(basis))), 'positive angular Gram diagonal')
        need(all(not gram[i][j] for i in range(len(basis)) for j in range(len(basis)) if i != j), 'orthogonal registered angular basis')
        bound = multiplicity*sum((abs(gram[i][j])*cc[i]*cc[j] for i in range(len(basis)) for j in range(len(basis))), F(0))
        block = observed['block_Gram_certificates'][name]
        need(gram == [[F(v) for v in row] for row in block['gram']], 'independent unreduced sphere Gram '+name)
        need(multiplicity == block['multiplicity'], 'positive Frobenius multiplicity '+name)
        need(bound == F(block['squared_norm_upper']), 'independent positive block norm '+name)
        total += bound
    need(total == F(observed['squared_upper']), 'independent total positive forcing norm')
    upper = isqrt(total.numerator//total.denominator)+1
    need(F(upper) == F(observed['F_H_upper']) and upper*upper >= total, 'independent rational square-root roundup')
    need(F(observed['retained_parent_F_H']) == F(171610068644000), 'same retained parent forcing comparison')
    need(observed['strictly_improves_parent_bulk_forcing_bound'] == (upper < 171610068644000), 'parent forcing improvement boolean')
    return str(upper), str(total)


def shell_adapters(sources):
    base = dict(leaves(sources['one_copy_outward_hessian']))
    for name, factor in (('two_copy_outward_hessian', 2), ('two_copy_inward_hessian', -2)):
        other = dict(leaves(sources[name]))
        need(set(other) == set(base), 'all physical shell families retained')
        for key, ends in base.items():
            need(other[key] == tuple(sorted(factor*v for v in ends)), 'exact shell orientation adapter '+name+'/'+key)
    forcing = dict(leaves(sources['coordinate_physical_forcing']))
    for key, ends in base.items():
        factor = -1 if key.startswith('scalar/') else -2
        need(forcing[key] == tuple(sorted(factor*v for v in ends)), 'exact coordinate shell forcing adapter '+key)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--optimized', type=Path)
    parser.add_argument('--maps', type=Path, default=CHECKPOINT/'projection/GEOMETRIC_COEFFICIENT_MAPS.json')
    parser.add_argument('--replay-selected', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse overwrite')
    raw = args.report.read_bytes()
    report = json.loads(raw)
    need(report['status'] == 'PASS_CONTINUOUS_ACTUAL_EVOLVED_GEOMETRIC_SOURCE_CELL_ENCLOSURES', 'actual source status')
    if args.optimized:
        need(raw == args.optimized.read_bytes(), 'normal optimized complete source report equality')
    maps = json.loads(args.maps.read_text())
    need(report['source_map']['sha256'] == sha(args.maps) and report['source_map']['bytes'] == args.maps.stat().st_size, 'source report exact final map identity')
    freeze_path = CHECKPOINT/'numeric/SOURCE_CODE_FREEZE.json'
    freeze = json.loads(freeze_path.read_text())
    for name, pin in freeze['files'].items():
        path = CHECKPOINT/'numeric'/name
        need(path.stat().st_size == pin['bytes'] and sha(path) == pin['sha256'], 'final numeric source-code freeze '+name)
    for name, digest in maps['engine_pins'].items():
        need(sha(CHECKPOINT/'geometry'/name) == digest, 'final source engine pin '+name)
    interval_review_path = CHECKPOINT/'review/SOURCE_INTERVALS_FINAL_NORMAL.json'
    interval_review = json.loads(interval_review_path.read_text())
    need(interval_review['status'] == 'PASS_INDEPENDENT_EXACT_SOURCE_INTERVAL_EVALUATION' and interval_review['maps_sha256'] == sha(args.maps), 'final independent source interval review')
    for name, digest in interval_review['input_sha256'].items():
        need(sha(CHECKPOINT/name) == digest, 'final independent interval review dependency '+name)
    source_schema = set()
    for name in ('QE', 'QEred', 'Qphi', 'QC', 'metric_forcing', 'scalar_forcing'):
        mapped = maps['maps'][name]
        if isinstance(mapped, dict):
            source_schema.update(name+'/'+field+'/'+str(i) for field, coefficients in mapped.items() for i in range(len(coefficients)))
        else:
            source_schema.update(name+'/'+str(i) for i in range(len(mapped)))
    cover = report['global_cover']
    d = cover['division_per_unit']
    need(isinstance(d, int) and d >= 1, 'nonempty rational cover scale')
    expected = [[i, k] for i in range(3*d) for k in range(8*d) if i+k < 8*d]
    rows = cover['cell_enclosures']
    need(len(rows) == cover['cells'] == len(expected), 'complete actual D cover cardinality')
    need([row['index'] for row in rows] == expected, 'all actual D cells including clipped diagonal')
    for row in rows:
        i, k = row['index']
        need(interval(row['x']) == (F(i, d), F(i+1, d)) and interval(row['r']) == (F(k, d), F(k+1, d)), 'exact actual D grid cell')
        need(row['intersection'] == 'x+r<=8', 'diagonal actual domain restriction')
    for key in ('global_cover', 'evolved_radial_profile'):
        part = report[key]
        caps, counts = summarize(part['cell_enclosures'])
        need(set(caps) == source_schema, 'every registered mapped source coefficient retained '+key)
        need(caps == {k: F(v) for k, v in part['summary']['absolute_coefficient_caps'].items()}, 'all stored interval cap maxima '+key)
        need(counts == part['summary']['strict_nonzero_cell_counts'], 'all stored strict-sign counts '+key)
    profile = report['evolved_radial_profile']
    n = profile['cells']
    need(n >= 8 and len(profile['cell_enclosures']) == n, 'nonempty complete evolved profile')
    for i, row in enumerate(profile['cell_enclosures']):
        need(row['index'] == i and interval(row['r']) == (1+F(i, n), 1+F(i+1, n)), 'gap-free closed evolved profile')
        need(row['x'] == profile['x'], 'whole time cell of evolved profile')
    shell = report['shell_whole_time_cover']
    need(shell['cells'] == 3*d == len(shell['cell_enclosures']), 'complete physical shell time cover')
    for i, row in enumerate(shell['cell_enclosures']):
        need(interval(row['x']) == (F(i, d), F(i+1, d)), 'gap-free shell time partition')
        shell_adapters(row['sources'])
    shell_adapters(report['shell_witness']['sources'])
    need(F(report['shell_witness']['sources']['one_copy_outward_hessian']['scalar'][0][0]) > 4*10**8, 'strict whole-cell positive physical scalar shell Hessian')
    for path, sign in report['bulk_witness']['strict_signs'].items():
        value = report['bulk_witness']['sources']
        for part in path.split('/'):
            value = value[int(part)] if isinstance(value, list) else value[part]
        lo, hi = interval(value)
        need((hi < 0) if sign == 'negative' else (lo > 0), 'reported whole-cell source coefficient sign '+path)
    need(set(dict(leaves(report['bulk_witness']['sources']))) == source_schema and set(dict(leaves(report['reflected_witness']['sources']))) == source_schema, 'complete two named witness source families')
    norm, square = audit_norm(cover['summary'])
    replays = []
    if args.replay_selected:
        sys.path.insert(0, str(CHECKPOINT/'numeric'))
        sys.path.insert(0, str(CHECKPOINT/'projection'))
        producer = load_module(CHECKPOINT/'numeric/evaluate_actual_sources.py', 'selected_source_replay')
        for row in (rows[0], rows[len(rows)//2], rows[-1]):
            values = producer.cell(producer.I(*interval(row['x'])), producer.I(*interval(row['r'])), args.maps, clipped=True)
            need(producer.encode(values) == row['sources'], 'selected actual complete-cell source receipt replay '+str(row['index']))
            replays.append(row['index'])
    paths = [CHECKPOINT/'numeric/evaluate_actual_sources.py', CHECKPOINT/'numeric/shell_source_values.py', CHECKPOINT/'numeric/fast_dyadic.py', CHECKPOINT/'projection/source_coefficients.py', CHECKPOINT/'projection/angular_generic.py', Path(__file__)]
    out = {'status': 'PASS_INDEPENDENT_COMPLETE_ACTUAL_SOURCE_REPORT',
           'producer_report_sha256': hashlib.sha256(raw).hexdigest(), 'source_map_sha256': sha(args.maps),
           'source_code_freeze_sha256': sha(freeze_path), 'independent_interval_review_sha256': sha(interval_review_path),
           'normal_optimized_bytes_equal': bool(args.optimized),
           'whole_causal_cells': len(rows), 'whole_evolved_profile_cells': n, 'whole_shell_time_cells': 3*d,
           'every_cell_cap_and_sign_count_independently_recomputed': True,
           'selected_complete_cell_numeric_replays': replays, 'independent_positive_forcing_H_upper': norm,
           'independent_positive_forcing_squared_upper': square,
           'input_sha256': {str(path.relative_to(CHECKPOINT)): sha(path) for path in paths},
           'scope': 'Whole-cell interval report audit and independent positive sphere norm conversion. Actual jet composition and universal geometry/map/interval-extension proofs are separately reviewed. No source quadrature, second response, nonlinear family, or invariant force.'}
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
    print(out['status'])


if __name__ == '__main__':
    main()
