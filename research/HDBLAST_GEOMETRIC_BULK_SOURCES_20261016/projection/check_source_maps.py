#!/usr/bin/env python3
"""Check exact collected-map joins, serialization, evaluation, and controls."""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'geometry'))
import angular_generic
import angular_project as projection
import source_coefficients
from formal_coefficients import P, formal_inputs
from check_complete_source_projection import load_parent, combine
from geometric_sources import evaluate_sources as geometric_sources

CONTROLS = {'gauge_gradient_sign': 'MAP_REDUCED_GAUGE_SIGN',
            'metric_forcing_factor_one': 'MAP_METRIC_FORCING_FACTOR',
            'tt_coefficients_deleted': 'MAP_TT_COMPLETE_JOIN',
            'qc_time_index_sign': 'MAP_QC_PARENT_ADAPTER'}


def require(condition, marker):
    if not condition: raise ValueError(marker)


def decode(terms):
    out = {}
    for row in terms:
        monomial = tuple(sorted(row['powers'].items()))
        require(monomial not in out, 'DUPLICATE_COEFFICIENT_MONOMIAL')
        out[monomial] = Q(row['coefficient'])
    return P(out)


def decode_maps(report):
    return {name: ({field: tuple(decode(row) for row in coefficients)
                    for field, coefficients in mapped.items()} if isinstance(mapped, dict)
                   else tuple(decode(row) for row in mapped))
            for name, mapped in report['maps'].items()}


def validate(path, mutation=None):
    raw = path.read_bytes(); report = json.loads(raw)
    require(report['status'] == 'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION', 'MAP_STATUS')
    mapped = decode_maps(report); m, s, b = formal_inputs(); parent = load_parent()
    first_symbols = set(m) | set(s)
    for mapped_value in mapped.values():
        coefficient_rows = mapped_value.values() if isinstance(mapped_value, dict) else [mapped_value]
        for coefficients in coefficient_rows:
            for coefficient in coefficients:
                require(all(sum(power for name, power in monomial if name in first_symbols) == 2
                            for monomial in coefficient.terms), 'MAP_QUADRATIC_FIRST_FIELD_HOMOGENEITY')
    if mutation == 'gauge_gradient_sign': mapped['QEred'] = combine(mapped['QE'], mapped['gauge_gradient'])
    elif mutation == 'metric_forcing_factor_one':
        mapped['metric_forcing'] = {name: tuple(value/2 for value in values)
                                    for name, values in mapped['metric_forcing'].items()}
    elif mutation == 'tt_coefficients_deleted':
        mapped['QEred']['ab'] = mapped['QEred']['ab'][:5] + (P(), P())
        expected = combine(mapped['QE'], mapped['gauge_gradient'], -1)
        require(mapped['QEred']['ab'] == expected['ab'], 'MAP_TT_COMPLETE_JOIN')
    elif mutation == 'qc_time_index_sign':
        mapped['QC']['0'] = tuple(-value for value in mapped['QC']['0'])
    require(mapped['QC'] == parent.gauge_source(m, b), 'MAP_QC_PARENT_ADAPTER')
    require(mapped['Qphi'] == parent.scalar_bulk_source(m, s, b), 'MAP_QPHI_PARENT_ADAPTER')
    require(mapped['QE_matter'] == parent.einstein_matter_source(m, s, b), 'MAP_MATTER_PARENT_ADAPTER')
    require(mapped['QE'] == combine(mapped['QRic'], mapped['QE_matter']), 'MAP_EINSTEIN_COMPLETE_JOIN')
    require(mapped['QEred'] == combine(mapped['QE'], mapped['gauge_gradient'], -1), 'MAP_REDUCED_GAUGE_SIGN')
    require(mapped['QEred_general'] == combine(mapped['QEred'], mapped['PC1']), 'MAP_GENERAL_PC1_JOIN')
    factor = b['ell']**2*b['rho']**2
    require(mapped['metric_forcing'] == {name: tuple(-2*factor*x for x in values)
                                        for name, values in mapped['QEred'].items()}, 'MAP_METRIC_FORCING_FACTOR')
    require(mapped['scalar_forcing'] == tuple(factor*x for x in mapped['Qphi']), 'MAP_SCALAR_FORCING_FACTOR')
    # Independent exact rational evaluation checks serialization and the public
    # evaluation path against complete indexed geometry and all-point projection.
    rm = {name: Q((i*7)%29-14, 11) for i, name in enumerate(sorted(m))}
    rs = {name: Q((i*3)%13-6, 7) for i, name in enumerate(sorted(s))}
    rb = {name: Q(i+1, 13) for i, name in enumerate(sorted(b))}
    rb.update(rho=Q(7,5), c=Q(5,4), theta=Q(3,5), ell=Q(1,10))
    direct = projection.project_geometric(geometric_sources(rm, rs, rb, angular_backend=angular_generic))
    evaluated = source_coefficients.evaluate_sources(rm, rs, rb, path=path)
    require(evaluated == direct, 'MAP_PUBLIC_EVALUATOR_EXACT_RATIONAL_GEOMETRY')
    forcing = source_coefficients.evaluate_bulk_forcing23(rm, rs, rb, path=path)
    require(forcing == projection.pack23({**direct['metric_forcing'], 'scalar': direct['scalar_forcing']}),
            'MAP_PUBLIC_FORCING23')
    integer_m = {name: i%5-2 for i, name in enumerate(sorted(m))}
    integer_s = {name: i%3-1 for i, name in enumerate(sorted(s))}
    integer_b = {name: 1 for name in b}; integer_b.update(rho=3, c=2, theta=0)
    integer_result = source_coefficients.evaluate_sources(integer_m, integer_s, integer_b, path=path)
    fraction_result = source_coefficients.evaluate_sources(
        {name: Q(value) for name, value in integer_m.items()},
        {name: Q(value) for name, value in integer_s.items()},
        {name: Q(value) for name, value in integer_b.items()}, path=path)
    require(integer_result == fraction_result, 'MAP_INTEGER_FRACTION_EXACT_EQUALITY')
    for mapped_value in integer_result.values():
        rows = mapped_value.values() if isinstance(mapped_value, dict) else [mapped_value]
        require(all(isinstance(value, (int,Q)) for row in rows for value in row),
                'MAP_INTEGER_INPUT_NO_FLOAT')
    return {'status': 'PASS_EXACT_SOURCE_MAP_JOINS_AND_EVALUATION',
            'maps_sha256': hashlib.sha256(raw).hexdigest(), 'maps_bytes': len(raw),
            'map_names': sorted(mapped),
            'independent_formal_jet_identities': True,
            'first_field_homogeneity_degree': 2,
            'direct_exact_rational_geometry_equal': True,
            'integer_fraction_public_api_exact_equal': True,
            'both_TT_maps_retained': [len(mapped['QEred']['ab'][i].terms) for i in (5,6)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maps', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--run-controls', action='store_true')
    parser.add_argument('--mutation', choices=tuple(CONTROLS))
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists(): raise ValueError('OUTPUT_MUST_BE_NEW')
    here = Path(__file__).resolve().parent
    if output == here or here in output.parents: raise ValueError('OUTPUT_MUST_BE_EXTERNAL')
    report = validate(args.maps.resolve(), args.mutation)
    controls = {}
    if args.run_controls:
        with tempfile.TemporaryDirectory(prefix='hdblast-source-map-controls-') as directory:
            for name, marker in CONTROLS.items():
                destination = Path(directory) / (name + '.json')
                command = [sys.executable, '-I', '-S']
                if sys.flags.optimize: command.append('-O')
                command += ['-B', str(Path(__file__).resolve()), '--maps', str(args.maps.resolve()),
                            '--output', str(destination), '--mutation', name]
                result = subprocess.run(command, text=True, capture_output=True, timeout=60)
                require(result.returncode != 0 and not destination.exists() and marker in result.stderr,
                        'CONTROL_NOT_REJECTED_' + name)
                controls[name] = {'rejected': True, 'required_stderr': marker, 'output_absent': True}
    report.update(controls=controls, control_count=len(controls))
    with output.open('x', encoding='utf-8') as handle:
        json.dump(report, handle, indent=2, sort_keys=True); handle.write('\n')
    print('PASS: exact source map joins, evaluator, and requested controls.')


if __name__ == '__main__': main()
