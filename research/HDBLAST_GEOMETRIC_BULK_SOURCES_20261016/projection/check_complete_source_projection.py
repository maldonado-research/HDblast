#!/usr/bin/env python3
"""Export universal exact ordinary-jet maps after all-point reconstruction.

All 42 metric two-jet entries and all 6 scalar two-jet entries are independent
Laurent-ring generators. No harmonic equations or sampled jet values are used
in this algebraic export. Actual-mode scalar/reduced formulas require the
actual C1=0 identity, as separately stated by the geometric engine.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import importlib.util
import json
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'geometry'))
import angular_generic
import angular_project as projection
from formal_coefficients import P, formal_inputs
from geometric_sources import evaluate_sources


def require(condition, marker):
    if not condition: raise ValueError(marker)


def load_parent():
    spec = importlib.util.spec_from_file_location('pinned_source_adapters', HERE / 'inputs' / 'source_adapters.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def combine(left, right, weight=1):
    return {name: tuple(x + weight*y for x, y in zip(left[name], right[name])) for name in left}


def encode_coefficients(mapped):
    if isinstance(mapped, dict):
        return {field: [P(value).encode() for value in values] for field, values in mapped.items()}
    return [P(value).encode() for value in mapped]


def build(mutation=None):
    m, s, b = formal_inputs()
    source = evaluate_sources(m, s, b, angular_backend=angular_generic)
    mapped = projection.project_geometric(source)
    parent = load_parent()
    require(mapped['QC'] == parent.gauge_source(m, b), 'FORMAL_QC_PARENT_ADAPTER')
    require(mapped['Qphi'] == parent.scalar_bulk_source(m, s, b), 'FORMAL_QPHI_PARENT_ADAPTER')
    require(mapped['QE_matter'] == parent.einstein_matter_source(m, s, b), 'FORMAL_MATTER_PARENT_ADAPTER')
    require(mapped['QE'] == combine(mapped['QRic'], mapped['QE_matter']), 'FORMAL_EINSTEIN_COMPLETE_JOIN')
    require(mapped['QEred'] == combine(mapped['QE'], mapped['gauge_gradient'], -1), 'FORMAL_REDUCED_GAUGE_SIGN')
    require(mapped['QEred_general'] == combine(mapped['QEred'], mapped['PC1']), 'FORMAL_GENERAL_PC1_JOIN')
    factor = b['ell']**2*b['rho']**2
    require(mapped['metric_forcing'] == {name: tuple(-2*factor*x for x in values)
                                        for name, values in mapped['QEred'].items()}, 'FORMAL_METRIC_FORCING_FACTOR')
    require(mapped['scalar_forcing'] == tuple(factor*x for x in mapped['Qphi']), 'FORMAL_SCALAR_FORCING_FACTOR')
    # Direct negative controls alter collected exact coefficients after all-point
    # projection and are checked against the retained complete join identities.
    if mutation == 'gauge_gradient_sign':
        broken = combine(mapped['QE'], mapped['gauge_gradient'])
        require(broken == mapped['QEred'], 'FORMAL_REDUCED_GAUGE_SIGN')
    elif mutation == 'metric_forcing_factor_one':
        broken = {name: tuple(-factor*x for x in values) for name, values in mapped['QEred'].items()}
        require(broken == mapped['metric_forcing'], 'FORMAL_METRIC_FORCING_FACTOR')
    elif mutation == 'tt_coefficients_deleted':
        broken = {**mapped['QEred'], 'ab': mapped['QEred']['ab'][:5] + (P(), P())}
        require(broken == mapped['QEred'], 'FORMAL_FULL_TT_RETENTION')
    return {'status': 'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION',
            'proof': 'independent Laurent jet symbols; exact coefficient witness; every polynomial remainder zero',
            'input_symbols': sorted([*m, *s, *b]),
            'input_units': {
                'derivatives': 'ordinary physical y,t derivatives; no factorial normalization',
                'frame': '0=y,1=orthonormal timelike,2..5=ambient S3 tangent; covariant slots',
                'rho': 'background warp factor', 'a': 'rho_y/rho', 'ay': 'ordinary y derivative of a',
                'c': 'cosh(t)', 'theta': 'tanh(t)', 'p': 'phi_y', 'py': 'phi_yy',
                'Vp': "V'(phi)", 'Vpp': "V''(phi)", 'Vppp': "V'''(phi)",
                'ell': 't=ell*x and dy=-rho*ell*dr'},
            'variation_convention': 'second parameter derivative d2/deps2 along affine g0+eps*q,phi0+eps*eta',
            'actual_mode_domain': 'QEred,Qphi require C1=0 as an identity, including its derivative; general alternatives retained',
            'gauge_convention': 'QC=second variation of full metric-lowered gauge covector; QC_background also retained',
            'certificate': projection.certificate(),
            'maps': {name: encode_coefficients(values) for name, values in mapped.items()},
            'exact_term_counts': {name: (sum(len(P(value).terms) for values in mapped_value.values() for value in values)
                                         if isinstance(mapped_value, dict) else sum(len(P(value).terms) for value in mapped_value))
                                  for name, mapped_value in mapped.items()},
            'parents': {'source_adapters.py': hashlib.sha256((HERE/'inputs'/'source_adapters.py').read_bytes()).hexdigest()},
            'engine_pins': {name: hashlib.sha256((HERE.parent/'geometry'/name).read_bytes()).hexdigest()
                            for name in ('geometric_sources.py', 'indexed_geometry.py')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mutation', choices=('gauge_gradient_sign', 'metric_forcing_factor_one', 'tt_coefficients_deleted'))
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists(): raise ValueError('OUTPUT_MUST_BE_NEW')
    if HERE == output or HERE in output.parents: raise ValueError('OUTPUT_MUST_BE_EXTERNAL')
    report = build(args.mutation)
    with output.open('x', encoding='utf-8') as handle:
        json.dump(report, handle, indent=2, sort_keys=True); handle.write('\n')
    print('PASS: complete formal geometric sources reconstructed at every angular point.')


if __name__ == '__main__': main()
