#!/usr/bin/env python3
"""Independent exact dependency-order audit of the complete 78-term TT row.

No frozen Python is imported. The exact JSON bytes are authenticated before
parsing. This audit establishes derivative requirements, not an error bound.
"""
import argparse
import hashlib
import json
from pathlib import Path

MAP_REL = Path('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json')
MAP_SHA256 = '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'
PROTECTED_LABELS = ('HDBLAST_COUPLED_MEAN_PILOT_20261025',
    'HDBLAST_INPUT_PRECISION_20261023', 'HDBLAST_COMMON_GAUGE_FORCING_20261016',
    'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016', 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
    'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023', 'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023',
    'HDBLAST_FULL_ANGULAR_RESPONSE_20261026', 'HDBLAST_SHELL_FEEDBACK_PRECISION_20261025')


def need(ok, message):
    if not ok:
        raise ValueError(message)


def guard_output(repo_root, output):
    need(not output.is_symlink(), 'OUTPUT_SYMLINK')
    need(not output.exists(), 'FRESH_OUTPUT')
    target = output.resolve()
    own_capsule = Path(__file__).resolve().parent.parent
    protected = [own_capsule]+[(Path(repo_root).resolve()/'research'/label).resolve() for label in PROTECTED_LABELS]
    need(all(target != root and root not in target.parents for root in protected), 'OUTPUT_INSIDE_PROTECTED_SOURCE_TREE')


def audit(repo_root):
    path = Path(repo_root) / MAP_REL
    need(not path.is_symlink(), 'MAP_SYMLINK')
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == MAP_SHA256, 'MAP_SHA256')
    data = json.loads(raw)
    need(data['status'] == 'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION', 'MAP_STATUS')
    rows = data['maps']['metric_forcing']['ab']
    terms = rows[5]
    need(len(rows) == 7 and rows[5] == rows[6] and len(terms) == 78, 'FULL_IDENTICAL_TT_ROWS')
    metric = set('ABCDELM')
    symbols = sorted({name for term in terms for name in term['powers']})
    jet_symbols = [name for name in symbols if name[0] in metric and name[1:] in ('', 'y', 't', 'yy', 'yt', 'tt')]
    background = sorted(set(symbols) - set(jet_symbols) - {'s'})
    second = sorted(name for name in jet_symbols if len(name) == 3)
    need(second == ['Mtt', 'Myt', 'Myy'], 'ONLY_M_SECOND_DERIVATIVES')
    need('s' in symbols and not any(name.startswith('s') and name != 's' for name in symbols), 'SCALAR_VALUE_ONLY')
    # Reconstruction differential orders are read from the registered formulas:
    # A..L depend on W and the first gauge derivatives; M=2d exactly.
    reconstruction_order = {name: (0 if name == 'M' else 1) for name in metric}
    requirements = {}
    for name in jet_symbols:
        order = reconstruction_order[name[0]] + len(name[1:]) + 1
        requirements[name] = {
            'physical_derivative_order_in_symbol': len(name[1:]),
            'reconstruction_master_gauge_order_upper': reconstruction_order[name[0]],
            'after_one_scaled_radial_derivative_order_upper': order,
        }
    requirements['s'] = {
        'physical_derivative_order_in_symbol': 0,
        'reconstruction_master_gauge_order_upper': 1,
        'after_one_scaled_radial_derivative_order_upper': 2,
    }
    maximum = max(row['after_one_scaled_radial_derivative_order_upper'] for row in requirements.values())
    need(maximum == 3, 'C3_SUFFICIENCY')
    return {
        'status': 'PASS_FULL_TT_RADIAL_SOURCE_REQUIRES_C3_ONLY',
        'map_relative_path': str(MAP_REL),
        'map_sha256': MAP_SHA256,
        'map_bytes': len(raw),
        'rows': 'metric_forcing/ab/5 = metric_forcing/ab/6',
        'all_terms_examined': len(terms),
        'symbols': symbols,
        'second_derivative_metric_symbols': second,
        'background_symbols': background,
        'symbol_derivative_requirements': requirements,
        'maximum_actual_master_gauge_derivative_order': maximum,
        'new_C4_error_certificate_required': False,
        'reason': 'The collected complete TT map contains second physical derivatives only of M=2d. Every A,B,C,D,E,L symbol has at most one physical derivative in the map. An additional scaled radial derivative therefore reaches total order three of actual master/gauge fields. This is an exact dependency audit of all 78 nonzero terms, not an assumed cancellation during interval evaluation.',
        'physical_conversion': 'partial_r=-ell*rho*partial_y; partial_y=-(ell*rho)^(-1)*partial_r. Differentiation also acts on rho and every background factor.',
        'background_requirement': 'Actual background ODE jets through total radial order three suffice for these selected symbols. The inherited ODE routine supplies order four.',
        'premises': [
            'Authenticated complete TT source coefficient map.',
            'Registered actual frame amplitudes A=-4u-2f_r/(ell*rho), B=f_x/(ell*rho)+b_r/ell, C=f/(rho*c)-c*d_r/ell, D=-2(u+a*f+b_x/ell), E=c*d_x/ell-b/c, L=2(u+a*f+H*b), M=2d; scalar s=v+p*f.',
            'Retained smooth actual master and common gauge solutions, and their authenticated C3 approximation errors on the backward-closed exterior restriction.',
        ],
        'scope': 'This verifies derivative sufficiency only; the actual radial source interval must be generated by a separate evaluator retaining all terms and all inherited nonzero errors.',
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    guard_output(args.repo_root,args.output)
    report = audit(args.repo_root)
    with args.output.open('x') as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(report['status'])


if __name__ == '__main__':
    main()
