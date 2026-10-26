#!/usr/bin/env python3
"""Bind the geometric coefficient implementation to the encoded map and bounds."""
import argparse
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'numeric'))
from output_safety import guard_output
from evaluate_encoded_map import Interval, need, sha
from evaluate_selected_magnetic import Polynomial, supplied_groups, encoded_polynomial, decimal_outward


def load(path):
    spec = importlib.util.spec_from_file_location('screen_coefficient_binding', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--map', type=Path, dest='map_path', required=True)
    parser.add_argument('--complete', type=Path, required=True)
    parser.add_argument('--inputs-sha256', required=True)
    parser.add_argument('--map-sha256', required=True)
    parser.add_argument('--complete-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    guard_output(args.output, source_file=__file__, input_paths=[args.inputs, args.map_path, args.complete])
    files = {'check_final_coefficient.py':Path(__file__),
             'theory/screen_coefficient.py':ROOT/'theory/screen_coefficient.py',
             'numeric/evaluate_selected_magnetic.py':ROOT/'numeric/evaluate_selected_magnetic.py',
             'numeric/evaluate_encoded_map.py':ROOT/'numeric/evaluate_encoded_map.py',
             'numeric/output_safety.py':ROOT/'numeric/output_safety.py'}
    before = {name:sha(path) for name,path in files.items()}
    for path,pin in ((args.inputs,args.inputs_sha256),(args.map_path,args.map_sha256),(args.complete,args.complete_sha256)):
        need(sha(path) == pin, 'EXTERNAL_INPUT_PIN:'+path.name)
    inputs, mapping, complete = [json.loads(p.read_text()) for p in (args.inputs,args.map_path,args.complete)]
    need(inputs['apex'] == complete['apex'] == {'x':'1/4','r':'5/4','t':'1/4000000'}, 'REGISTERED_APEX')
    theory = load(files['theory/screen_coefficient.py'])
    names = {name for terms in mapping['components'].values() for term in terms for name in term['powers']}
    symbolic = {name:Polynomial.symbol(name) for name in names}
    derived = theory.coefficients(symbolic, symbolic, symbolic['u_x'])
    rename = {'quadratic':'quadratic','linear':'linear','Tplus_total':'plus','Tminus_total':'minus'}
    for name,target in rename.items():
        need(derived[target].terms == encoded_polynomial(mapping['components'][name]).terms,
             'THEORY_IMPLEMENTATION_EQUALS_ENCODED_MAP:'+name)
    env = {name:Interval.cv(row) for name,row in inputs['symbols'].items()}
    actual = theory.coefficients(env,env,env['u_x'])
    grouped = supplied_groups(env)
    for name,target in rename.items():
        need(actual[target].row() == grouped[name].row(), 'THEORY_INTERVAL_GROUPING:'+name)
    certified = {name:Interval(*row['exact']) for name,row in complete['certified_intersections'].items()}
    for name in rename:
        need(grouped[name].lo <= certified[name].lo <= certified[name].hi <= grouped[name].hi,
             'CERTIFIED_INTERVAL_IS_WITHIN_DERIVED_ROUTE:'+name)
    linear, quadratic, total = certified['linear'],certified['quadratic'],certified['Tplus_total']
    need(linear.hi < 0 and total.hi < 0, 'STRICT_NEGATIVE_LINEAR_AND_TOTAL')
    qcap = max(abs(quadratic.lo),abs(quadratic.hi))
    ratio = qcap/(-linear.hi)
    need(qcap < F(2,10**18) and ratio < F(5,10**16), 'QUADRATIC_CANNOT_CANCEL_LINEAR')
    norm = Interval(total.hi**2/3,total.lo**2/3)
    need(norm.lo > 0, 'POSITIVE_MEAN_SQUARED_TT_DIAGNOSTIC')
    need({name:sha(path) for name,path in files.items()} == before, 'IMPLEMENTATIONS_CHANGED')
    result = {'schema':'HDBLAST_REGISTERED_SCREEN_MAGNETIC_FINAL_BINDING_V1',
      'status':'CONDITIONAL_REGISTERED_SCREEN_MAGNETIC_COEFFICIENT_EXCLUDES_ZERO',
      'apex':inputs['apex'], 'input_sha256':args.inputs_sha256,
      'map_sha256':args.map_sha256, 'complete_report_sha256':args.complete_sha256,
      'source_sha256':before, 'exact_theory_to_encoded_Laurent_identity':True,
      'quadratic_terms':len(mapping['components']['quadratic']),
      'Tplus_total':complete['certified_intersections']['Tplus_total'],
      'Tminus_total':complete['certified_intersections']['Tminus_total'],
      'quadratic':complete['certified_intersections']['quadratic'],
      'absolute_quadratic_upper_bound_exact':str(qcap),
      'absolute_quadratic_strictly_less_than':'2/1000000000000000000',
      'quadratic_to_linear_ratio_upper_bound_exact':str(ratio),
      'quadratic_to_linear_ratio_strictly_less_than':'5/10000000000000000',
      'mean_squared_TT_projection':{'exact':norm.row(),'outward_decimal':decimal_outward(norm,12),
                                  'formula':'mean |Pi_TT H2|^2=h^2/3; mean1=1'},
      'angular_pointwise_norm':'|h(Tplus-Tminus)|^2=2*h^2*S*(1-S); S=n0^2+n1^2',
      'angular_zero_set':'S=0 or S=1; two great circles, so not everywhere nonzero',
      'limits':['Conditional formal registered response; U is the unit future normal of t slices and n is their positive spatial y normal.',
                'The angular screen diagnostic is a projection of five-dimensional Weyl, not its entire observer magnetic tensor.',
                'The observer/clock/point prescription has not been geometrically anchored; no second-order gauge-invariance claim.',
                'Full compatible second-order Einstein–matter boundary solution and finite-amplitude nonlinear family remain unfinished.',
                'No gravitational radiation, detector-force, energy or cosmological-origin conclusion.']}
    with args.output.open('x') as handle:
        json.dump(result,handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps({'status':result['status'],'quadratic_terms':14,
                      'Tplus_outward_decimal':result['Tplus_total']['outward_decimal']}))


if __name__ == '__main__':
    main()
