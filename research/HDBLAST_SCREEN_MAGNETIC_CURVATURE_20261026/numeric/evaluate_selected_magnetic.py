#!/usr/bin/env python3
"""Grouped and encoded evaluation of the explicitly supplied selected map.

No actual input oracle calls. The geometric formula is the parent-provided
consensus expression; this module verifies its exact transcription and gives
independent interval groupings plus exact synthetic arithmetic controls.
"""
import argparse
import hashlib
import json
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from evaluate_encoded_map import Interval, evaluate_map, need, sha
from output_safety import guard_output, source_snapshot, verify_source_snapshot, source_labels
FIELDS = ('A','B','C','D','E','L','M')


class Polynomial:
    """Exact finite Laurent polynomials for transcription identity checks."""
    def __init__(self, terms=None):
        self.terms = {key:F(val) for key,val in (terms or {}).items() if val}

    @staticmethod
    def cv(value):
        if isinstance(value,Polynomial):
            return value
        return Polynomial({():F(value)})

    @staticmethod
    def symbol(name):
        return Polynomial({((name,1),):F(1)})

    def __add__(self,other):
        out = dict(self.terms)
        for key,val in self.cv(other).terms.items():
            out[key] = out.get(key,F(0))+val
        return Polynomial(out)
    __radd__ = __add__

    def __neg__(self):
        return Polynomial({key:-val for key,val in self.terms.items()})

    def __sub__(self,other):
        return self+-self.cv(other)

    def __rsub__(self,other):
        return self.cv(other)+-self

    def __mul__(self,other):
        out = {}
        for a,av in self.terms.items():
            for b,bv in self.cv(other).terms.items():
                powers = dict(a)
                for name,power in b:
                    powers[name] = powers.get(name,0)+power
                key = tuple(sorted((name,power) for name,power in powers.items() if power))
                out[key] = out.get(key,F(0))+av*bv
        return Polynomial(out)
    __rmul__ = __mul__

    def __pow__(self,power):
        need(type(power) is int,'POLYNOMIAL_INTEGER_POWER')
        if power < 0:
            need(len(self.terms) == 1,'LAURENT_MONOMIAL_DIVISOR')
            key,val = next(iter(self.terms.items()))
            return Polynomial({tuple((name,power*exponent) for name,exponent in key):val**power})
        out = Polynomial.cv(1)
        for _ in range(power):
            out *= self
        return out

    def __truediv__(self,other):
        return self*self.cv(other)**(-1)

    def __rtruediv__(self,other):
        return self.cv(other)*self**(-1)


def supplied_groups(v):
    """The exact formula supplied by two theory reviewers via the parent."""
    rho,c,theta,ell = (v[key] for key in ('rho','c','theta','ell'))
    B,C,D,E,L,M,Ct,Ey,Lt,My,Mt,ux = (v[key] for key in
        ('B','C','D','E','L','M','Ct','Ey','Lt','My','Mt','u_x'))
    R = rho*c
    g1 = ((D+4*M-L)*Mt-2*M*Lt)/(4*rho)
    g2 = -E*(D+L)/(2*R)
    g3 = (C/R-My/2)*(R*(Ct-theta*C)/(2*rho)-R*Ey/2+B/2)
    force = g1+g2+g3
    quadratic = 4*force/R
    linear = 3*ux/(2*ell*rho**2*c)
    return {'F_group1':g1,'F_group2':g2,'F_group3':g3,'F':force,
            'linear':linear,'quadratic':quadratic,
            'Tplus_total':linear+quadratic,'Tminus_total':-(linear+quadratic)}


def encoded_polynomial(terms):
    out = Polynomial()
    for term in terms:
        powers = tuple(sorted((name,power) for name,power in term['powers'].items() if power))
        out += Polynomial({powers:F(term['coefficient'])})
    return out


def verify_transcription(encoded):
    names = {'rho','c','theta','ell','u_x','B','C','D','E','L','M','Ct','Ey','Lt','My','Mt'}
    symbolic = supplied_groups({name:Polynomial.symbol(name) for name in names})
    for name,terms in encoded['components'].items():
        need(name in symbolic,'UNEXPECTED_MAP_COMPONENT:'+name)
        need(encoded_polynomial(terms).terms == symbolic[name].terms,'SYMBOLIC_TRANSCRIPTION:'+name)
    need(set(encoded['components']) == {'linear','quadratic','Tplus_total','Tminus_total'},'FOUR_REQUIRED_COMPONENTS')
    need(len(symbolic['quadratic'].terms) == 14,'FOURTEEN_QUADRATIC_TERMS')
    for monomial in symbolic['quadratic'].terms:
        degree = sum(power for name,power in monomial if name[0] in FIELDS)
        need(degree == 2,'HOMOGENEOUS_QUADRATIC_FIRST_FIELD')
    return {'all_four_map_components_equal_supplied_formula':True,
            'quadratic_laurent_terms':14,'homogeneous_first_field_degree':2,
            'arithmetic':'Exact Fraction Laurent-polynomial coefficient comparison.',
            'scope':'Transcription equality to supplied formula; geometric proof is separate.'}


def interval_intersection(a,b,label):
    lo,hi = max(a.lo,b.lo),min(a.hi,b.hi)
    need(lo <= hi,'DISJOINT_EQUIVALENT_INTERVAL_ROUTES:'+label)
    return Interval(lo,hi)


def decimal_outward(value,places):
    scale = 10**places
    def render(k):
        sign = '-' if k < 0 else ''
        whole,frac = divmod(abs(k),scale)
        return f'{sign}{whole}.{frac:0{places}d}'
    return [render((value.lo*scale).__floor__()),render((value.hi*scale).__ceil__())]


def controls(encoded):
    names = {symbol for terms in encoded['components'].values() for term in terms for symbol in term['powers']}
    blank = {name:Interval(0) for name in names}
    blank.update({'rho':Interval(2),'c':Interval(3),'ell':Interval('1/1000000')})
    results = []
    fixtures = [
      ('zero_fields',{},F(0)),
      ('E_D_connection_term',{'D':5,'E':7},-F(35,18)),
      ('C_B_product',{'C':7,'B':5},F(35,18)),
      ('C_Ey_derivative',{'C':7,'Ey':5},-F(35,3)),
      ('M_Lt_shift',{'M':5,'Lt':7},-F(35,6)),
      ('Mt_D_L_difference',{'D':5,'L':7,'Mt':11},-F(11,6)),
      ('My_Ct_cross',{'My':5,'Ct':7},-F(35,2)),
    ]
    for label,values,expected in fixtures:
        env = {**blank,**{key:Interval(val) for key,val in values.items()}}
        direct = supplied_groups(env)['quadratic']
        expanded = evaluate_map(encoded,env)['quadratic']
        need(direct.lo == direct.hi == expected and direct.row() == expanded.row(),'SINGLETON_FIXTURE:'+label)
        results.append({'name':label,'passed':True,'exact':[str(expected),str(expected)]})
    env = {**blank,'u_x':Interval(-1)}
    actual = evaluate_map(encoded,env)
    need(actual['quadratic'].row() == ['0','0'] and actual['linear'].lo < 0 and
         actual['Tplus_total'].row() == actual['linear'].row(),'LINEAR_ONLY_FIXTURE')
    results.append({'name':'zero_first_field_retains_linear_response','passed':True})
    sample = {**blank,**{key:Interval(F(index+1,7)) for index,key in enumerate(sorted(names-{'rho','c','ell','u_x'}))}}
    doubled = dict(sample)
    for key in sample:
        if key[0] in FIELDS:
            doubled[key] = 2*sample[key]
    q1,q2 = supplied_groups(sample)['quadratic'],supplied_groups(doubled)['quadratic']
    need(q2.row() == (4*q1).row(),'QUADRATIC_SCALING_FIXTURE')
    results.append({'name':'first_field_scaling_degree_two','passed':True})
    # The omitted 4*M*Mt term is a deliberately rejected transcription mutant.
    mutant = json.loads(json.dumps(encoded))
    mutant['components']['quadratic'] = [term for term in mutant['components']['quadratic']
      if not (term['coefficient'] == '4' and term['powers'].get('M') == 1 and term['powers'].get('Mt') == 1)]
    try:
        verify_transcription(mutant)
    except ValueError as error:
        need(str(error) == 'SYMBOLIC_TRANSCRIPTION:quadratic','WRONG_MUTANT_REJECTION')
    else:
        raise ValueError('OMITTED_TERM_MUTANT_ACCEPTED')
    results.append({'name':'omitted_M_Mt_term_mutant','passed':True,'expected_rejection':'SYMBOLIC_TRANSCRIPTION:quadratic'})
    return {'all_passed':True,'controls':results,'actual_solution_calls':0,
            'scope':'Exact formula-transcription and arithmetic controls; independent geometric controls are separate.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--map',type=Path,dest='map_path',required=True)
    parser.add_argument('--require-input-sha256',required=True)
    parser.add_argument('--require-map-sha256',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    guard_output(args.output,source_file=__file__,input_paths=[args.inputs,args.map_path])
    sources_before = source_snapshot([__file__,HERE/'evaluate_encoded_map.py',HERE/'output_safety.py'])
    own_pin = sources_before[str(Path(__file__).resolve())]
    need(sha(args.inputs) == args.require_input_sha256,'INPUT_EXTERNAL_PIN')
    need(sha(args.map_path) == args.require_map_sha256,'MAP_EXTERNAL_PIN')
    inputs,encoded = json.loads(args.inputs.read_text()),json.loads(args.map_path.read_text())
    transcription = verify_transcription(encoded)
    control_results = controls(encoded)
    env = {key:Interval.cv(value) for key,value in inputs['symbols'].items()}
    grouped,expanded = supplied_groups(env),evaluate_map(encoded,env)
    intersections = {name:interval_intersection(grouped[name],value,name) for name,value in expanded.items()}
    need(intersections['Tminus_total'].row() == (-intersections['Tplus_total']).row(),'OPPOSITE_TENSOR_COEFFICIENTS')
    need(sha(args.inputs) == args.require_input_sha256 and sha(args.map_path) == args.require_map_sha256 and sha(__file__) == own_pin,
         'INPUT_MAP_OR_EVALUATOR_CHANGED')
    verify_source_snapshot(sources_before)
    report = {'schema':'HDBLAST_ACTUAL_SELECTED_COMPLETE_MAGNETIC_INTERVAL_V1',
      'status':'CONDITIONAL_COMPLETE_SELECTED_MAGNETIC_COEFFICIENT_EXCLUDES_ZERO' if intersections['Tplus_total'].sign() != 'inconclusive'
               else 'CONDITIONAL_COMPLETE_SELECTED_MAGNETIC_COEFFICIENT_INCONCLUSIVE',
      'apex':inputs['apex'],'input_sha256':args.require_input_sha256,'map_sha256':args.require_map_sha256,
      'evaluator_sha256':own_pin,'generic_evaluator_sha256':sha(HERE/'evaluate_encoded_map.py'),
      'source_integrity_before_and_after':source_labels(sources_before),
      'actual_input_oracle_calls_total':inputs['actual_apex_calls'],'additional_actual_solution_calls':0,
      'campaign_runs':0,'supplied_formula':encoded['supplied_grouped_formula'],
      'transcription_audit':transcription,'controls':control_results,
      'grouped_intermediates':{name:val.row() for name,val in grouped.items()},
      'expanded_intervals':{name:val.row() for name,val in expanded.items()},
      'certified_intersections':{name:{'exact':val.row(),'sign':val.sign(),
            'outward_decimal':decimal_outward(val,24 if name == 'quadratic' else 12),
            'display_only':[float(val.lo),float(val.hi)]} for name,val in intersections.items()},
      'all_authenticated_actual_first_field_errors_retained':True,
      'upstream_ten_capsule_authentication_before_and_after':inputs['all_ten_complete_capsules_authenticated_before_project_import_and_after_run'],
      'limits':['Conditional registered gauge and chosen slicing observer U with spatial y normal.',
                'Geometric formula validity and angular projection rely on the separate theory and audit artifacts.',
                'No observer-invariant, gauge-invariant, detector-force or nonlinear-family claim.'],
      'arithmetic':'Exact rational intervals, with exact symbolic equality between grouped and expanded expressions. Their intersection is valid for the same retained actual fields.'}
    with args.output.open('x') as handle:
        json.dump(report,handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps({'output':str(args.output),'sha256':sha(args.output),'status':report['status'],
      'intervals':{name:val['outward_decimal'] for name,val in report['certified_intersections'].items()}}))


if __name__ == '__main__':
    main()
