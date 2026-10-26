#!/usr/bin/env python3
"""Stdlib exact interval evaluator for an explicitly supplied encoded map.

Required map shape:
  {"schema":"...", "components":{
     "component-name":[{"coefficient":"2/3", "powers":{"A":1,"By":1}}]
  }}
All coefficients/endpoints are rational strings or integers. No curvature
formula is embedded here; an authenticated new map is a required CLI input.
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
from output_safety import guard_output, source_snapshot, verify_source_snapshot, source_labels


def need(ok, message):
    if not ok:
        raise ValueError(message)


def rational(value):
    need(type(value) in (str,int), 'EXACT_RATIONAL_REQUIRED')
    return F(value)


class Interval:
    """Exact rational endpoint operations, with no floating point decisions."""
    __slots__ = ('lo','hi')

    def __init__(self, lo=0, hi=None):
        self.lo = lo if isinstance(lo,F) else rational(lo)
        self.hi = self.lo if hi is None else (hi if isinstance(hi,F) else rational(hi))
        need(self.lo <= self.hi, 'INTERVAL_ORDER')

    @staticmethod
    def cv(value):
        if isinstance(value,Interval):
            return value
        if isinstance(value,F):
            return Interval(value)
        if isinstance(value,(list,tuple)):
            need(len(value) == 2, 'INTERVAL_ENDPOINT_PAIR')
            return Interval(*value)
        return Interval(value)

    def __add__(self, other):
        other = self.cv(other)
        return Interval(self.lo+other.lo,self.hi+other.hi)
    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi,-self.lo)

    def __sub__(self, other):
        return self+-self.cv(other)

    def __rsub__(self, other):
        return self.cv(other)+-self

    def __mul__(self, other):
        other = self.cv(other)
        vals = [a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Interval(min(vals),max(vals))
    __rmul__ = __mul__

    def inv(self):
        need(self.lo > 0 or self.hi < 0, 'ZERO_DENOMINATOR')
        return Interval(1/self.hi,1/self.lo)

    def __truediv__(self, other):
        return self*self.cv(other).inv()

    def __rtruediv__(self, other):
        return self.cv(other)*self.inv()

    def __pow__(self, power):
        need(type(power) is int, 'INTEGER_POWER_REQUIRED')
        if power < 0:
            return self.inv()**(-power)
        if power == 0:
            return Interval(1)
        if power % 2:
            return Interval(self.lo**power,self.hi**power)
        return Interval(0 if self.lo <= 0 <= self.hi else min(self.lo**power,self.hi**power),
                        max(self.lo**power,self.hi**power))

    def row(self):
        return [str(self.lo),str(self.hi)]

    def mag(self):
        return max(abs(self.lo),abs(self.hi))

    def sign(self):
        return 'positive' if self.lo > 0 else ('negative' if self.hi < 0 else 'inconclusive')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def evaluate_encoded(terms, symbols):
    """Evaluate one explicit list of rational coefficient/powers terms."""
    need(isinstance(terms,list), 'ENCODED_TERM_LIST_REQUIRED')
    out, cache = Interval(), {}
    for term in terms:
        need(isinstance(term,dict) and set(term) == {'coefficient','powers'}, 'ENCODED_TERM_SHAPE')
        product = Interval(rational(term['coefficient']))
        need(isinstance(term['powers'],dict), 'POWERS_MAPPING_REQUIRED')
        for name,power in sorted(term['powers'].items()):
            need(isinstance(name,str) and name in symbols, 'UNKNOWN_SYMBOL:'+str(name))
            need(type(power) is int, 'INTEGER_POWER_REQUIRED:'+name)
            key = (name,power)
            if key not in cache:
                cache[key] = Interval.cv(symbols[name])**power
            product *= cache[key]
        out += product
    return out


def evaluate_map(encoded, symbols):
    need(isinstance(encoded,dict) and isinstance(encoded.get('components'),dict), 'COMPONENT_MAP_REQUIRED')
    return {name:evaluate_encoded(terms,symbols) for name,terms in encoded['components'].items()}


def controls():
    """Exact independent arithmetic controls; no actual solution calls."""
    env = {'A':Interval(-2,3),'By':Interval(4,5),'rho':Interval(2,4)}
    tests = []
    def equal(name, value, endpoints):
        need(value.row() == endpoints, 'CONTROL:'+name)
        tests.append({'name':name,'passed':True,'exact':value.row()})
    equal('cross_product',evaluate_encoded([{'coefficient':'2/3','powers':{'A':1,'By':1}}],env),['-20/3','10'])
    equal('square_crossing_zero',evaluate_encoded([{'coefficient':'-1/2','powers':{'A':2}}],env),['-9/2','0'])
    equal('negative_power',evaluate_encoded([{'coefficient':'3','powers':{'rho':-2}}],env),['3/16','3/4'])
    equal('empty_map_is_zero',evaluate_encoded([],env),['0','0'])
    for name, terms, symbols, message in [
      ('missing_symbol',[{'coefficient':'1','powers':{'missing':1}}],env,'UNKNOWN_SYMBOL:missing'),
      ('zero_denominator',[{'coefficient':'1','powers':{'A':-1}}],env,'ZERO_DENOMINATOR'),
      ('float_coefficient',[{'coefficient':0.5,'powers':{}}],env,'EXACT_RATIONAL_REQUIRED'),
      ('noninteger_power',[{'coefficient':'1','powers':{'A':'2'}}],env,'INTEGER_POWER_REQUIRED:A')]:
        try:
            evaluate_encoded(terms,symbols)
        except ValueError as error:
            need(str(error) == message, 'CONTROL_WRONG_REJECTION:'+name)
            tests.append({'name':name,'passed':True,'expected_rejection':message})
        else:
            raise ValueError('CONTROL_MISSING_REJECTION:'+name)
    return {'schema':'HDBLAST_EXACT_ENCODED_MAP_ARITHMETIC_CONTROLS_V1','all_passed':True,'controls':tests,
            'actual_solution_calls':0,'scope':'Arithmetic controls only; no magnetic Weyl formula is certified.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs',type=Path)
    parser.add_argument('--map',type=Path,dest='map_path')
    parser.add_argument('--require-map-sha256')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--controls',action='store_true')
    args = parser.parse_args()
    guard_output(args.output,source_file=__file__,input_paths=[args.inputs,args.map_path])
    sources_before = source_snapshot([__file__,HERE/'output_safety.py'])
    if args.controls:
        result = controls()
    else:
        need(args.inputs is not None and args.map_path is not None and args.require_map_sha256 is not None,
             'EXPLICIT_INPUTS_AND_EXTERNALLY_PINNED_MAP_REQUIRED')
        need(sha(args.map_path) == args.require_map_sha256,'MAP_EXTERNAL_PIN')
        encoded = json.loads(args.map_path.read_text())
        inputs = json.loads(args.inputs.read_text())
        outputs = evaluate_map(encoded,inputs['symbols'])
        result = {'schema':'HDBLAST_SCREEN_MAGNETIC_SUPPLIED_MAP_EVALUATION_V1','status':'EVALUATED_SUPPLIED_MAP_ONLY',
          'input_sha256':sha(args.inputs),'map_sha256':sha(args.map_path),'evaluator_sha256':sources_before[str(Path(__file__).resolve())],
          'map_schema':encoded.get('schema'),
          'arithmetic':'Exact rational interval endpoints; no binary floating point decisions.',
          'components':{name:{'exact':value.row(),'sign':value.sign(),
            'term_count':len(encoded['components'][name]),
            'symbols':sorted({symbol for term in encoded['components'][name] for symbol in term['powers']}),
            'display_only':[float(value.lo),float(value.hi)]} for name,value in outputs.items()},
          'formula_validity':'The supplied map requires independent geometric derivation and controls; evaluation does not establish it.'}
        need(sha(args.map_path) == args.require_map_sha256,'MAP_CHANGED_DURING_EVALUATION')
    verify_source_snapshot(sources_before)
    result['source_integrity_before_and_after'] = source_labels(sources_before)
    with args.output.open('x') as handle:
        json.dump(result,handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps({'output':str(args.output),'sha256':sha(args.output),'all_passed':result.get('all_passed')}))


if __name__ == '__main__':
    main()
