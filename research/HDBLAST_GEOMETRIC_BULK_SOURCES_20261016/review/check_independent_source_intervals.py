#!/usr/bin/env python3
"""Independent exact endpoint interval evaluator for exported source maps.

It reads the actual-jet oracle as an authenticated analytic premise, but does
not use its arithmetic class or the producer coefficient evaluator for the
independent source enclosure calculation. No numerical midpoint is substituted.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, sys
ROOT=Path(__file__).resolve().parents[1]
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'numeric'));sys.path.insert(0,str(ROOT/'projection'))
import actual_first_jets as actual
import source_coefficients as producer

class Box:
    def __init__(self,lo=0,hi=None):
        self.lo=F(lo);self.hi=self.lo if hi is None else F(hi)
        if self.lo>self.hi:raise ValueError('REVERSED_INDEPENDENT_BOX')
    @staticmethod
    def cv(x):return x if isinstance(x,Box) else Box(x)
    def __add__(self,other):
        other=Box.cv(other);return Box(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __mul__(self,other):
        other=Box.cv(other);ends=[a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Box(min(ends),max(ends))
    __rmul__=__mul__
    def __pow__(self,n):
        if n==0:return Box(1)
        if n<0:
            if self.lo<=0<=self.hi:raise ValueError('INDEPENDENT_DENOMINATOR_CROSSES_ZERO')
            return Box(1/self.hi,1/self.lo)**(-n)
        if n%2:return Box(self.lo**n,self.hi**n)
        return Box(0 if self.lo<=0<=self.hi else min(self.lo**n,self.hi**n),max(self.lo**n,self.hi**n))
    def row(self):return [str(self.lo),str(self.hi)]

def one_coefficient(terms,values):
    out=Box()
    for term in terms:
        value=Box(F(term['coefficient']))
        for name,power in sorted(term['powers'].items()):value=value*values[name]**power
        out=out+value
    return out

def verify_maps(report):
    first={name+suffix for name in 'ABCDELM' for suffix in ('','y','t','yy','yt','tt')}|{'s','sy','st','syy','syt','stt'}
    count=0
    for name,fields in report['maps'].items():
        rows=[row for values in fields.values() for row in values] if isinstance(fields,dict) else fields
        for terms in rows:
            for term in terms:
                powers=term['powers']
                if any(powers.get(k,0)<0 for k in first):raise ValueError('NEGATIVE_FIRST_JET_POWER')
                degree=sum(powers.get(k,0) for k in first)
                if degree!=2:raise ValueError('SOURCE_NOT_QUADRATIC_'+name)
                count+=1
    return count

def verify_cell(report,path,label,x,r):
    m,s,b=actual.actual_jet(actual.I(*x),actual.I(*r))
    raw={**m,**s,**b}
    values={key:Box(value.lo,value.hi) if isinstance(value,actual.I) else Box(value) for key,value in raw.items()}
    names=('QRic','QE','QEred','QC','Qphi','metric_forcing','scalar_forcing')
    prod=producer.evaluate_sources(m,s,b,names,path)
    selected={};ncoeff=0
    for name in names:
        fields=report['maps'][name]
        mapping=fields if isinstance(fields,dict) else {'scalar':fields}
        for field,rows in mapping.items():
            for slot,terms in enumerate(rows):
                independent=one_coefficient(terms,values)
                observed=prod[name][field][slot] if isinstance(fields,dict) else prod[name][slot]
                if not(observed.lo<=independent.lo<=independent.hi<=observed.hi):
                    raise ValueError(f'PRODUCER_INTERVAL_NOT_OUTWARD_{name}_{field}_{slot}')
                ncoeff+=1
                if slot==0 and ((name in ('QE','QEred','metric_forcing') and field=='yy') or name in ('Qphi','scalar_forcing')):
                    selected[f'{name}:{field}:mean']={'exact_independent_interval':independent.row(),'display_only':[float(independent.lo),float(independent.hi)]}
    if label=='evolved_exterior_packet':
        if not F(selected['QE:yy:mean']['exact_independent_interval'][1])<0:raise ValueError('ACTUAL_NATURAL_EINSTEIN_SIGN_NOT_SEPARATED')
        if not F(selected['QEred:yy:mean']['exact_independent_interval'][1])<0:raise ValueError('ACTUAL_REDUCED_EINSTEIN_SIGN_NOT_SEPARATED')
    return {'label':label,'x':list(map(str,x)),'r':list(map(str,r)),'coefficients_checked':ncoeff,'selected':selected,
            'scope':'Continuous closed cell and whole S3 through exact projected coefficients; inherited actual C3 error premise, not comparator values alone.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--maps',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    if args.output.exists():raise ValueError('OUTPUT_ALREADY_EXISTS')
    raw=args.maps.read_bytes();report=json.loads(raw)
    if report['status']!='PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION':raise ValueError('MAP_STATUS')
    for name,digest in report['engine_pins'].items():
        if hashlib.sha256((ROOT/'geometry'/name).read_bytes()).hexdigest()!=digest:raise ValueError('MAP_ENGINE_PIN_'+name)
    terms=verify_maps(report)
    h=F(1,4096)
    cells=[verify_cell(report,args.maps,'evolved_exterior_packet',(F(1,4),F(1,4)+h),(F(5,4),F(5,4)+h)),
           verify_cell(report,args.maps,'reflected_packet',(F(3,2),F(3,2)+h),(F(1,4),F(1,4)+h))]
    inputs=[Path(__file__),ROOT/'numeric/actual_first_jets.py',ROOT/'numeric/exterior_error_caps.py',ROOT/'numeric/interval_core.py',ROOT/'projection/source_coefficients.py']
    out={'status':'PASS_INDEPENDENT_EXACT_SOURCE_INTERVAL_EVALUATION','maps_sha256':hashlib.sha256(raw).hexdigest(),
         'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
         'quadratic_monomials_checked':terms,'cells':cells,
         'arithmetic':'Independent exact Fraction endpoint interval arithmetic with no rounding; producer rounded intervals contain every independently evaluated result.',
         'not_claimed':'No new full-domain cap, actual background selection, gauge error reduction, second-response solve, or nonlinear remainder.'}
    args.output.write_text(json.dumps(out,indent=2)+'\n');print(out['status'])

if __name__=='__main__':main()
