"""Regression for compact certified error transport; no physical builder calls."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import sys
sys.dont_write_bytecode=True
BASE=Path(__file__).resolve().parent
ROOT=BASE.parent/'implementation'
sys.path[:0]=[str(ROOT/'engine'),str(ROOT/'source')]
from endpoint_engine import EndpointFlow,geometry,ENDPOINTS,SCALE
from later_source import outward_error

def require(ok,msg):
    if not ok:raise RuntimeError(msg)
tail=Q(1,15*2**90)
raws=[tail+Q(1,(16384-(2*j-63)**2)**201) for j in range(64)]
rounded=[outward_error(q) for q in raws]
for before,after in zip(raws,rounded):
    require(before<=after<before+Q(1,2**512),'outward dyadic error not within one unit')
    require((after*2**512).denominator==1,'source error is not dyadic512')
require(outward_error(Q(0))==0 and outward_error(Q(3,2**512))==Q(3,2**512),'exact grid rounding changed')
for bad in (Q(-1),0.0):
    try:outward_error(bad)
    except ValueError:pass
    else:raise RuntimeError('invalid error type/sign accepted')
rows=[dict(g,coefficients=(Q(0),)*25,uniform_error=q) for g,q in zip(geometry(),rounded)]
cases=[]
for endpoint in ENDPOINTS:
    flow=EndpointFlow(rows,endpoint);budget=flow.budget();blob=json.dumps(budget,sort_keys=True)
    active=[(g,q) for g,q in zip(geometry(),raws) if g['right']<=endpoint]
    expected_w=sum((2*g['half']*q for g,q in active),Q(0))
    expected_u=sum((2*g['half']*(endpoint-g['center'])*q for g,q in active),Q(0))
    result=flow.evaluate(Q(0),(Q(0),Q(0)),(Q(0),Q(0)))
    for part,exact in (('U',expected_u),('W',expected_w)):
        lo,hi=result[part]['real']
        require(Q(lo,SCALE)<=-exact<=exact<=Q(hi,SCALE),'attained real-source uncertainty escaped')
    duration=endpoint+Q(9,2)
    require(Q(0)<=flow.source_error_W-expected_w<duration/Q(2**512),'W enlargement counted incorrectly')
    require(Q(0)<=flow.source_error_U-expected_u<duration*duration/Q(2**513),'U enlargement counted incorrectly')
    cases.append({'endpoint':str(endpoint),'budget_serialized_bytes':len(blob),
                  'original_error_denominator_bits':expected_w.denominator.bit_length(),
                  'rounded_error_denominator_bits':flow.source_error_W.denominator.bit_length()})
receipt={'status':'PASS_INDEPENDENT_SOURCE_ERROR_ROUNDING','cases':cases,'cells':64,
         'pure_error_rounding_calls':68,'physical_source_calls':0,'retained_decodes':0,
         'source_pins':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in (ROOT/'source/later_source.py',ROOT/'engine/endpoint_engine.py')}}
with (BASE/'SOURCE_ERROR_ROUNDING_CONTROL.json').open('x') as out:json.dump(receipt,out,sort_keys=True,indent=2);out.write('\n')
print(json.dumps(receipt,sort_keys=True))
