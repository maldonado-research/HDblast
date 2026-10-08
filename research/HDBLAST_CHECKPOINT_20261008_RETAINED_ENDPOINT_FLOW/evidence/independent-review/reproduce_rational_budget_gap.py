"""Fabricated source errors expose many-denominator serialization growth."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import sys
sys.dont_write_bytecode=True
BASE=Path(__file__).resolve().parent
ROOT=BASE.parent/'implementation'
sys.path.insert(0,str(ROOT/'engine'))
from endpoint_engine import EndpointFlow,geometry,ENDPOINTS

rows=[dict(g,coefficients=(Q(0),)*25,uniform_error=Q(1,(16384-(2*j-63)**2)**201))
      for j,g in enumerate(geometry())]
results=[]
for endpoint in ENDPOINTS:
    flow=EndpointFlow(rows,endpoint)
    entry={'endpoint':str(endpoint),'source_error_W_denominator_bits':flow.source_error_W.denominator.bit_length(),
           'source_error_U_denominator_bits':flow.source_error_U.denominator.bit_length()}
    try:flow.budget()
    except ValueError as error:entry.update(budget_serialization_passed=False,error=str(error))
    else:entry['budget_serialization_passed']=True
    results.append(entry)
receipt={'status':'INITIAL_MANY_DENOMINATOR_BUDGET_CONTROL','cases':results,
         'engine_sha256':hashlib.sha256((ROOT/'engine/endpoint_engine.py').read_bytes()).hexdigest(),
         'retained_decodes':0,'physical_source_calls':0}
with (BASE/'INITIAL_RATIONAL_BUDGET_GAP.json').open('x') as out:json.dump(receipt,out,sort_keys=True,indent=2);out.write('\n')
print(json.dumps(receipt,sort_keys=True))
