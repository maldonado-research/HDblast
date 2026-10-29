#!/usr/bin/env python3
"""Compare independent local-quadrature and cone-norm candidate enclosures."""
from fractions import Fraction as F
from pathlib import Path
import argparse,json

def need(c,m):
    if not c:raise RuntimeError(m)
def interval(v):
    lo,hi=map(F,v);need(lo<=hi,'interval order');return lo,hi
p=argparse.ArgumentParser();p.add_argument('--complete',required=True,type=Path);p.add_argument('--coarse',required=True,type=Path);a=p.parse_args()
x=json.loads(a.complete.read_text());y=json.loads(a.coarse.read_text())
need(x['status']==y['status']=='PASS' and x['point']==y['point']==['1/1000','1499/1000'] and x['component']==y['component']=='w:0','same successful coefficient point')
need(len(x['full23_point_table'])==23,'complete candidate table')
checks=0
for k in ('value','dx','dr'):
    d=interval(x['candidate_intervals'][k]);c=interval(y['candidate_component_intervals'][k])
    need(c[0]<=d[0]<=d[1]<=c[1],'detailed candidate outside coarse bound: '+k);checks+=1
    cap=F(y['B1_value_H_cap'] if k=='value' else y['B1_derivative_H_cap']);b=interval(x['B1'][k])
    need(-cap<=b[0]<=b[1]<=cap,'detailed feedback outside cone norm cap: '+k);checks+=1
    u=interval(x['conditional_actual_formal_response_intervals'][k]);v=interval(y['actual_formal_component_intervals'][k])
    need(max(u[0],v[0])<=min(u[1],v[1]),'disjoint independently transferred intervals: '+k);checks+=1
need(F(x['C1_error_exact_recomputed_upper'])==F(y['global_C1_error']),'independent transfer arithmetic agrees')
need(interval(x['conditional_actual_formal_response_intervals']['value'])[0]>0 and interval(y['actual_formal_component_intervals']['value'])[0]>0,'positive registered scalar mean')
need(interval(x['conditional_actual_formal_response_intervals']['dr'])[1]<0 and interval(y['actual_formal_component_intervals']['dr'])[1]<0,'negative registered radial derivative')
for data,key in ((x,'conditional_actual_formal_response_intervals'),(y,'actual_formal_component_intervals')):
    lo,hi=interval(data[key]['dx']);need(lo<0<hi,'time derivative remains undecided')
print(json.dumps({'status':'PASS','cross_method_interval_checks':checks,'full23_rows':23,'same_exact_transfer_bound':True,'primary_formal_signs':'value positive, radial derivative negative, time derivative undecided','scope':'one registered formal modal point; not a C2 error, Einstein residual or physical observable'},sort_keys=True))
