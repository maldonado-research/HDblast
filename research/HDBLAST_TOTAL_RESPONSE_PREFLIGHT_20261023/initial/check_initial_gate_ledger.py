#!/usr/bin/env python3
"""Portable exact report arithmetic check; does not replay numerical producer."""
import argparse,json,hashlib,math
from fractions import Fraction as F
from pathlib import Path

def need(ok,msg):
 if not ok:raise ValueError(msg)
def run(path):
 raw=path.read_bytes();r=json.loads(raw);need(r['status']=='FAIL_PROVISIONAL_INITIAL_FREE_WAVE_ACCURACY_GATE_BOUND_INSUFFICIENT','status')
 need(not r['actual_error_lower_bound_proved'] and not r['mathematical_impossibility_claimed'],'honest failure scope')
 e=r['error_budget'];caps=[[F(x) for x in rr] for rr in e['component_error_caps']];acc=[[F(0)]*5 for _ in range(3)]
 rows=r['partition']['rows'];edge=F(1)
 for row in rows:
  lo,hi=map(F,row['r']);need(lo==edge and hi-lo==F(r['partition']['step']),'whole footprint cover');edge=hi
  er=[[F(x) for x in rr] for rr in row['component_error_caps']]
  need(all(x>=0 for rr in er for x in rr),'nonnegative error caps')
  acc=[[max(x,y) for x,y in zip(cr,rr)] for cr,rr in zip(acc,er)]
 need(edge==F(7,4) and acc==caps,'exact whole footprint accumulation')
 delta=list(map(F,e['H_error_caps']))
 for ds,cc in zip(delta,caps):need(ds>=0 and ds*ds>=sum(w*c*c for w,c in zip((1,2,1,3,1),cc)),'weighted H norm')
 ef=max((delta[0]+F(3,2)*delta[2])/3,delta[1]+delta[2]);need(F(e['free_wave_X_error_cap'])==ef,'whole cone C1 free wave estimate')
 need(not e['passes'] and ef>10**28,'failed provisional bound gate')
 for p in r['point_preflight']:
  aa=p['actual_initial_intervals'];radius0=(F(aa[0][-1][1])-F(aa[0][-1][0]))/2;radiusr=(F(aa[1][-1][1])-F(aa[1][-1][0]))/2
  need(radius0==F(p['scalar_value_oracle_radius']) and radiusr==F(p['scalar_radial_oracle_radius']),'point oracle radius')
 need(r['reauthenticated_after_evaluation'],'postevaluation authentication')
 return {'status':'PASS_PORTABLE_INITIAL_GATE_LEDGER','producer_report_sha256':hashlib.sha256(raw).hexdigest(),'whole_cells':len(rows),'free_wave_X_error_cap':str(ef),'scope':'Exact arithmetic and report consistency only; numerical-input and functional enclosure replay requires certify_initial_representation.py and authenticated dependencies.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--report',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run(a.report);a.output.write_text(json.dumps(r,indent=2)+'\n');print(r['status'])
