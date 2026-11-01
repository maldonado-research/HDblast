#!/usr/bin/env python3
"""Independent aggregation reviewer controls; mutates data, never source."""
from pathlib import Path
import argparse,copy,importlib.util,json
spec=importlib.util.spec_from_file_location('independent_review',Path(__file__).resolve().parent/'review_trace.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--raw-report',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 raw=json.loads(a.raw_report.read_text());pin=json.loads((Path(__file__).resolve().parent/'SOURCE_PINS.json').read_text())['raw_box_projection_sha256'];results=[]
 def reject(name,change,guard=False,expected=None):
  r=copy.deepcopy(raw);change(r)
  try:
   if guard:m.guard(r,pin)
   else:m.reconstruct(r)
  except m.AuditError as e:
   if expected is not None:m.need(expected==str(e),'CONTROL_WRONG_REJECTION:'+name+':'+str(e))
   results.append({'name':name,'status':'REJECTED','reason':str(e)});return
  raise RuntimeError('BAD_MUTATION_ACCEPTED:'+name)
 def multiply(r,key,n):
  row=r['target_rows'][-1];row[key]=[str(m.F(x)*n) for x in row[key]]
 reject('wrong_kernel_power',lambda r:r['target_rows'][-1]['forcing_u_terms'][0].update({'exact_weight':r['target_rows'][-1]['forcing_ux_terms'][0]['exact_weight']}),expected='ENCLOSURE:weight8:1')
 reject('point_rectangle_instead_of_integral_weight',lambda r:r['target_rows'][-1]['forcing_u_terms'][0].update({'exact_weight':[str((m.F(5,4)-m.F(r['history_rows'][0]['scaled_x'][0]))*(m.F(r['history_rows'][0]['scaled_x'][1])-m.F(r['history_rows'][0]['scaled_x'][0])))]*2}),expected='ENCLOSURE:weight8:1')
 reject('missing_half_in_double_integrated_forcing',lambda r:multiply(r,'forcing_integrated_u',2),expected='ENCLOSURE:target8:forcing_integrated_u')
 reject('missing_rho_in_proper_time_derivative',lambda r:multiply(r,'u_prime',130),expected='ENCLOSURE:target8:u_prime')
 reject('missing_ell_in_coordinate_derivative',lambda r:multiply(r,'recovered_u_t',m.F(1,10**6)),expected='ENCLOSURE:target8:recovered_u_t')
 reject('second_derivative_factor_two',lambda r:multiply(r,'curvature_defect_second',m.F(1,2)),expected='ENCLOSURE:target8:curvature_defect_second')
 reject('drop_feedback_acceleration_error',lambda r:r.update({'feedback_acceleration_error_cap':'0'}),expected='UPPER:feedback_acceleration_error_cap')
 reject('history_gap',lambda r:r['history_rows'][1].update({'scaled_x':['1/10',r['history_rows'][1]['scaled_x'][1]]}),expected='HISTORY_COVER')
 reject('drop_complete_actual_error_premise',lambda r:r['scope'].update({'actual_scalar_source_with_complete_inherited_errors':False}),guard=True,expected='ACTUAL_SOURCE_ERROR_GUARD:actual_scalar_source_with_complete_inherited_errors')
 reject('precontact_zero_shortcut',lambda r:r['scope'].update({'finite_propagation_zero_before_one_used':True}),guard=True,expected='NO_PRECONTACT_ZERO_SHORTCUT')
 reject('shrink_raw_source_box_to_comparison',lambda r:r['history_rows'][0]['actual_scalar'].update({'s':['0','0']}),guard=True,expected='RAW_BOX_PROJECTION_PIN')
 # Independent exact polynomial controls, including truncated and not-yet-active cells.
 known=[(m.F(3),m.F(0),m.F(1),0,m.F(1)),(m.F(3),m.F(0),m.F(1),1,m.F(5,2)),(m.F(3),m.F(0),m.F(1),2,m.F(19,3)),(m.F(1,2),m.F(0),m.F(1),2,m.F(1,24)),(m.F(0),m.F(1),m.F(2),2,m.F(0))]
 for x,lo,hi,k,value in known:m.need(m.weight(x,lo,hi,k)==value,'EXACT_KERNEL_FIXTURE')
 m.need(m.I(-2,3).sq().lo==0 and m.I(-2,3).sq().hi==9,'ZERO_CROSSING_SQUARE')
 report={'status':'PASS','expected_rejections':results,'exact_kernel_fixtures':len(known),'zero_crossing_square_verified':True,'oracle_campaigns':0,'producer_imports':False}
 m.need(not a.output.exists(),'FRESH_OUTPUT');a.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');print(json.dumps({'status':'PASS','expected_rejections':len(results),'output_sha256':m.sha(a.output)},sort_keys=True))
if __name__=='__main__':main()
