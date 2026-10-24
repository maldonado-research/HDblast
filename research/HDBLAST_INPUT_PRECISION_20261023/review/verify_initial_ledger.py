#!/usr/bin/env python3
"""Independent exact ledger and coverage audit, not a numeric-profile replay."""
import argparse,hashlib,json,math
from fractions import Fraction as F
from pathlib import Path
checks=[]
def need(ok,msg):
 if not ok:raise ValueError(msg)
 checks.append(msg)
def sqrtup(x):
 grid=1<<112;return F(math.isqrt((x*grid*grid).__ceil__())+1,grid)
def hnorm(cs):return sqrtup(sum(w*c*c for w,c in zip((1,2,1,3,1),cs)))
def free(cs):return max((cs[0]+F(3,2)*cs[2])/3,cs[1]+cs[2])
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--report',type=Path,required=True);ap.add_argument('--branch',type=Path,default=Path('/tmp/hdblast-input-precision-20261023/branch/BRANCH_COMPARISON_REPORT.json'));ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 r=json.loads(args.report.read_text());b=json.loads(args.branch.read_text());bh=hashlib.sha256(args.branch.read_bytes()).hexdigest()
 need(bh=='024acc3720696322e9fa2f864de4c75b17c88d81ec1a6b80c020e928cefd9cd0','independently_reviewed_branch_certificate_hash')
 need(r['case']=='branch' and r['branch_input_certificate']['sha256']==bh,'actual_branch_join')
 need(r['means']==['K0','L0','N0','J0','w0'] and r['gram']==[1,2,1,3,1],'registered_mean_order_and_gram')
 need(F(r['ell'])==F(1,1000000),'physical_scale')
 beta=list(map(F,r['beta_interval']));exact=list(map(F,b['beta']))
 need(beta[0]<=exact[0]<=exact[1]<=beta[1] and beta[1]-beta[0]<F(1,10**16),'outward_certified_beta_interval')
 need(r['partition_indices']==[128,224] and r['whole_cells']==96 and r['full_footprint_covered'],'complete_footprint_metadata')
 need(len(r['rows'])==96,'complete_row_count')
 caps=[[[F(0) for _ in range(5)] for _ in range(3)] for _ in range(3)]
 for idx,row in enumerate(r['rows'],128):
  need(row['index']==idx and list(map(F,row['r']))==[F(idx,128),F(idx+1,128)],'whole_cell_'+str(idx))
  vals=row['channel_component_error_caps'];need(len(vals)==3 and all(len(v)==3 and all(len(q)==5 for q in v) for v in vals),'channel_shape_'+str(idx))
  for c in range(3):
   for k in range(3):
    for m in range(5):
     v=F(vals[c][k][m]);need(v>=0,'nonnegative_cap_'+str(idx)+'_'+str(c)+'_'+str(k)+'_'+str(m));caps[c][k][m]=max(caps[c][k][m],v)
 need(caps==[[list(map(F,row)) for row in chan] for chan in r['channel_component_error_caps']],'reported_channel_component_maxima')
 hs=[[hnorm(row) for row in chan] for chan in caps]
 need(hs==[list(map(F,row)) for row in r['channel_H_caps']],'weighted_H_caps_exact')
 total=[sum(hs[c][k] for c in range(3)) for k in range(3)]
 need(total==list(map(F,r['total_H_caps'])),'channel_triangle_total_exact')
 cap=free(total);need(cap==F(r['free_wave_X_cap']),'whole_cone_free_wave_formula_exact')
 pass_actual=cap<=10**28
 need(r['passes_1e28']==pass_actual,'actual_gate_boolean')
 need(r['reauthenticated_after'] and r['actual_error_lower_bound_proved'] is False,'reauthentication_and_scope_flags')
 coefficients=r['representation_background_coefficients_ascending'];need(set(coefficients)=={'rho','eta','p','a','alpha'} and all(len(v)==5 for v in coefficients.values()),'one_degree_four_nominal_background')
 sh=hashlib.sha256(json.dumps(coefficients,sort_keys=True,separators=(',',':')).encode()).hexdigest();need(sh==r['representation_background_coefficients_sha256'],'nominal_coefficient_identity')
 need(F(coefficients['alpha'][0])==1,'alpha_initial_normalization')
 # The design table is sufficient but not necessary; report failed rows
 # separately from the exact free-wave gate to avoid silently reclassifying it.
 allocations=[[F(5*10**27),F(2*10**27),F(5*10**26)],[F(10**28),F(4*10**27),F(10**27)],[F(5*10**27),F(2*10**27),F(5*10**26)]]
 failures=[{'channel':r['channels'][c],'quantity':['delta0_H','deltar_H','deltav_H'][k],'actual_cap':str(hs[c][k]),'actual_cap_display':float(hs[c][k]),'design_allocation':str(allocations[c][k])} for c in range(3) for k in range(3) if hs[c][k]>allocations[c][k]]
 total_allocation=[F(2*10**28),F(8*10**27),F(2*10**27)]
 totalfail=[{'quantity':['delta0_H','deltar_H','deltav_H'][k],'actual_cap':str(total[k]),'design_allocation':str(total_allocation[k])} for k in range(3) if total[k]>total_allocation[k]]
 out={'status':'PASS_INDEPENDENT_INITIAL_LEDGER_AUDIT','initial_accuracy_gate_passes':pass_actual,'candidate_repaired':r['candidate_changed'],'checks_passed':len(checks),'report_sha256':hashlib.sha256(args.report.read_bytes()).hexdigest(),'branch_sha256':bh,'candidate_coefficients_sha256':sh,'whole_cells':96,'channel_H_caps_exact':[[str(v) for v in row] for row in hs],'total_H_caps_exact':list(map(str,total)),'free_wave_X_cap_exact':str(cap),'free_wave_X_cap_display':float(cap),'failed_proposed_design_rows':failures,'failed_proposed_total_rows':totalfail,'review_scope':'Exact per-cell coverage, maxima, Gram norms, triangle sum, free-wave bound, actual-beta dependency join, candidate identity and gate flags; does not replay profile interval arithmetic. Separate source/proof review is needed.','limits':['Initial/free-wave bound only.','Illustrative design component limits may fail although the exact free-wave gate passes.','No total U, B1, observable sign or error lower bound inferred.']}
 args.output.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],out['initial_accuracy_gate_passes'],float(cap))
if __name__=='__main__':main()
