#!/usr/bin/env python3
"""Read-only failure/provenance audit; never calls scientific producers."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path

def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--replay',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    root=args.replay.resolve();require(not args.output.exists(),'Refusing overwrite')
    state=read(root/'VALIDATION.json');records=read(root/'EXECUTION.json');plan=read(root/'COMMAND_PLAN.json')
    partial=read(root/'fresh/primary/partial_results.json');failure=read(root/'fresh/primary/failure.json');start=read(root/'fresh/primary/started.json')
    reg=read(root/'checkpoint/FULL_REGISTRATION.json');experiment=read(root/'checkpoint/EXPERIMENT.json')
    require(state['status']=='FAIL' and state['source_verification_before']==state['source_verification_after']=='PASS','Original classification/source verification differs')
    require(len(plan)==32 and len(records)==27 and len([x for x in records if x['status']=='PASS'])==26,'Actual command coverage differs')
    require(all(not x['physical_producer'] and x['exit_code']==0 for x in records[:-1]),'Unexpected physical or failed early command')
    require(records[-1]['check']=='primary' and records[-1]['physical_producer'] and records[-1]['status']=='FAIL' and records[-1]['exit_code']==1 and records[-1]['timed_out'] is False,'Primary stop differs')
    remaining=[c['check'] for c in plan[len(records):]]
    require(remaining==['independent_modes','validation','validation_optimized','summary','figures'],'Unattempted stage coverage differs')
    require('IntegrationWarning' in failure['exception'] and 'roundoff' in failure['exception'],'Failure type differs')
    require(len(partial['rows'])==1 and partial['active']=={'source':'positive_B','eta':-4.5},'Preserved partial progress differs')
    row=partial['rows'][0]
    require(row['source']=='positive_B' and row['eta']==-5.5 and [x['K'] for x in row['finite_k']]==[64,128,256],'Completed zero-row coverage differs')
    response=['q','q_prime','q_second','rho','p','current']
    for point in [row['continuum'],*row['finite_k']]:require(all(point['values'][x]==0 for x in response),'Completed row is not source-free zero response')
    require(not (root/'fresh/primary/results.json').exists() and not (root/'fresh/independent').exists(),'Unexpected completed physical data')
    source_pins=[]
    for name,pin in reg['files'].items():
        q=root/'checkpoint'/name;require(q.is_file() and sha(q)==pin,'Frozen source mismatch: '+name);source_pins.append(name)
    for name,pin in state['fresh_output_sha256'].items():require(sha(root/'fresh'/name)==pin,'Captured output hash differs: '+name)
    evidence=[root/'VALIDATION.json',root/'EXECUTION.json',root/'COMMAND_PLAN.json',root/'INPUT_MANIFEST.json',root/'fresh/primary/partial_results.json',root/'fresh/primary/failure.json',root/'fresh/primary/started.json',root/'logs/primary.stderr.log',root/'logs/primary.stdout.log',root/'logs/primary.exit.json',root/'checkpoint/FULL_REGISTRATION.json',root/'checkpoint/EXPERIMENT.json',root/'checkpoint/FREEZE_RECEIPT.json']
    result={'audit_status':'PASS_READ_ONLY_FAILURE_CLASSIFICATION','original_registered_calibration_status':'FAIL','recorded_utc':datetime.now(timezone.utc).isoformat(),'source_replay_directory':str(root),'new_physical_evaluations':0,'public_freeze_commit':start['public_freeze_commit'],'registration_sha256':start['registration_sha256'],'independent_manifest_sha256':read(root/'checkpoint/FREEZE_RECEIPT.json')['independent_manifest_sha256'],'planned_commands':len(plan),'attempted_commands':len(records),'pure_commands_completed_pass':len(records)-1,'physical_producers_attempted':1,'physical_producers_completed':0,'failed_command':'primary','failed_primary_exit_code':1,'timed_out':False,'warning':failure['exception'],'primary_process_elapsed_seconds':records[-1]['elapsed_seconds'],'producer_failure_report_elapsed_seconds':failure['elapsed_seconds'],'completed_rows':1,'completed_nonzero_response_rows':0,'completed_source_free_row':{'source':row['source'],'eta':row['eta'],'finite_cutoffs':[x['K'] for x in row['finite_k']],'response_fields_zero':response},'active_row_when_failed':partial['active'],'continuum_derivative_order_that_triggered_warning':'not recorded; cannot infer from preserved traceback','not_attempted_commands':remaining,'frozen_sources_reverified':len(source_pins),'captured_output_hashes_reverified':len(state['fresh_output_sha256']),'source_verification_before_and_after':'PASS','primary_completed_results_exist':False,'independent_mode_results_exist':False,'validators_executed':False,'physical_instability_established':False,'model_rejected_by_physical_evidence':False,'coupled_stability_or_heating_established':False,'original_primary_configuration':experiment['primary'],'original_gates':experiment['gates'],'evidence_files':[{'path':str(p.relative_to(root)),'sha256':sha(p),'bytes':p.stat().st_size} for p in evidence],'audit_source_sha256':sha(Path(__file__))}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('audit_status','original_registered_calibration_status','planned_commands','attempted_commands','pure_commands_completed_pass','physical_producers_completed','frozen_sources_reverified')}))

if __name__=='__main__':main()
