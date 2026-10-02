#!/usr/bin/env python3
"""Read-only provenance and coverage audit of two preserved failed rounds."""
import argparse
from datetime import datetime,timezone
import hashlib,json
from pathlib import Path

def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def read(p):return json.loads(p.read_text())
def sha(p):
    digest=hashlib.sha256()
    with p.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()
def verify_frozen(root):
    reg=read(root/'checkpoint/FULL_REGISTRATION.json')
    for path,pin in reg['files'].items():require(sha(root/'checkpoint'/path)==pin,'Frozen source differs: '+path)
    return len(reg['files'])

def main():
    p=argparse.ArgumentParser();p.add_argument('--replay',type=Path,required=True);p.add_argument('--original',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    require(not a.output.exists(),'Refusing overwrite')
    root=a.replay.resolve();old=a.original.resolve()
    state=read(root/'VALIDATION.json');records=read(root/'EXECUTION.json');plan=read(root/'COMMAND_PLAN.json')
    primary=read(root/'fresh/primary/results.json');modes=read(root/'fresh/independent/results.json');receipt=read(root/'checkpoint/FREEZE_RECEIPT.json')
    require(state['status']=='FAIL' and state['source_verification_before']==state['source_verification_after']=='PASS','Followup status/source verification differs')
    require(len(plan)==32 and len(records)==28 and len([r for r in records if r['status']=='PASS'])==27,'Command coverage differs')
    require(all(r['status']=='PASS' and not r['physical_producer'] for r in records[:26]),'Pure stage coverage differs')
    require(records[26]['check']=='primary' and records[26]['status']=='PASS' and records[26]['exit_code']==0,'Primary command differs')
    require(records[27]['check']=='independent_modes' and records[27]['status']=='FAIL' and records[27]['exit_code']==1 and not records[27]['timed_out'],'Independent command failure differs')
    require(primary['status']=='completed' and len(primary['rows'])==12 and sum(len(r['finite_k']) for r in primary['rows'])==36,'Primary completion coverage differs')
    require(primary['continuum_algorithm']['decimal_precisions']==[50,70],'Continuum precision schedule differs')
    require(modes['status']=='failed' and [(r['source'],r['setting']) for r in modes['runs']]==[('positive_B','coarse'),('positive_B','fine')],'Independent completed-run coverage differs')
    require(all(len(r['rows'])==6 and sum(len(p['finite_k']) for p in r['rows'])==18 for r in modes['runs']),'Independent observation coverage differs')
    require('rows' not in modes and modes['active_run'] is None,'Unexpected completed combination or active run')
    require(modes['resources']['peak_rss_kib']==306088 and modes['configuration']['rss_budget_kib']==262144,'Memory resource capture differs')
    require(modes['exception']['type']=='RuntimeError' and modes['exception']['message']=='Frozen 256-MiB peak-memory budget exceeded','Failure reason differs')
    remaining=[r['check'] for r in plan[len(records):]]
    require(remaining==['validation','validation_optimized','summary','figures'],'Unattempted stage coverage differs')
    require(not any((root/'fresh/independent').glob('*signed_uB*')),'Unexpected signed-pulse mode archive')
    nnew=verify_frozen(root);nold=verify_frozen(old);require(nnew==331 and nold==265,'Frozen input counts differ')
    oldstate=read(old/'VALIDATION.json');require(oldstate['status']=='FAIL','Original must remain failed')
    for name,pin in state['fresh_output_sha256'].items():require(sha(root/'fresh'/name)==pin,'Fresh output differs: '+name)
    archives=[]
    for run in modes['runs']:
        q=root/'fresh/independent'/run['archive']['path'];require(sha(q)==run['archive']['sha256'],'Saved raw archive differs')
        archives.append({'path':str(q.relative_to(root)),'sha256':sha(q),'bytes':q.stat().st_size,'source':run['source'],'setting':run['setting'],'observation_rows':len(run['rows']),'finite_K_points':18})
    evidence=['VALIDATION.json','EXECUTION.json','COMMAND_PLAN.json','INPUT_MANIFEST.json','fresh/primary/results.json','fresh/primary/EXECUTION.json','fresh/primary/started.json','fresh/independent/results.json','logs/primary.exit.json','logs/independent_modes.stderr.log','logs/independent_modes.stdout.log','logs/independent_modes.exit.json','checkpoint/FULL_REGISTRATION.json','checkpoint/EXPERIMENT.json','checkpoint/FREEZE_RECEIPT.json','checkpoint/independent/forced_metric.py']
    result={'audit_status':'PASS_READ_ONLY_FAILURE_CLASSIFICATION','followup_registered_calibration_status':'FAIL','original_registered_calibration_status':'FAIL','recorded_utc':datetime.now(timezone.utc).isoformat(),'new_physical_evaluations':0,'followup_public_freeze_commit':receipt['public_freeze_commit'],'followup_registration_sha256':receipt['registration_sha256'],'followup_independent_manifest_sha256':receipt['independent_manifest_sha256'],'original_freeze_receipt':read(old/'checkpoint/FREEZE_RECEIPT.json'),'planned_commands':32,'attempted_commands':28,'successful_commands':27,'pure_commands_completed_pass':26,'physical_producers_attempted':2,'physical_producers_completed':1,'primary_completed_rows':12,'primary_completed_finite_K_points':36,'primary_producer_elapsed_seconds':primary['elapsed_seconds'],'primary_process_elapsed_seconds':records[26]['elapsed_seconds'],'independent_complete_runs':[{'source':r['source'],'setting':r['setting'],'observation_rows':6,'finite_K_points':18} for r in modes['runs']],'independent_combined_rows_produced':False,'independent_combined_refinement_and_Ward_gates_executed':False,'independent_signed_uB_physical_modes_executed':False,'independent_producer_elapsed_seconds':modes['resources']['elapsed_seconds'],'independent_process_elapsed_seconds':records[27]['elapsed_seconds'],'peak_rss_kib':306088,'frozen_peak_rss_budget_kib':262144,'excess_peak_rss_kib':43944,'failure_type':modes['exception']['type'],'failure_message':modes['exception']['message'],'independent_exit_code':1,'timed_out':False,'not_attempted_commands':remaining,'followup_frozen_files_reverified':nnew,'original_frozen_files_reverified':nold,'fresh_output_hashes_reverified':len(state['fresh_output_sha256']),'archives':archives,'allocation_cause_measured':False,'source_based_memory_contributor':'sha(path) reads the entire archive into bytes while evolve retains the archive dictionary; called after its post-write resource check and before the main resource check that failed. Capture does not isolate allocation causality.','cross_validators_executed':False,'physical_instability_established':False,'coupled_stability_or_heating_established':False,'evidence_files':[{'path':x,'sha256':sha(root/x),'bytes':(root/x).stat().st_size} for x in evidence],'audit_source_sha256':sha(Path(__file__))}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('audit_status','followup_registered_calibration_status','planned_commands','attempted_commands','pure_commands_completed_pass','physical_producers_completed','followup_frozen_files_reverified','original_frozen_files_reverified','fresh_output_hashes_reverified')}))

if __name__=='__main__':main()
