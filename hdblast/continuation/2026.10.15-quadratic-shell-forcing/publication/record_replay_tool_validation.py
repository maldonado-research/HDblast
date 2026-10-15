"""Bind completed copied-fixture replay, all artifacts and reviewed helper inputs."""
from pathlib import Path
from hashlib import sha256
import json
from datetime import datetime,timezone
P=Path(__file__).resolve().parent;ROOT=P.parent
FIXTURE=P/'runner-scientific-fixture-final'
CAP='HDBLAST_QUADRATIC_SHELL_FORCING_20261015'

def require(ok,reason):
    if not ok:raise RuntimeError(reason)
def item(path,base=ROOT):
    require(path.is_file() and not path.is_symlink(),'Real evidence file required')
    raw=path.read_bytes();return {'path':path.relative_to(base).as_posix(),'bytes':len(raw),'sha256':sha256(raw).hexdigest()}
def authenticate(path,record):
    actual=item(path, path.parent)
    require(actual['bytes']==record['bytes'] and actual['sha256']==record['sha256'],'Evidence bytes differ: '+str(path))

replay=json.loads((FIXTURE/'replay/REPRODUCTION_RECEIPT.json').read_bytes())
require(replay['status']=='PASS','Complete fixture replay required')
counts={k:replay[k] for k in ['successful_checker_runs','external_corruption_processes_rejected','embedded_wrong_formula_controls_rejected','finite_positive_check_groups_each_mode']}
require(counts=={'successful_checker_runs':12,'external_corruption_processes_rejected':70,'embedded_wrong_formula_controls_rejected':16,'finite_positive_check_groups_each_mode':111},'Reviewed final roster counts')
source=FIXTURE/'frozen'/CAP;manifest=json.loads((source/'MANIFEST.json').read_bytes());require(len(manifest['files'])==33,'Minimal copied scientific fixture source inventory')
for row in manifest['files']:authenticate(source/row['path'],row)
for name in ['reproduce.py','REPLAY_ROSTER.json','verify_manifest.py']:require((source/name).read_bytes()==(P/'replay-tools'/name).read_bytes(),'Accepted tool bytes changed')
build=json.loads((P/'replay-tools/RUNNER_BUILD_RECEIPT.json').read_bytes())
for row in build['source_input_pins']:authenticate(ROOT/row['path'],row)
require(replay['authentication']['manifest_sha256']==sha256((source/'MANIFEST.json').read_bytes()).hexdigest(),'Fixture manifest pin')
require(replay['embedded_configuration_sha256']==sha256((P/'RUNNER_CONFIG.json').read_bytes()).hexdigest(),'Embedded approved config pin')
artifacts={}
for row in replay['steps']:
    for record in [row['log'],row['output']]+[x['output'] for x in row['additional_outputs']]:
        require(record['path'] not in artifacts,'Distinct artifact record required');authenticate(FIXTURE/'replay'/record['path'],record);artifacts[record['path']]=record
    if row['kind']=='external_corruption':
        data=json.loads((FIXTURE/'replay'/row['output']['path']).read_bytes());log=(FIXTURE/'replay'/row['log']['path']).read_text()
        require(data['status']=='FAIL' and data['failure']==row['expected_reason']==row['reported_failure'] and data['complete_diagnostic']==log and data['exit_code']==row['exit_code']==1 and data['passing_output_absent'] is True,'Named failure artifact consistency')
require(len(artifacts)==166 and len(replay['steps'])==82,'Complete retained artifact inventory')
synthetic=json.loads((P/'RUNNER_PLUMBING_VALIDATION.json').read_bytes())
for name,pin in synthetic['helper_sha256'].items():require(sha256((P/name).read_bytes()).hexdigest()==pin,'Tested helper source pin')
require(len(synthetic['rejected_integrity_and_contract_cases'])==24 and all(synthetic['rejected_integrity_and_contract_cases'].values()),'Complete synthetic hardening inventory')
receipt={'status':'PASS','observed_utc':datetime.now(timezone.utc).isoformat(),'capsule':CAP,'runner':item(P/'replay-tools/reproduce.py'),'roster':item(P/'replay-tools/REPLAY_ROSTER.json'),'verifier':item(P/'replay-tools/verify_manifest.py'),'approved_build_config':item(P/'RUNNER_CONFIG.json'),'build_receipt':item(P/'replay-tools/RUNNER_BUILD_RECEIPT.json'),'synthetic_plumbing_validation':item(P/'RUNNER_PLUMBING_VALIDATION.json'),'independent_tooling_review':item(ROOT/'referee/RUNNER_PLUMBING_REVIEW.json'),'fixture_freeze':item(FIXTURE/'frozen/SOURCE_FREEZE_RECEIPT.json'),'fixture_replay':item(FIXTURE/'replay/REPRODUCTION_RECEIPT.json'),'fixture_source_files':33,'fixture_archive_members':34,'fresh_processes':82,'fresh_artifact_files_excluding_receipt':166,'fresh_artifact_files_including_receipt':167,'replay_counts':counts,'source_input_pins_rechecked_unchanged':len(build['source_input_pins']),'artifact_records':[artifacts[n] for n in sorted(artifacts)],'historical_diagnostic_contract_fix':'The earlier copied fixture stopped at the first archive fault because its exact REJECTED prefix was not yet configured. The revised explicit prefix contract, accepted independently, rejects unexpected leading/trailing bytes and wrong reasons; this completed fixture uses the revised runner. No scientific source was changed.','scope':'Complete minimal copied-source fixture replay and publication-tool validation. This is distinct from the root final full-source frozen-ZIP replay and does not certify remote CI, merge, deployment, a full U2 response or analytic PDE premises.'}
(P/'REPLAY_TOOL_VALIDATION.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':'PASS','receipt_sha256':sha256((P/'REPLAY_TOOL_VALIDATION.json').read_bytes()).hexdigest(),'counts':counts,'artifacts':len(artifacts)},indent=2))
