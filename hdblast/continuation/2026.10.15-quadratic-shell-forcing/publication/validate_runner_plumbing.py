"""Synthetic replay tests only; never execute HDBLAST scientific code."""
from pathlib import Path
import hashlib,importlib.util,json,shutil,subprocess,sys,tempfile,copy
HERE=Path(__file__).resolve().parent
CAP='HDBLAST_QUADRATIC_SHELL_FORCING_20261015'
def require(ok,reason):
    if not ok:raise RuntimeError(reason)
def run(arguments,expected=0):
    p=subprocess.run([sys.executable,'-I','-B',*map(str,arguments)],capture_output=True,text=True,timeout=60)
    require(p.returncode==expected,'Unexpected subprocess exit: '+str(arguments)+'\n'+p.stdout+p.stderr)
    return p
fixture=Path(tempfile.mkdtemp(prefix='runner-synthetic-',dir=HERE));source=fixture/'source';source.mkdir()
checker='''import argparse,json,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--fault',default='none');a=p.parse_args()
reason=None
if a.fault=='sign':reason='synthetic wrong sign'
if Path('coefficient.txt').read_text()!='1':reason='synthetic coefficient bytes'
if json.loads(Path('input.json').read_text())['coefficient']!='1':reason='synthetic coefficient value'
if a.fault=='traceback':
 print('REJECTED: synthetic named exception',file=sys.stderr)
 raise ValueError('synthetic named exception')
if a.fault=='json':print(json.dumps({'status':'FAIL','failure':'synthetic named JSON failure'}));sys.exit(1)
if reason:print('FAIL: '+reason);sys.exit(1)
print(json.dumps({'status':'PASS','checks':['alpha','beta'],'count':2,'embedded':{'wrong_alpha':True},'embedded_count':1,'runtime':{'optimization':sys.flags.optimize},'value':1},indent=2)+'\\n',end='')
'''
generator='''import json
from pathlib import Path
p=Path('generated');p.mkdir();(p/'output.json').write_text(json.dumps({'status':'PASS','cells':2},indent=2)+'\\n');(p/'output.csv').write_text('a,b\\n1,2\\n')
'''
(source/'checker.py').write_text(checker);(source/'generator.py').write_text(generator);(source/'coefficient.txt').write_text('1');(source/'input.json').write_text('{"coefficient":"1"}\n')
for mode,optimization in [('normal',0),('optimized',1)]:
 (source/(mode+'.json')).write_text(json.dumps({'status':'PASS','checks':['alpha','beta'],'count':2,'embedded':{'wrong_alpha':True},'embedded_count':1,'runtime':{'optimization':optimization},'value':1},indent=2)+'\n')
(source/'generated.json').write_text(json.dumps({'status':'PASS','cells':2},indent=2)+'\n');(source/'generated.csv').write_text('a,b\n1,2\n')
config={'schema':'hdblast-build-time-replay-config-v1','status':'APPROVED_FOR_FROZEN_REPLAY','capsule':CAP,'python_version':sys.version.split()[0],'dependency_versions':{},'checkers':[
 {'name':'synthetic','script':'checker.py','dependencies':['coefficient.txt','input.json'],'golden':'normal.json','optimized_golden':'optimized.json','status':'PASS','checks':{'path':['checks'],'kind':'names'},'check_count_path':['count'],'embedded_checks':{'path':['embedded'],'kind':'true_map'},'embedded_count_path':['embedded_count'],'mode_fields':[{'path':['runtime','optimization'],'normal':0,'optimized':1}], 'required_fields':[{'path':['value'],'equals':1}]},
 {'name':'generator','script':'generator.py','golden':'generated.json','output':'generated/output.json','additional_outputs':[['generated/output.csv','generated.csv']],'status':'PASS','byte_equal':True,'required_fields':[{'path':['cells'],'equals':2}]}], 'controls':[]}
for name,reason,fmt,extra in [
 ('sign','synthetic wrong sign','fail_line',{'arguments':['--fault','sign']}),
 ('json_failure','synthetic named JSON failure','json',{'arguments':['--fault','json'],'failure_path':['failure']}),
 ('exception','synthetic named exception','traceback',{'arguments':['--fault','traceback'],'exception':'ValueError','diagnostic_prefix':'REJECTED: synthetic named exception\n'}),
 ('append_input','synthetic coefficient bytes','fail_line',{'target':'coefficient.txt','append':'x'}),
 ('replace_input','synthetic coefficient bytes','fail_line',{'target':'coefficient.txt','before':'1','after':'-1'}),
 ('replace_json','synthetic coefficient value','fail_line',{'target':'input.json','json_path':['coefficient'],'value':'-1'}),
 ('scale_json','synthetic coefficient value','fail_line',{'target':'input.json','json_path':['coefficient'],'factor':'-1'})]:
 config['controls'].append({'checker':'synthetic','name':name,'reason':reason,'failure_format':fmt,'mechanism':'synthetic copied input mutation' if 'target' in extra else 'synthetic production formula injection',**extra})
configpath=fixture/'CONFIG.json';configpath.write_text(json.dumps(config,indent=2)+'\n')
run([HERE/'build_replay_runner.py','--config',configpath,'--source-root',source,'--output-dir',fixture/'built'])
for name in ['reproduce.py','REPLAY_ROSTER.json']:shutil.copyfile(fixture/'built'/name,source/name)
shutil.copyfile(HERE/'replay-tools/verify_manifest.py',source/'verify_manifest.py');(source/'requirements.txt').write_text('')
run([HERE/'freeze_source.py','--source',source,'--output-dir',fixture/'frozen','--capsule',CAP])
frozen=fixture/'frozen'/CAP;pin=hashlib.sha256((frozen/'MANIFEST.json').read_bytes()).hexdigest()
run([frozen/'reproduce.py','--manifest-sha256',pin,'--output-dir',fixture/'replay'])
receipt=json.loads((fixture/'replay/REPRODUCTION_RECEIPT.json').read_bytes());require(receipt['successful_checker_runs']==4 and receipt['external_corruption_processes_rejected']==14 and receipt['finite_positive_check_groups_each_mode']==2 and receipt['embedded_wrong_formula_controls_rejected']==2,'Synthetic complete process counts')
rejections={}
for name,payloadpin in [('wrong_external_pin','0'*64)]:
 output=fixture/name;result=run([frozen/'reproduce.py','--manifest-sha256',payloadpin,'--output-dir',output],1);require(not output.exists(),'Bad external pin created output');rejections[name]=True
(frozen/'coefficient.txt').write_text('changed');output=fixture/'changed_source';run([frozen/'reproduce.py','--manifest-sha256',pin,'--output-dir',output],1);require(not output.exists(),'Changed source created output');rejections['changed_scientific_source']=True;(frozen/'coefficient.txt').write_text('1')
spec=importlib.util.spec_from_file_location('runner',frozen/'reproduce.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
roster=json.loads((frozen/'REPLAY_ROSTER.json').read_bytes())
for name,change in {
 'same_counts_wrong_checker':lambda x:x['positive_processes'][0].update(name='other'),
 'same_counts_wrong_script':lambda x:x['positive_processes'][0].update(script='other.py'),
 'same_counts_wrong_group':lambda x:x['positive_processes'][0]['finite_positive_check_names'].__setitem__(0,'other'),
 'same_counts_wrong_embedded':lambda x:x['positive_processes'][0]['embedded_control_names'].__setitem__(0,'other'),
 'same_counts_wrong_reason':lambda x:x['external_controls'][0].update(expected_reason='other'),
 'same_counts_wrong_mechanism':lambda x:x['external_controls'][0].update(mechanism='other')}.items():
 candidate=copy.deepcopy(roster);change(candidate)
 try:module.validate_declared_roster(candidate)
 except RuntimeError:rejections[name]=True
 else:raise RuntimeError('Wrong roster accepted: '+name)
base=json.loads((source/'normal.json').read_bytes());checker=module.CHECKERS[0]
for name,change in {'wrong_mode_value':lambda d:d['runtime'].update(optimization=1),'wrong_mode_type':lambda d:d['runtime'].update(optimization=False),'wrong_named_inventory':lambda d:d['checks'].__setitem__(0,'other'),'false_embedded_control':lambda d:d['embedded'].update(wrong_alpha=False),'wrong_declared_count':lambda d:d.update(count=3)}.items():
 candidate=copy.deepcopy(base);change(candidate)
 try:module.positive_inventory(checker,candidate,'normal')
 except RuntimeError:rejections[name]=True
 else:raise RuntimeError('Wrong inventory accepted: '+name)
require(not module.json_equal({'x':True},{'x':1}) and not module.json_equal({'x':1},{'x':1.0}),'Strict JSON scalar types');rejections['boolean_integer_json_alias']=True;rejections['integer_float_json_alias']=True
# A producer-created symlink cannot count as a fresh JSON/CSV file.
external=fixture/'outside.json';external.write_text('{"status":"PASS"}')
symlinkwork=fixture/'symlink-output';symlinkwork.mkdir();(symlinkwork/'result.json').symlink_to(external)
try:module.read_generated(symlinkwork/'result.json',symlinkwork)
except RuntimeError:rejections['producer_symlink_result']=True
else:raise RuntimeError('Producer symlink result accepted')
real=fixture/'external-directory';real.mkdir();(real/'output.csv').write_text('x\n')
(symlinkwork/'nested').symlink_to(real,target_is_directory=True)
try:module.read_generated(symlinkwork/'nested/output.csv',symlinkwork)
except RuntimeError:rejections['producer_symlink_parent_for_additional_output']=True
else:raise RuntimeError('Producer symlink output parent accepted')
strictchecker=copy.deepcopy(checker);strictchecker['required_fields']=[{'path':['runtime'],'equals':{'optimization':False}}]
try:module.positive_inventory(strictchecker,base,'normal')
except RuntimeError:rejections['nested_required_field_boolean_integer_alias']=True
else:raise RuntimeError('Nested required-field type alias accepted')
trace_control=next(c for c in module.CONTROLS if c['name']=='exception')
trace_log=(fixture/'replay/normal/synthetic/exception/execution.log').read_text()
for name,diagnostic in {
    'traceback_unexpected_leading_bytes':'unexpected\n'+trace_log,
    'traceback_unexpected_trailing_bytes':trace_log+'unexpected\n',
    'traceback_trailing_blank_line':trace_log+'\n',
    'traceback_wrong_final_reason':trace_log.rsplit('ValueError: synthetic named exception',1)[0]+'ValueError: unrelated\n',
    'traceback_wrong_prefix_reason':trace_log.replace('REJECTED: synthetic named exception','REJECTED: unrelated',1),
    'traceback_missing_exact_prefix':trace_log.split('\n',1)[1],
}.items():
    try:module.validate_failure_diagnostic(trace_control,diagnostic)
    except RuntimeError:rejections[name]=True
    else:raise RuntimeError('Unrelated traceback contract accepted: '+name)
pins={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ['runner_template.py','build_replay_runner.py','validate_runner_plumbing.py']}
result={'status':'PASS_SYNTHETIC_REPLAY_PLUMBING','helper_sha256':pins,'fixture':str(fixture),'synthetic_positive_processes':4,'synthetic_named_external_failures':14,'synthetic_embedded_controls':2,'synthetic_groups_per_mode':2,'rejected_integrity_and_contract_cases':rejections,'complete_json_csv_generated_in_both_modes':True,'http_requests':0,'git_commands':0,'scientific_HDBLAST_processes_executed':0,'scope':'Synthetic plumbing and corruption tests only. New scientific source and roster have not been accepted or replayed.'}
(HERE/'RUNNER_PLUMBING_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
