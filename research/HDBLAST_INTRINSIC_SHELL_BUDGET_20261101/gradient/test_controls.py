#!/usr/bin/env python3
"""Strict isolated subprocess negatives and separate relocation positives."""
import argparse,copy,hashlib,json,shutil,subprocess,sys
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
ESSENTIAL=['producer.py','audit_aggregation.py','SOURCE_PINS.json','FREEZE.json','FROZEN_REFERENCE.json','PILOT_ONLY.json']
class Failure(Exception):pass
def need(b,s):
    if not b:raise Failure(s)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x):return (json.dumps(x,sort_keys=True,indent=2)+'\n').encode()
def command(script,optimized,args):return [sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])+[str(script),*map(str,args)]
def run(script,optimized,args,expected_witness,positive):
    p=subprocess.run(command(script,optimized,args),capture_output=True,text=True,timeout=60)
    lines=p.stdout.splitlines()
    need(len(lines)==1 and not p.stderr,'EXACT_SINGLE_FINAL_WITNESS')
    final=json.loads(lines[0]);need(final.get('status')==('PASS' if positive else 'FAIL') and final.get('witness')==expected_witness,'INTENDED_FINAL_WITNESS')
    need(p.returncode==(0 if positive else 2),'EXPECTED_RETURN_CODE')
    if not positive:need(final=={'status':'FAIL','witness':expected_witness},'EXACT_FAILURE_OBJECT')
    return final
def scale_pair(value,factor):return [str(F(x)*factor) for x in value]
def mutate(name,base):
    r=copy.deepcopy(base);row=r['rows'][0]
    if name=='wrong_second_derivative_factor':
        for key in ['mean_K2','mean_S2','mean_Q2','direct_mean_Q2']:row[key]=scale_pair(row[key],F(1,2))
    elif name=='wrong_angular_eigenvalue':row['mean_S2']=scale_pair(row['mean_S2'],F(7,8))
    elif name=='bulk_radius_substitution':row['scaled_r']=['1/8','1/8']
    elif name=='wrong_ell_scale':row['actual_background']['ell']='1/1000'
    elif name=='wrong_background_metric':row['actual_background']['rho']=['1','1']
    elif name=='metric_instead_of_scalar_jet':
        row['actual_scalar']['s']=row['actual_metric']['L'];row['actual_scalar']['st']=row['actual_metric']['Lt']
    elif name=='corrupted_raw_metric_data':row['actual_metric']['A']=['0','0']
    elif name=='report_integrity_extra_field':r['unregistered_control_field']=True
    else:raise Failure('UNKNOWN_MUTATION')
    return r
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output-dir',required=True,type=Path);a=ap.parse_args()
    need(sys.version_info[:3]==(3,12,14),'PYTHON_VERSION')
    repo=a.repo_root.resolve(strict=True);out=a.output_dir.absolute()
    need(not out.exists() and all(not p.is_symlink() for p in [out,*out.parents]),'FRESH_OUTPUT')
    need(all(x!=out and x not in out.parents for x in [HERE,repo]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
    out.mkdir(parents=True,exist_ok=False);records=[];baseline=None;report_bytes=[];audit_bytes=[]
    # Positive controls are kept distinct from the rejection cases.
    relocated=out/'relocated';code=relocated/'component';relorepo=relocated/'source';code.mkdir(parents=True)
    for name in ESSENTIAL:shutil.copy2(HERE/name,code/name)
    pins=json.loads((HERE/'SOURCE_PINS.json').read_bytes())
    for label in pins['capsules']:
        destination=relorepo/'research'/label;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copytree(repo/'research'/label,destination)
    for label,component,source in [('original',HERE,repo),('relocated',code,relorepo)]:
        for optimized in [False,True]:
            mode='optimized' if optimized else 'normal';tag=label+'_'+mode
            target=out/(tag+'_producer');audit=out/(tag+'_audit')
            produced=run(component/'producer.py',optimized,['--repo-root',source,'--output-dir',target],'FROZEN_NINE_CALL_CAMPAIGN',True)
            report=target/'REPORT.json';report_bytes.append(report.read_bytes())
            checked=run(component/'audit_aggregation.py',optimized,['--repo-root',source,'--report',report,'--output-dir',audit],'INDEPENDENT_EXACT_RATIONAL_AGGREGATION',True)
            audit_bytes.append((audit/'AUDIT.json').read_bytes())
            records.append({'kind':'positive_relocation_or_baseline','location':label,'mode':mode,'producer_final_witness':produced,'auditor_final_witness':checked})
            if baseline is None:baseline=json.loads(report.read_bytes())
    need(all(x==report_bytes[0] for x in report_bytes) and all(x==audit_bytes[0] for x in audit_bytes),'BYTE_IDENTICAL_NORMAL_OPTIMIZED_AND_RELOCATION')
    cases=[('wrong_second_derivative_factor','SECOND_DERIVATIVE_FACTOR'),('wrong_angular_eigenvalue','ANGULAR_GRADIENT_FACTOR'),
           ('bulk_radius_substitution','PHYSICAL_SHELL'),('wrong_ell_scale','ELL_UNITS'),('wrong_background_metric','BACKGROUND_JET_REFERENCE'),
           ('metric_instead_of_scalar_jet','SCALAR_JET_REFERENCE'),('corrupted_raw_metric_data','METRIC_JET_REFERENCE'),
           ('report_integrity_extra_field','REPORT_INTEGRITY')]
    negative_root=out/'negative_reports';negative_root.mkdir()
    for name,witness in cases:
        report=negative_root/(name+'.json');report.write_bytes(canonical(mutate(name,baseline)))
        for optimized in [False,True]:
            mode='optimized' if optimized else 'normal';target=out/('negative_'+name+'_'+mode)
            final=run(HERE/'audit_aggregation.py',optimized,['--repo-root',repo,'--report',report,'--output-dir',target],witness,False)
            need(not target.exists(),'FAILURE_DID_NOT_WRITE_OUTPUT')
            records.append({'kind':'negative','case':name,'mode':mode,'final_witness':final,'output_directory_absent':True})
    # Mutate only a disposable relocated PUBLIC source capsule, then reject
    # before any scientific import. The source checkout itself stays untouched.
    altered=relorepo/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py'
    altered.write_bytes(altered.read_bytes()+b'\n# deliberate disposable source-control mutation\n')
    for optimized in [False,True]:
        mode='optimized' if optimized else 'normal';target=out/('negative_source_capsule_'+mode)
        final=run(HERE/'producer.py',optimized,['--repo-root',relorepo,'--output-dir',target],'SOURCE_MEMBER_BYTES',False)
        need(not target.exists(),'FAILURE_DID_NOT_WRITE_OUTPUT')
        records.append({'kind':'negative','case':'source_capsule_member_corruption','mode':mode,'final_witness':final,'output_directory_absent':True})
    result={'schema':'HDBLAST_PORTABLE_SHELL_GRADIENT_STRICT_SUBPROCESS_CONTROLS_V1','status':'PASS_STRICT_SUBPROCESS_CONTROLS',
        'test_source_sha256':sha(Path(__file__)),'python_version':'3.12.14','isolation_flags':['-I','-S','-B'],
        'optimization_modes':['normal','optimized'],'negative_cases_per_mode':9,'negative_subprocesses':18,
        'positive_producer_subprocesses':4,'positive_auditor_subprocesses':4,'positive_relocation_separate_from_negatives':True,
        'report_byte_identical_across_modes_and_relocation':True,'audit_byte_identical_across_modes_and_relocation':True,
        'producer_report_sha256':hashlib.sha256(report_bytes[0]).hexdigest(),'audit_report_sha256':hashlib.sha256(audit_bytes[0]).hexdigest(),
        'every_negative_checks_exact_final_object_return_code_empty_stderr_and_absent_output':True,'records':records}
    target=out/'TEST_REPORT.json';target.write_bytes(canonical(result))
    print(json.dumps({'status':'PASS','witness':'STRICT_SUBPROCESS_AND_RELOCATION_CONTROLS','test_report_sha256':sha(target)},sort_keys=True))
if __name__=='__main__':
    try:main()
    except Failure as exc:
        print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True));sys.exit(2)
