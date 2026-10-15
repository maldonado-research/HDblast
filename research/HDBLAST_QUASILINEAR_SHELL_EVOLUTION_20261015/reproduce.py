"""Authenticated fresh replay of exact finite checks, not computer certification of PDE proofs."""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

CHECKERS = [
    {'name':'energy','script':'check_variable_energy.py','golden':'VARIABLE_ENERGY_CHECKS_NORMAL.json',
     'controls':[('--corrupt','cross_factor','variable_scalar_energy_divergence_identity'),
                 ('--corrupt','energy_sign','variable_scalar_energy_divergence_identity'),
                 ('--corrupt','boundary_shift','four_D_twelve_N_zero_principal_boundary_flux'),
                 ('--corrupt','mixed_flux','four_D_twelve_N_zero_principal_boundary_flux')]},
    {'name':'astra','script':'astra/check_variable_energy.py','golden':'astra/CHECK_NORMAL_none.json',
     'controls':[('--mutation','reverse_flux','local_energy_identity'),
                 ('--mutation','drop_shift','local_energy_identity'),
                 ('--mutation','mix_DN','homogeneous_split_flux_zero')]},
    {'name':'gauge','script':'gauge-parameter/verify_gauge_subsidiary.py','golden':'gauge-parameter/GAUGE_SUBSIDIARY_CHECK_RESULTS.json','controls':[]},
]

def record(path, relative_to=None):
    raw=path.read_bytes()
    return {'path':path.relative_to(relative_to).as_posix() if relative_to else str(path),'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}

def require(condition,message):
    if not condition:
        raise RuntimeError(message)

def scientific(name,data):
    return {k:v for k,v in data.items() if k!='optimized_python'} if name=='gauge' else data

def execute(root,output,checker,mode,control=None):
    work=output/mode/checker['name']/(control[1] if control else 'positive')
    work.mkdir(parents=True)
    source=root/checker['script']; script=work/source.name
    shutil.copy2(source,script)
    require(source.read_bytes()==script.read_bytes(),'Fresh executable bytes differ')
    flags=['-I','-B']+(['-O'] if mode=='optimized' else [])
    arguments=['fresh.json'] if checker['name']=='gauge' else ['--output','fresh.json']
    if control:
        arguments+=list(control[:2])
    started=time.monotonic()
    with (work/'execution.log').open('w') as log:
        completed=subprocess.run([sys.executable,*flags,str(script),*arguments],cwd=work,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    row={'kind':'external_corruption' if control else 'positive','checker':checker['name'],'mode':mode,'script':checker['script'],'source_sha256':record(source)['sha256'],'flags':flags,'arguments':arguments,'exit_code':completed.returncode,'wall_seconds':round(time.monotonic()-started,3),'copied_script':record(script,output),'log':record(work/'execution.log',output)}
    result=work/'fresh.json'
    if result.is_file():
        data=json.loads(result.read_text()); row['output']=record(result,output)
    if control:
        require(completed.returncode==1,'Corruption did not reject with exit code one')
        diagnostic=(work/'execution.log').read_text()
        if checker['name']=='astra':
            require(result.is_file() and data['passed'] is False and data['checks'][control[2]] is False,'Astra mutation rejected for wrong reason')
        else:
            require(control[2] in diagnostic and not result.exists(),'Energy mutation rejected for wrong reason')
        row.update({'control':control[1],'expected_reason':control[2],'expected_reason_detected':True})
    else:
        require(completed.returncode==0 and result.is_file(),'Fresh checker failed')
        require(data.get('status')=='PASS' or data.get('passed') is True,'Checker did not report success')
        golden=root/checker['golden']; expected=json.loads(golden.read_text())
        require(scientific(checker['name'],data)==scientific(checker['name'],expected),'Scientific output differs from frozen reference')
        if checker['name']=='gauge':
            require(data['optimized_python'] is (mode=='optimized'),'Python optimization mode mismatch')
            require(data['rejected_count']==6 and len(data['wrong_formula_controls_rejected'])==6,'Missing embedded gauge controls')
            row['embedded_wrong_formula_rejections']=6
        row.update({'golden':record(golden,root),'scientific_output_matches':True,'ignored_top_level_fields':['optimized_python'] if checker['name']=='gauge' else []})
    return row

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256',required=True)
    parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    spec=importlib.util.spec_from_file_location('quasilinear_integrity',root/'verify_manifest.py')
    verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)
    auth=verifier.verify(root,args.manifest_sha256)
    requested=args.output_dir.absolute()
    require(not requested.is_symlink(),'Real output directory required')
    output=requested.resolve()
    require(not output.is_relative_to(root),'Output must be outside frozen capsule')
    require(not output.exists() or (output.is_dir() and not any(output.iterdir())),'New or empty output directory required')
    require(importlib.metadata.version('sympy')=='1.14.0','Pinned sympy 1.14.0 required')
    output.mkdir(parents=True,exist_ok=True)
    rows=[]
    for mode in ('normal','optimized'):
        for checker in CHECKERS:
            rows.append(execute(root,output,checker,mode))
        for checker in CHECKERS:
            for control in checker['controls']:
                rows.append(execute(root,output,checker,mode,control))
    positives=sum(row['kind']=='positive' for row in rows)
    external=sum(row['kind']=='external_corruption' for row in rows)
    embedded=sum(row.get('embedded_wrong_formula_rejections',0) for row in rows)
    require((positives,external,embedded)==(6,14,12),'Incomplete replay')
    receipt={'schema':'hdblast-quasilinear-shell-fresh-replay-v1','status':'PASS','capsule':'HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015','authentication':auth,'runner':record(Path(__file__).resolve(),root),'integrity_verifier':record(root/'verify_manifest.py',root),'requirements':record(root/'requirements.txt',root),'python':{'version':sys.version,'executable':sys.executable},'sympy_version':importlib.metadata.version('sympy'),'successful_checker_runs':positives,'external_corruption_processes_rejected':external,'embedded_gauge_wrong_formula_controls_rejected':embedded,'finite_positive_check_groups_each_mode':19,'scope':'Exact finite algebra and index checks only. Analytic existence and parameter proofs are separately written and internally reviewed; this replay is not a proof-assistant certificate or external peer review. External processes and embedded controls are reported separately.','steps':rows}
    (output/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','successful_checker_runs':positives,'external_corruption_processes_rejected':external,'embedded_gauge_controls_rejected':embedded,'output':str(output)},indent=2))

if __name__=='__main__':
    main()
