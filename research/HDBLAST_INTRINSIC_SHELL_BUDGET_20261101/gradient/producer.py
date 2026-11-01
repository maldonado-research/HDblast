#!/usr/bin/env python3
"""Portable frozen nine-call actual shell-gradient campaign; Python 3.12.14."""
import argparse
from fractions import Fraction as F
import hashlib, importlib.util, json, sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
PINS_SHA='9de076a120a0eb5b84afc45e626f145f129a0c1bddb87a1a0921ae9c3441155e'
FREEZE_SHA='67f33c24285cf252b7dde4f40da5b71b1fafde4c41bc1836effba7b6541ed3bc'
REFERENCE_SHA='71d200264d7ba2bc7cfed87e8c1abd6978ab8cd92571c10e1e34a55670cef3f5'
class Failure(Exception):pass
def need(ok,witness):
    if not ok:raise Failure(witness)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def regular(p):return p.is_file() and not p.is_symlink() and all(not x.is_symlink() for x in p.absolute().parents)
def read(p):return json.loads(p.read_bytes())
def canonical(x):return (json.dumps(x,sort_keys=True,indent=2)+'\n').encode()
def encode(x):
    if hasattr(x,'lo') and hasattr(x,'hi'):return [str(x.lo),str(x.hi)]
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    return x
def authenticate(repo,pins):
    need(repo.is_dir() and not repo.is_symlink(),'SOURCE_ROOT')
    answer={}
    for label,pin in sorted(pins['capsules'].items()):
        root=repo/'research'/label;manifest=root/'MANIFEST.json'
        need(root.is_dir() and not root.is_symlink(),'SOURCE_CAPSULE_ROOT')
        need(regular(manifest) and sha(manifest)==pin,'SOURCE_MANIFEST_BYTES')
        seen={'MANIFEST.json'}
        for row in read(manifest)['files']:
            rel=Path(row['path']);need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'SOURCE_SAFE_MEMBER')
            q=root/rel
            need(regular(q) and q.stat().st_size==row['bytes'] and sha(q)==row['sha256'],'SOURCE_MEMBER_BYTES')
            need(pins['all_manifest_member_pins'][label+'/'+str(rel)]=={'bytes':row['bytes'],'sha256':row['sha256']},'SOURCE_FLAT_PIN')
            seen.add(str(rel))
        actual=set()
        for q in root.rglob('*'):
            need(not q.is_symlink(),'SOURCE_SYMLINK')
            if q.is_file():actual.add(str(q.relative_to(root)))
            else:need(q.is_dir(),'SOURCE_SPECIAL_FILE')
        need(actual==seen,'SOURCE_CLOSED_INVENTORY')
        answer[label]={'manifest_sha256':pin,'members_including_manifest':len(seen)}
    return answer
def outward(v):return [str(v.lo.__floor__()),str(-(-v.hi).__floor__())]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
    need(sys.version_info[:3]==(3,12,14),'PYTHON_VERSION')
    repo=a.repo_root.absolute();out=a.output_dir.absolute()
    need(all(not p.is_symlink() for p in [repo,*repo.parents,out,*out.parents]),'PATH_SYMLINK')
    repo=repo.resolve(strict=True);out=out.resolve()
    need(not out.exists(),'FRESH_OUTPUT_DIRECTORY')
    need(all(base!=out and base not in out.parents for base in [repo,HERE]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
    own=sha(Path(__file__));local={}
    for name,pin in [('SOURCE_PINS.json',PINS_SHA),('FREEZE.json',FREEZE_SHA),('FROZEN_REFERENCE.json',REFERENCE_SHA)]:
        p=HERE/name;need(regular(p) and sha(p)==pin,'LOCAL_'+name);local[name]=pin
    pins,freeze,reference=map(read,[HERE/'SOURCE_PINS.json',HERE/'FREEZE.json',HERE/'FROZEN_REFERENCE.json'])
    pilot=HERE/'PILOT_ONLY.json';need(regular(pilot) and sha(pilot)==freeze['exploratory_pilot_sha256'],'SEPARATED_PILOT')
    need(freeze['source_pins_sha256']==PINS_SHA and freeze['time_cells']==8 and freeze['actual_oracle_calls']==9
         and freeze['selected_target']=={'scaled_x':'5/4','scaled_r':'0'}
         and freeze['continuous_postcontact_interval']==['31/25','63/50'] and freeze['no_adaptive_refinement'] is True,'FROZEN_CAMPAIGN')
    before=authenticate(repo,pins)
    # All scientific dependencies are in these authenticated closed capsules.
    path=repo/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py'
    spec=importlib.util.spec_from_file_location('authenticated_actual_shell_gradient_oracle',path)
    jets=importlib.util.module_from_spec(spec);spec.loader.exec_module(jets)
    I=jets.I;need(jets.ELL==F(1,10**6) and jets.core.GRID==1<<112,'ARITHMETIC_UNITS')
    _,errors,errorpins=jets.inputs()
    need(errors['status']=='PASS_ACTUAL_MASTER_AND_GAUGE_CONTINUOUS_SCALED_C3_APPROXIMATION_ERRORS','ACTUAL_ERROR_STATUS')
    lo,hi=F(31,25),F(63,50);N=8
    need(1<lo<hi<F(3,2),'POSTCONTACT_EXACT_FAMILY_DOMAIN')
    specs=[('continuous_cell',[lo+(hi-lo)*i/N,lo+(hi-lo)*(i+1)/N]) for i in range(N)]+[('selected_target',[F(5,4),F(5,4)])]
    rows=[]
    for kind,xb in specs:
        m,s,b=jets.actual_jet(I(*xb),I(0),mode='actual')
        K2=s['st']**2/(3*b['rho']**2)
        S2=8*s['s']**2/(3*b['rho']**2*b['c']**2)
        Q2=S2-K2;ratio=S2/K2
        direct=(-s['st']**2+8*s['s']**2/b['c']**2)/(3*b['rho']**2)
        need(K2.lo>0 and S2.lo>0 and Q2.hi<0 and ratio.hi<1 and direct.hi<0,'STRICT_SCIENTIFIC_SIGN')
        rows.append({'kind':kind,'scaled_x':xb,'scaled_r':[F(0),F(0)],'actual_metric':m,'actual_scalar':s,'actual_background':b,
                     'mean_K2':K2,'mean_S2':S2,'mean_Q2':Q2,'direct_mean_Q2':direct,'spatial_to_temporal_ratio':ratio,'outward_integer_Q2':outward(Q2)})
    cells=rows[:N]
    need(cells[0]['scaled_x'][0]==lo and cells[-1]['scaled_x'][1]==hi and all(cells[i]['scaled_x'][1]==cells[i+1]['scaled_x'][0] for i in range(N-1)),'CONTINUOUS_COVER')
    summary={key:I(min(row[key].lo for row in cells),max(row[key].hi for row in cells)) for key in ['mean_K2','mean_S2','mean_Q2','spatial_to_temporal_ratio']}
    scientific=encode({'rows':rows,'selected_target':rows[-1],'continuous_interval_hulls':summary,
        'continuous_outward_integer_Q2':outward(summary['mean_Q2']),'actual_oracle_calls':len(rows),
        'first_derivative_zero':True,'derivative_convention':reference['derivative_convention'],'scope':reference['scope']})
    need(all(scientific[k]==reference[k] for k in reference['scientific_reference_keys']),'UNCHANGED_FROZEN_SCIENTIFIC_VALUES')
    after=authenticate(repo,pins);need(after==before,'SOURCE_POSTCHECK')
    need(sha(Path(__file__))==own and all(sha(HERE/name)==pin for name,pin in local.items()) and sha(pilot)==freeze['exploratory_pilot_sha256'],'LOCAL_POSTCHECK')
    result={**scientific,'schema':'HDBLAST_PORTABLE_INTRINSIC_SHELL_GRADIENT_SECOND_COEFFICIENT_V1',
        'status':'PASS_CONDITIONAL_CONTINUOUS_POSTCONTACT_NEGATIVE_MEAN_SHELL_GRADIENT_Q2',
        'producer_sha256':own,'source_pins_sha256':PINS_SHA,'freeze_sha256':FREEZE_SHA,'scientific_reference_sha256':REFERENCE_SHA,
        'original_final_report_sha256':reference['original_final_report_sha256'],'unchanged_frozen_scientific_values':True,
        'source_authentication':before,'all_sources_authenticated_before_import_and_after':True,'frozen_campaign':freeze,
        'python_version':'3.12.14','arithmetic_grid_denominator':str(1<<112),
        'actual_error_input_pins':{k:{'bytes':v['bytes'],'sha256':v['sha256']} for k,v in sorted(errorpins.items())}}
    out.mkdir(parents=True,exist_ok=False);report=out/'REPORT.json';report.write_bytes(canonical(result))
    print(json.dumps({'status':'PASS','witness':'FROZEN_NINE_CALL_CAMPAIGN','report_sha256':sha(report)},sort_keys=True))
if __name__=='__main__':
    try:main()
    except Failure as exc:
        print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True));sys.exit(2)
