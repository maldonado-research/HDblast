#!/usr/bin/env python3
"""Check narrow path-only review changes and the final portable bindings."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
ORIGINAL_RECEIPT='2a3f7029ecb6e8661cb9f07750a4ffac39c3a710f6112a354657ea3d556e9076'
OLD_GEOMETRY_CHECKER='e6dd2175df2cdc03cb3ab842b74177a45ccd1cb3f15e8b048422784fd861feef'
PATCHED=('review/check_frame_against_coordinate.py','review/check_producer_against_coordinate.py')
SMOKES=('independent_coordinate_oracle.py','check_producer_against_coordinate.py',
        'check_frame_against_coordinate.py','check_independent_source_intervals.py',
        'check_independent_shell_zeros.py','check_portable_review.py')

def need(ok,label):
    if not ok:raise ValueError(label)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def read(name):return json.loads((ROOT/name).read_text())
def pin(name):
    raw=(ROOT/name).read_bytes();return {'sha256':sha(raw),'bytes':len(raw)}

def build():
    original_raw=(ROOT/'review/REVIEW_RECEIPT.json').read_bytes()
    need(sha(original_raw)==ORIGINAL_RECEIPT,'PRESERVED_ORIGINAL_REVIEW_RECEIPT')
    original=json.loads(original_raw);current={};patches={}
    needle=b'sys.path.insert(0,str(Path(__file__).resolve().parent))\n'
    for name,old in original['accepted_files'].items():
        raw=(ROOT/name).read_bytes();current[name]=pin(name)
        if name in PATCHED:
            need(raw.count(needle)==1,'EXACT_ONE_REVIEW_BOOTSTRAP_'+name)
            need(sha(raw.replace(needle,b'',1))==old['sha256'],'ONLY_REVIEW_BOOTSTRAP_CHANGED_'+name)
            patches[name]={'before_sha256':old['sha256'],'after_sha256':sha(raw),'inverse_patch_matches_old_bytes':True}
        else:need(current[name]==old,'UNCHANGED_SCIENTIFIC_INPUT_'+name)
    geom='geometry/check_geometric_sources.py';raw=(ROOT/geom).read_bytes()
    before=raw.replace(b'sys.path.insert(0,str(ROOT))\n',b'',1)
    before=before.replace(b"[sys.executable,'-I','-S','-B']",b'[sys.executable]',1)
    need(sha(before)==OLD_GEOMETRY_CHECKER,'ONLY_GEOMETRY_IMPORT_AND_CHILD_FLAGS_CHANGED')
    patches[geom]={'before_sha256':OLD_GEOMETRY_CHECKER,'after_sha256':sha(raw),'inverse_patch_matches_old_bytes':True}
    freeze=read('geometry/GEOMETRY_FREEZE.json')
    for name,expected in freeze['files'].items():need(pin('geometry/'+name)==expected,'FINAL_GEOMETRY_FREEZE_'+name)
    for name,expected in freeze['dependencies'].items():need(pin('geometry/'+name)==expected,'FINAL_GEOMETRY_DEPENDENCY_'+name)
    need(freeze['normal_optimized_equal'],'GEOMETRY_NORMAL_OPTIMIZED_PREMISE')
    for mode in ('NORMAL','OPTIMIZED'):
        report=read('geometry/GEOMETRY_'+mode+'.json')
        need(report['status']=='pass' and len(report['controls'])==6,'GEOMETRY_SIX_CONTROLS_'+mode)
        need(all(row['rejected'] for row in report['controls']),'GEOMETRY_ALL_CONTROLS_REJECTED_'+mode)
    rebindings=[]
    for old,new,key,checker in (
        ('FRAME_COORDINATE_FINAL_NORMAL.json','PORTABLE_FRAME_EXPECTED.json','sha256',PATCHED[0]),
        ('PRODUCER_CORE_NORMAL.json','PORTABLE_CORE_EXPECTED.json','input_sha256',PATCHED[1])):
        expected=read('review/'+old);expected[key][checker]=pin(checker)['sha256']
        need(expected==read('review/'+new),'METADATA_ONLY_PORTABLE_EXPECTATION_'+new)
        rebindings.append({'actual_prior_run':'review/'+old,'portable_expected':'review/'+new,'changed_field':key+'.'+checker,'scientific_payload_unchanged':True})
    for name in SMOKES:
        command=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])+[str(ROOT/'review'/name),'--help']
        result=subprocess.run(command,text=True,capture_output=True)
        need(result.returncode==0,'ISOLATED_IMPORT_SMOKE_'+name)
    additions=['review/REVIEW_RECEIPT.json','review/REVIEW_ADDENDUM.md','review/check_portable_review.py',
               'review/PORTABLE_FRAME_EXPECTED.json','review/PORTABLE_CORE_EXPECTED.json','review/REPLAY_ROSTER.json',
               'review/PUBLIC_FILES.json','geometry/GEOMETRY_FREEZE.json','geometry/check_geometric_sources.py',
               'geometry/GEOMETRY_NORMAL.json','geometry/GEOMETRY_OPTIMIZED.json']
    current.update({name:pin(name) for name in additions})
    return {'status':'ACCEPTED_FINAL_INTERNAL_GEOMETRIC_REVIEW_WITH_ISOLATED_REPLAY',
       'blocking_findings':[],'preserved_original_receipt_sha256':ORIGINAL_RECEIPT,
       'scope':'Narrow import-path and isolated-child invocation changes only; previous scientific acceptance unchanged.',
       'verified_path_only_changes':patches,'portable_expected_output_rebindings':rebindings,
       'isolated_import_smokes':list(SMOKES),'negative_controls':[],
       'accepted_current_files':current,
       'scientific_replay_provenance':'Portable expected outputs are explicit metadata-only rebindings of prior actual normal/optimized evidence. Full fresh capsule replay separately executes current checkers. No new mathematical run is claimed by an import smoke.',
       'numeric_arithmetic_scope':'New fast_dyadic implementation and full-domain numerical certificate are bound by numeric-review, not inferred from this addendum.'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    need(not args.output.exists(),'OUTPUT_ALREADY_EXISTS')
    report=build();args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(report['status'])

if __name__=='__main__':main()
