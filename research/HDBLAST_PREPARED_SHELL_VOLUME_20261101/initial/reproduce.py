#!/usr/bin/env python3
"""Portable replay of the unchanged, frozen actual initial mean campaign."""
import argparse,hashlib,importlib.util,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
CAP_REL=Path('research/HDBLAST_PREPARED_SHELL_VOLUME_20261101')
NUM_REL=Path('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric')
sys.dont_write_bytecode=True

class Reject(Exception):
 def __init__(self,code):self.code=code

def require(ok,code):
 if not ok:raise Reject(code)

def digest(path):
 b=path.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}

def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'))

def auth_local():
 pins=json.loads((ROOT/'LOCAL_PINS.json').read_text())
 for rel,pin in pins.items():
  p=ROOT/rel;require(p.is_file() and digest(p)==pin,'LOCAL_AUTHENTICATION')
 return pins

def auth_source(repo,pins):
 result={}
 for rel,pin in pins.items():
  p=repo/rel;require(p.is_file() and digest(p)==pin,'SOURCE_AUTHENTICATION');result[rel]=pin
 return result

def auth_gauge(gauge,pins):
 result={}
 for name,expected in pins.items():
  p=gauge/name;require(p.is_file() and digest(p)['sha256']==expected,'GAUGE_AUTHENTICATION');result[name]=digest(p)
 return result

def load(path,name):
 sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

def execute(args):
 # Every local scientific file and every transitive scientific input is
 # authenticated before any scientific module is imported.
 local=auth_local();freeze=json.loads((ROOT/'CAMPAIGN_FREEZE.json').read_text());source_pins=json.loads((ROOT/'SOURCE_PINS.json').read_text())
 repo=args.repo_root.resolve();gauge=(args.gauge_api or repo/CAP_REL/'gauge').resolve();out=args.output_dir.resolve()
 require(not out.exists(),'OUTPUT_ALREADY_EXISTS')
 for source_root in (repo,ROOT,gauge):require(not out.is_relative_to(source_root),'OUTPUT_INSIDE_SOURCE')
 before=auth_source(repo,source_pins);gauge_before=auth_gauge(gauge,freeze['gauge_api'])
 for rel,pin in freeze['reference_science'].items():require(digest(ROOT/rel)==pin,'REFERENCE_AUTHENTICATION')
 sci=load(ROOT/'scientific.py','prepared_initial_science');m=load(repo/NUM_REL/'interval_core.py','prepared_actual_initial_interval_core')
 m.mean_adapter=load(gauge/'mean_initial.py','prepared_original_mean_gauge');g=load(gauge/'conditioned_mean_initial.py','prepared_conditioned_mean_gauge')
 data,_=m.load()
 calls={'initial_background':0,'evolved_actual_field':0};original_background=m.background
 def counted_background(*a,**kw):calls['initial_background']+=1;return original_background(*a,**kw)
 m.background=counted_background
 raw=sci.generate(m,data)
 static_path=repo/NUM_REL/'inputs/ACTUAL_ANNULUS_STATIC_JETS.json';static=json.loads(static_path.read_text())
 conditioned=sci.conditioned(m,g,raw,static,digest(static_path)['sha256'])
 require(calls['initial_background']==freeze['source_background_calls_per_generation'],'INITIAL_BACKGROUND_CALL_COUNT')
 generated={'RAW_BOXES.json':raw,'CONDITIONED_RAW_BOXES.json':conditioned};comparisons={}
 for name,value in generated.items():
  ref=json.loads((ROOT/'references'/name).read_text());require(value==ref,'SCIENTIFIC_REFERENCE_MISMATCH');comparisons[name]={'every_scientific_field_equal':True,'canonical_sha256':hashlib.sha256(canonical(value).encode()).hexdigest()}
 after=auth_source(repo,source_pins);gauge_after=auth_gauge(gauge,freeze['gauge_api']);require(auth_local()==local,'LOCAL_AUTHENTICATION_AFTER')
 require(before==after and gauge_before==gauge_after,'SOURCE_CUSTODY')
 receipt={'status':'PASS_PORTABLE_FROZEN_INITIAL_GENERATION','campaign':{k:freeze[k] for k in ('N','targets','cells','tail_width','grid_denominator')},'source_capsule_before':before,'source_capsule_after':after,'gauge_before':gauge_before,'gauge_after':gauge_after,'essential_local_hashes_before_scientific_imports':local,'source_counts':calls,'comparisons':comparisons,'representation':'Same exact rational scientific rows as pinned path-free original references; no new targets, grid or source profile tuning.','scope':'Actual initial seed/background only; no evolved actual-field oracle or shell-sign calculation.'}
 out.mkdir(parents=True)
 for name,value in generated.items():(out/name).write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')
 (out/'RECEIPT.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
 return {'status':'PASS','witness':'FROZEN_INITIAL_SCIENCE_MATCH','initial_background_calls':calls['initial_background'],'evolved_actual_field_calls':0,'scientific_outputs':list(generated)}

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--gauge-api',type=Path,help='Directory containing mean_initial.py and conditioned_mean_initial.py; default CAP/gauge');p.add_argument('--output-dir',type=Path,required=True);args=p.parse_args()
 try:final=execute(args);code=0
 except Reject as e:final={'status':'REJECTED','witness':e.code};code=2
 print(json.dumps({'FINAL':final},sort_keys=True));return code
if __name__=='__main__':raise SystemExit(main())
