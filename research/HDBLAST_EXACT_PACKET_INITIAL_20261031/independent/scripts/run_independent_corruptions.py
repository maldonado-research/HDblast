#!/usr/bin/env python3
"""Exercise actual on-disk corruptions with absent successful artifacts."""
import argparse,copy,hashlib,json,subprocess,sys
from pathlib import Path
sys.dont_write_bytecode=True
def need(ok,why):
    if not ok:raise RuntimeError(why)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--report',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--source-pins',type=Path,default=Path(__file__).resolve().parent.parent/'source-pins.json');a=ap.parse_args()
    need(not a.output_dir.exists(),'FRESH_DIRECTORY_REQUIRED');a.output_dir.mkdir(parents=True)
    script=Path(__file__).resolve();checker=script.with_name('check_independent_report.py');helper=script.with_name('check_independent_exact.py');rawfixture=script.with_name('bad_raw_quadratic_fixture.py')
    before={str(p.name):digest(p) for p in (script,checker,helper,rawfixture,a.source_pins)};base=json.loads(a.report.read_bytes());results=[];python=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])
    mutations=[]
    b=copy.deepcopy(base);del b['independent_full_displacement_coverage']['leaves'][10];b['independent_full_displacement_coverage']['leaves_count']-=1;mutations.append(('skipped_interior_leaf',b,'FULL_FOOTPRINT_COVERAGE'))
    b=copy.deepcopy(base);del b['independent_full_displacement_coverage']['leaves'][0];b['independent_full_displacement_coverage']['leaves_count']-=1;mutations.append(('skipped_flat_tail',b,'FOOTPRINT_ENDPOINTS'))
    b=copy.deepcopy(base);b['actual_target_H_interval']=['3249000000000000000000000000000000']*2;mutations.append(('global_denominator_used_locally',b,'LOCAL_DENOMINATOR_FROM_LOCAL_ALL23_COMPONENTS'))
    b=copy.deepcopy(base)
    for row in b['local_components']:row['paired_error_rr']=['0','0'];row['weighted_error_square_upper']='0'
    mutations.append(('zero_background_difference',b,'LOCAL_ERROR_RECONSTRUCTION'))
    b=copy.deepcopy(base);b['independent_full_displacement_coverage']['initial_gate_thresholds'][0]='0';mutations.append(('changed_initial_threshold',b,'FROZEN_INITIAL_EXACT_THRESHOLDS'))
    for name,bad,expected in mutations:
        case=a.output_dir/name;case.mkdir();inputpath=case/'CORRUPTED.json';inputpath.write_text(json.dumps(bad,sort_keys=True,indent=2)+'\n');out=case/'SUCCESS.json'
        cmd=python+[str(checker),'--source',str(a.source),'--report',str(inputpath),'--output',str(out),'--source-pins',str(a.source_pins)]
        p=subprocess.run(cmd,capture_output=True,text=True);need(p.returncode!=0 and expected in p.stderr and not out.exists(),'CORRUPTION_NOT_REJECTED:'+name)
        (case/'stderr.txt').write_text(p.stderr);results.append({'case':name,'returncode':p.returncode,'failure':p.stderr.splitlines()[-1],'success_artifact_absent':not out.exists(),'corrupted_input_sha256':digest(inputpath)})
    case=a.output_dir/'source_pin_map';case.mkdir();badpins=json.loads(a.source_pins.read_text());badpins[next(iter(badpins))]='0'*64;pinpath=case/'CORRUPTED_SOURCE_PINS.json';pinpath.write_text(json.dumps(badpins,sort_keys=True)+'\n');out=case/'SUCCESS.json'
    p=subprocess.run(python+[str(checker),'--source',str(a.source),'--report',str(a.report),'--output',str(out),'--source-pins',str(pinpath)],capture_output=True,text=True)
    need(p.returncode!=0 and 'PINNED_SOURCE_MAP' in p.stderr and not out.exists(),'SOURCE_MAP_CORRUPTION_NOT_REJECTED');(case/'stderr.txt').write_text(p.stderr);results.append({'case':'source_pin_map','returncode':p.returncode,'failure':p.stderr.splitlines()[-1],'success_artifact_absent':True,'corrupted_input_sha256':digest(pinpath)})
    case=a.output_dir/'raw_packet_quadratic';case.mkdir();out=case/'SUCCESS.json'
    p=subprocess.run(python+[str(rawfixture),'--source',str(a.source),'--helper',str(helper),'--helper-sha256',digest(helper),'--source-pins',str(a.source_pins),'--output',str(out)],capture_output=True,text=True)
    need(p.returncode!=0 and 'REQUIRED_RAW_QUADRATIC_SCALAR_VALUE_MATCH_FAILED' in p.stderr and not out.exists(),'RAW_QUADRATIC_CORRUPTION_NOT_REJECTED');(case/'stderr.txt').write_text(p.stderr);results.append({'case':'raw_packet_quadratic','returncode':p.returncode,'failure':p.stderr.splitlines()[-1],'success_artifact_absent':True,'fixture_sha256':digest(rawfixture)})
    need(before=={str(p.name):digest(p) for p in (script,checker,helper,rawfixture,a.source_pins)},'EXECUTED_HELPERS_MAP_UNCHANGED')
    result={'status':'PASS_REAL_INDEPENDENT_DISK_CORRUPTIONS','cases':results,'executed_helpers_and_map_sha256_before_after':before,'report_sha256':digest(a.report)}
    (a.output_dir/'CONTROLS.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(result['status'])
if __name__=='__main__':main()
