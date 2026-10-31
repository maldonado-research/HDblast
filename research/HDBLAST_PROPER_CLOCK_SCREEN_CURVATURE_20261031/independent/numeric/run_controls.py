#!/usr/bin/env python3
"""Real subprocess rejection controls for the independent numeric reviewer."""
import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def need(ok,msg):
    if not ok:raise ValueError(msg)


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def final_witness(proc):
    lines=[line.strip() for line in proc.stderr.splitlines()if line.strip()]
    line=lines[-1]if lines else ''
    prefix='ValueError: '
    return line[len(prefix):]if line.startswith(prefix)else None


def has_success_stdout(proc):
    for line in proc.stdout.splitlines():
        try:value=json.loads(line)
        except(ValueError,TypeError):continue
        if isinstance(value,dict)and isinstance(value.get('status'),str)and value['status'].startswith(('PASS','CERTIFIED')):
            return True
    return False


def rejection_matches(proc,out,expected):
    return(proc.returncode!=0 and final_witness(proc)==expected and not out.exists()
           and not has_success_stdout(proc))


def main():
    ap=argparse.ArgumentParser()
    for flag in ('repo-root','producer-root','report','background-root','tt-root','formula-root','output'):
        ap.add_argument('--'+flag,type=Path,required=True)
    a=ap.parse_args();need(not a.output.exists(),'FRESH_CONTROL_OUTPUT')
    checker=Path(__file__).resolve().with_name('review_numeric.py')
    original=json.loads(a.report.read_bytes())
    jobs=[]
    def scalar(key,value):
        def change(r):r[key]=value
        return change
    def cell(key,value):
        def change(r):r['cells'][7][key]=value
        return change
    def point_invariant(r):r['point_weyl_invariants']['Ay']=['0','0']
    def raw_cell(r):r['cells'][9]['metric_m1']['D']=['0','0']
    def raw_apex(r):r['apex_metric_m2']['Myy']=['0','0']
    def time_gap(r):r['cells'][11]['x'][0]='0'
    def omit_cell(r):r['cells'].pop()
    jobs=[
      ('wrong_SP',cell('SP',['0','0']),'INTERVAL_RECONSTRUCTION:SP:7'),
      ('wrong_SQ',cell('SQ',['0','0']),'INTERVAL_RECONSTRUCTION:SQ:7'),
      ('halve_proper_weight',cell('proper_length_upper',str(F(original['cells'][7]['proper_length_upper'])/2)),'RATIONAL_RECONSTRUCTION:PROPER_LENGTH:7'),
      ('zero_M_cell',cell('M_cell','0'),'RATIONAL_RECONSTRUCTION:M_cell:7'),
      ('zero_K_cell',cell('K_cell','0'),'RATIONAL_RECONSTRUCTION:K_cell:7'),
      ('zero_source_integral',scalar('source_integral','0'),'RATIONAL_RECONSTRUCTION:source_integral'),
      ('drop_drift_denominator',scalar('den','1'),'RATIONAL_RECONSTRUCTION:den'),
      ('drop_X_denominator',scalar('den_X','1'),'RATIONAL_RECONSTRUCTION:den_X'),
      ('zero_position_X',scalar('position_X_cap','0'),'RATIONAL_RECONSTRUCTION:position_X_cap'),
      ('zero_position_Y',scalar('position_Y_cap','0'),'RATIONAL_RECONSTRUCTION:position_Y_cap'),
      ('zero_apex_Ay',point_invariant,'INTERVAL_RECONSTRUCTION:WEYL_Ay'),
      ('drop_radial_quadratic',scalar('Fy',['0','0']),'INTERVAL_RECONSTRUCTION:APEX_FY'),
      ('drop_full_tilt',scalar('full_geodesic_tilt_interval',['0','0']),'INTERVAL_RECONSTRUCTION:full_geodesic_tilt_interval'),
      ('drop_TT_term',scalar('metric_linear_TT',['0','0']),'INTERVAL_RECONSTRUCTION:metric_linear_TT'),
      ('fabricate_negative_full_interval',scalar('full_coefficient',['-1','-1/10']),'INTERVAL_RECONSTRUCTION:full_coefficient'),
      ('narrow_decimal_without_proof',scalar('outward_decimal',['-0.404312161091','-0.02']),'OUTWARD_DECIMAL_EXACT'),
      ('extend_to_invalid_aux_clock',scalar('restricted_tau_upper','1/5000'),'RATIONAL_RECONSTRUCTION:RESTRICTED_CLOCK'),
      ('change_raw_oracle_cell',raw_cell,'CANONICAL_RAW_ORACLE_CELL:9'),
      ('change_raw_oracle_apex',raw_apex,'CANONICAL_RAW_ORACLE_APEX'),
      ('missing_source_cell',omit_cell,'CANDIDATE_CELL_COUNT'),
      ('time_cell_gap',time_gap,'SOURCE_CELL_TIME_COVERAGE:11'),
    ]
    results=[]
    with tempfile.TemporaryDirectory(prefix='proper-clock-independent-controls-')as tmp:
        tmp=Path(tmp)
        # A token in an earlier traceback line must not validate a different
        # final failure. This is a matcher self-test, separate from the 48
        # scientific report/source corruptions below.
        token='INTENDED_WITNESS_ONLY_IN_EARLIER_LINE'
        fixture=subprocess.run([sys.executable,'-I','-S','-B','-c',
            'import sys;sys.stderr.write("INTENDED_WITNESS_ONLY_IN_EARLIER_LINE\\nValueError: DIFFERENT_FINAL_FAILURE\\n");sys.exit(1)'],
            capture_output=True,text=True,timeout=10)
        fixture_out=tmp/'witness-contract-must-not-exist.json'
        need(not rejection_matches(fixture,fixture_out,token)
             and final_witness(fixture)=='DIFFERENT_FINAL_FAILURE',
             'FINAL_WITNESS_CONTRACT_SELFTEST')
        common=['--repo-root',str(a.repo_root),'--producer-root',str(a.producer_root),'--report',str(a.report),
                '--background-root',str(a.background_root),'--tt-root',str(a.tt_root),'--formula-root',str(a.formula_root)]
        def reject(name,extra,expected,replace=None):
            for optimized in (False,True):
                argv=[sys.executable,'-B']+(['-O']if optimized else [])+[str(checker)]+list(common)
                if replace:
                    flag,value=replace;argv[argv.index(flag)+1]=str(value)
                out=tmp/(name+('-O'if optimized else '')+'.out.json')
                argv+=extra+['--output',str(out)]
                proc=subprocess.run(argv,capture_output=True,text=True,timeout=20)
                need(rejection_matches(proc,out,expected),'CONTROL_NOT_REJECTED:'+name)
                results.append({'name':name,'optimized':optimized,'expected_rejection':expected,
                    'expected_witness':expected,'observed_witness':final_witness(proc),
                    'subprocess_exit':proc.returncode,'output_absent':True,
                    'success_artifact_absent':True,'success_status_absent_from_stdout':True})
        for name,change,expected in jobs:
            altered=deepcopy(original);change(altered);path=tmp/(name+'.json');path.write_text(json.dumps(altered)+'\n')
            reject(name,['--candidate-report',str(path)],expected)
        changed=tmp/'badcanonical.json';changed.write_text(json.dumps({**original,'source_integral':'0'})+'\n')
        reject('canonical_report_pin',[],'CANONICAL_REPORT_PIN',('--report',changed))
        copied=tmp/'producer';copied.mkdir()
        for f in ('check_proper_clock.py','SOURCE_PINS.json','FREEZE.json'):shutil.copyfile(a.producer_root/f,copied/f)
        with (copied/'check_proper_clock.py').open('a')as f:f.write('\nraise RuntimeError("MUST_NOT_IMPORT_MUTATED_PRODUCER")\n')
        reject('producer_source_pin',[],'FROZEN_PRODUCER_FILE:check_proper_clock.py',('--producer-root',copied))
        formula=tmp/'formula';formula.mkdir()
        for f in ('proper_clock_coefficient.py','FINAL_FORMULA_NORMAL.json'):shutil.copyfile(a.formula_root/f,formula/f)
        with(formula/'proper_clock_coefficient.py').open('a')as f:f.write('\nraise RuntimeError("MUST_NOT_IMPORT_MUTATED_FORMULA")\n')
        reject('formula_source_pin',[],'COMPONENT_PIN:FORMULA_API',('--formula-root',formula))
    report={'status':'PASS_INDEPENDENT_NUMERIC_REPORT_AND_SOURCE_CORRUPTION_CONTROLS',
            'checker_sha256':sha(checker),'controls_source_sha256':sha(__file__),
            'distinct_controls':len(jobs)+3,'subprocess_rejections':len(results),
            'both_normal_and_optimized':True,'numerical_controls_bypass_only_candidate_hash_with_raw_oracle_boxes_still_canonically_bound':True,
            'witness_contract_selftest':{'synthetic_subprocesses':1,
                'earlier_expected_token_with_different_final_failure_rejected':True,
                'counted_as_scientific_rejection_control':False},
            'results':results}
    with a.output.open('x')as f:json.dump(report,f,sort_keys=True,indent=2);f.write('\n')
    print(json.dumps({'status':report['status'],'subprocess_rejections':len(results),'sha256':sha(a.output)}))


if __name__=='__main__':main()
