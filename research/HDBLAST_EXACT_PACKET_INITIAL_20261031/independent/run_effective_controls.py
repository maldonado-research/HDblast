#!/usr/bin/env python3
"""Real copied-script, selected-source and rehashed-report rejections."""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,sort_keys=True,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--repo-root',type=Path,required=True)
    p.add_argument('--initial',type=Path,required=True)
    p.add_argument('--producer-response',type=Path,required=True)
    p.add_argument('--initial-source',type=Path,required=True)
    p.add_argument('--initial-freeze',type=Path,required=True)
    p.add_argument('--response-source',type=Path,required=True)
    p.add_argument('--checker',type=Path,required=True)
    p.add_argument('--review',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    need(not a.output.exists(),'fresh absent controls output')
    fixtures=a.output.parent/(a.output.stem+'-fixtures')
    need(not fixtures.exists(),'fresh absent fixtures directory')
    fixtures.mkdir(parents=True)
    pre={x:sha(x) for x in (Path(__file__).resolve(),a.checker,a.initial,a.producer_response,a.initial_source,a.initial_freeze,a.response_source,a.review)}
    initial=json.loads(a.initial.read_text()); response=json.loads(a.producer_response.read_text())
    review=json.loads(a.review.read_text()); terms=review['named_nonnegative_terms']
    cases=[]
    def report_case(name, which, mutate, expected, witness):
        c=copy.deepcopy(initial if which=='initial' else response)
        mutate(c)
        path=fixtures/(name+'.json');write(path,c)
        cases.append(dict(name=name,initial=path if which=='initial' else a.initial,
                          response=path if which=='response' else a.producer_response,
                          checker=a.checker,root=a.repo_root,expected=expected,witness=str(witness),fixture=path))
    report_case('omit_velocity_error_charge','initial',lambda d:d['campaigns'][0]['initial_error_H_caps'].update(v0_r='0'),
                'unchanged complete velocity error charges',-Q(initial['campaigns'][0]['initial_error_H_caps']['v0_r']))
    report_case('absolute_error_jet_interchange','initial',lambda d:d['campaigns'][0]['candidate_initial_absolute_H_caps'].update(u_rr=d['campaigns'][0]['initial_error_H_caps']['e0_rr']),
                'fresh absolute displacement caps aggregated from full leaves',
                Q(initial['campaigns'][0]['initial_error_H_caps']['e0_rr'])-Q(initial['campaigns'][0]['candidate_initial_absolute_H_caps']['u_rr']))
    report_case('zero_background_channel','initial',lambda d:d['campaigns'][0]['paired_channel_H_caps'][0].__setitem__(2,'0'),
                'global paired channel maxima',-Q(initial['campaigns'][0]['paired_channel_H_caps'][0][2]))
    report_case('skip_near_edge_tail_cell','initial',lambda d:(d['campaigns'][0]['leaves'].pop(1),d['campaigns'][0].update(leaf_count=d['campaigns'][0]['leaf_count']-1)),
                'complete ordered leaf cover including tail',Q(1,64))
    report_case('global_denominator_used_locally','initial',lambda d:d['campaigns'][0]['local']['actual_H_norm_interval'].__setitem__(0,'3249000000000000000000000000000000'),
                'independently aggregate local actual lower and new error upper',
                Q(3249*10**30)-Q(initial['campaigns'][0]['local']['actual_H_norm_interval'][0]))
    report_case('filtered_Fx_report_substitution','response',lambda d:d['constants'].update(Fx='66768050521240320'),
                'every reported operator source contraction and commutator constant',Q(66768050521240320)-Q(905405545954284))
    report_case('copied_old_point_values','response',lambda d:d.update(full23_point_table={'inherited':'old shell point (9/8,0)'}),
                'no copied old point or quadrature table in ledger-only report','forbidden full23_point_table inserted')
    report_case('incorrect_response_numeric_cap','response',lambda d:d['ledger'].update(r0_initial_H_cap=str(Q(d['ledger']['r0_initial_H_cap'])+1)),
                'independent exact equation mismatch: r0_initial_H_cap',Q(1))
    report_case('different_affine_shell_law','response',lambda d:d['boundary_error'].update(g_x='1'),
                'same original candidate source and affine error boundary law',Q(1))
    script=a.checker.read_text()
    for name,old,new,expected,witness in [
        ('omit_frozen_initial_acceleration','initial_acceleration=P*fzero,','initial_acceleration=Q_ZERO,','W0_minus_Bs_time_derivative_X_cap',-Q(terms['W0_minus_Bs_time_derivative_X_cap']['initial_acceleration'])),
        ('omit_bulk_initial_acceleration','initial_acceleration=P*rzero,','initial_acceleration=Q_ZERO,','correction_time_derivative_X_cap',-Q(terms['correction_time_derivative_X_cap']['initial_acceleration'])),
        ('omit_original_source_derivative',"original_source_derivative=k['Fx'])",'original_source_derivative=Q_ZERO)', 'r0_time_derivative_H_cap',-Q(905405545954284)),
        ('omit_velocity_gradient_in_initial_J',"e['v0_r']+w)",'Q_ZERO+w)', 'differentiated_initial_free_wave_X_cap',-Q(initial['campaigns'][0]['initial_error_H_caps']['v0_r'])),
    ]:
        need(script.count(old)==1,'unique mutation anchor: '+name)
        changed=script.replace(old,new).replace('import argparse\n','import argparse\nQ_ZERO = 0\n',1)
        path=fixtures/(name+'.py');path.write_text(changed)
        cases.append(dict(name=name,initial=a.initial,response=a.producer_response,checker=path,root=a.repo_root,
                          expected='independent exact equation mismatch: '+expected,witness=str(witness),fixture=path))
    # Copy only the selected source closure needed by the checker. Originals are
    # never edited; each source case has a complete independent physical copy.
    selected=set(initial['authentication_before_after'])
    selected.update(['research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TRANSFER.json',
        'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TIME.json',
        'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json',
        'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/CANDIDATE.json',
        'research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json'])
    for name,relative,mutator,expected,witness in [
        ('registered_formula_source_corruption','research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/full23_initial.py',
         lambda path:path.write_text(path.read_text()+'\n# nonzero deliberate source corruption\n'),
         'initial selected source closure hash','appended nonempty byte string changes pinned source SHA256'),
        ('actual_Fx_input_corruption','research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json',
         lambda path:write(path,dict(json.loads(path.read_text()),genuine_scaled_F_x_H_cap='905405545954285')),
         'retained hash mismatch: original_Fx',Q(1))]:
        root=fixtures/name/'repo'
        for rel in selected:
            source=a.repo_root/rel; destination=root/rel
            destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,destination)
        mutator(root/relative)
        cases.append(dict(name=name,initial=a.initial,response=a.producer_response,checker=a.checker,root=root,
                          expected=expected,witness=str(witness),fixture=root/relative))
    receipts=[]
    for case in cases:
        runs=[]
        for mode in ('normal','optimized'):
            out=fixtures/(case['name']+'-'+mode+'-MUST_NOT_EXIST.json')
            command=[sys.executable,'-I','-S','-B']+(['-O'] if mode=='optimized' else [])+[str(case['checker']),
                '--repo-root',str(case['root']),'--initial',str(case['initial']),'--initial-sha256',sha(case['initial']),
                '--producer-response',str(case['response']),'--producer-response-sha256',sha(case['response']),
                '--initial-source',str(a.initial_source),'--initial-freeze',str(a.initial_freeze),
                '--response-source',str(a.response_source),'--output',str(out)]
            result=subprocess.run(command,capture_output=True,text=True)
            need(result.returncode!=0 and case['expected'] in result.stderr and not out.exists(),
                 'effective rejected child execution: '+case['name']+' '+mode+'\n'+result.stderr)
            runs.append(dict(mode=mode,returncode=result.returncode,expected_failure=case['expected'],
                             successful_output_absent=True,stderr_last_line=result.stderr.strip().splitlines()[-1]))
        receipts.append(dict(name=case['name'],exact_nonzero_witness_or_failed_identity=case['witness'],
                             fixture_sha256=sha(case['fixture']),runs=runs))
    need(all(sha(x)==value for x,value in pre.items()),'all original selected inputs and controls helper unchanged')
    write(a.output,dict(status='PASS_EFFECTIVE_RESPONSE_INDEPENDENT_CONTROLS',case_count=len(cases),
          child_execution_count=2*len(cases),cases=receipts,original_files_unchanged=True,
          input_sha256={role:sha(path) for role,path in dict(controls_self=Path(__file__).resolve(),checker=a.checker,
              initial_report=a.initial,producer_response=a.producer_response,initial_source=a.initial_source,
              initial_freeze=a.initial_freeze,response_source=a.response_source,review=a.review).items()},
          scope='Each real copied source/script or rehashed report corruption executes normal and optimized checker and fails before any success output.'))


if __name__=='__main__':
    main()
