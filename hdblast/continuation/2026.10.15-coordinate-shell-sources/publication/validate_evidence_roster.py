"""Network-free synthetic tests for the dynamic frozen-roster evidence contract."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent


def main():
    path=ROOT/'build_final_evidence.py';spec=importlib.util.spec_from_file_location('builder',path)
    builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)
    counts={'successful_checker_runs':2,'external_corruption_processes_rejected':2,
            'embedded_wrong_formula_controls_rejected':2,'finite_positive_check_groups_each_mode':3}
    roster={'schema':'hdblast-replay-roster-v1','status':'APPROVED_FOR_FROZEN_REPLAY',
            'capsule':'HDBLAST_COORDINATE_SHELL_SOURCES_20261015','modes':['normal','optimized'],'counts':counts,
            'positive_processes':[{'name':'synthetic','script':'synthetic.py','finite_positive_groups_each_mode':3,'embedded_controls_each_mode':1,'finite_positive_check_names':['a','b','c'],'embedded_control_names':['wrong']} ],
            'external_controls':[{'checker':'synthetic','name':'sign','expected_reason':'wrong sign','mechanism':'synthetic formula injection'}]}
    steps=[]
    for mode in roster['modes']:
        steps.extend([{'mode':mode,'kind':'positive','checker':'synthetic','script':'synthetic.py','exit_code':0,'scientific_output_matches':True,
                       'finite_positive_check_groups':3,'embedded_wrong_formula_controls':1,'finite_positive_check_names':['a','b','c'],'embedded_control_names':['wrong']},
                      {'mode':mode,'kind':'external_corruption','checker':'synthetic','script':'synthetic.py','control':'sign','reported_failure':'wrong sign',
                       'exit_code':1,'expected_failure_exit_code':1,'expected_reason_detected':True,'expected_reason':'wrong sign','control_mechanism':'synthetic formula injection'}])
    replay={'status':'PASS','capsule':roster['capsule'],**counts,'steps':steps}
    if builder.validate_roster(roster,replay)!=counts:raise RuntimeError('Valid synthetic new counts rejected')
    rejected={}
    controls={
      'old_checkpoint_counts':lambda r,p:p.update(finite_positive_check_groups_each_mode=counts['finite_positive_check_groups_each_mode']+1),
      'unapproved_roster':lambda r,p:r.update(status='PENDING'),
      'missing_optimized_process':lambda r,p:p['steps'].pop(),
      'duplicate_process':lambda r,p:p['steps'].append(copy.deepcopy(p['steps'][0])),
      'unknown_mode':lambda r,p:p['steps'][0].update(mode='unchecked'),
      'missing_expected_rejection':lambda r,p:p['steps'][1].update(expected_reason_detected=False),
      'unrelated_failure_exit':lambda r,p:p['steps'][1].update(exit_code=2),
      'mismatched_finite_groups':lambda r,p:p['steps'][2].update(finite_positive_check_groups=4),
      'different_mode_control':lambda r,p:p['steps'][3].update(control='another'),
      'positive_output_not_reproduced':lambda r,p:p['steps'][0].update(scientific_output_matches=False),
      'same_counts_wrong_script':lambda r,p:[row.update(script='wrong.py') for row in p['steps']],
      'same_counts_wrong_checker':lambda r,p:[row.update(checker='different') for row in p['steps']],
      'same_counts_wrong_reason':lambda r,p:[row.update(expected_reason='unrelated',reported_failure='unrelated') for row in p['steps'] if row['kind']=='external_corruption'],
      'same_counts_wrong_mechanism':lambda r,p:[row.update(control_mechanism='unrelated mechanism') for row in p['steps'] if row['kind']=='external_corruption'],
      'changed_roster_mechanism':lambda r,p:r['external_controls'][0].update(mechanism='unapproved mechanism'),
      'same_counts_wrong_group_name':lambda r,p:p['steps'][0]['finite_positive_check_names'].__setitem__(0,'unapproved'),
      'same_counts_wrong_embedded_name':lambda r,p:p['steps'][0]['embedded_control_names'].__setitem__(0,'unapproved'),
      'changed_roster_script':lambda r,p:r['positive_processes'][0].update(script='another.py'),
    }
    for name,change in controls.items():
        r,p=copy.deepcopy(roster),copy.deepcopy(replay);change(r,p)
        try:builder.validate_roster(r,p)
        except RuntimeError:rejected[name]=True
        else:raise RuntimeError('Invalid synthetic roster accepted: '+name)
    helpers=['build_final_evidence.py','readback_evidence.py','verify_publication_scope.py']
    pins={}
    for name in helpers:
        raw=(ROOT/name).read_bytes();compile(raw,str(ROOT/name),'exec')
        result=subprocess.run([sys.executable,'-I','-B',str(ROOT/name),'--help'],capture_output=True)
        if result.returncode:raise RuntimeError('Help failed: '+name)
        pins[name]=hashlib.sha256(raw).hexdigest()
    receipt={'status':'PASS_SYNTHETIC_DYNAMIC_ROSTER_GATES','positive_schema_cases':1,
             'rejections':rejected,'helper_sha256':pins,'http_requests':0,'git_commands':0,
             'scope':'Synthetic schema and CLI checks only. No scientific replay, frozen source, PR, merge or publication has been validated.'}
    (ROOT/'EVIDENCE_ROSTER_TOOL_VALIDATION.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
