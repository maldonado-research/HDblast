#!/usr/bin/env python3
"""Reason-specific output/source guard controls, with zero actual calls."""
import argparse
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from output_safety import (guard_output, need, packaged_ancestors, require_hash,
    sha, source_labels, source_snapshot, verify_source_snapshot)


def load_exporter():
    spec = importlib.util.spec_from_file_location('guard_test_exporter',HERE/'export_actual_inputs.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    guard_output(args.output,source_file=__file__)
    sources_before = source_snapshot([HERE/name for name in ('check_numeric_guards.py','output_safety.py',
        'export_actual_inputs.py','evaluate_encoded_map.py','evaluate_selected_magnetic.py')])
    exporter = load_exporter()
    results = []
    def reject(label,reason,fn):
        try:
            fn()
        except ValueError as error:
            need(str(error) == reason,'WRONG_REJECTION:'+label+':'+str(error))
        else:
            raise ValueError('MISSING_REJECTION:'+label)
        results.append({'name':label,'passed':True,'expected_reason':reason})
    def accept(label,fn):
        fn()
        results.append({'name':label,'passed':True})
    with tempfile.TemporaryDirectory(prefix='numeric-guard-fixtures-',dir=args.output.parent) as temporary:
        root = Path(temporary)
        outside = root/'scratch'
        outside.mkdir()
        source_file = HERE/'export_actual_inputs.py'
        guard = lambda output,**kw:guard_output(output,source_file=source_file,**kw)
        accept('fresh_external_path',lambda:guard(outside/'fresh.json'))
        existing = outside/'existing.json';existing.write_text('{}\n')
        reject('existing_file','OUTPUT_EXISTS',lambda:guard(existing))
        reject('existing_directory','OUTPUT_EXISTS',lambda:guard(outside))
        dangling = outside/'dangling.json';dangling.symlink_to(outside/'missing')
        reject('dangling_leaf_symlink','OUTPUT_SYMLINK',lambda:guard(dangling))
        linked = outside/'linked.json';linked.symlink_to(existing)
        reject('existing_leaf_symlink','OUTPUT_SYMLINK',lambda:guard(linked))
        reject('source_tree_child','OUTPUT_INSIDE_SOURCE_TREE',lambda:guard(HERE/'never-written.json'))
        reject('source_tree_nested','OUTPUT_INSIDE_SOURCE_TREE',lambda:guard(HERE/'absent'/'never-written.json'))
        reject('source_tree_dotdot','OUTPUT_INSIDE_SOURCE_TREE',lambda:guard(HERE/'absent'/'..'/'never-written.json'))
        alias = outside/'source-alias';alias.symlink_to(HERE,target_is_directory=True)
        reject('symlink_parent_into_source','OUTPUT_INSIDE_SOURCE_TREE',lambda:guard(alias/'never-written.json'))
        sibling = outside/'numeric-shared-prefix';sibling.mkdir()
        accept('shared_name_prefix_sibling',lambda:guard(sibling/'allowed.json'))
        repo = root/'fake-repo'
        deps = [('DEPENDENCY:'+key,repo/'research'/label) for key,(label,pin) in exporter.CAPSULES.items()]
        for reason,protected in deps:
            reject('dependency_'+reason.split(':',1)[1],'OUTPUT_INSIDE_'+reason,
                   lambda p=protected:guard(p/'absent'/'never-written.json',protected_roots=deps))
        capsule = root/'research'/'HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026'
        packaged_source = capsule/'numeric'/'code.py'
        reject('own_packaged_capsule_ancestor','OUTPUT_INSIDE_PACKAGED_CAPSULE',lambda:
            guard_output(capsule/'other'/'never-written.json',source_file=packaged_source,protected_roots=deps))
        input_capsule = root/'external-repo'/'research'/'HDBLAST_EXTERNAL_INPUTS'
        reject('input_capsule_ancestor','OUTPUT_INSIDE_INPUT_CAPSULE',lambda:
            guard(input_capsule/'other'/'never-written.json',input_paths=[input_capsule/'numeric'/'inputs.json']))
        map_capsule = root/'external-repo'/'research'/'HDBLAST_EXTERNAL_MAP'
        reject('map_capsule_ancestor','OUTPUT_INSIDE_INPUT_CAPSULE',lambda:
            guard(map_capsule/'other'/'never-written.json',input_paths=[None,map_capsule/'projection'/'map.json']))
        reject('missing_parent','OUTPUT_PARENT_MISSING',lambda:guard(outside/'missing-parent'/'report.json'))
        fixture_source = outside/'fixture_source.py';fixture_source.write_text('original\n')
        snapshot = source_snapshot([fixture_source])
        accept('unchanged_source_snapshot',lambda:verify_source_snapshot(snapshot))
        need(list(source_labels(snapshot)) == ['fixture_source.py'],'STABLE_BASENAME_SOURCE_LABEL')
        results.append({'name':'stable_source_label','passed':True})
        fixture_source.write_text('changed\n')
        reject('source_tamper','SOURCE_CHANGED:fixture_source.py',lambda:verify_source_snapshot(snapshot))
        reject('oracle_source_pin','ORACLE_EXTERNAL_PIN',lambda:require_hash(fixture_source,'0'*64,'ORACLE_EXTERNAL_PIN'))
        # Production CLI rejection ordering proves no repository authentication,
        # source import, or actual oracle evaluation is reached for bad outputs.
        def cli_reject(label,script,arguments,reason):
            command = [sys.executable,'-I','-S','-B']
            if sys.flags.optimize:
                command.append('-O')
            command.extend([str(HERE/script),*map(str,arguments)])
            result = subprocess.run(command,text=True,capture_output=True)
            need(result.returncode != 0 and result.stderr.rstrip().endswith('ValueError: '+reason),
                 'CLI_GUARD_REJECTION:'+label)
            results.append({'name':label,'passed':True,'expected_reason':reason})
        cli_reject('export_guard_before_repo_auth','export_actual_inputs.py',
            ['--repo-root',repo,'--output',HERE/'never-written.json'],'OUTPUT_INSIDE_SOURCE_TREE')
        cli_reject('export_dependency_guard_before_repo_auth','export_actual_inputs.py',
            ['--repo-root',repo,'--output',deps[0][1]/'never-written.json'],'OUTPUT_INSIDE_DEPENDENCY:exterior')
        cli_reject('selected_input_pin_before_json','evaluate_selected_magnetic.py',
            ['--inputs',existing,'--map',existing,'--require-input-sha256','0'*64,
             '--require-map-sha256',sha(existing),'--output',outside/'bad-input-output.json'],'INPUT_EXTERNAL_PIN')
        cli_reject('selected_map_pin_before_json','evaluate_selected_magnetic.py',
            ['--inputs',existing,'--map',existing,'--require-input-sha256',sha(existing),
             '--require-map-sha256','0'*64,'--output',outside/'bad-map-output.json'],'MAP_EXTERNAL_PIN')
        cli_reject('generic_map_pin_before_json','evaluate_encoded_map.py',
            ['--inputs',existing,'--map',existing,'--require-map-sha256','0'*64,
             '--output',outside/'bad-generic-output.json'],'MAP_EXTERNAL_PIN')
    verify_source_snapshot(sources_before)
    report = {'schema':'HDBLAST_NUMERIC_PORTABLE_OUTPUT_SOURCE_GUARDS_V1','all_passed':True,
      'actual_solution_calls':0,'frozen_repository_imports':0,'controls':results,
      'source_integrity_before_and_after':source_labels(sources_before)}
    with args.output.open('x') as handle:
        json.dump(report,handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps({'output':str(args.output),'sha256':sha(args.output),'controls':len(results),'all_passed':True}))


if __name__ == '__main__':
    main()
