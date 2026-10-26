#!/usr/bin/env python3
"""Reason-specific tests of fresh external output guards, with no frozen writes."""
import argparse
import hashlib
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SOURCE_SHA = 'd36f8aca0b257ed8bc9de6b2cbacbdcde24ebf3b447823e54d8af006638e08cb'


def need(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--repo-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    path = HERE/'actual_third_jets.py'
    need(hashlib.sha256(path.read_bytes()).hexdigest() == SOURCE_SHA,'SOURCE_PIN')
    spec = importlib.util.spec_from_file_location('tested_actual_source_guards',path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.guard_output(a.repo_root,a.output)
    checks = []

    def rejects(label, target, expected):
        try:
            module.guard_output(a.repo_root,target)
        except ValueError as error:
            reason = str(error)
            need(reason == expected, 'WRONG_GUARD_REASON:'+label+':'+reason)
            checks.append({'case':label,'expected_rejection':expected,'actual_rejection':reason})
        else:
            raise ValueError('UNREJECTED_OUTPUT:'+label)

    with tempfile.TemporaryDirectory(prefix='output_guard_control_') as temp:
        temp = Path(temp)
        existing = temp/'existing'
        existing.write_text('retained\n')
        rejects('existing_file',existing,'FRESH_OUTPUT_REQUIRED')
        link = temp/'dangling_link'
        link.symlink_to(temp/'missing_target')
        rejects('dangling_output_symlink',link,'OUTPUT_SYMLINK')
        for label in module.PROTECTED_LABELS:
            target = a.repo_root/'research'/label/'REJECTED_NEW_OUTPUT_GUARD_CONTROL.json'
            need(not target.exists() and not target.is_symlink(),'GUARD_CONTROL_COLLISION')
            rejects('frozen_'+label,target,'OUTPUT_INSIDE_FROZEN_CAPSULE:'+label)
            need(not target.exists() and not target.is_symlink(),'GUARD_CREATED_FROZEN_OUTPUT')
        own_file = module.__file__
        packaged = temp/'research'/'NEW_CAPSULE'
        module.__file__ = str(packaged/'inputs'/'actual_third_jets.py')
        rejects('own_packaged_capsule',packaged/'new.json','OUTPUT_INSIDE_PACKAGED_CAPSULE')
        module.guard_output(a.repo_root,temp/'permitted_external.json')
        checks.append({'case':'fresh_external_output','accepted':True})
        module.__file__ = own_file
    report = {'status':'PASS_ACTUAL_SOURCE_OUTPUT_GUARDS', 'source_sha256':SOURCE_SHA,
              'checks':checks, 'all_nine_frozen_upstream_roots_protected':True,
              'packaged_capsule_protected':True, 'frozen_files_written':False}
    with a.output.open('x') as handle:
        json.dump(report,handle,indent=2,sort_keys=True)
        handle.write('\n')
    print(report['status'])


if __name__ == '__main__':
    main()
