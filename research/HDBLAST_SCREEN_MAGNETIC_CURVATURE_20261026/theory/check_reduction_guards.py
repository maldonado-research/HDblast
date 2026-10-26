#!/usr/bin/env python3
"""Small non-scientific controls for reduction output and pre-import guards."""
import argparse
import importlib.util
import json
import sys
import tempfile
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    spec=importlib.util.spec_from_file_location('screen_reduction_guards_target',HERE/'check_reduction.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.validate_output(args.output)
    checks=[]
    def reject(name,call,reason):
        try:call()
        except ValueError as exc:
            if str(exc)!=reason:raise ValueError('WRONG_REJECTION_'+name) from exc
            checks.append(name)
        else:raise ValueError('ACCEPTED_'+name)
    with tempfile.TemporaryDirectory(prefix='hdblast-screen-theory-guards-') as directory:
        root=Path(directory)
        if mod.validate_output(root/'fresh.json')!=root/'fresh.json':raise ValueError('SCRATCH_TARGET')
        checks.append('external_fresh_scratch_accepted')
        existing=root/'existing.json';existing.write_text('preserve\n')
        reject('existing_file_rejected',lambda:mod.validate_output(existing),'FRESH_NONSYMLINK_OUTPUT_REQUIRED')
        link=root/'link.json';link.symlink_to(existing)
        reject('existing_symlink_rejected',lambda:mod.validate_output(link),'FRESH_NONSYMLINK_OUTPUT_REQUIRED')
        dangling=root/'dangling.json';dangling.symlink_to(root/'absent.json')
        reject('dangling_symlink_rejected',lambda:mod.validate_output(dangling),'FRESH_NONSYMLINK_OUTPUT_REQUIRED')
        reject('fresh_source_child_rejected',lambda:mod.validate_output(HERE/'__guard_fresh__.json'),'OUTPUT_INSIDE_SOURCE_OR_RESEARCH_CAPSULE')
        source_alias=root/'source-alias';source_alias.symlink_to(HERE,target_is_directory=True)
        reject('resolved_source_alias_rejected',lambda:mod.validate_output(source_alias/'__guard_fresh__.json'),'OUTPUT_INSIDE_SOURCE_OR_RESEARCH_CAPSULE')
        capsule=root/'research'/'HDBLAST_GUARD_FIXTURE';nested=capsule/'theory'
        nested.mkdir(parents=True)
        original_here=mod.HERE
        try:
            mod.HERE=nested
            reject('packaged_capsule_sibling_rejected',lambda:mod.validate_output(capsule/'fresh.json'),'OUTPUT_INSIDE_SOURCE_OR_RESEARCH_CAPSULE')
            if mod.validate_output(root/'outside-capsule.json')!=root/'outside-capsule.json':raise ValueError('PACKAGED_SCRATCH_TARGET')
            checks.append('packaged_external_scratch_accepted')
        finally:mod.HERE=original_here
        # The marker would prove execution if digest authentication were bypassed.
        marker=root/'executed.txt';mutant=root/'mutated-angular.py'
        mutant.write_text('from pathlib import Path\nPath('+repr(str(marker))+').write_text("executed")\n')
        reject('angular_mutant_rejected_before_execution',lambda:mod.load_angular(mutant),'ANGULAR_SOURCE_SHA256_MISMATCH')
        if marker.exists():raise ValueError('ANGULAR_MUTANT_EXECUTED')
        mod.authenticated_angular_bytes(mod.DEFAULT)
        checks.append('immutable_angular_pin_accepts_copy')
        if existing.read_text()!='preserve\n':raise ValueError('EXISTING_OUTPUT_MODIFIED')
    report={'status':'PASS_SCREEN_REDUCTION_OUTPUT_AND_PREIMPORT_GUARDS','checks':checks,'scientific_replay':False,'angular_source_sha256':mod.ANGULAR_SHA256}
    with mod.validate_output(args.output).open('x',encoding='utf-8') as stream:stream.write(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':main()
