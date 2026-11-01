#!/usr/bin/env python3
"""Standard-library rejection controls for the portable checker entrypoint."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    repo_root = args.repo_root.resolve(strict=True)
    if args.output is not None:
        output_path = args.output.resolve()
        if output_path.is_relative_to(ROOT) or output_path.is_relative_to(repo_root):
            raise ValueError('--output must be outside the code and source roots')
    if args.output is not None and args.output.exists():
        raise FileExistsError('Refusing to replace existing control report')
    pins = json.loads((ROOT/'SOURCE_PINS.json').read_text())
    source_bytes = {}
    for entry in pins['source_files']:
        raw = (repo_root/entry['path']).read_bytes()
        if len(raw) != entry['bytes'] or hashlib.sha256(raw).hexdigest() != entry['sha256']:
            raise ValueError('Changed control input: '+entry['path'])
        source_bytes[entry['path']] = raw
    controls = []
    flags = ['-I','-B'] + (['-O'] if sys.flags.optimize else [])

    def rejected(name, command, expected_final, artifact, preserved_bytes=None):
        result = subprocess.run(command, text=True, capture_output=True, check=False)
        lines = result.stderr.splitlines()
        observed_final = lines[-1] if lines else ''
        if result.returncode != 1 or observed_final != expected_final or result.stdout != '':
            raise ValueError(name+' did not reject as expected: '+result.stdout+result.stderr)
        if preserved_bytes is None:
            if artifact.exists() or artifact.is_symlink():
                raise ValueError(name+' created a forbidden output artifact')
            artifact_result = {'path':str(artifact.relative_to(work)),
                               'expected_state':'absent','observed_state':'absent'}
        else:
            if not artifact.is_file() or artifact.read_bytes() != preserved_bytes:
                raise ValueError(name+' changed the existing output artifact')
            digest = hashlib.sha256(preserved_bytes).hexdigest()
            artifact_result = {'path':str(artifact.relative_to(work)),
                               'expected_state':'preserved','observed_state':'preserved',
                               'expected_sha256':digest,'observed_sha256':digest}
        # Normalize dynamic paths only after requiring the complete raw final
        # exception line to equal its case-specific expected witness.
        controls.append({'case':name,'expected_final_witness':expected_final.replace(str(work),'<TEMP>'),
                         'observed_final_witness':observed_final.replace(str(work),'<TEMP>'),
                         'returncode':result.returncode,'stdout_empty':True,
                         'artifact':artifact_result})

    with tempfile.TemporaryDirectory(prefix='weighted-mean-io-') as temp:
        work = Path(temp)
        copied_repo = work/'repo'
        for path,raw in source_bytes.items():
            target = copied_repo/path
            target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes(raw)
        package = work/'package'
        (package/'checks').mkdir(parents=True)
        script = package/'checks/check_weighted_mean.py'
        shutil.copyfile(ROOT/'checks/check_weighted_mean.py',script)
        pin_path = package/'SOURCE_PINS.json'
        pin_bytes = (ROOT/'SOURCE_PINS.json').read_bytes()
        pin_path.write_bytes(pin_bytes)
        selected = next(iter(source_bytes))
        target = copied_repo/selected
        target.write_bytes(source_bytes[selected]+b'\n')
        output = work/'before-import.json'
        rejected('source_tamper_before_SymPy_import',
                 [sys.executable,*flags,'-S',str(script),'--repo-root',str(copied_repo),'--output',str(output)],
                 'ValueError: Changed source: '+selected,output)
        target.write_bytes(source_bytes[selected])
        report = work/'existing.json'
        report.write_text('preserved\n')
        rejected('existing_report_preserved',
                 [sys.executable,*flags,'-S',str(script),'--repo-root',str(copied_repo),'--output',str(report)],
                 'FileExistsError: Refusing to replace existing report: '+str(report),report,b'preserved\n')
        for name,output in [('fresh_output_inside_code_root',package/'fresh.json'),
                            ('fresh_output_inside_source_root',copied_repo/'fresh.json')]:
            rejected(name,
                     [sys.executable,*flags,'-S',str(script),'--repo-root',str(copied_repo),'--output',str(output)],
                     'ValueError: --output must be outside the code and source roots',output)
        driver = '''import importlib.util,sys
from pathlib import Path
spec=importlib.util.spec_from_file_location('producer',sys.argv[1])
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
repo=Path(sys.argv[2]);target=Path(sys.argv[3]);output=Path(sys.argv[4])
def mutate(repo_root,pins):
    target.write_bytes(target.read_bytes()+b'\\n')
    return {}
module.evaluate=mutate
sys.argv=['producer','--repo-root',str(repo),'--output',str(output)]
module.main()
'''
        output = work/'after-evaluation.json'
        rejected('source_tamper_after_evaluation',
                 [sys.executable,*flags,'-S','-c',driver,str(script),str(copied_repo),str(target),str(output)],
                 'ValueError: Changed source: '+selected,output)
        target.write_bytes(source_bytes[selected])
        output = work/'pin-tamper.json'
        rejected('pin_file_tamper_after_evaluation',
                 [sys.executable,*flags,'-S','-c',driver,str(script),str(copied_repo),str(pin_path),str(output)],
                 'ValueError: SOURCE_PINS.json changed during execution',output)
        pin_path.write_bytes(pin_bytes)
        original_script = script.read_text()
        needle = '+ 3*a*D(Ar)/2'
        if original_script.count(needle) != 1:
            raise ValueError('Mathematical mutation target is not unique')
        script.write_text(original_script.replace(needle,'- 3*a*D(Ar)/2'))
        output = work/'radial-sign-mutation.json'
        rejected('radial_metric_sign_mutation',
                 [sys.executable,*flags,str(script),'--repo-root',str(copied_repo),'--output',str(output)],
                 'ValueError: radial_metric_operator: 3*a(y)*Derivative(A(y), y)',output)
    # Verify the actual repository inputs have remained unchanged.
    for path,raw in source_bytes.items():
        if (repo_root/path).read_bytes() != raw:
            raise ValueError('Source changed while controls ran: '+path)
    result = {'status':'PASS','control_count':len(controls),'controls':controls,
              'scope':'Packaging and fail-closed controls; separate from 34 scientific checks and 8 formula controls.'}
    raw = json.dumps(result,indent=2)+'\n'
    if args.output is None:
        print(raw,end='')
    else:
        with args.output.open('x',encoding='utf-8') as stream:
            stream.write(raw)
        print(json.dumps({'status':result['status'],'control_count':len(controls)}))


if __name__ == '__main__':
    main()
