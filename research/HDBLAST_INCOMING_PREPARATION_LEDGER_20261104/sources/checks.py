#!/usr/bin/env python3
"""Two retained canonical children and eight strict rejection subprocesses."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def need(ok, reason):
    if not ok:
        raise ValueError(reason)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def regular(path):
    need(not any(p.is_symlink() for p in (path,*path.parents)), 'SYMLINK_PATH')
    need(path.is_file(), 'REGULAR_INPUT_REQUIRED')
    return path.read_bytes()


def auth():
    pins = json.loads(regular(HERE/'LOCAL_PINS.json'))
    for rel, pin in pins['files'].items():
        raw = regular(HERE/rel)
        need(digest(raw) == pin['sha256'] and len(raw) == pin['bytes'],
             'LOCAL_SOURCE_AUTHENTICATION:' + rel)
    return pins


def launch(source, repo, output, optimized):
    flags = ['-I','-S','-B'] + (['-O'] if optimized else [])
    command = [sys.executable,*flags,str(source),'--repo-root',str(repo),'--output',str(output)]
    process = subprocess.run(command,capture_output=True,text=True,check=False)
    need(process.stderr == '', 'CHILD_STDERR_NOT_EMPTY')
    # A complete single JSON object is mandatory; json.loads rejects any
    # additional stdout record. Trailing newline is allowed and expected.
    need(len(process.stdout.splitlines()) == 1, 'CHILD_STDOUT_SINGLE_LINE_REQUIRED')
    answer = json.loads(process.stdout)
    need(isinstance(answer,dict), 'CHILD_STDOUT_JSON_OBJECT_REQUIRED')
    return process.returncode, answer, flags


def execute(args):
    before = auth()
    need(args.repo_root.is_dir(), 'REPOSITORY_REQUIRED')
    need(not any(p.is_symlink() for p in (args.repo_root,*args.repo_root.parents)), 'SYMLINK_PATH')
    repo = args.repo_root.resolve()
    need(args.output.is_absolute(), 'ABSOLUTE_OUTPUT_REQUIRED')
    need(not any(p.is_symlink() for p in (args.output,*args.output.parents)), 'SYMLINK_PATH')
    output = args.output.resolve()
    artifacts = output.with_name(output.stem+'.artifacts')
    need(not output.exists() and not artifacts.exists(), 'OUTPUT_EXISTS')
    need(not output.is_relative_to(HERE) and not output.is_relative_to(repo)
         and not artifacts.is_relative_to(HERE) and not artifacts.is_relative_to(repo),
         'OUTPUT_INSIDE_SOURCE')
    need(output.parent.is_dir(), 'OUTPUT_PARENT_MISSING')
    artifacts.mkdir()
    rows = []
    canonical = []
    for optimized, label in ((False,'NORMAL'),(True,'OPTIMIZED')):
        target = artifacts/(label+'.json')
        code, answer, flags = launch(HERE/'source_budget.py',repo,target,optimized)
        raw = regular(target)
        expected = regular(HERE/'references'/f'{label}.json')
        need(code == 0 and answer == {'status':'PASS_EXACT_ZERO_ORACLE_SOURCE_LEDGER',
                                    'output_sha256':digest(raw)}, 'CANONICAL_CHILD_PROTOCOL')
        need(raw == expected, 'CANONICAL_REFERENCE_MISMATCH')
        canonical.append(raw)
        rows.append({'mode':label,'kind':'positive_subprocess','flags':flags,'exit_code':code,
                     'single_json_stdout':True,'stderr_empty':True,
                     'retained_report':label+'.json','sha256':digest(raw)})

        for case, witness in (
                ('tampered_source_pin','LOCAL_PIN_AUTHENTICATION'),
                ('altered_design_scale','DESIGN_AUTHENTICATION'),
                ('occupied_output','OUTPUT_EXISTS'),
                ('relative_output','ABSOLUTE_OUTPUT_REQUIRED')):
            tag = label+'-'+case
            source = HERE/'source_budget.py'
            target = artifacts/(tag+'.json')
            preserved = None
            if case in ('tampered_source_pin','altered_design_scale'):
                mutant = artifacts/(tag+'.component')
                mutant.mkdir()
                # This is a real altered component read by a separate child.
                for name in ('source_budget.py','SOURCE_PINS.json','SOURCE_CLOSURE_INPUTS.json',
                             'PACKET_SOURCE_JOIN.json','FROZEN_DESIGN.json'):
                    shutil.copyfile(HERE/name,mutant/name)
                if case == 'tampered_source_pin':
                    path = mutant/'SOURCE_PINS.json'
                    data = json.loads(path.read_text())
                    first = next(iter(data['files']))
                    data['files'][first]['sha256'] = '0'*64
                else:
                    path = mutant/'FROZEN_DESIGN.json'
                    data = json.loads(path.read_text())
                    data['ell'] = '1/100000'
                path.write_text(json.dumps(data,sort_keys=True,indent=2)+'\n')
                source = mutant/'source_budget.py'
            elif case == 'occupied_output':
                preserved = b'occupied output must remain byte-identical\n'
                target.write_bytes(preserved)
            else:
                target = Path('sources-relative-output-rejection.json')
                need(not target.exists(), 'RELATIVE_CONTROL_TARGET_ALREADY_EXISTS')

            code, answer, flags = launch(source,repo,target,optimized)
            need(code == 2 and answer == {'status':'REJECTED','witness':witness},
                 'NEGATIVE_CHILD_PROTOCOL:' + tag)
            if preserved is None:
                need(not target.exists(), 'REJECTED_CHILD_CREATED_REPORT:' + tag)
                semantics = 'no_output_artifact'
            else:
                need(target.read_bytes() == preserved, 'OCCUPIED_OUTPUT_CHANGED:' + tag)
                semantics = 'existing_output_preserved_exactly'
            rows.append({'mode':label,'kind':'negative_subprocess','case':case,'flags':flags,
                         'exit_code':code,'single_json_stdout':True,'stderr_empty':True,
                         'witness':witness,'artifact_semantics':semantics})

    need(canonical[0] == canonical[1], 'NORMAL_OPTIMIZED_BYTE_IDENTITY')
    need(auth() == before, 'LOCAL_AUTHENTICATION_AFTER')
    report = {
        'schema':'HDBLAST_INCOMING_SOURCES_STRICT_CONTROLS_V1',
        'status':'PASS_TWO_CANONICAL_AND_EIGHT_REJECTED_CHILDREN',
        'canonical_report_pairs_reproduced':1,
        'positive_control_subprocesses_executed':2,
        'actual_rejected_subprocesses':8,
        'inprocess_mutant_evaluations':0,
        'canonical_sha256':digest(canonical[0]),
        'retained_canonical_artifacts':['NORMAL.json','OPTIMIZED.json'],
        'artifact_directory_contract':'Sibling directory OUTPUT_STEM.artifacts',
        'scientific_activity':{'project_imports':0,'source_oracle_calls':0,'background_queries':0,
                               'bump_queries':0,'PDE_solves':0},
        'rows':rows,
        'scope':'Two actual exact-postprocessor children plus eight actual protocol/input rejection children; inherited source and PDE proofs are not rerun.'
    }
    raw = (json.dumps(report,sort_keys=True,indent=2)+'\n').encode()
    with output.open('xb') as stream:
        stream.write(raw)
    print(json.dumps({'status':report['status'],'output_sha256':digest(raw)},sort_keys=True))


def main():
    class Parser(argparse.ArgumentParser):
        def error(self,message):
            raise ValueError('CLI_ARGUMENTS')
    parser = Parser(add_help=False)
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    execute(parser.parse_args())


if __name__ == '__main__':
    try:
        main()
    except (ValueError,OSError,ArithmeticError,KeyError,TypeError,IndexError,AttributeError,UnicodeError) as exc:
        witness = str(exc) if isinstance(exc,ValueError) else type(exc).__name__
        print(json.dumps({'status':'REJECTED','witness':witness},sort_keys=True))
        raise SystemExit(2)
