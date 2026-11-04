#!/usr/bin/env python3
"""Portable isolated replay and exact rejection-witness controls."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

COMPONENT = Path(__file__).resolve().parent
SOURCE_PATHS = (
    Path('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json'),
    Path('research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json'),
    Path('research/HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/numeric/references/ORIGINAL_REPORT.json'),
)
CONTRACT = {'success_subprocesses': 2, 'actual_negative_subprocesses': 12, 'inprocess_checks': 2, 'total_controls': 16}

class Rejected(Exception):
    pass

class Parser(argparse.ArgumentParser):
    def error(self, message):
        raise Rejected('invalid arguments')

def need(ok, message):
    if not ok:
        raise Rejected(message)

def under(path, parent):
    return path == parent or parent in path.parents

def validate_output(out, repo):
    need(out.is_absolute(), 'absolute output path')
    need(not any(part.is_symlink() for part in (out, *out.parents)), 'output path contains symlink')
    need(out == out.resolve(), 'normalized output path')
    need(out.parent.is_dir(), 'output parent exists')
    need(not out.exists() and not out.is_symlink(), 'output absent')
    need(not under(out, COMPONENT) and not under(out, repo.resolve()), 'output outside code and repository')

def main():
    parser = Parser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    need(args.repo_root.is_absolute(), 'absolute repository path')
    validate_output(args.output, args.repo_root)
    baseline = (COMPONENT / 'references/NORMAL.json').read_bytes()
    optimized = (COMPONENT / 'references/OPTIMIZED.json').read_bytes()
    need(baseline == optimized, 'stored normal optimized equality')
    base_json = json.loads(baseline)
    need(base_json['verification']['finite_checks'] == 38, 'frozen 38 scientific checks')
    inprocess = [
        {'name': 'stored_normal_optimized_identical', 'pass': True},
        {'name': 'frozen_scientific_check_count_38', 'pass': True},
    ]
    successes, negatives = [], []
    with tempfile.TemporaryDirectory(prefix='hdblast-bounds-controls-') as directory:
        tmp = Path(directory)
        def child(name, out, witness=None, options=(), repo=None, optimized_mode=False, existing=None):
            root = args.repo_root if repo is None else repo
            command = [sys.executable, '-I', '-S', '-B']
            if optimized_mode:
                command.append('-O')
            command.extend([str(COMPONENT/'produce.py'), '--repo-root', str(root), '--output', str(out), *map(str, options)])
            result = subprocess.run(command, capture_output=True, text=True, cwd=tmp)
            need(result.stderr == '', name + ': empty stderr')
            need(len(result.stdout.splitlines()) == 1, name + ': single JSON stdout')
            try:
                summary = json.loads(result.stdout)
            except ValueError:
                raise Rejected(name + ': JSON stdout') from None
            need(isinstance(summary, dict), name + ': JSON object stdout')
            if witness is None:
                need(result.returncode == 0, name + ': success exit')
                need(summary.get('status') == 'PASS_CONDITIONAL_PHYSICAL_INITIAL_PROFILE_CLASS_ROBUSTNESS', name + ': success status')
                need(summary.get('finite_checks') == 38, name + ': 38 scientific checks')
                need(out.read_bytes() == baseline, name + ': exact baseline bytes')
                need(summary.get('output_sha256') == hashlib.sha256(baseline).hexdigest(), name + ': output hash')
                successes.append({'name': name, 'pass': True, 'exit_code': 0, 'finite_checks': 38,
                                  'interpreter_flags': ['-I', '-S', '-B'] + (['-O'] if optimized_mode else [])})
            else:
                need(result.returncode == 2, name + ': rejection exit 2')
                need(summary == {'status': 'REJECTED_INCONCLUSIVE', 'reason': witness}, name + ': exact rejection witness')
                if existing is None:
                    need(not out.exists() and not out.is_symlink(), name + ': no output artifact')
                else:
                    need(out.read_bytes() == existing, name + ': existing bytes preserved')
                negatives.append({'name': name, 'pass': True, 'exit_code': 2, 'reason': witness,
                                  'single_JSON_stdout': True, 'empty_stderr': True,
                                  'output_preserved' if existing is not None else 'no_output_artifact': True,
                                  'interpreter_flags': ['-I', '-S', '-B'] + (['-O'] if optimized_mode else [])})
        child('normal_exact_replay', tmp/'normal.json')
        child('optimized_exact_replay', tmp/'optimized.json', optimized_mode=True)
        marker = b'preserved existing output\n'
        prior = tmp/'existing.json'
        prior.write_bytes(marker)
        child('existing_output_rejected', prior, 'output absent', existing=marker)
        child('relative_output_rejected', Path('relative.json'), 'absolute output path')
        child('parent_escape_rejected', tmp/'absent'/'..'/'escape.json', 'normalized output path')
        real = tmp/'real'
        real.mkdir()
        link = tmp/'linked'
        link.symlink_to(real, target_is_directory=True)
        child('symlink_parent_rejected', link/'report.json', 'output path contains symlink')
        child('inside_component_rejected', COMPONENT/'control-never-created.json', 'output outside code and repository')
        child('inside_repository_rejected', args.repo_root/'control-never-created.json', 'output outside code and repository')
        child('missing_output_parent_rejected', tmp/'missing'/'report.json', 'output parent exists')
        child('missing_source_root_rejected', tmp/'missing-root-out.json', 'source unreadable: static', repo=tmp/'missing-root')
        bad_design = tmp/'bad-design.json'
        bad_design.write_bytes((COMPONENT/'FROZEN_DESIGN.json').read_bytes()+b' ')
        child('mutated_design_rejected', tmp/'design-out.json', 'frozen design byte pin', options=('--design',bad_design))
        child('optimized_mutated_design_rejected', tmp/'design-O-out.json', 'frozen design byte pin', options=('--design',bad_design), optimized_mode=True)
        bad_follow = tmp/'bad-follow.json'
        bad_follow.write_bytes(b'{}\n')
        child('mutated_local_reference_rejected', tmp/'follow-out.json', 'follow-on byte pin', options=('--follow-on',bad_follow))
        bad_root = tmp/'bad-root'
        for relative in SOURCE_PATHS:
            target = bad_root/relative
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(args.repo_root/relative,target)
        with (bad_root/SOURCE_PATHS[0]).open('ab') as stream:
            stream.write(b' ')
        child('mutated_static_source_rejected', tmp/'static-out.json', 'static byte pin', repo=bad_root)
    counts = {'success_subprocesses':len(successes), 'actual_negative_subprocesses':len(negatives),
              'inprocess_checks':len(inprocess), 'total_controls':len(successes)+len(negatives)+len(inprocess)}
    need(counts == CONTRACT, 'fixed control count contract')
    report = {'status':'PASS_PORTABLE_ISOLATED_REPLAY_AND_REJECTION_CONTROLS', 'count_contract':CONTRACT,
              'counts':counts, 'successes':successes, 'actual_negatives':negatives, 'inprocess':inprocess,
              'normal_optimized_identical':True, 'report_sha256':hashlib.sha256(baseline).hexdigest(),
              'oracle_calls':0,'PDE_solves':0,'new_bump_queries':0}
    payload = (json.dumps(report,indent=2,sort_keys=True)+'\n').encode()
    with args.output.open('xb') as stream:
        stream.write(payload)
    print(json.dumps({'status':report['status'],'counts':counts,'output_sha256':hashlib.sha256(payload).hexdigest()},sort_keys=True))

if __name__ == '__main__':
    try:
        main()
    except Rejected as exc:
        print(json.dumps({'status':'REJECTED_INCONCLUSIVE','reason':str(exc)},sort_keys=True))
        sys.exit(2)
    except (OSError,ValueError,KeyError,TypeError):
        print(json.dumps({'status':'REJECTED_INCONCLUSIVE','reason':'control input or filesystem error'},sort_keys=True))
        sys.exit(2)
