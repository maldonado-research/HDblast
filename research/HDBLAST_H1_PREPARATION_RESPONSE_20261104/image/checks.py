#!/usr/bin/env python3
"""One canonical normal/optimized pair and four exact CLI rejections per mode."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
OUTPUT_BASENAME = 'IMAGE_CONTROLS.json'
ARTIFACT_DIRECTORY = 'IMAGE_CONTROLS.artifacts'
CERTIFICATE_SHA = 'af2e1c9e54a04e603a6b721a71c773200bbbf67f6dd0f32fedd924613672dd75'


def need(ok, witness):
    if not ok:
        raise ValueError(witness)


class Parser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError('invalid arguments')


def main():
    parser = Parser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output
    need(out.is_absolute() and out == out.resolve(), 'absolute normalized output path')
    need(not any(p.is_symlink() for p in [out, *out.parents]), 'output path contains symlink')
    need(out.parent.is_dir() and not out.exists() and not out.is_symlink(), 'fresh output with existing parent')
    need(not out.is_relative_to(HERE), 'output outside component')
    need(out.name == OUTPUT_BASENAME, 'fixed output basename IMAGE_CONTROLS.json')
    artifacts = out.parent / ARTIFACT_DIRECTORY
    need(not artifacts.exists() and not artifacts.is_symlink(), 'fresh canonical artifact directory')
    expected = (HERE / 'references/NORMAL.json').read_bytes()
    need(hashlib.sha256(expected).hexdigest() == CERTIFICATE_SHA and
         expected == (HERE / 'references/OPTIMIZED.json').read_bytes(), 'frozen canonical certificate bytes')
    artifacts.mkdir()
    successes, negatives = [], []
    with tempfile.TemporaryDirectory(prefix='hdblast-image-controls-') as folder:
        tmp = Path(folder)

        def execute(case, target, optimized, witness=None, prior=None):
            flags = ['-I', '-S', '-B'] + (['-O'] if optimized else [])
            command = [sys.executable, *flags, str(HERE / 'verify_algebra.py'), '--output', str(target)]
            child = subprocess.run(command, cwd=tmp, capture_output=True, text=True)
            need(child.stderr == '' and len(child.stdout.splitlines()) == 1, case + ': exact stdout/stderr framing')
            document = json.loads(child.stdout)
            actual_path = target if target.is_absolute() else tmp / target
            if witness is None:
                need(child.returncode == 0 and document == {
                    'status': 'PASS_EXACT_FORMAL_ALGEBRA_AND_CONTINUOUS_BOUND_CONSTANTS',
                    'sha256': CERTIFICATE_SHA}, case + ': exact success witness')
                need(actual_path.read_bytes() == expected, case + ': exact canonical bytes')
                successes.append({'case': case, 'optimized': optimized, 'interpreter_flags': flags,
                    'exit_code': 0, 'stdout_complete_JSON': document, 'stderr_empty': True,
                    'artifact': ARTIFACT_DIRECTORY + '/' + target.name, 'sha256': CERTIFICATE_SHA})
            else:
                need(child.returncode == 2 and document == {'status': 'REJECTED', 'witness': witness},
                     case + ': exact rejection witness')
                if prior is None:
                    need(not actual_path.exists() and not actual_path.is_symlink(), case + ': no output artifact')
                else:
                    need(actual_path.read_bytes() == prior, case + ': existing output preserved')
                negatives.append({'case': case, 'optimized': optimized, 'interpreter_flags': flags,
                    'exit_code': 2, 'stdout_complete_JSON': document, 'stderr_empty': True,
                    'no_new_output_artifact': True, 'existing_bytes_preserved': prior is not None})

        real = tmp / 'real'
        real.mkdir()
        link = tmp / 'linked'
        link.symlink_to(real, target_is_directory=True)
        for optimized, mode in [(False, 'NORMAL'), (True, 'OPTIMIZED')]:
            execute('canonical', artifacts / (mode + '.json'), optimized)
            execute('relative-output', Path('relative-' + mode + '.json'), optimized, 'absolute output path')
            prior_path = tmp / ('existing-' + mode + '.json')
            prior = b'preexisting image control bytes\n'
            prior_path.write_bytes(prior)
            execute('existing-output', prior_path, optimized, 'output absent', prior)
            execute('inside-component', HERE / ('never-created-' + mode + '.json'), optimized, 'output outside component')
            execute('symlink-parent', link / ('never-created-' + mode + '.json'), optimized, 'output path contains symlink')
    need(len(successes) == 2 and len(negatives) == 8 and
         (artifacts / 'NORMAL.json').read_bytes() == (artifacts / 'OPTIMIZED.json').read_bytes(),
         'exact canonical and negative roster')
    report = {'status': 'PASS_PORTABLE_ACTUAL_IMAGE_ALGEBRA_CONTROLS',
        'canonical_artifact_directory': ARTIFACT_DIRECTORY,
        'canonical_child_processes': 2, 'canonical_pairs_byte_identical': True,
        'actual_CLI_negative_control_count': 8, 'negative_cases': 4, 'total_child_processes': 10,
        'canonical_children': successes, 'actual_CLI_negative_controls': negatives,
        'scope': {'actual_background_queries': 0, 'new_oracle_queries': 0, 'new_PDE_solves': 0,
                  'numerical_bump_queries': 0, 'new_observable_sign_claim': False}}
    payload = (json.dumps(report, indent=2, sort_keys=True) + '\n').encode()
    with out.open('xb') as stream:
        stream.write(payload)
    print(json.dumps({'status': report['status'], 'actual_CLI_negative_control_count': 8,
                      'output_sha256': hashlib.sha256(payload).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(json.dumps({'status': 'REJECTED', 'witness': str(exc)}, sort_keys=True))
        sys.exit(2)
