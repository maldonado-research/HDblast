#!/usr/bin/env python3
"""Authenticate and replay the fixed-budget shell precision capsule."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tempfile
import zipfile

sys.dont_write_bytecode = True


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def authenticate(root, pin):
    need(not root.is_symlink(), 'ROOT_SYMLINK')
    manifest = root/'MANIFEST.json'
    need(not manifest.is_symlink() and sha(manifest) == pin, 'MANIFEST_PIN')
    data = json.loads(manifest.read_text())
    need(data['schema'] == 'HDBLAST_SHELL_PRECISION_MANIFEST_V1', 'MANIFEST_SCHEMA')
    expected = {'MANIFEST.json'}
    for row in data['files']:
        rel = row['path']
        path = PurePosixPath(rel)
        need(rel == path.as_posix() and rel not in expected and not path.is_absolute()
             and '..' not in path.parts and '.' not in path.parts and '\\' not in rel,
             'MANIFEST_MEMBER_PATH')
        expected.add(rel)
        member = root/rel
        need(member.is_file() and not member.is_symlink(), 'MANIFEST_MEMBER_TYPE')
        need(member.stat().st_size == row['bytes'] and sha(member) == row['sha256'],
             'MANIFEST_MEMBER_BYTES')
    actual = {'MANIFEST.json'}
    for path in root.rglob('*'):
        need(not path.is_symlink(), 'TREE_SYMLINK')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
    need(actual == expected, 'MANIFEST_COMPLETE_CLOSURE')
    return {'manifest_sha256': pin, 'files_including_manifest': len(expected)}


def archive_check(root, archive, pin):
    need(pin and sha(archive) == pin, 'ARCHIVE_PIN')
    with zipfile.ZipFile(archive) as packed:
        need(packed.testzip() is None, 'ARCHIVE_CRC')
        names = packed.namelist()
        prefix = root.name+'/'
        expected = {prefix+p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
        need(len(names) == len(set(names)) and set(names) == expected, 'ARCHIVE_EXACT_MEMBERS')
        for name in names:
            need(packed.read(name) == (root/name[len(prefix):]).read_bytes(), 'ARCHIVE_MEMBER_BYTES')
    return {'sha256': pin, 'members': len(names), 'all_members_match_source': True}


def corruption_controls(root, pin, work):
    """Mutate a disposable copy, never any source or dependency capsule."""
    copy = work/'corruption-copy'
    shutil.copytree(root, copy)
    manifest = copy/'MANIFEST.json'
    original = manifest.read_bytes()
    rows = json.loads(original)['files']
    target = copy/rows[0]['path']
    saved = target.read_bytes()
    passed = []
    def reject(name, supplied_pin=pin):
        try:
            authenticate(copy, supplied_pin)
        except (ValueError, OSError, KeyError):
            passed.append(name)
        else:
            raise ValueError('CORRUPTION_ACCEPTED:' + name)
    reject('wrong_external_manifest_pin', '0'*64)
    target.write_bytes(saved+b'!')
    reject('changed_member_bytes')
    target.write_bytes(saved)
    extra = copy/'undeclared-file'
    extra.write_bytes(b'undeclared')
    reject('undeclared_member')
    extra.unlink()
    target.unlink()
    target.symlink_to(root/rows[0]['path'])
    reject('member_symlink')
    target.unlink()
    target.write_bytes(saved)
    for name, change in (
        ('duplicate_manifest_member', lambda d: d['files'].append(d['files'][0].copy())),
        ('parent_traversal', lambda d: d['files'][0].update(path='../outside')),
        ('absolute_member_path', lambda d: d['files'][0].update(path='/outside')),
    ):
        changed = json.loads(original)
        change(changed)
        manifest.write_text(json.dumps(changed))
        reject(name, sha(manifest))
        manifest.write_bytes(original)
    need(authenticate(copy, pin) == authenticate(root, pin), 'RESTORED_COPY_AUTHENTICATION')
    shutil.rmtree(copy)
    return passed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--archive-sha256')
    parser.add_argument('--certificates-only', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    out = args.output_dir.resolve()
    need(not out.exists() and root not in out.parents, 'FRESH_EXTERNAL_OUTPUT_DIRECTORY')
    before = authenticate(root, args.manifest_sha256)  # Before every project import or execution.
    out.mkdir(parents=True)
    archive = archive_check(root, args.archive, args.archive_sha256) if args.archive else None
    controls = corruption_controls(root, args.manifest_sha256, out)
    records = []
    def execute(script, extra, destination, reference, optimized=False):
        command = [sys.executable, '-I', '-S', '-B']+(['-O'] if optimized else [])
        command += [str(root/script), '--repo-root', str(args.repo_root.resolve())]
        command += [str(v) for v in extra]+['--output', str(out/destination)]
        result = subprocess.run(command, capture_output=True, text=True)
        log = out/(destination+'.log')
        log.write_text(result.stdout+result.stderr)
        need(result.returncode == 0, 'REPLAY_FAILED:' + script + ':' + result.stderr[-2500:])
        need((out/destination).read_bytes() == (root/reference).read_bytes(), 'GOLDEN_BYTES:' + script)
        records.append({'script': script, 'optimized': optimized, 'output': destination,
                        'sha256': sha(out/destination), 'golden_match': True})
        print(json.dumps(records[-1], sort_keys=True), flush=True)
    for optimized in (False, True):
        suffix = 'optimized' if optimized else 'normal'
        if not args.certificates_only:
            for script, golden, label in (
                ('numeric/component_ledger.py', 'numeric/COMPONENT_LEDGER_TEN_INTEGRATED_NORMAL.json', 'component'),
                ('numeric/second_feedback_ledger.py', 'numeric/SECOND_FEEDBACK_NORMAL.json', 'second-feedback'),
                ('balanced/exact_boundary.py', 'balanced/EXACT_BOUNDARY_NORMAL.json', 'boundary'),
                ('balanced/b2/exact_b2.py', 'balanced/b2/EXACT_B2_NORMAL.json', 'b2'),
            ):
                execute(script, [], label+'-'+suffix+'.json', golden, optimized)
        execute('balanced/check_boundary_report.py', [], 'boundary-review-'+suffix+'.json',
                'balanced/INDEPENDENT_BOUNDARY_REVIEW.json', optimized)
        review = out/('b2-review-'+suffix+'.json')
        execute('balanced/b2/check_b2.py', ['--review-output', review], 'b2-controls-'+suffix+'.json',
                'balanced/b2/B2_POLYNOMIAL_CONTROLS_FINAL.json', optimized)
        need(review.read_bytes() == (root/'balanced/b2/B2_INDEPENDENT_REVIEW.json').read_bytes(),
             'B2_INDEPENDENT_REVIEW_GOLDEN')
        execute('numeric/check_second_feedback.py', [], 'feedback-review-'+suffix+'.json',
                'numeric/SECOND_FEEDBACK_FINAL_REVIEW_NORMAL.json', optimized)
        execute('analysis/join_results.py', [], 'joined-'+suffix+'.json', 'results/JOINED_NORMAL.json', optimized)
        execute('review-volume/check_volume_join.py', [], 'volume-review-'+suffix+'.json',
                'review-volume/FINAL_VOLUME_JOIN_REVIEW_NORMAL.json', optimized)
    need(authenticate(root, args.manifest_sha256) == before, 'CAPSULE_CHANGED_DURING_REPLAY')
    receipt = {'schema': 'HDBLAST_SHELL_PRECISION_REPRODUCTION_V1', 'status': 'PASS',
               'python_version': sys.version, 'authentication': before, 'archive': archive,
               'full_bounded_numerics': not args.certificates_only,
               'corruption_controls_rejected': controls, 'runs': records,
               'source_unchanged_after_replay': True,
               'limits': ['Internal computational verification; not proof-assistant formalization or external peer review.',
                          'Conditional coefficient/kinematical-jet results, not a nonlinear physical solution or origin evidence.']}
    (out/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt, sort_keys=True, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'runs': len(records), 'corruption_controls': len(controls)}))


if __name__ == '__main__':
    main()
