#!/usr/bin/env python3
"""Reproduce exact symmetry goldens and rejection receipts in a fresh output.

No project or source campaign imports. The requested directory must be fresh,
external to this source and any packaged screen-curvature capsule, and have a
real existing parent. Scientific reports contain no host-specific paths.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
CHECKER = HERE / 'check_scalar_magnetic_symmetry.py'

def need(ok, label):
    if not ok:
        raise ValueError(label)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(output, optimized=False, mutation=None):
    command = [sys.executable, '-B'] + (['-O'] if optimized else [])
    command += [str(CHECKER), '--output', str(output)]
    if mutation:
        command += ['--mutation', mutation]
    return subprocess.run(command, capture_output=True, text=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('independent_symmetry_guard', CHECKER)
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    destination = checker.guard_output(args.output_dir)
    source_before = {p.name: sha(p) for p in (Path(__file__), CHECKER)}
    destination.mkdir() # exclusive creation after preflight, no parent creation
    normal = destination / 'SYMMETRY_NORMAL.json'
    optimized = destination / 'SYMMETRY_OPTIMIZED.json'
    for output, is_optimized in ((normal, False), (optimized, True)):
        result = run(output, is_optimized)
        need(result.returncode == 0, 'POSITIVE_SUBPROCESS')
    need(normal.read_bytes() == optimized.read_bytes(), 'NORMAL_OPTIMIZED_IDENTITY')
    rejections = []
    scientific = {'magnetic_sign':'COMPLETE_MAGNETIC_TT_IDENTITY',
                  'electric_half':'ELECTRIC_TT_PROJECTION',
                  'reflection_slots':'REFLECTION_EXCHANGE_INCLUDING_SLOTS'}
    for is_optimized in (False, True):
        for mutation, marker in scientific.items():
            output = destination / ('reject_' + mutation + ('_optimized' if is_optimized else '') + '.json')
            result = run(output, is_optimized, mutation)
            rejected = result.returncode != 0 and marker in result.stderr
            need(rejected and not output.exists(), 'SCIENTIFIC_MUTATION_REJECTION')
            rejections.append({'optimized':is_optimized, 'mutation':mutation,
                               'expected_marker':marker, 'rejected':True, 'output_absent':True})
    output_guards = []
    with tempfile.TemporaryDirectory(prefix='screen-symmetry-guards-') as temporary:
        scratch = Path(temporary)
        sentinel = scratch / 'sentinel.txt'
        sentinel.write_text('retained sentinel\n', encoding='utf-8')
        symlink = scratch / 'linked.json'
        symlink.symlink_to(sentinel)
        linked_parent = scratch / 'linked-parent'
        linked_parent.symlink_to(scratch, target_is_directory=True)
        package = scratch / 'research' / checker.CAPSULE
        package.mkdir(parents=True)
        own_tree = HERE / 'forbidden_guard_probe.json'
        need(not own_tree.exists(), 'FRESH_SOURCE_GUARD_PROBE')
        probes = [('existing', sentinel, 'OUTPUT_EXISTS'),
                  ('leaf_symlink', symlink, 'OUTPUT_SYMLINK'),
                  ('ancestor_symlink', linked_parent / 'new.json', 'OUTPUT_ANCESTOR_SYMLINK'),
                  ('own_source', own_tree, 'OUTPUT_INSIDE_SOURCE_TREE'),
                  ('packaged_capsule', package / 'new.json', 'OUTPUT_INSIDE_PACKAGED_CAPSULE'),
                  ('missing_parent', scratch / 'missing' / 'new.json', 'OUTPUT_PARENT_MISSING')]
        for is_optimized in (False, True):
            for label, output, marker in probes:
                # A scientific fault would fail later; observing the output
                # marker proves the guard ran before the proof calculation.
                result = run(output, is_optimized, 'magnetic_sign')
                need(result.returncode != 0 and marker in result.stderr, 'EARLY_OUTPUT_GUARD')
                need(sentinel.read_text(encoding='utf-8') == 'retained sentinel\n', 'SENTINEL_PRESERVED')
                if label not in ('existing', 'leaf_symlink'):
                    need(not output.exists(), 'GUARD_OUTPUT_ABSENT')
                output_guards.append({'optimized':is_optimized, 'case':label,
                                      'expected_marker':marker, 'rejected_before_calculation':True})
    need(source_before == {p.name:sha(p) for p in (Path(__file__),CHECKER)}, 'SOURCE_INTEGRITY')
    report = {'status':'PASS_SYMMETRY_SCIENTIFIC_AND_OUTPUT_CONTROLS',
              'normal_optimized_identical':True, 'report_sha256':sha(normal),
              'source_sha256':source_before, 'rejections':rejections,
              'output_guards':output_guards, 'source_unchanged':True,
              'exclusive_report_creation':True}
    with (destination / 'SYMMETRY_CONTROLS.json').open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status':report['status'], 'output_directory_name':destination.name}))

if __name__ == '__main__':
    main()
