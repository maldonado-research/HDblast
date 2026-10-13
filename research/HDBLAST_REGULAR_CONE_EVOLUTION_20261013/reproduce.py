#!/usr/bin/env python3
"""Authenticate and replay the reviewed regular-cone-evolution capsule.

The immutable Git revision and CI bind the current manifest externally.
--manifest-sha256 enforces that external pin before executing any checker.
--verify-only uses only the standard library and never writes a file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import stat
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
PRIOR_MANIFEST_SHA = '6021915456815936a80c864d5e503b05725325a2cb15cbd882051fade0890194'
MANIFEST_STATUS = 'FROZEN_REVIEWED_REGULAR_CONE_EVOLUTION_CAPSULE'
CHECK_NAMES = ('reconstruction', 'independent', 'canonical', 'boundary')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'Duplicate JSON key: '+key)
        result[key] = value
    return result


def read_json(path):
    return json.loads(path.read_text(), object_pairs_hook=unique_object)


def relative_path(name):
    require(isinstance(name, str) and name and '\\' not in name,
            'Expected a portable nonempty path')
    path = PurePosixPath(name)
    require(not path.is_absolute() and '..' not in path.parts
            and name == path.as_posix() and name != '.',
            'Unsafe or noncanonical relative path: '+name)
    return Path(*path.parts)


def snapshot(root):
    require(root.is_dir() and not root.is_symlink(), 'Real capsule directory required')
    entries = {}
    for path in sorted(root.rglob('*')):
        mode = path.lstat().st_mode
        require(stat.S_ISDIR(mode) or stat.S_ISREG(mode),
                'Only real directories and regular files belong in a capsule')
        if stat.S_ISREG(mode):
            entries[path.relative_to(root).as_posix()] = sha(path)
    return entries


def inventory(data, item):
    positive = data['positive_checks']
    controls = data['wrong_formula_controls']
    require(isinstance(positive, dict) and isinstance(controls, dict),
            'Dictionary inventories required: '+item['name'])
    require(all(value == 'passed' for value in positive.values()),
            'Positive identity did not pass: '+item['name'])
    require(all(isinstance(value, str) and value not in ('', '0')
                for value in controls.values()),
            'Wrong-formula residual was not nonzero: '+item['name'])
    require(sorted(positive) == item['positive_check_names']
            and sorted(controls) == item['wrong_formula_control_names'],
            'Exact check inventory differs: '+item['name'])
    actual = (data['positive_count'], data['wrong_formula_count'])
    require(actual == (len(positive), len(controls))
            and actual == (item['positive_checks'], item['wrong_formula_controls']),
            'Check counts differ: '+item['name'])
    if item['name'] == 'independent':
        require(data['passed'] is True and data['mutation'] is None,
                'Independent checker did not pass its unmodified equations')
        rows = data['checks']
        names = [row['name'] for row in rows]
        require(len(names) == len(set(names)) and sorted(names) == sorted([*positive, *controls])
                and data['check_count'] == len(names)
                and all(row['passed'] is True for row in rows),
                'Independent detailed inventory differs')
    return actual


def contracts(root):
    contract = read_json(root/'CHECK_CONTRACTS.json')
    items = contract['checks']
    require([item['name'] for item in items] == list(CHECK_NAMES),
            'Four reviewed checker contracts required in canonical order')
    for item in items:
        name = item['name']
        require(relative_path(item['script']).parts[0] == name,
                'Checker must live in its contracted directory')
        require(item['expected_receipt'] == 'results/'+name+'.json',
                'Unexpected reviewed receipt path')
        require(item['positive_check_names'] == sorted(set(item['positive_check_names']))
                and item['wrong_formula_control_names'] == sorted(set(item['wrong_formula_control_names'])),
                'Contract inventory must be sorted and unique')
        require(not set(item['positive_check_names']) & set(item['wrong_formula_control_names']),
                'Identity and wrong-formula inventories overlap')
        require(item['positive_checks'] == len(item['positive_check_names'])
                and item['wrong_formula_controls'] == len(item['wrong_formula_control_names']),
                'Contract count and inventory differ')
        require((root/relative_path(item['script'])).is_file(), 'Missing checker: '+name)
        inventory(read_json(root/relative_path(item['expected_receipt'])), item)
    require(contract['positive_checks_per_mode'] == sum(item['positive_checks'] for item in items)
            and contract['wrong_formula_controls_per_mode'] == sum(item['wrong_formula_controls'] for item in items),
            'Contract totals differ')
    return contract


def verify(root, manifest_pin=None):
    manifest_path = root/'MANIFEST.json'
    if manifest_pin is not None:
        require(re.fullmatch(r'[0-9a-f]{64}', manifest_pin) is not None,
                'External manifest pin must be 64 lowercase hexadecimal characters')
        require(sha(manifest_path) == manifest_pin, 'Manifest differs from external pin')
    manifest = read_json(manifest_path)
    require(manifest['status'] == MANIFEST_STATUS, 'Unexpected manifest status')
    actual = snapshot(root)
    names = set()
    total = 0
    for row in manifest['files']:
        name = row['path']
        part = relative_path(name)
        require(name not in names and name != 'MANIFEST.json', 'Duplicate/self manifest path')
        require(type(row['bytes']) is int and row['bytes'] >= 0,
                'Manifest byte count must be a nonnegative integer')
        require(re.fullmatch(r'[0-9a-f]{64}', row['sha256']) is not None,
                'Manifest digest must be SHA256')
        names.add(name)
        path = root/part
        require(actual.get(name) == row['sha256'] and path.stat().st_size == row['bytes'],
                'Capsule file differs: '+name)
        total += row['bytes']
    require(names == set(actual)-{'MANIFEST.json'}, 'Manifest coverage differs')
    require(len(names) == manifest['file_count_excluding_this_manifest']
            and total == manifest['bytes_excluding_this_manifest'], 'Manifest totals differ')

    old = root/'basis'/'PRIOR_MANIFEST.json'
    require(sha(old) == PRIOR_MANIFEST_SHA, 'Prior manifest differs')
    old_entries = {}
    for row in read_json(old)['files']:
        relative_path(row['path'])
        require(row['path'] not in old_entries, 'Duplicate prior manifest path')
        old_entries[row['path']] = row
    source = read_json(root/'basis'/'SOURCE_INPUTS.json')
    require(source['manifest_sha256'] == PRIOR_MANIFEST_SHA,
            'Prior-input declaration does not name the pinned manifest')
    rows = source['copied_files']
    require(bool(rows), 'At least one copied prior input is required')
    copied_names = set()
    for row in rows:
        name = row['path']
        part = relative_path(name)
        require(name not in copied_names and row == old_entries.get(name),
                'Prior input declaration differs: '+name)
        copied_names.add(name)
        path = root/'basis'/'prior-checkpoint'/part
        require(sha(path) == row['sha256'] and path.stat().st_size == row['bytes'],
                'Prior input differs: '+name)
    require(set(snapshot(root/'basis'/'prior-checkpoint')) == copied_names,
            'Copied prior-input coverage differs')
    contracts(root)
    return actual


def command(script, *arguments, optimized=False):
    return [sys.executable, '-I', '-B']+(['-O'] if optimized else [])+[str(script), *map(str, arguments)]


def run(argv, log, expected_failure=False, diagnostic=None):
    process = subprocess.run(argv, cwd=log.parent,
                             env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'),
                             capture_output=True, text=True, timeout=300)
    output = process.stdout+process.stderr
    log.write_text(output)
    require((process.returncode != 0) if expected_failure else (process.returncode == 0),
            'Unexpected exit; inspect '+str(log))
    if diagnostic is not None:
        require(diagnostic in output, 'Expected rejection reason absent; inspect '+str(log))
    return process


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--manifest-sha256')
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    before = verify(ROOT, args.manifest_sha256)
    if args.verify_only:
        print('PASS_CAPSULE_AND_PRIOR_INPUT_INTEGRITY')
        return
    require(args.output is not None, 'Provide a fresh --output directory')
    target = args.output.resolve()
    require(target != ROOT and ROOT not in target.parents, 'Output must be outside capsule')
    require(not target.exists(), 'Output exists; choose a fresh directory')
    import sympy
    import mpmath
    require(platform.python_version() == '3.12.14', 'Recorded replay uses Python 3.12.14')
    require(sympy.__version__ == '1.14.0' and mpmath.__version__ == '1.3.0',
            'Install the recorded pinned dependencies')
    target.mkdir(parents=True)
    work = target/'relocated-capsule'
    shutil.copytree(ROOT, work)
    require(verify(work, args.manifest_sha256) == before, 'Relocated capsule differs')
    contract = contracts(work)
    results = []
    for mode in ('normal', 'optimized'):
        stage = target/mode
        stage.mkdir()
        for item in contract['checks']:
            out = stage/(item['name']+'.json')
            run(command(work/item['script'], '--output', out, optimized=mode == 'optimized'),
                stage/(item['name']+'.log'))
            require(out.read_bytes() == (work/item['expected_receipt']).read_bytes(),
                    'Reviewed receipt differs: '+item['name']+'/'+mode)
            actual_counts = inventory(read_json(out), item)
            results.append({'mode':mode, 'name':item['name'],
                            'positive_checks':actual_counts[0],
                            'wrong_formula_controls':actual_counts[1], 'sha256':sha(out)})
    for item in contract['checks']:
        name = item['name']+'.json'
        require((target/'normal'/name).read_bytes() == (target/'optimized'/name).read_bytes(),
                'Modes differ: '+item['name'])

    corrupt = target/'corrupted-input-copy'
    shutil.copytree(work, corrupt)
    input_name = read_json(corrupt/'basis'/'SOURCE_INPUTS.json')['copied_files'][0]['path']
    input_file = corrupt/'basis'/'prior-checkpoint'/relative_path(input_name)
    input_file.write_bytes(input_file.read_bytes()+b'\nDeliberate input-integrity control.\n')
    # Altered current manifest controls bind to the externally recorded digest.
    changed = target/'altered-manifest-copy'
    shutil.copytree(work, changed)
    original_pin = sha(changed/'MANIFEST.json')
    with (changed/'MANIFEST.json').open('ab') as stream:
        stream.write(b'\n')
    # Re-seal only the outer manifest in a negative-test copy. This reaches the
    # independent fixed prior pin instead of failing at outer file coverage.
    prior_corrupt = target/'altered-prior-manifest-copy'
    shutil.copytree(work, prior_corrupt)
    prior_path = prior_corrupt/'basis'/'PRIOR_MANIFEST.json'
    with prior_path.open('ab') as stream:
        stream.write(b'\n')
    resealed = read_json(prior_corrupt/'MANIFEST.json')
    for row in resealed['files']:
        if row['path'] == 'basis/PRIOR_MANIFEST.json':
            row['bytes'] = prior_path.stat().st_size
            row['sha256'] = sha(prior_path)
    resealed['bytes_excluding_this_manifest'] = sum(row['bytes'] for row in resealed['files'])
    (prior_corrupt/'MANIFEST.json').write_text(json.dumps(resealed, indent=2, sort_keys=True)+'\n')
    for optimized in (False, True):
        mode = 'optimized' if optimized else 'normal'
        run(command(corrupt/'reproduce.py', '--verify-only', optimized=optimized),
            target/('corrupted-input-'+mode+'.log'), True,
            diagnostic='Capsule file differs: basis/prior-checkpoint/'+input_name)
        run(command(changed/'reproduce.py', '--verify-only', '--manifest-sha256', original_pin,
                    optimized=optimized), target/('altered-manifest-'+mode+'.log'), True,
            diagnostic='Manifest differs from external pin')
        run(command(prior_corrupt/'reproduce.py', '--verify-only', optimized=optimized),
            target/('altered-prior-manifest-'+mode+'.log'), True,
            diagnostic='Prior manifest differs')
        for item in contract['checks']:
            out = target/('MUST_NOT_EXIST_'+item['name']+'_'+mode+'.json')
            run(command(work/item['script'], '--inject-failure', '--output', out, optimized=optimized),
                target/('injected-failure-'+item['name']+'-'+mode+'.log'), True)
            require(not out.exists(), 'Injected failure wrote a success receipt: '+item['name'])
        for mutation in ('master-potential', 'scalar-source', 'null-jet'):
            out = target/('MUST_NOT_EXIST_mutation_'+mutation+'_'+mode+'.json')
            run(command(work/'independent'/'check_regular_master_independently.py',
                        '--mutation', mutation, '--output', out, optimized=optimized),
                target/('mutation-'+mutation+'-'+mode+'.log'), True)
            require(not out.exists(), 'Mutation wrote a success receipt: '+mutation)
    require(snapshot(ROOT) == before and snapshot(work) == before, 'Replay changed a capsule')
    receipt = {
        'status':'PASS_REVIEWED_REGULAR_CONE_EVOLUTION_REPLAY',
        'scope':'Exact local algebra and jets, authenticated copied inputs, explicit boundary/evolution controls, and replay integrity. Analytic PDE existence and boundary/evolution hypotheses remain in the reviewed notes; these machine checks are not an analytic theorem proof.',
        'python':platform.python_version(), 'sympy':sympy.__version__, 'mpmath':mpmath.__version__,
        'source_manifest_sha256':sha(ROOT/'MANIFEST.json'),
        'positive_checks_per_mode':contract['positive_checks_per_mode'],
        'wrong_formula_controls_per_mode':contract['wrong_formula_controls_per_mode'],
        'exact_identity_and_control_inventories_verified':True,
        'normal_and_optimized_receipts_identical':True,
        'corrupted_prior_input_rejected_both_modes':True,
        'externally_pinned_altered_manifest_rejected_both_modes':True,
        'altered_prior_manifest_rejected_after_outer_reseal_both_modes':True,
        'injected_failure_rejected_without_success_receipt_all_checkers_both_modes':True,
        'independent_equation_mutations_rejected_without_success_receipt_both_modes':True,
        'source_and_relocated_capsules_unchanged':True,
        'prior_scientific_capsule_replayed':False,
        'checks':results,
    }
    (target/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status':receipt['status'],
                      'positive_checks_per_mode':receipt['positive_checks_per_mode'],
                      'wrong_formula_controls_per_mode':receipt['wrong_formula_controls_per_mode']}, sort_keys=True))


if __name__ == '__main__':
    main()
