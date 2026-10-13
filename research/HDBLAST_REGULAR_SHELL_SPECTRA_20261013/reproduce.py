#!/usr/bin/env python3
"""Verify and replay the frozen regular-shell spectral research capsule."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
STATUS = 'PASS_REVIEWED_REGULAR_SHELL_SPECTRA_20261013'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def snapshot(root):
    result = {}
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(), 'Symbolic links are not accepted')
        if path.is_file():
            result[path.relative_to(root).as_posix()] = digest(path)
    return result


def verify(root):
    manifest = json.loads((root / 'MANIFEST.json').read_text())
    require(manifest['status'] == 'FROZEN_REVIEWED_RESEARCH_CAPSULE', 'Manifest is not frozen')
    actual = snapshot(root)
    names = set()
    for row in manifest['files']:
        name = row['path']; part = Path(name)
        require(not part.is_absolute() and '..' not in part.parts and name not in names,
                'Unsafe or duplicate manifest path')
        names.add(name)
        path = root / part
        require(actual.get(name) == row['sha256'] and path.stat().st_size == row['bytes'],
                'Manifest file mismatch: ' + name)
    require(names == set(actual) - {'MANIFEST.json'}, 'Manifest coverage mismatch')
    return actual


def command(argv, cwd, log, expect_failure=False):
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    p = subprocess.run(argv, cwd=cwd, env=environment, capture_output=True, text=True)
    log.write_text(p.stdout + p.stderr)
    require((p.returncode != 0) if expect_failure else (p.returncode == 0),
            'Unexpected command outcome; inspect ' + str(log))
    return p


def py(script, *args, optimized=False):
    return [sys.executable, '-I', '-B'] + (['-O'] if optimized else []) + [str(script), *map(str, args)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    target = args.output.resolve()
    require(target != ROOT and ROOT not in target.parents, 'Output must be outside the capsule')
    require(not target.exists(), 'Output already exists; choose a fresh directory')
    import sympy, mpmath
    require(platform.python_version() == '3.12.14', 'Pinned receipt replay requires Python 3.12.14')
    require(sympy.__version__ == '1.14.0' and mpmath.__version__ == '1.3.0', 'Install pinned requirements')
    before = verify(ROOT)
    target.mkdir(parents=True)
    work = target / 'capsule'
    shutil.copytree(ROOT, work)
    require(snapshot(work) == before, 'Copied capsule differs')
    contracts = json.loads((work / 'CHECK_CONTRACTS.json').read_text())
    results = []
    for mode in ('normal', 'optimized'):
        stage = target / mode; stage.mkdir()
        for item in contracts['checks']:
            out = stage / (item['name'] + '.json')
            command(py(work / item['script'], '--output', out, optimized=mode=='optimized'),
                    target, stage / (item['name'] + '.log'))
            require(out.read_bytes() == (work / item['receipt']).read_bytes(),
                    'Replay differs from reviewed receipt: ' + item['name'] + '/' + mode)
            results.append({'mode':mode, 'name':item['name'], 'sha256':digest(out),
                            'positive_checks':item['positive_checks'],
                            'wrong_formula_controls':item['wrong_formula_controls']})
    for item in contracts['checks']:
        name = item['name'] + '.json'
        require((target/'normal'/name).read_bytes() == (target/'optimized'/name).read_bytes(),
                'Optimized execution differs: ' + name)
    # The original branch arithmetic is a dependency of the spectral claims.
    command(py(work/'basis/static/reproduce.py', '--output', target/'static-replay'),
            target, target/'static-replay.log')
    baseline = json.loads((target/'static-replay/REPRODUCTION_RECEIPT.json').read_text())
    require(baseline['status'] == 'PASS_SELF_CONTAINED_EXTERNAL_REPLAY' and not baseline['provisional'],
            'Static basis replay did not pass')
    # A changed mathematical basis must be rejected before a success receipt exists.
    bad_basis = target / 'corrupted-basis'; shutil.copytree(work/'basis', bad_basis)
    proof = bad_basis/'bridge/PROSPECTIVE_ZERO_MODE_STRATEGY.md'
    proof.write_bytes(proof.read_bytes() + b'\nDeliberate corrupt-input control.\n')
    bad_out = target/'MUST_NOT_EXIST.json'
    for mode in (False, True):
        command(py(work/'independent/check_scalar.py', '--basis', bad_basis, '--output', bad_out,
                   optimized=mode), target, target/('corrupted-basis-' + str(mode) + '.log'), True)
        require(not bad_out.exists(), 'Corrupt basis produced a success receipt')
    # Verify manifest coverage catches corruption as well, separately from pinned basis.
    bad_copy = target/'manifest-control'; bad_copy.mkdir()
    (bad_copy/'MANIFEST.json').write_text(json.dumps({'status':'FROZEN_REVIEWED_RESEARCH_CAPSULE',
         'files':[{'path':'proof.md','bytes':1,'sha256':hashlib.sha256(b'x').hexdigest()}]}))
    (bad_copy/'proof.md').write_bytes(b'y')
    rejected = False
    try:
        verify(bad_copy)
    except RuntimeError as exc:
        rejected = 'Manifest file mismatch' in str(exc)
    require(rejected, 'Manifest failed to reject altered content')
    require(snapshot(ROOT) == before and snapshot(work) == before, 'Frozen capsule changed during replay')
    report = {'status':STATUS, 'source_manifest_sha256':digest(ROOT/'MANIFEST.json'),
              'source_and_relocated_copy_unchanged':True, 'normal_and_optimized_receipts_identical':True,
              'python':platform.python_version(),'sympy':sympy.__version__,'mpmath':mpmath.__version__,
              'new_positive_checks_per_mode':sum(c['positive_checks'] for c in contracts['checks']),
              'new_wrong_formula_controls_per_mode':sum(c['wrong_formula_controls'] for c in contracts['checks']),
              'checks':results, 'static_basis_replay':baseline['status'],
              'corrupted_basis_rejected_without_output_both_modes':True,
              'altered_manifest_content_rejected':True,
              'scope':'Specified algebra, exact sign, and source-integrity checks; analytic and physical hypotheses remain as stated in the reviewed proofs.'}
    (target/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:report[k] for k in ('status','new_positive_checks_per_mode',
           'new_wrong_formula_controls_per_mode','normal_and_optimized_receipts_identical')},indent=2))


if __name__ == '__main__':
    main()
