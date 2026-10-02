#!/usr/bin/env python3
"""Synthetic inert payloads only: test real64MiB split and standalone assembly."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import build_package as B
import reassemble_package as R

HERE = Path(__file__).resolve().parent


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def fixture(parent, large=False):
    root = parent/B.CHECKPOINT_NAME
    (root/'code').mkdir(parents=True); (root/'independent').mkdir()
    for name in ('build_package.py','reassemble_package.py'):
        shutil.copyfile(HERE/name, root/'code'/name)
    (root/'independent/forced_metric.py').write_text('# Inert synthetic source, never executed.\n')
    write(root/'independent/MANIFEST.json', {'files':[{'path':'forced_metric.py',
          'sha256':R.file_sha(root/'independent/forced_metric.py')}]})
    runtime = {'python_major_minor':[3,12], 'physical_producers_optimized':False,
               'numpy':'2.2.6','scipy':'1.15.3','sympy':'1.14.0','mpmath':'1.3.0','matplotlib':'3.10.1'}
    experiment = {'gates':{'synthetic':1},'runtime':runtime}
    write(root/'EXPERIMENT.json', experiment)
    files = {p.relative_to(root).as_posix():R.file_sha(p) for p in sorted(root.rglob('*')) if p.is_file()}
    registration = {'files':files, 'independent_manifest_sha256':R.file_sha(root/'independent/MANIFEST.json'),
                    'frozen_configuration':experiment, 'frozen_gates':experiment['gates']}
    write(root/'FULL_REGISTRATION.json', registration)
    write(root/'FREEZE_RECEIPT.json', {'public_freeze_commit':'1'*40,
        'registration_sha256':R.file_sha(root/'FULL_REGISTRATION.json'),
        'independent_manifest_sha256':registration['independent_manifest_sha256']})
    (root/'README.md').write_text('Inert synthetic package; no scientific results.\n')
    if large:
        with (root/'synthetic_incompressible.bin').open('wb') as stream:
            for count in range(65):
                stream.write(hashlib.shake_256(('synthetic-payload-'+str(count)).encode()).digest(1024*1024))
    return root


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replay-driver', type=Path, required=True)
    parser.add_argument('--small-only', action='store_true')
    args = parser.parse_args()
    R.require(not args.output.exists(), 'Fresh synthetic output directory required')
    args.output.mkdir(parents=True)
    checks = []
    def record(name, condition):
        R.require(condition, 'Synthetic package check failed: '+name); checks.append(name)
    spec = importlib.util.spec_from_file_location('frozen_input_verifier_only', args.replay_driver)
    replay = importlib.util.module_from_spec(spec); spec.loader.exec_module(replay)
    sizes = (False,) if args.small_only else (False, True)
    results = []
    for large in sizes:
        kind = 'split' if large else 'single'
        root = fixture(args.output/kind, large)
        before = B.payload_files(root)
        first = B.build(root, args.output/(kind+'_external_first'))
        second = B.build(root, args.output/(kind+'_external_second'))
        record(kind+' complete deterministic builds byte-identical',
               B.equal_files(Path(first['zip_path']), Path(second['zip_path'])) and
               first['zip_sha256'] == second['zip_sha256'])
        record(kind+' all original payload unchanged', B.payload_files(root) == before)
        record(kind+' canonical exclusions accepted by unchanged32command driver',
               replay.verify_inputs(root)['package_manifest_sha256'] == R.file_sha(root/'MANIFEST.json'))
        ledger = json.loads((root/R.PARTS_MANIFEST_NAME).read_text())
        record(kind+' all tracked files64MiB maximum',
               all(p['bytes'] <= R.PART_BYTES for p in ledger['parts']))
        record(kind+' correct storage policy', (first['tracked_storage'] == kind and
               first['tracked_zip_present'] is (not large)))
        helper = args.output/(kind+'_standalone_helper.py')
        shutil.copyfile(HERE/'reassemble_package.py', helper)
        reconstructed = args.output/(kind+'_standalone_reconstructed.zip')
        process = subprocess.run([sys.executable, str(helper), '--parts-manifest',
            str(root/R.PARTS_MANIFEST_NAME), '--output', str(reconstructed),
            '--expected-sha256', first['zip_sha256']], text=True, capture_output=True)
        R.require(process.returncode == 0, process.stderr)
        record(kind+' standalone fresh ZIP byte-identical and fully verified',
               B.equal_files(reconstructed, Path(first['zip_path'])) and
               json.loads(process.stdout)['status'] == 'PASS_REASSEMBLY_ONLY')
        results.append(first)
        if not large:
            continue
        original_ledger = (root/R.PARTS_MANIFEST_NAME).read_bytes()
        def reject_ledger(name, change):
            document = json.loads(original_ledger); change(document)
            write(root/R.PARTS_MANIFEST_NAME, document)
            try:
                R.reassemble(root/R.PARTS_MANIFEST_NAME, args.output/(name+'.zip'))
            except (RuntimeError, KeyError, ValueError, TypeError):
                checks.append(name+' rejected')
            else:
                raise RuntimeError('Assembly mutation survived: '+name)
            finally:
                (root/R.PARTS_MANIFEST_NAME).write_bytes(original_ledger)
        reject_ledger('wrong_logical_SHA', lambda d:d.update(zip_sha256='0'*64))
        reject_ledger('wrong_part_SHA', lambda d:d['parts'][0].update(sha256='0'*64))
        reject_ledger('wrong_part_order', lambda d:d['parts'].reverse())
        reject_ledger('wrong_manifest_SHA', lambda d:d.update(package_manifest_sha256='0'*64))
        reject_ledger('wrong_public_freeze', lambda d:d.update(public_freeze_commit='0'*40))
        reject_ledger('unsafe_part_path', lambda d:d['parts'][0].update(path='../escape'))
        reject_ledger('oversized_part', lambda d:d['parts'][0].update(bytes=R.PART_BYTES+1))
        reject_ledger('missing_part', lambda d:d['parts'].pop())
        part = root/ledger['parts'][0]['path']
        with part.open('r+b') as stream:
            byte = stream.read(1); stream.seek(0); stream.write(bytes([byte[0]^1]))
        try:
            R.reassemble(root/R.PARTS_MANIFEST_NAME, args.output/'corrupt_part.zip')
        except RuntimeError:
            checks.append('corrupt physical part rejected')
        else:
            raise RuntimeError('Corrupt part accepted')
        finally:
            with part.open('r+b') as stream:
                stream.write(byte)
    def reject_builder(name, mutate):
        root = fixture(args.output/name); mutate(root)
        try:
            B.build(root, args.output/(name+'_external'))
        except (RuntimeError, KeyError, ValueError, TypeError):
            checks.append(name+' rejected')
        else:
            raise RuntimeError('Builder mutation survived: '+name)
    reject_builder('missing_receipt', lambda p:(p/'FREEZE_RECEIPT.json').unlink())
    reject_builder('mutated_frozen_source', lambda p:(p/'independent/forced_metric.py').write_text('changed'))
    reject_builder('payload_symlink', lambda p:(p/'alias').symlink_to(p/'README.md'))
    reject_builder('payload_cache', lambda p:(p/'bad.pyc').write_bytes(b'cache'))
    report = {'status':'PASS_SYNTHETIC_PACKAGING_ONLY','physical_evaluations':0,'child_science_producers':0,
              'check_count':len(checks),'checks':checks,'python_optimization':sys.flags.optimize,
              'large_fixture_executed':not args.small_only,'results':results,
              'source_pins':{name:R.file_sha(HERE/name) for name in
                             ('build_package.py','reassemble_package.py','test_package_parts.py')},
              'replay_driver_sha256':R.file_sha(args.replay_driver)}
    write(args.output/'CHECKS.json', report)
    print(json.dumps({k:report[k] for k in ('status','check_count','physical_evaluations','large_fixture_executed')}))


if __name__ == '__main__':
    main()
