#!/usr/bin/env python3
"""Authenticate the separate attribution package and run four tiny checks.

Standard library only. No core scientific replay, network or Git mutation.
The externally supplied pin is mandatory before any arithmetic program runs.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys

PACKAGE = 'HDBLAST_RESPONSE_ATTRIBUTION_20261016'
AUTHORIZATION = '3cdeb704cc8b74319912ba307ead9dfe04248ac6e5443d20146fa02cd23902f8'
JOBS = [('producer', 'followon-response/check_budget.py',
         'followon-response/BUDGET_NORMAL.json', 'followon-response/BUDGET_OPTIMIZED.json',
         'PASS_PROSPECTIVE_BUDGET_ARITHMETIC_ONLY'),
        ('independent', 'followon-response-review/check_independent_budget.py',
         'followon-response-review/BUDGET_REVIEW_NORMAL.json',
         'followon-response-review/BUDGET_REVIEW_OPTIMIZED.json',
         'PASS_INDEPENDENT_PROSPECTIVE_RESPONSE_BUDGET')]


def need(ok, reason):
    if not ok:
        raise ValueError(reason)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def load(raw):
    def pairs(rows):
        result = {}
        for key, value in rows:
            need(key not in result, 'duplicate_JSON_key')
            result[key] = value
        return result
    def nonfinite(value):
        raise ValueError('nonfinite_JSON_' + value)
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=nonfinite)


def typed_equal(a, b):
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(typed_equal(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(typed_equal(x, y) for x, y in zip(a, b))
    return a == b


def authenticate(root, pins, pins_raw):
    need(pins['package'] == PACKAGE and pins['destination_prefix'] == 'research/' + PACKAGE and
         pins['authorization_sha256'] == AUTHORIZATION, 'reviewed_package_identity')
    need(pins['source_files_including_manifest'] == 24 and pins['approved_scientific_assets'] == 21 and
         pins['replay_processes'] == 4 and pins['mutation_processes'] == 0,
         'separate_24_file_four_process_contract')
    listed = pins['files']
    need(isinstance(listed, dict) and len(listed) == 23 and 'PUBLICATION_PINS.json' not in listed,
         'complete_nonself_pin_inventory')
    for name, pin in listed.items():
        rel = PurePosixPath(name)
        need(bool(name) and name != '.' and chr(92) not in name and not rel.is_absolute() and
             '..' not in rel.parts and rel.as_posix() == name, 'canonical_package_member')
        need(type(pin['bytes']) is int and pin['bytes'] >= 0 and
             re.fullmatch('[0-9a-f]{64}', pin['sha256']) is not None, 'complete_member_pin')
        target = root / name
        need(target.is_file() and not any(p.is_symlink() for p in (target, *target.parents)),
             'regular_package_member_' + name)
        raw = target.read_bytes()
        need(len(raw) == pin['bytes'] and sha(raw) == pin['sha256'], 'member_pin_differs_' + name)
    walked = list(root.rglob('*'))
    need(all(not p.is_symlink() and (p.is_file() or p.is_dir()) for p in walked), 'special_package_entry')
    actual = {p.relative_to(root).as_posix() for p in walked if p.is_file()}
    need(actual == set(listed) | {'PUBLICATION_PINS.json'}, 'complete_package_inventory_differs')
    need((root / 'PUBLICATION_PINS.json').read_bytes() == pins_raw, 'publication_pins_changed')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pins-sha256', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    need(re.fullmatch('[0-9a-f]{64}', args.pins_sha256) is not None and args.pins_sha256 != '0' * 64,
         'concrete_external_publication_pin_required')
    root = Path(__file__).resolve().parent
    pins_raw = (root / 'PUBLICATION_PINS.json').read_bytes()
    need(sha(pins_raw) == args.pins_sha256, 'external_publication_pin_differs')
    pins = load(pins_raw)
    authenticate(root, pins, pins_raw)
    out = args.output_dir.resolve()
    need(not out.exists() and out != root and root not in out.parents, 'new_external_output_directory_required')
    out.mkdir(parents=True, exist_ok=False)
    records = []
    for name, script, normal, optimized, status in JOBS:
        for mode, flags, golden in [('normal', [], normal), ('optimized', ['-O'], optimized)]:
            output = out / (name + '.' + mode + '.json')
            command = [sys.executable, '-I', '-S', '-B', *flags, str(root / script), '--output', str(output)]
            result = subprocess.run(command, cwd=root, capture_output=True, timeout=30)
            (out / (name + '.' + mode + '.stdout')).write_bytes(result.stdout)
            (out / (name + '.' + mode + '.stderr')).write_bytes(result.stderr)
            need(result.returncode == 0 and output.is_file(), 'arithmetic_process_failed_' + name + '_' + mode)
            raw, expected = output.read_bytes(), (root / golden).read_bytes()
            parsed, accepted = load(raw), load(expected)
            need(typed_equal(parsed, accepted) and raw == expected and parsed['status'] == status,
                 'complete_typed_and_raw_golden_disagreement_' + name + '_' + mode)
            records.append({'name': name, 'mode': mode, 'script': script, 'golden': golden,
                            'report_sha256': sha(raw), 'complete_typed_JSON_matches': True,
                            'raw_bytes_match': True})
    authenticate(root, pins, pins_raw)
    report = {'status': 'PASS_PINNED_FOLLOWON_FOUR_ARITHMETIC_PROCESSES', 'package': PACKAGE,
              'publication_pins_sha256': args.pins_sha256, 'source_files_including_manifest': 24,
              'arithmetic_processes': len(records), 'mutation_processes': 0,
              'core_scientific_replays': 0, 'source_unchanged': True, 'records': records,
              'scope': 'Conditional auxiliary increment/prospective budget arithmetic only; no computed waveform, candidate, residual, nonzero increment or separate physical constraint solution.'}
    (out / 'VERIFICATION_RECORD.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
