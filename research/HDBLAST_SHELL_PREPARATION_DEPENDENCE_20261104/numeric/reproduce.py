#!/usr/bin/env python3
"""Portable frozen preparation-response replay, with standard-library guards."""
import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path, PurePosixPath

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
BASE_REL = Path('research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/FROZEN_SCIENTIFIC_REFERENCE.json')
INJECTIONS = {
    'drop_conditioned_scalar': 'M2_slope',
    'skip_initial_norm': 'N_h',
    'skip_resolvent_feedback': 'E_h',
    'add_unchanged_first_forcing': 'E_h',
    'skip_target_feedback': 'M2_slope',
    'skip_kernel_comparison': 'volume_history_error',
    'drop_proper_time_factors': 'volume_slope',
    'reverse_lambda_sign': 'finite_lambda_witness',
}


class Reject(Exception):
    def __init__(self, witness):
        self.witness = witness


def require(ok, witness):
    if not ok:
        raise Reject(witness)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def digest(path):
    raw = path.read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def regular(path):
    return path.is_file() and not path.is_symlink() and all(
        not p.is_symlink() for p in path.absolute().parents)


def relative(name):
    path = PurePosixPath(name)
    require(isinstance(name, str) and path.as_posix() == name and
            not path.is_absolute() and '..' not in path.parts and name != '.',
            'CANONICAL_RELATIVE_PATH')
    return Path(*path.parts)


def authenticate(repo):
    local_path = HERE/'LOCAL_PINS.json'
    source_path = HERE/'SOURCE_PINS.json'
    require(regular(local_path) and regular(source_path), 'PIN_FILES_REGULAR')
    local = json.loads(local_path.read_text())
    source = json.loads(source_path.read_text())
    for root, rows, witness in ((HERE, local, 'LOCAL_MEMBER_BYTES'),
                                (repo, source, 'SOURCE_MEMBER_BYTES')):
        for name, expected in rows.items():
            path = root/relative(name)
            require(regular(path) and digest(path) == expected, witness+':'+name)
    freeze = json.loads((HERE/'CAMPAIGN_FREEZE.json').read_text())
    require(freeze['direction']['P'] == '1' and freeze['target'] == '5/4' and
            freeze['partition']['N'] == 1024, 'FROZEN_REGISTERED_CAMPAIGN')
    original = json.loads((HERE/'references/ORIGINAL_REPORT.json').read_text())
    reference = json.loads((HERE/'references/SCIENTIFIC_REFERENCE.json').read_text())
    require(reference == {k: v for k, v in original.items()
                          if k != 'source_authentication'}, 'ORIGINAL_SCIENCE_PROJECTION')
    return {'local': local, 'source': source,
            'local_pins_sha256': digest(local_path)['sha256'],
            'source_pins_sha256': digest(source_path)['sha256']}


def load_science():
    spec = importlib.util.spec_from_file_location('portable_preparation_science', HERE/'scientific.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_result(result, reference, inject='none'):
    # A control first checks the field that its real arithmetic mutation changes.
    if inject != 'none':
        field = INJECTIONS[inject]
        require(result[field] == reference[field], 'FROZEN_SCIENCE_MISMATCH:'+field)
        raise Reject('CONTROL_FAILED_TO_CHANGE_REGISTERED_FIELD')
    require(result.keys() == reference.keys(), 'FROZEN_SCIENTIFIC_FIELD_SET')
    for field in sorted(reference):
        require(result[field] == reference[field], 'FROZEN_SCIENCE_MISMATCH:'+field)


def output_guard(output, repo):
    require(all(not p.is_symlink() for p in (output, *output.parents)), 'OUTPUT_SYMLINK')
    require(not output.exists(), 'OUTPUT_ALREADY_EXISTS')
    resolved = output.resolve()
    require(all(resolved != root and root not in resolved.parents
                for root in (repo.resolve(), HERE.resolve())), 'OUTPUT_INSIDE_SOURCE')


def execute(args):
    repo = args.repo_root.resolve(strict=True)
    output = args.output_dir.absolute()
    output_guard(output, repo)
    before = authenticate(repo)
    reference = json.loads((HERE/'references/SCIENTIFIC_REFERENCE.json').read_text())
    science = load_science()  # All scientific bytes authenticated above.
    if args.use_frozen_raw_for_control:
        require(args.inject != 'none', 'FROZEN_RAW_MODE_REQUIRES_CONTROL')
        raw = json.loads((HERE/'references/RAW_WITNESS.json').read_text())
    else:
        raw = science.generate(repo, digest(HERE/'CAMPAIGN_FREEZE.json')['sha256'])
    require(canonical(raw) == canonical(json.loads((HERE/'references/RAW_WITNESS.json').read_text())),
            'FROZEN_RAW_WITNESS_MATCH')
    base = json.loads((repo/BASE_REL).read_text())
    result = science.aggregate(raw, base, inject=args.inject)
    check_result(result, reference, args.inject)
    after = authenticate(repo)
    require(before == after, 'PRE_POST_AUTHENTICATION')
    original = json.loads((HERE/'references/ORIGINAL_REPORT.json').read_text())
    report = {**result, 'source_authentication': before,
              'portable_provenance': {
                  'schema': 'PREPARATION_PORTABLE_PROVENANCE_V1',
                  'original_report': digest(HERE/'references/ORIGINAL_REPORT.json'),
                  'original_raw_witness': digest(HERE/'references/RAW_WITNESS.json'),
                  'original_receipt': digest(HERE/'references/ORIGINAL_RECEIPT.json'),
                  'every_original_scientific_field_equal': True,
                  'scientific_fields_checked': len(reference),
                  'changed_fields': ['source_authentication', 'portable_provenance'],
                  'source_authentication_scope': 'Direct byte pins; enclosing capsule authenticates the full retained closure.',
                  'same_original_status': result['status'] == original['status']}}
    output.mkdir(parents=True)
    (output/'RAW_WITNESS.json').write_text(json.dumps(raw, indent=2, sort_keys=True)+'\n')
    (output/'REPORT.json').write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    receipt = {'status': 'PASS_PORTABLE_FROZEN_PREPARATION_SCIENCE',
               'report': digest(output/'REPORT.json'), 'raw_witness': digest(output/'RAW_WITNESS.json'),
               'original_scientific_fields_equal': len(reference),
               'raw_witness_byte_identical_to_original': (output/'RAW_WITNESS.json').read_bytes() ==
                    (HERE/'references/RAW_WITNESS.json').read_bytes(),
               'pre_post_authentication_identical': True,
               'oracle_calls': raw['calls'],
               'provenance_version': 'PREPARATION_PORTABLE_PROVENANCE_V1'}
    (output/'RECEIPT.json').write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    return {'status': 'PASS', 'witness': 'FROZEN_PREPARATION_SCIENCE_AND_RAW_MATCH',
            'original_scientific_fields_equal': len(reference), 'oracle_calls': raw['calls']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--inject', choices=('none', *INJECTIONS), default='none')
    parser.add_argument('--use-frozen-raw-for-control', action='store_true',
                        help='Controls only: reaggregate authenticated raw, generating no new bump boxes.')
    args = parser.parse_args()
    try:
        final = execute(args)
        code = 0
    except Reject as error:
        final = {'status': 'REJECTED', 'witness': error.witness}
        code = 2
    except (OSError, ValueError, KeyError, TypeError, StopIteration):
        final = {'status': 'REJECTED', 'witness': 'INPUT_OR_SCIENTIFIC_FAILURE'}
        code = 2
    print(json.dumps(final, sort_keys=True))
    return code


if __name__ == '__main__':
    raise SystemExit(main())
