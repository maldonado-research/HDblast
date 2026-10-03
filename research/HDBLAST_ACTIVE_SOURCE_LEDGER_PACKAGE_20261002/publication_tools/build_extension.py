#!/usr/bin/env python3
"""Build a locally pinned v25 extension; no network, secret, publish or science run."""
import argparse
import copy
from datetime import date
import hashlib
import json
from pathlib import Path
import re
import shutil
import stat
import sys
import zipfile

sys.path.insert(0, str(Path(__file__).parent / 'helpers'))
from extension_inventory import (checked_file, read_json, require, sha256,
    validate_extension_inventory, authenticated_target, validate_target_receipt)

ROLES = ('active_source_complete_package', 'active_source_final_result',
         'active_source_external_fresh_replay')
NAMES = ('HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER.zip',
         'HDBLAST_ACTIVE_SOURCE_LEDGER_REPORT_20261003.md',
         'HDBLAST_ACTIVE_SOURCE_LEDGER_FRESH_ZIP_REPLAY_20261003.json')
CHECKPOINT_NAME = 'HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER'


def json_write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def input_file(path):
    path = Path(path).absolute()
    require(path.is_file() and not path.is_symlink()
            and all(not parent.is_symlink() for parent in path.parents), 'Missing or symlink input file')
    return path


def record(path, filename, role):
    require(re.fullmatch('[A-Za-z0-9_.-]+', filename) is not None,
            'Unsafe attachment filename')
    h, m = hashlib.sha256(), hashlib.md5()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block); m.update(block)
    size = path.stat().st_size
    require(size > 0, 'Empty attachment')
    return {'filename': filename, 'role': role, 'bytes': size,
            'sha256': h.hexdigest(), 'md5': m.hexdigest(),
            'source': {'kind': 'candidate', 'path': 'attachments/' + filename,
                       'bytes': size, 'sha256': h.hexdigest()}}


def package_metadata(path, expected_manifest):
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        manifest_name = CHECKPOINT_NAME + '/MANIFEST.json'
        require(len(names) == len(set(names)) and manifest_name in names,
                'Duplicate or absent ZIP manifest')
        for item in archive.infolist():
            require(item.filename and not Path(item.filename).is_absolute()
                    and '\\' not in item.filename and '\x00' not in item.orig_filename
                    and item.filename.startswith(CHECKPOINT_NAME + '/')
                    and all(part not in ('', '.', '..') for part in item.filename.split('/'))
                    and stat.S_ISREG(item.external_attr >> 16)
                    and not item.flag_bits & 1, 'Unsafe or nonregular ZIP member')
        require(archive.getinfo(manifest_name).file_size <= 8 * 1024 * 1024,
                'Oversize package manifest')
        raw = archive.read(manifest_name)
        require(hashlib.sha256(raw).hexdigest() == expected_manifest, 'ZIP manifest pin differs')
        def pairs(items):
            result = {}
            for key, value in items:
                require(key not in result, 'Duplicate package JSON key')
                result[key] = value
            return result
        manifest = json.loads(raw, object_pairs_hook=pairs)
        require(type(manifest) is dict and type(manifest.get('files')) is dict,
                'Malformed ZIP manifest')
        table = manifest['files']
        require(set(names) == {manifest_name} | {CHECKPOINT_NAME + '/' + name for name in table},
                'ZIP payload membership differs from manifest')
        for name, digest in table.items():
            require(type(digest) is str and re.fullmatch('[0-9a-f]{64}', digest),
                    'Invalid ZIP payload digest')
            h = hashlib.sha256()
            with archive.open(CHECKPOINT_NAME + '/' + name) as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b''):
                    h.update(block)
            require(h.hexdigest() == digest, 'ZIP payload bytes differ from manifest')
    return len(names) - 1


def completion_receipt(receipt, evidence, zip_record):
    require(type(receipt) is dict and receipt.get('status') == 'PASS_FRESH_STANDALONE_ZIP_REPLAY'
            and receipt.get('classification') == evidence.get('classification')
            and receipt.get('old_metric_status') == 'FAIL', 'Fresh receipt is incomplete or disagrees')
    require(receipt.get('zip_sha256') == zip_record['sha256']
            and type(receipt.get('zip_bytes')) is int and receipt['zip_bytes'] == zip_record['bytes'],
            'Fresh receipt describes different package bytes')
    replay = receipt.get('replay')
    require(type(replay) is dict and replay.get('status') == 'COMPLETED_REPLAY'
            and replay.get('plan_only') is False
            and type(replay.get('physical_routes_completed')) is int
            and replay['physical_routes_completed'] == 2
            and type(replay.get('successful_commands')) is int
            and replay['successful_commands'] == 4
            and replay.get('source_verification_before') == 'PASS'
            and replay.get('source_verification_after') == 'PASS'
            and replay.get('classification') == evidence.get('classification')
            and replay.get('old_metric_status') == 'FAIL', 'Fresh route/validator execution is incomplete')
    checks = receipt.get('checks')
    require(type(checks) is dict and {'primary_all_scientific_records_equal',
            'independent_all_scientific_records_equal'} <= set(checks)
            and all(value is True for value in checks.values()),
            'Full original/fresh scientific comparison is absent or unequal')
    return receipt


def build(old_candidate, expected_old_manifest, expected_old_metadata, active_paths,
          evidence_path, publication_date, output, published_get_path=None, draft_get_path=None):
    old_candidate = Path(old_candidate).absolute()
    manifest_path = checked_file(old_candidate, 'FILE_MANIFEST.json')
    metadata_path = checked_file(old_candidate, 'METADATA.json')
    require(sha256(manifest_path) == expected_old_manifest, 'Original candidate manifest pin differs')
    require(sha256(metadata_path) == expected_old_metadata, 'Original metadata pin differs')
    old = read_json(manifest_path); old_metadata = read_json(metadata_path)
    require(old['metadata_request']['sha256'] == expected_old_metadata,
            'Original candidate metadata differs internally')
    require(old.get('concept_doi') == '10.5281/zenodo.17088132'
            and old.get('latest_published_record') == 22347452
            and old.get('expected_final_files_count') == 19
            and len(old.get('inherited_files', [])) == 10
            and len(old.get('new_files', [])) == 9,
            'Wrong original main-family candidate')
    require(published_get_path is not None and draft_get_path is not None,
            'Authenticated reconciled latest/new-draft GET evidence is required')
    published_path, draft_path = input_file(published_get_path), input_file(draft_get_path)
    target_receipt = {'status': 'VERIFIED_RECONCILED_MAIN_TARGET',
        'published': authenticated_target(read_json(published_path), sha256(published_path), old['inherited_files'], True),
        'draft': authenticated_target(read_json(draft_path), sha256(draft_path), old['inherited_files'], False)}
    validate_target_receipt(target_receipt, old['inherited_files'], old_metadata['metadata'])
    parsed_date = date.fromisoformat(publication_date)
    evidence = read_json(input_file(evidence_path))
    paths = [input_file(path) for path in active_paths]
    rows = [record(path, name, role) for path, name, role in zip(paths, NAMES, ROLES)]
    require(evidence.get('zip_sha256') == rows[0]['sha256']
            and evidence.get('result_report_sha256') == rows[1]['sha256']
            and evidence.get('fresh_receipt_sha256') == rows[2]['sha256'],
            'Actual attachments differ from reviewed completion evidence')
    payload_count = package_metadata(paths[0], evidence.get('package_manifest_sha256'))
    # These are byte/metadata checks. Scientific acceptance is supplied by the
    # externally reviewed completion receipt and its own source pins.
    completion_receipt(read_json(paths[2]), evidence, rows[0])
    result_report = paths[1].read_text()
    require(evidence.get('classification') in result_report and 'FAIL' in result_report,
            'Final report does not disclose the selected outcome and old FAIL')
    output = Path(output).absolute()
    require(not output.exists() and not output.is_symlink()
            and all(not parent.is_symlink() for parent in output.parents),
            'Derived candidate output must be fresh with no symlink ancestors')
    require(not output.resolve().is_relative_to(old_candidate.resolve()),
            'Never write inside original candidate')
    derived = copy.deepcopy(old)
    derived.update(candidate_version=parsed_date.strftime('%Y.%m.%d') + '-v25',
        expected_final_files_count=22, active_main_draft=23114217,
        latest_published_record=23112891, latest_published_version=None,
        newversion_url='https://zenodo.org/api/deposit/depositions/23112891/actions/newversion',
        target_reconciliation=target_receipt,
        automatic_github_archiving='ALL_EIGHT_OFF',
        active_source_classification=evidence.get('classification'),
        active_source_evidence=evidence, source_artifact_commit=evidence.get('science_commit'),
        parent_candidate={'manifest_path': 'source_free_snapshot/FILE_MANIFEST.json',
            'manifest_sha256': expected_old_manifest, 'metadata_path': 'source_free_snapshot/METADATA.json',
            'metadata_sha256': expected_old_metadata},
        status='PREPARED_DERIVED_CANDIDATE_NOT_PROOF_OF_PUBLICATION')
    derived['new_files'].extend(rows)
    request = copy.deepcopy(old_metadata); metadata = request['metadata']
    metadata.update(version=derived['candidate_version'], publication_date=parsed_date.isoformat(),
        title='HDBLAST Zenodo v25: Registered Metric-Response Failures and Source-Free/Active-Source Ledger Diagnostics')
    metadata['description'] += (
        '<p>The completed active-source follow-up is included as a separate self-contained package '
        'with the final measured-results report and an external fresh-ZIP replay receipt. Its '
        'registered classification is ' + str(evidence.get('classification')) +
        '. The report preserves the first execution-schema failure, the publicly disclosed '
        'repair chronology, all fixed cases and controls, and every earlier scientific FAIL. '
        'Original and fresh registered full scientific frames match exactly according to the '
        'reviewed external receipt. The native source/phase arithmetic and empirical quadrature '
        'controls are distinguished from MP80/100 final reductions. This finite-cutoff prescribed '
        'calculation does not prove a higher-dimensional Big Bang origin, coupled backreaction, '
        'continuum-momentum certification or a new law of nature.</p>')
    metadata['notes'] += (
        '<p>Start with HDBLAST_ACTIVE_SOURCE_LEDGER_REPORT_20261003.md for the completed active-source '
        'round. The preserved HDBLAST_V25_PUBLICATION_OVERVIEW_20261002.md describes the preceding '
        'source-free nine-attachment candidate snapshot. The complete active-source ZIP SHA256 is '
        + rows[0]['sha256'] + ', size ' + str(rows[0]['bytes']) + ' bytes, with '
        + str(payload_count) + ' payload files plus MANIFEST.json. Its external fresh replay receipt '
        'remains outside the ZIP it verifies. Automatic GitHub archiving is OFF for all eight public '
        'repositories. The published current record23112891 contains only the ten inherited files '
        'and has no semantic-version field. The intended unpublished v25 update targets the '
        'reconciled main-family draft23114217; no test '
        'release or new standalone main record was created.</p>')
    for key in ('freeze_commit', 'science_commit'):
        relationship = {'identifier':
            'https://github.com/maldonado-research/HDblast/commit/' + str(evidence.get(key)),
            'relation': 'references', 'resource_type': 'software', 'scheme': 'url'}
        if relationship not in metadata['related_identifiers']:
            metadata['related_identifiers'].append(relationship)
    output.mkdir(parents=True)
    try:
        (output / 'source_free_snapshot').mkdir(); (output / 'attachments').mkdir()
        shutil.copyfile(manifest_path, output / 'source_free_snapshot/FILE_MANIFEST.json')
        shutil.copyfile(metadata_path, output / 'source_free_snapshot/METADATA.json')
        for row in old['new_files']:
            if row['source']['kind'] == 'candidate':
                original = checked_file(old_candidate, row['source']['path'])
                require(row['source']['bytes'] == row['bytes']
                        and row['source']['sha256'] == row['sha256'],
                        'Original candidate-owned source pin disagrees with attachment')
                original_record = record(original, row['filename'], row['role'])
                require(all(original_record[key] == row[key] for key in ('bytes', 'sha256', 'md5')),
                        'Original candidate-owned attachment bytes changed')
                destination = output / row['source']['path']; destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(original, destination)
                copied_record = record(destination, row['filename'], row['role'])
                require(all(copied_record[key] == row[key] for key in ('bytes', 'sha256', 'md5')),
                        'Copied previous candidate attachment differs')
        for source, row in zip(paths, rows):
            shutil.copyfile(source, output / row['source']['path'])
            require(sha256(output / row['source']['path']) == row['sha256'], 'Copied active attachment changed')
        for name in ('prepare_candidate.py', 'verify_saved_record.py', 'extension_inventory.py'):
            shutil.copyfile(Path(__file__).parent / 'helpers' / name, output / name)
        json_write(output / 'METADATA.json', request)
        derived['metadata_request'] = {'path': 'METADATA.json', 'sha256': sha256(output / 'METADATA.json')}
        derived['publication_requirements'] = [
            'Reuse reconciled main-family draft23114217 and preserve all ten inherited files.',
            'Saved metadata must exactly match METADATA.json and contain nonempty required fields.',
            'All twenty-two remote filenames, byte sizes and MD5 pins must match the candidate.',
            'Require submitted=true, state=done and the intended main-family record DOI after publication.',
            'No unchanged retry of the known metadata HTTP500; use diagnosed corrected transport or existing-draft browser fallback.']
        json_write(output / 'FILE_MANIFEST.json', derived)
        validate_extension_inventory(derived, output)
        json_write(output / 'BUILD_RECEIPT.json', {'status': 'PASS_LOCAL_DERIVED_CANDIDATE_BUILD',
            'concept_doi': derived['concept_doi'], 'draft': 23114217, 'candidate_version': derived['candidate_version'],
            'manifest_sha256': sha256(output / 'FILE_MANIFEST.json'),
            'metadata_sha256': derived['metadata_request']['sha256'], 'inherited_files': 10,
            'preserved_previous_additions': 9, 'active_source_additions': 3, 'expected_final_files': 22,
            'remote_mutations': 0, 'scientific_commands': 0,
            'scope': 'Local attachment byte checks, ZIP metadata and externally reviewed completion evidence only.'})
        shutil.copyfile(Path(__file__).parent / 'CANDIDATE_README.md', output / 'README.md')
        return derived
    except Exception as error:
        json_write(output / 'BUILD_FAILED.json', {'status': 'FAILED_LOCAL_BUILD',
            'error_type': type(error).__name__, 'remote_mutations': 0})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old-candidate', type=Path, required=True)
    parser.add_argument('--expected-old-manifest-sha256', required=True)
    parser.add_argument('--expected-old-metadata-sha256', required=True)
    parser.add_argument('--active-zip', type=Path, required=True)
    parser.add_argument('--active-report', type=Path, required=True)
    parser.add_argument('--fresh-receipt', type=Path, required=True)
    parser.add_argument('--completion-evidence', type=Path, required=True)
    parser.add_argument('--published-get', type=Path, required=True)
    parser.add_argument('--draft-get', type=Path, required=True)
    parser.add_argument('--publication-date', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = build(args.old_candidate, args.expected_old_manifest_sha256,
        args.expected_old_metadata_sha256, [args.active_zip, args.active_report, args.fresh_receipt],
        args.completion_evidence, args.publication_date, args.output, args.published_get, args.draft_get)
    print(json.dumps({'status': result['status'], 'expected_files': 22,
                      'manifest_sha256': sha256(args.output / 'FILE_MANIFEST.json')}))


if __name__ == '__main__':
    main()
