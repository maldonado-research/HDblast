#!/usr/bin/env python3
"""Prepare verified publication files locally; no network, secrets or science run."""
import argparse
import hashlib
import json
import re
import shutil
from extension_inventory import validate_extension_inventory
from pathlib import Path

class CandidateError(RuntimeError):
    pass

def require(condition, message):
    if not condition:
        raise CandidateError(message)

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def safe_relative(value):
    require(isinstance(value, str) and value != '', 'Empty or invalid relative path')
    p = Path(value)
    require(not p.is_absolute() and all(x not in ('', '.', '..') for x in value.split('/')), 'Unsafe relative path')
    require('\\' not in value, 'Backslash path is not portable')
    return p

def safe_file(root, relative):
    p = root / safe_relative(relative)
    q = root
    require(root.is_dir() and not root.is_symlink(), 'Source root must be a real directory')
    for part in safe_relative(relative).parts:
        q = q / part
        require(not q.is_symlink(), 'Symlink input is not accepted')
    require(p.is_file(), 'Source file is absent or not a regular file: ' + relative)
    require(p.resolve().is_relative_to(root.resolve()), 'Source escapes root')
    return p

def check_pin(path, record):
    require(type(record.get('bytes')) is int and record['bytes'] > 0, 'Invalid positive byte count')
    require(isinstance(record.get('sha256'), str) and re.fullmatch('[0-9a-f]{64}', record['sha256']), 'Invalid SHA256')
    require(path.stat().st_size == record['bytes'], 'Source byte count changed: ' + path.name)
    require(sha256(path) == record['sha256'], 'Source hash changed: ' + path.name)

def plan_sources(manifest, repository, candidate, prepared_previous=None):
    validate_extension_inventory(manifest, candidate)
    files = manifest.get('new_files')
    require(isinstance(files, list) and len(files) == 12, 'Expected exactly twelve new candidate files')
    require(manifest.get('concept_doi') == '10.5281/zenodo.17088132', 'Wrong main DOI family')
    require(manifest.get('latest_published_record') == 23112891, 'Wrong reconciled latest published source record')
    require(isinstance(manifest.get('candidate_version'), str) and re.fullmatch(r'[0-9]{4}\.[0-9]{2}\.[0-9]{2}-v25', manifest['candidate_version']), 'Wrong candidate version')
    require(manifest.get('old_metric_scientific_status') == 'FAIL', 'Historical scientific FAIL must be retained')
    require(manifest.get('ledger_classification') == 'LEDGER_ERROR_DEMONSTRATED', 'Ledger classification changed')
    inherited = manifest.get('inherited_files')
    require(isinstance(inherited, list) and len(inherited) == 10, 'Expected ten preserved historical files')
    names = set()
    source_sets = []
    inherited_names = {x['filename'] for x in inherited}
    for historical in inherited:
        require(isinstance(historical['filename'], str) and re.fullmatch('[A-Za-z0-9_.-]+', historical['filename']) and historical['filename'] not in ('.', '..'), 'Unsafe inherited filename')
        require(type(historical.get('bytes')) is int and historical['bytes'] > 0, 'Invalid inherited byte count')
        require(isinstance(historical.get('md5'), str) and re.fullmatch('[0-9a-f]{32}', historical['md5']), 'Invalid inherited MD5')
        require(historical.get('policy') == 'RETAIN_UNCHANGED_HISTORICAL_V24', 'Individual inherited file policy changed')
    require(len(inherited_names) == 10, 'Duplicate inherited file names')
    require(manifest.get('inherited_file_policy') == 'RETAIN_ALL_TEN_UNCHANGED_LABEL_HISTORICAL_V24', 'Historical file policy changed')
    require(type(manifest.get('expected_final_files_count')) is int and manifest.get('expected_final_files_count') == 22, 'Wrong final file inventory')
    for index, f in enumerate(files):
        require(type(f.get('bytes')) is int and f['bytes'] > 0, 'Invalid positive output byte count')
        require(isinstance(f.get('sha256'), str) and re.fullmatch('[0-9a-f]{64}', f['sha256']), 'Invalid output SHA256')
        name = f.get('filename')
        require(isinstance(name, str) and re.fullmatch('[A-Za-z0-9_.-]+', name) and name not in ('.', '..'), 'Unsafe output filename')
        require(name not in {'METADATA_REQUEST.json', 'PREPARATION_RECEIPT.json', 'PREPARATION_FAILED.json'} and not name.endswith('.partial'), 'Filename collides with helper sidecar or temporary file')
        require(name not in names and name not in inherited_names, 'Duplicate or inherited-colliding filename')
        names.add(name)
        require(isinstance(f.get('md5'), str) and re.fullmatch('[0-9a-f]{32}', f['md5']), 'Invalid output MD5')
        if prepared_previous is not None and index < 9:
            path = safe_file(prepared_previous, name)
            check_pin(path, f)
            source_sets.append((f, [path]))
            continue
        source = f.get('source', {})
        if source.get('kind') == 'concatenate':
            m = source['parts_manifest']
            mp = safe_file(repository, m['path'])
            require(sha256(mp) == m['sha256'], 'Parts manifest changed')
            parts_doc = json.loads(mp.read_text())
            require(parts_doc['zip_bytes'] == f['bytes'] and parts_doc['zip_sha256'] == f['sha256'], 'Parts describe a different ZIP')
            require(isinstance(source.get('parts'), list) and len(source['parts']) > 0, 'Concatenation parts must not be empty')
            require(isinstance(parts_doc.get('parts'), list) and len(parts_doc['parts']) > 0, 'Parts manifest must not be empty')
            require(len(parts_doc['parts']) == len(source['parts']), 'Parts count mismatch')
            sources = []
            for p, old in zip(source['parts'], parts_doc['parts']):
                require(p.get('kind') == 'repository', 'Invalid part source kind')
                require(p['path'] == str(Path(m['path']).parent / old['path']) and p['bytes'] == old['bytes'] and p['sha256'] == old['sha256'], 'Part order or identity mismatch')
                path = safe_file(repository, p['path'])
                check_pin(path, p)
                sources.append(path)
            require(sum(p['bytes'] for p in source['parts']) == f['bytes'], 'Part sizes do not sum to ZIP size')
        else:
            kind = source.get('kind')
            require(kind in ('repository', 'candidate'), 'Unknown input source kind')
            path = safe_file(repository if kind == 'repository' else candidate, source['path'])
            check_pin(path, source)
            require(source['bytes'] == f['bytes'] and source['sha256'] == f['sha256'], 'Direct source pin disagrees with output')
            sources = [path]
        source_sets.append((f, sources))
    return source_sets

def prepare(repository, candidate, manifest_path, expected_manifest_sha256, output, prepared_previous=None):
    require(re.fullmatch('[0-9a-f]{64}', expected_manifest_sha256) is not None, 'Expected manifest SHA256 is required')
    require(manifest_path == safe_file(candidate, manifest_path.relative_to(candidate).as_posix()), 'Manifest must be inside candidate directory')
    require(sha256(manifest_path) == expected_manifest_sha256, 'Candidate manifest hash changed')
    manifest = json.loads(manifest_path.read_text())
    meta = manifest['metadata_request']
    metadata_path = safe_file(candidate, meta['path'])
    require(sha256(metadata_path) == meta['sha256'], 'Metadata request hash changed')
    request = json.loads(metadata_path.read_text())
    require(set(request) == {'metadata'}, 'Expected documented metadata wrapper')
    md = request['metadata']
    require(md.get('version') == manifest['candidate_version'] and md.get('upload_type') == 'software' and md.get('license') == 'cc-by-4.0', 'Metadata version/type/license changed')
    require(md.get('title') and md.get('description') and md.get('creators'), 'Metadata required fields are empty')
    sources = plan_sources(manifest, repository, candidate, prepared_previous)
    output = output.absolute()
    require(not output.exists() and not output.is_symlink(), 'Output directory must not exist')
    output = output.parent.resolve(strict=True) / output.name
    require(not output.is_relative_to(repository.resolve()), 'Place prepared files outside the repository')
    output.mkdir(parents=False)
    outputs = []
    try:
        for record, inputs in sources:
            target = output / record['filename']
            partial = output / (record['filename'] + '.partial')
            h = hashlib.sha256()
            m = hashlib.md5()
            size = 0
            with partial.open('xb') as writer:
                for source in inputs:
                    with source.open('rb') as reader:
                        for block in iter(lambda: reader.read(1024 * 1024), b''):
                            writer.write(block)
                            h.update(block)
                            m.update(block)
                            size += len(block)
            require(size == record['bytes'] and h.hexdigest() == record['sha256'] and m.hexdigest() == record['md5'], 'Prepared bytes differ: ' + record['filename'])
            partial.rename(target)
            outputs.append({k: record[k] for k in ('filename', 'role', 'bytes', 'sha256', 'md5')})
        # This is a request file for the API/UI, not a tenth Zenodo attachment.
        shutil.copyfile(metadata_path, output / 'METADATA_REQUEST.json')
        require(sha256(output / 'METADATA_REQUEST.json') == meta['sha256'], 'Copied metadata changed')
        receipt = {'status': 'PASS_LOCAL_FILE_PREPARATION', 'scope': 'Twelve additions prepared and hash-verified locally; ten historical files preserved; no network request or publication.', 'concept_doi': manifest['concept_doi'], 'candidate_version': manifest['candidate_version'], 'manifest_sha256': expected_manifest_sha256, 'metadata_request_sha256': meta['sha256'], 'outputs': outputs, 'inherited_files': manifest['inherited_files'], 'remote_mutations': 0}
        (output / 'PREPARATION_RECEIPT.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
        return receipt
    except Exception as e:
        (output / 'PREPARATION_FAILED.json').write_text(json.dumps({'status': 'FAILED_LOCAL_PREPARATION', 'message': str(e), 'completed_files': outputs, 'remote_mutations': 0}, indent=2) + '\n')
        raise

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository-root', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--expected-manifest-sha256', required=True)
    parser.add_argument('--prepared-previous-dir', type=Path,
                        help='Reuse nine byte-pinned previously prepared files instead of rebuilding their multipart ZIPs')
    args = parser.parse_args()
    candidate = Path(__file__).absolute().parent
    receipt = prepare(args.repository_root.absolute(), candidate, candidate / 'FILE_MANIFEST.json', args.expected_manifest_sha256, args.output_dir,
                      args.prepared_previous_dir.absolute() if args.prepared_previous_dir else None)
    print(json.dumps({'status': receipt['status'], 'new_files': len(receipt['outputs']), 'remote_mutations': 0}))

if __name__ == '__main__':
    main()
