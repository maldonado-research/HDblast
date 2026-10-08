#!/usr/bin/env python3
"""Mutation controls for authenticated saved witnesses; no network/source work."""
from pathlib import Path
import argparse
import copy
import hashlib
import importlib.util
import json
import re
import sys

sys.dont_write_bytecode = True

def require(value, message):
    if not value:
        raise ValueError(message)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verifier', required=True)
    parser.add_argument('--verifier-sha256', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--witness-manifest', required=True)
    parser.add_argument('--witness-manifest-sha256', required=True)
    parser.add_argument('--expected-record-id', required=True)
    parser.add_argument('--inventory-sha256', required=True)
    parser.add_argument('--metadata-sha256', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    path = Path(args.verifier).absolute()
    require(path.is_file() and all(not p.is_symlink() for p in (path, *path.parents)),
            'Real verifier source required')
    require(re.fullmatch('[0-9a-f]{64}', args.verifier_sha256 or '')
            and hashlib.sha256(path.read_bytes()).hexdigest() == args.verifier_sha256,
            'External verifier pin differs before import')
    spec = importlib.util.spec_from_file_location('authenticated_saved_witness_verifier', path)
    v = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(v)
    original = v.verify_directory(args.root, args.witness_manifest, args.witness_manifest_sha256,
                                  args.expected_record_id, args.inventory_sha256, args.metadata_sha256)
    require(original['status'] == 'PASS_OFFLINE_27_FILE_PUBLICATION_WITNESSES',
            'Original saved witnesses must pass before mutation controls')
    root = Path(args.root).absolute()
    manifest = v.load(Path(args.witness_manifest).read_bytes())
    roles = manifest['witnesses']
    data = {k: v.load(v.safe_file(root, name).read_bytes())
            for k, name in roles.items() if k != 'controller_journal'}
    journal = [v.load(line) for line in v.safe_file(root, roles['controller_journal']).read_bytes().splitlines()]
    preserved = {rid: {key: v.load(v.safe_file(root, name).read_bytes())
                      for key, name in group.items()}
                 for rid, group in manifest['preserved_records'].items()}
    baseline = [data['inventory'], data['metadata'], data['new_public_record'], data['new_public_files'],
                data['latest_public_record'], data['verified_public'], data['prior_full_stream'],
                journal, preserved, args.expected_record_id, args.inventory_sha256, args.metadata_sha256]
    rejections = []
    def rejected(name, mutate):
        values = copy.deepcopy(baseline)
        mutate(values)
        try:
            v.verify_values(*values)
        except (ValueError, TypeError, KeyError):
            rejections.append(name)
        else:
            raise ValueError('Accepted corrupted saved-witness mutation: ' + name)
    def first_file(document):
        return next(iter(v.file_map(document).values()))
    def first_stream(values):
        url = 'https://zenodo.org/api/records/' + args.expected_record_id + '/files'
        markers = [i for i, event in enumerate(values[7]) if event.get('kind') == 'GET'
                   and event.get('url') == url and event.get('status') == 200
                   and event.get('authenticated_request') is False]
        return next(event for event in values[7][markers[-1] + 1:]
                    if event.get('kind') == 'FULL_CONTENT_STREAM_VERIFIED')
    rejected('wrong externally expected new identity', lambda x: x.__setitem__(9, '1'))
    rejected('stale latest version identity', lambda x: x[4].__setitem__('id', '23114217'))
    rejected('new public record marked unpublished', lambda x: x[2].__setitem__('is_published', False))
    rejected('new public record marked draft', lambda x: x[2].__setitem__('is_draft', True))
    rejected('new public record wrong family', lambda x: x[2]['parent'].__setitem__('id', '22922927'))
    rejected('new public record wrong DOI', lambda x: x[2]['pids']['doi'].__setitem__('identifier', '10.5281/zenodo.1'))
    rejected('new public embedded files deleted', lambda x: x[2].__setitem__('files', []))
    rejected('new public listing deleted', lambda x: x.__setitem__(3, {'entries': []}))
    rejected('new public size corrupted', lambda x: first_file(x[3]).__setitem__('size', 1))
    rejected('new public size type corrupted', lambda x: first_file(x[3]).__setitem__('size', True))
    rejected('new public MD5 corrupted', lambda x: first_file(x[3]).__setitem__('checksum', 'md5:'+'0'*32))
    rejected('new public file marked pending', lambda x: first_file(x[3]).__setitem__('status', 'pending'))
    rejected('complete metadata title changed', lambda x: x[2]['metadata'].__setitem__('title', 'Altered'))
    rejected('complete metadata references deleted', lambda x: x[2]['metadata'].pop('references'))
    rejected('complete metadata description altered', lambda x: x[2]['metadata'].__setitem__('description', 'Proof of a Big Bang cause'))
    rejected('complete metadata unexpected field', lambda x: x[2]['metadata'].__setitem__('unregistered', True))
    rejected('complete metadata custom fields deleted', lambda x: x[2].__setitem__('custom_fields', {}))
    rejected('complete metadata access altered', lambda x: x[2]['access'].__setitem__('files', 'restricted'))
    rejected('verification receipt marked draft', lambda x: x[5].__setitem__('published', False))
    rejected('verification receipt inventory hash corrupted', lambda x: x[5].__setitem__('inventory_sha256', '0'*64))
    rejected('verification receipt metadata hash corrupted', lambda x: x[5].__setitem__('metadata_sha256', '0'*64))
    rejected('verification receipt count corrupted', lambda x: x[5].__setitem__('file_count', 26))
    rejected('verification receipt total corrupted', lambda x: x[5].__setitem__('total_bytes', 1))
    rejected('SHA256 basis row deleted', lambda x: x[5]['content'].pop())
    rejected('SHA256 basis row duplicated', lambda x: x[5]['content'].append(copy.deepcopy(x[5]['content'][0])))
    rejected('SHA256 basis hash corrupted', lambda x: x[5]['content'][0].__setitem__('sha256', '0'*64))
    rejected('SHA256 basis unrecognized', lambda x: x[5]['content'][0].__setitem__('basis', 'SELF_ATTESTED'))
    rejected('addition falsely reuses inherited stream', lambda x: x[5]['content'][-1].__setitem__('basis', 'SAME_IMMUTABLE_FILE_AND_CONTENT_VERSION_AS_PRIOR_FULL_STREAM'))
    rejected('fresh public stream bytes truncated', lambda x: first_stream(x)['observed'].__setitem__('bytes', 1))
    rejected('fresh public stream MD5 corrupted', lambda x: first_stream(x)['observed'].__setitem__('md5', '0'*32))
    rejected('fresh public stream SHA256 corrupted', lambda x: first_stream(x)['observed'].__setitem__('sha256', '0'*64))
    rejected('complete stream journal deleted', lambda x: x.__setitem__(7, []))
    rejected('prior full-stream checksum corrupted', lambda x: x[6]['rows'][0]['observed'].__setitem__('sha256', '0'*64))
    rejected('prior full-stream row deleted', lambda x: x[6]['rows'].pop())
    rejected('preserved prior record missing', lambda x: x[8].pop('22347452'))
    for rid in ('23114217', '22347452', '23111008'):
        rejected('preserved '+rid+' embedded inventory deleted',
                 lambda x, rid=rid: x[8][rid]['record'].__setitem__('files', []))
        rejected('preserved '+rid+' file listing deleted',
                 lambda x, rid=rid: x[8][rid].__setitem__('files', {'entries': []}))
        rejected('preserved '+rid+' metadata altered',
                 lambda x, rid=rid: x[8][rid]['record']['metadata'].__setitem__('title', 'Altered'))
    rejected('unsealed upload inventory', lambda x: x[0].__setitem__('sealed_upload_manifest', False))
    rejected('upload inventory byte total changed', lambda x: x[0].__setitem__('total_bytes_final', 1))
    rejected('upload inventory inherited row deleted', lambda x: x[0]['inherited_files'].pop())
    rejected('upload inventory new row duplicated', lambda x: x[0]['new_files'].append(copy.deepcopy(x[0]['new_files'][0])))
    require(v.verify_values(*baseline)['status'] == original['status'], 'Controls mutated originals')
    output = Path(args.output).absolute()
    require(not output.exists() and all(not p.is_symlink() for p in (output, *output.parents)), 'Fresh control receipt required')
    output.parent.mkdir(parents=True, exist_ok=True)
    result = {'status': 'PASS_OFFLINE_AUTHENTICATED_SAVED_WITNESS_MUTATION_CONTROLS',
              'scope': 'IN_MEMORY_CORRUPTIONS_OF_EXTERNALLY_PINNED_SAVED_WITNESSES_ONLY',
              'positive_checks': 2, 'rejected_mutants': len(rejections), 'rejections': rejections,
              'record_id': args.expected_record_id, 'witness_manifest_sha256': args.witness_manifest_sha256,
              'verifier_sha256': args.verifier_sha256, 'controls_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'network_calls': 0, 'controller_imports': 0, 'numeric_imports': 0,
              'source_callbacks': 0, 'array_decodes': 0, 'remote_writes': 0}
    with output.open('x') as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({key:result[key] for key in ('status','positive_checks','rejected_mutants','network_calls')}))

if __name__ == '__main__':
    main()
