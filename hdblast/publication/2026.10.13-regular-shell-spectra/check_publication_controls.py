#!/usr/bin/env python3
"""Eight focused negative controls for actual saved publication witnesses; offline."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile

spec = importlib.util.spec_from_file_location('publication_verifier',Path(__file__).with_name('verify_publication.py'))
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)

def run(root, pin):
    root = Path(root).resolve()
    v.verify(root,pin)
    checks = []
    # This positive control documents the narrowly permitted server decoration.
    expected = v.read(root/'FINAL_METADATA_MODERN.json')
    decorated = copy.deepcopy(expected)
    decorated['metadata']['related_identifiers'][-1]['relation_type']['title'] = {'en':'Server display wording'}
    decorated['access']['status'] = 'open'
    v.compare_metadata(decorated,expected)
    checks.append({'name':'recognized_display_decoration_allowed','outcome':'PASS'})
    public_index = v.read(root/'PUBLIC_STREAM_OBSERVATIONS.json')['preceding_public_file_list_event']['source_journal_index']
    cases = [
        ('old_33_file_count_rejected','PUBLIC_INVENTORY.json','FILE_COUNT_OR_DUPLICATE_DIFFERS',
         lambda o: o.update(file_count=33)),
        ('scientific_text_change_rejected','NEW_PUBLIC_RECORD.json','COMPLETE_METADATA_DIFFERS',
         lambda o: o['metadata'].update(description=o['metadata']['description']+' Altered scientific claim.')),
        ('authenticated_stream_round_rejected','PUBLIC_STREAM_OBSERVATIONS.json','PUBLIC_STREAM_ROUND_SCOPE_DIFFERS',
         lambda o: o['preceding_public_file_list_event'].update(authenticated_request=True)),
        ('public_round_before_publish_rejected','PUBLIC_STREAM_OBSERVATIONS.json','PUBLIC_STREAM_ROUND_PRECEDES_PUBLISH_INTENT',
         lambda o: o['preceding_public_file_list_event'].update(source_journal_index=0)),
        ('write_readback_index_collision_rejected','WRITE_SUMMARY.json','WRITE_AND_PUBLIC_READBACK_EVENT_INDEX_COLLISION',
         lambda o: o['write_responses'][0].update(source_journal_index=public_index)),
        ('complete_stream_sha256_change_rejected','PUBLIC_STREAM_OBSERVATIONS.json','COMPLETE_STREAM_OBSERVATION_DIFFERS_',
         lambda o: o['files'][0]['observed'].update(sha256='0'*64)),
        ('source_leaf_omission_rejected','IMMUTABLE_PUBLIC_SOURCE_READBACK.json','SOURCE_FILE_LIST_COUNT_DIFFERS',
         lambda o: o['files'].pop()),
        ('prior_record_change_rejected','preserved/after_RECORD_23294979.json','PRESERVED_RECORD_CHANGED_23294979',
         lambda o: o['metadata'].update(title=o['metadata']['title']+' altered')),
    ]
    for name, relative, error, mutate in cases:
        with tempfile.TemporaryDirectory(prefix='publication-control-') as directory:
            candidate = Path(directory)/'proof'
            shutil.copytree(root,candidate)
            path = candidate/relative
            obj = v.read(path)
            mutate(obj)
            path.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
            manifest = v.read(candidate/'WITNESS_MANIFEST.json')
            manifest['files'][relative] = {'bytes':path.stat().st_size,'sha256':v.sha(path.read_bytes())}
            manifest_path = candidate/'WITNESS_MANIFEST.json'
            manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
            try:
                v.verify(candidate,v.sha(manifest_path.read_bytes()))
            except ValueError as rejected:
                v.require(str(rejected).startswith(error),'CONTROL_FAILED_AT_UNEXPECTED_GUARD_'+name+': '+str(rejected))
            else:
                raise ValueError('INVALID_WITNESS_ACCEPTED_'+name)
        checks.append({'name':name,'outcome':'EXPECTED_REJECTION','guard':error})
    return {'status':'PASS_PUBLICATION_WITNESS_CONTROLS','valid_packet_checks':1,'allowed_decoration_controls':1,
            'negative_controls':len(cases),'network_calls':0,'manifest_sha256':pin,'checks':checks,
            'scope':'Actual saved packet passes first; each mutated packet is re-manifested to exercise content guards beyond outer hashes.'}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',required=True)
    p.add_argument('--manifest-sha256',required=True)
    p.add_argument('--output')
    args = p.parse_args()
    result = run(args.root,args.manifest_sha256)
    encoded = json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:
        output = Path(args.output).resolve()
        v.require(not output.is_relative_to(Path(args.root).resolve()), 'OUTPUT_MUST_BE_OUTSIDE_WITNESS_ROOT')
        output.write_text(encoded)
    print(encoded,end='')

if __name__ == '__main__':
    main()
