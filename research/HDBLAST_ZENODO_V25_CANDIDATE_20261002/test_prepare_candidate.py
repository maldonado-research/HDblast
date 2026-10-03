#!/usr/bin/env python3
"""Local synthetic review of prepare_candidate.py; no scientific or network calls."""
import argparse
import copy
import hashlib
import importlib.util
import json
import shutil
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True


def digest(data, kind='sha256'):
    return hashlib.new(kind, data).hexdigest()


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')
    return digest(path.read_bytes())


def fixture(template, directory):
    root = Path(tempfile.mkdtemp(prefix='case-', dir=directory))
    repo, candidate, output = root/'repository', root/'candidate', root/'prepared'
    repo.mkdir()
    candidate.mkdir()
    manifest = copy.deepcopy(template)
    inherited = {}
    for i, row in enumerate(manifest['inherited_files']):
        data = f'synthetic inherited file {i}\n'.encode()
        path = repo/'historical'/row['filename']
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(data)
        inherited[path] = data
    expected = {}
    for i, row in enumerate(manifest['new_files']):
        data = f'synthetic candidate {i}: {row["role"]}\n'.encode()
        expected[row['filename']] = data
        row.update(bytes=len(data), sha256=digest(data), md5=digest(data, 'md5'))
        source = row['source']
        if source['kind'] == 'concatenate':
            subdir = repo/f'parts-{i}'
            subdir.mkdir()
            chunks = [data[:7], data[7:17], data[17:]]
            parts, declared = [], []
            for j, chunk in enumerate(chunks):
                name = f'payload.part{j:02}'
                path = subdir/name
                path.write_bytes(chunk)
                part = {'path':name,'bytes':len(chunk),'sha256':digest(chunk)}
                declared.append(part)
                parts.append(dict(part, kind='repository', path=f'parts-{i}/{name}'))
            parts_doc = {'zip_bytes':len(data),'zip_sha256':digest(data),'parts':declared}
            pin = write_json(subdir/'PARTS.json', parts_doc)
            row['source'] = {'kind':'concatenate','parts':parts,'parts_manifest':{'path':f'parts-{i}/PARTS.json','sha256':pin}}
        else:
            name = f'input-{i}.bin'
            (candidate if source['kind'] == 'candidate' else repo).joinpath(name).write_bytes(data)
            row['source'] = {'kind':source['kind'],'path':name,'bytes':len(data),'sha256':digest(data)}
    metadata = {'metadata':{'version':'2026.10.02-v25','upload_type':'software','license':'cc-by-4.0','title':'Synthetic local fixture','description':'No remote publication.','creators':[{'name':'Synthetic fixture'}]}}
    metadata_path = candidate/'METADATA.json'
    manifest['metadata_request']['sha256'] = write_json(metadata_path, metadata)
    manifest_path = candidate/'FILE_MANIFEST.json'
    pin = write_json(manifest_path, manifest)
    return {'root':root,'repo':repo,'candidate':candidate,'output':output,'manifest':manifest,'manifest_path':manifest_path,'pin':pin,'metadata_path':metadata_path,'metadata':metadata,'expected':expected,'inherited':inherited}


def repin(f):
    f['pin'] = write_json(f['manifest_path'], f['manifest'])


def repin_metadata(f):
    f['manifest']['metadata_request']['sha256'] = write_json(f['metadata_path'], f['metadata'])
    repin(f)


def source_path(f, index=0):
    source = f['manifest']['new_files'][index]['source']
    return (f['candidate'] if source['kind']=='candidate' else f['repo'])/source['path']


def mutate_file(f):
    p = source_path(f)
    p.write_bytes(p.read_bytes().replace(b'synthetic', b'Synthetic', 1))


def mutate_metadata(f):
    f['metadata_path'].write_bytes(f['metadata_path'].read_bytes()+b' ')


def symlink_file(f):
    p = source_path(f)
    outside = f['root']/'external.bin'
    p.rename(outside)
    p.symlink_to(outside)


def symlink_parent(f):
    p = source_path(f)
    directory = f['root']/'external-dir'
    directory.mkdir()
    p.rename(directory/p.name)
    (f['repo']/'linked-dir').symlink_to(directory, target_is_directory=True)
    f['manifest']['new_files'][0]['source']['path'] = f'linked-dir/{p.name}'
    repin(f)


def symlink_repo(f):
    link = f['root']/'linked-repository'
    link.symlink_to(f['repo'], target_is_directory=True)
    f['repo'] = link


def symlink_metadata(f):
    external = f['root']/'external-metadata.json'
    f['metadata_path'].rename(external)
    f['metadata_path'].symlink_to(external)


def symlink_manifest(f):
    external = f['root']/'external-manifest.json'
    f['manifest_path'].rename(external)
    f['manifest_path'].symlink_to(external)


def swap_parts(f):
    parts = f['manifest']['new_files'][2]['source']['parts']
    parts[0], parts[1] = parts[1], parts[0]
    repin(f)


def alter_part(f):
    p = f['repo']/f['manifest']['new_files'][2]['source']['parts'][0]['path']
    p.write_bytes(p.read_bytes().upper())


def alter_parts_manifest(f):
    p = f['repo']/f['manifest']['new_files'][2]['source']['parts_manifest']['path']
    p.write_bytes(p.read_bytes()+b' ')


def existing_output(f):
    f['output'].mkdir()
    (f['output']/'existing.txt').write_text('existing output must be preserved\n')


def output_symlink(f):
    other = f['root']/'external-output'
    other.mkdir()
    f['output'].symlink_to(other, target_is_directory=True)


def broken_output_symlink(f):
    f['output'].symlink_to(f['root']/'missing-output', target_is_directory=True)


def output_inside_repo(f):
    f['output'] = f['repo']/'prepared'


def output_inside_repo_via_dotdot(f):
    directory = f['root']/'normal-parent'
    directory.mkdir()
    f['output'] = directory/'..'/'repository'/'prepared'


def output_inside_repo_via_symlink(f):
    alias = f['root']/'repository-alias'
    alias.symlink_to(f['repo'], target_is_directory=True)
    f['output'] = alias/'prepared'


def change_source_path(f, path):
    f['manifest']['new_files'][0]['source']['path'] = path
    repin(f)


def sidecar_collision(f, filename):
    first = f['manifest']['new_files'][0]
    data = f['expected'].pop(first['filename'])
    first['filename'] = filename
    f['expected'][filename] = data
    repin(f)


def mutate_inherited(f, field, value):
    f['manifest']['inherited_files'][0][field] = value
    repin(f)


def invalid_output_pin(f, field, value):
    f['manifest']['new_files'][2][field] = value
    repin(f)


def empty_concatenation(f):
    row = f['manifest']['new_files'][2]
    row.update(bytes=0, sha256=digest(b''), md5=digest(b'', 'md5'))
    row['source']['parts'] = []
    parts_manifest = f['repo']/row['source']['parts_manifest']['path']
    row['source']['parts_manifest']['sha256'] = write_json(parts_manifest, {'zip_bytes':0,'zip_sha256':digest(b''),'parts':[]})
    f['expected'][row['filename']] = b''
    repin(f)


def run_case(module, f, name, mutator, expectation):
    if mutator:
        mutator(f)
    output_before = f['output'].exists() or f['output'].is_symlink()
    existing_output_contents = {p.relative_to(f['output']).as_posix():p.read_bytes() for p in f['output'].rglob('*') if p.is_file()} if f['output'].is_dir() else {}
    error, receipt = None, None
    try:
        receipt = module.prepare(f['repo'],f['candidate'],f['manifest_path'],f['pin'],f['output'])
    except Exception as exc:
        error = f'{type(exc).__name__}: {exc}'
    unchanged = all(path.read_bytes()==data for path, data in f['inherited'].items())
    existing_output_preserved = all((f['output']/filename).read_bytes()==data for filename,data in existing_output_contents.items())
    outputs_match = bool(receipt) and all((f['output']/filename).read_bytes()==data for filename,data in f['expected'].items())
    if expectation == 'pass':
        expected_receipt_rows = [{key:row[key] for key in ('filename','role','bytes','sha256','md5')} for row in f['manifest']['new_files']]
        expected_directory = set(f['expected']) | {'METADATA_REQUEST.json','PREPARATION_RECEIPT.json'}
        exact_inventory = {p.name for p in f['output'].iterdir()} == expected_directory if f['output'].is_dir() else False
        passed = error is None and outputs_match and unchanged and exact_inventory and receipt['outputs']==expected_receipt_rows and receipt['inherited_files']==f['manifest']['inherited_files'] and len(receipt['outputs'])==9 and len(receipt['inherited_files'])==10 and receipt['remote_mutations']==0 and receipt['status']=='PASS_LOCAL_FILE_PREPARATION' and 'no network request or publication' in receipt['scope']
    elif expectation == 'reject_preflight':
        passed = error is not None and (output_before or not f['output'].exists()) and unchanged and existing_output_preserved
    else:
        raise RuntimeError(expectation)
    row = {'case':name,'expectation':expectation,'passed':passed,'exception':error,'output_created':not output_before and f['output'].exists(),'outputs_match':outputs_match,'inherited_unchanged':unchanged,'existing_output_preserved':existing_output_preserved}
    if receipt:
        row.update(status=receipt['status'],remote_mutations=receipt['remote_mutations'],prepared_count=len(receipt['outputs']))
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--helper',type=Path,required=True)
    parser.add_argument('--template',type=Path,required=True)
    parser.add_argument('--results',type=Path,required=True)
    parser.add_argument('--fixtures-root',type=Path,required=True)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('candidate_helper_under_review', args.helper)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    template = json.loads(args.template.read_text())
    args.fixtures_root.mkdir(parents=True,exist_ok=True)
    cases = [
        ('normal_preparation',None,'pass'),
        ('wrong_expected_manifest_hash',lambda f:f.update(pin='0'*64),'reject_preflight'),
        ('manifest_changed_without_repin',lambda f:f['manifest_path'].write_bytes(f['manifest_path'].read_bytes()+b' '),'reject_preflight'),
        ('metadata_changed_without_repin',mutate_metadata,'reject_preflight'),
        ('metadata_invalid_version',lambda f:(f['metadata']['metadata'].update(version='invalid'),repin_metadata(f)),'reject_preflight'),
        ('metadata_empty_creators',lambda f:(f['metadata']['metadata'].update(creators=[]),repin_metadata(f)),'reject_preflight'),
        ('source_hash_changed',mutate_file,'reject_preflight'),
        ('source_missing',lambda f:source_path(f).unlink(),'reject_preflight'),
        ('source_symlink',symlink_file,'reject_preflight'),
        ('source_parent_symlink',symlink_parent,'reject_preflight'),
        ('repository_root_symlink',symlink_repo,'reject_preflight'),
        ('metadata_symlink',symlink_metadata,'reject_preflight'),
        ('manifest_symlink',symlink_manifest,'reject_preflight'),
        ('parts_order_changed',swap_parts,'reject_preflight'),
        ('part_hash_changed',alter_part,'reject_preflight'),
        ('parts_manifest_changed',alter_parts_manifest,'reject_preflight'),
        ('source_dotdot_path',lambda f:change_source_path(f,'../external.bin'),'reject_preflight'),
        ('source_absolute_path',lambda f:change_source_path(f,str(source_path(f))),'reject_preflight'),
        ('source_backslash_path',lambda f:change_source_path(f,'bad\\name.bin'),'reject_preflight'),
        ('output_already_exists',existing_output,'reject_preflight'),
        ('output_symlink_exists',output_symlink,'reject_preflight'),
        ('output_broken_symlink_exists',broken_output_symlink,'reject_preflight'),
        ('output_inside_repository',output_inside_repo,'reject_preflight'),
        ('output_inside_repository_via_dotdot',output_inside_repo_via_dotdot,'reject_preflight'),
        ('output_inside_repository_via_parent_symlink',output_inside_repo_via_symlink,'reject_preflight'),
        ('output_metadata_sidecar_collision',lambda f:sidecar_collision(f,'METADATA_REQUEST.json'),'reject_preflight'),
        ('output_receipt_sidecar_collision',lambda f:sidecar_collision(f,'PREPARATION_RECEIPT.json'),'reject_preflight'),
        ('inherited_policy_changed',lambda f:mutate_inherited(f,'policy','DELETE'),'reject_preflight'),
        ('inherited_md5_invalid',lambda f:mutate_inherited(f,'md5','invalid'),'reject_preflight'),
        ('inherited_filename_unsafe',lambda f:mutate_inherited(f,'filename','../historical.bin'),'reject_preflight'),
        ('concatenation_zero_bytes',empty_concatenation,'reject_preflight'),
        ('output_concat_sha256_invalid',lambda f:invalid_output_pin(f,'sha256','invalid'),'reject_preflight'),
        ('output_concat_md5_invalid',lambda f:invalid_output_pin(f,'md5','invalid'),'reject_preflight'),
    ]
    results = [run_case(module,fixture(template,args.fixtures_root),name,mutator,expected) for name,mutator,expected in cases]
    report = {'helper':str(args.helper),'helper_sha256':digest(args.helper.read_bytes()),'optimization':sys.flags.optimize,'count':len(results),'passed':sum(row['passed'] for row in results),'failures':[row['case'] for row in results if not row['passed']],'results':results,'scope':'Synthetic byte copying and concatenation only; no scientific command or network call.'}
    write_json(args.results,report)
    print(json.dumps({key:report[key] for key in ('helper_sha256','optimization','count','passed','failures')}))
    return 0 if not report['failures'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
