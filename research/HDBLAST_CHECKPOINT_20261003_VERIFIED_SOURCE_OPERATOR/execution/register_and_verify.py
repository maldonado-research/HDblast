"""Custodian-only registration/full public blob readback; no source evaluation."""
import argparse
import base64
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

from registration_guard import load, require, sha

REPOSITORY = 'maldonado-research/HDblast'


def write(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, sort_keys=True, indent=2)
        handle.write('\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('operation', choices=('register', 'verify'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--repository-root', type=Path, required=True)
    parser.add_argument('--branch', required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repository = args.repository_root.resolve()
    relative = root.relative_to(repository).as_posix()
    def run(*command):
        return subprocess.check_output(command, cwd=repository)
    def api(path):
        return load(run('gh', 'api', 'repos/'+REPOSITORY+('/'+path if path else '')))
    if args.operation == 'register':
        review = load((root/'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json').read_bytes())
        require(review['status'] == 'GO_FOR_PROSPECTIVE_PUBLIC_FREEZE', 'Independent GO required')
        files = {}
        for path in sorted(root.rglob('*')):
            require(not path.is_symlink(), 'Checkpoint symlink forbidden')
            if path.is_file():
                require(path.name != 'FULL_REGISTRATION.json' and
                        path.suffix not in ('.pyc','.pyo','.so') and
                        '__pycache__' not in path.parts, 'Unexpected existing registration/cache')
                raw = path.read_bytes()
                files[path.relative_to(root).as_posix()] = {'sha256':sha(raw),'bytes':len(raw)}
        for name, pin in review['reviewed_source_sha256'].items():
            require(name in files and files[name]['sha256'] == pin, 'Reviewed source changed')
        contract = load((root/'REGISTRATION_CONTRACT.json').read_bytes())
        registration = {'schema_version':1,
            'created_utc':datetime.now(timezone.utc).isoformat(),
            'scientific_outcome_before_freeze':'UNCOMPUTED_SOURCE_OPERATOR',
            'new_target_evaluations_before_freeze':{'primary':0,'independent':0},
            'files':files,'frozen_configuration':contract}
        write(root/'FULL_REGISTRATION.json', registration)
        print(json.dumps({'files':len(files),'registration_sha256':sha(
            (root/'FULL_REGISTRATION.json').read_bytes())}), flush=True)
        return
    commit = run('git','rev-parse','HEAD').decode().strip()
    require(api('')['private'] is False, 'Public registration repository required')
    require(api('git/ref/heads/'+args.branch)['object']['sha'] == commit,
            'Public branch does not bind local freeze commit')
    tree = api('git/trees/'+commit+'?recursive=1')
    require(tree.get('truncated') is False, 'Complete remote Git tree required')
    mapping = {v['path']:v for v in tree['tree']}
    regraw = (root/'FULL_REGISTRATION.json').read_bytes()
    registration = load(regraw)
    files = dict(registration['files'])
    files['FULL_REGISTRATION.json'] = {'sha256':sha(regraw),'bytes':len(regraw)}
    remote_names = {name[len(relative)+1:] for name,item in mapping.items()
                    if name.startswith(relative+'/') and item['type'] == 'blob'}
    require(remote_names == set(files), 'Unregistered or missing remote checkpoint member')
    def check(item):
        name, pin = item
        raw = (root/name).read_bytes()
        require(sha(raw) == pin['sha256'] and len(raw) == pin['bytes'], 'Local registered member changed')
        blob = hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        obj = mapping[relative+'/'+name]
        require(obj['type'] == 'blob' and obj['mode'] == '100644' and obj['sha'] == blob,
                'Remote Git blob/path differs')
        fetched = api('git/blobs/'+blob)
        require(fetched['encoding'] == 'base64', 'Unexpected API blob encoding')
        decoded = base64.b64decode(fetched['content'])
        require(sha(decoded) == pin['sha256'] and len(decoded) == pin['bytes'],
                'Full remote byte readback differs')
        return name, {**pin,'git_blob':blob,
            'independent_public_blob_download_sha256_verified':True}
    checked = {}
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = [pool.submit(check, item) for item in files.items()]
        for future in as_completed(futures):
            name, result = future.result()
            checked[name] = result
            if len(checked)%20 == 0:
                print(json.dumps({'verified_files':len(checked),'total':len(files)}),flush=True)
    receipt = {'schema_version':1,'status':'PASS_REMOTE_REGISTERED_SOURCE_GO',
        'public_repository':'https://github.com/'+REPOSITORY,
        'public_branch':args.branch,'public_path':relative,
        'freeze_commit':commit,'registration_sha256':sha(regraw),
        'verified_utc':datetime.now(timezone.utc).isoformat(),
        'input_scope':'NO_RETAINED_ARRAYS_SOURCE_OPERATOR_ONLY',
        'new_target_evaluations_before_public_verification':{'primary':0,'independent':0},
        'prior_saved_data_and_failures_known':True,
        'all_registered_public_blobs_independently_downloaded':True,
        'remote_tree_sha':tree['sha'],'remote_files_verified':len(checked),
        'remote_registered_files':checked}
    write(root/'FREEZE_RECEIPT.json',receipt)
    print(json.dumps({'status':receipt['status'],'freeze_commit':commit,
        'registration_sha256':sha(regraw),'remote_files_verified':len(checked),
        'receipt_sha256':sha((root/'FREEZE_RECEIPT.json').read_bytes())}),flush=True)


if __name__ == '__main__':
    main()
