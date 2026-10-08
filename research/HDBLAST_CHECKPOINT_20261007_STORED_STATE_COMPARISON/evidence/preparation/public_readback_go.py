"""Independently download registered public Git bytes before real calculation."""
import argparse
import base64
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import urllib.request

PREFIX = 'research/HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON'
API = 'https://api.github.com/repos/maldonado-research/HDblast/'

def require(value, message):
    if not value:
        raise ValueError(message)

class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError('Git API redirect refused')

def get(path, anonymous=False):
    token = os.environ.get('GH_TOKEN') or os.environ.get('GITHUB_TOKEN')
    headers = {'Accept': 'application/vnd.github+json', 'User-Agent': 'HDBLAST-independent-public-byte-review',
               'X-GitHub-Api-Version': '2022-11-28'}
    if token and not anonymous:
        headers['Authorization'] = 'Bearer ' + token
    request = urllib.request.Request(API + path if path else API.rstrip('/'), headers=headers)
    with urllib.request.build_opener(NoRedirect()).open(request, timeout=60) as response:
        require(response.status == 200, 'Git API response not 200')
        raw = response.read(32 * 1024 * 1024 + 1)
        require(len(raw) <= 32 * 1024 * 1024, 'Git API response oversized')
        return json.loads(raw)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--candidate', required=True)
    p.add_argument('--commit', required=True)
    p.add_argument('--registration-sha256', required=True)
    p.add_argument('--guard-sha256', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    require(re.fullmatch('[0-9a-f]{40}', a.commit), 'Full immutable commit required')
    root = Path(a.candidate).resolve()
    registration_raw = (root / 'FULL_REGISTRATION.json').read_bytes()
    require(hashlib.sha256(registration_raw).hexdigest() == a.registration_sha256, 'External registration pin differs')
    registration = json.loads(registration_raw)
    files = registration['files']
    guard_path = root / 'execution/registration_guard.py'
    require(hashlib.sha256(guard_path.read_bytes()).hexdigest() == a.guard_sha256,
            'Independent guard pin differs')
    spec = importlib.util.spec_from_file_location('reviewed_stored_registration_guard', guard_path)
    guard = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(guard)
    guard.authenticate_local(root, a.registration_sha256)
    review_name = 'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json'
    require(review_name in files, 'Registered final review required')
    review_raw = (root / review_name).read_bytes()
    require({'bytes': len(review_raw), 'sha256': hashlib.sha256(review_raw).hexdigest()} == files[review_name],
            'Pre-freeze review changed after authentication')
    review = guard.load(review_raw)
    require(review['status'] == 'GO_FOR_PROSPECTIVE_PUBLIC_FREEZE', 'Final review not GO')
    zeros = dict.fromkeys(('registered_source_callbacks', 'retained_array_decodes', 'stored_state_comparisons',
                          'physical_trajectories', 'likelihood_evaluations'), 0)
    require(guard.zero_counts(review['new_target_evaluations_before_freeze']), 'Premature or ill-typed physical counters')
    reviewed = review['reviewed_source_sha256']
    require(set(reviewed) == {name for name in files if name.endswith('.py')}, 'Complete reviewed Python closure required')
    require(all(files[name]['sha256'] == pin for name, pin in reviewed.items()), 'Reviewed Python source changed')
    repository = get('')
    require(repository['full_name'] == 'maldonado-research/HDblast' and repository['private'] is False,
            'Repository is not publicly visible')
    anonymous = get('contents/' + PREFIX + '/README.md?ref=' + a.commit, anonymous=True)
    require(anonymous['encoding'] == 'base64' and
            base64.b64decode(anonymous['content']) == (root / 'README.md').read_bytes(),
            'Anonymous immutable public README readback failed')
    names = ['FULL_REGISTRATION.json', *sorted(files)]
    def download(name):
        local = (root / name).read_bytes()
        meta = get('contents/' + PREFIX + '/' + name + '?ref=' + a.commit)
        require(meta['type'] == 'file' and meta['path'] == PREFIX + '/' + name, 'Public path identity differs')
        git_sha = hashlib.sha1(b'blob ' + str(len(local)).encode() + b'\0' + local).hexdigest()
        require(meta['sha'] == git_sha and meta['size'] == len(local), 'Public Git path/blob differs')
        blob = get('git/blobs/' + git_sha)
        require(blob['sha'] == git_sha and blob['encoding'] == 'base64', 'Public blob encoding/identity differs')
        raw = base64.b64decode(blob['content'], validate=False)
        require(blob['size'] == len(raw) and raw == local, 'Downloaded public bytes differ')
        item = {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest(), 'git_blob': git_sha,
                'immutable_contents_url': API + 'contents/' + PREFIX + '/' + name + '?ref=' + a.commit,
                'independent_public_blob_download_sha256_verified': True}
        if name != 'FULL_REGISTRATION.json':
            require({'bytes': item['bytes'], 'sha256': item['sha256']} == files[name], 'Registration payload pin differs')
        if name == review_name:
            require(raw == review_raw, 'Downloaded review differs from authenticated review decision')
            downloaded_review = guard.load(raw)
            require(downloaded_review['status'] == 'GO_FOR_PROSPECTIVE_PUBLIC_FREEZE' and
                    guard.zero_counts(downloaded_review['new_target_evaluations_before_freeze']) and
                    downloaded_review['reviewed_source_sha256'] == reviewed,
                    'Downloaded pre-freeze review does not support GO')
        return name, item
    downloaded = {}
    with ThreadPoolExecutor(max_workers=6) as executor:
        for index, (name, item) in enumerate(executor.map(download, names), 1):
            downloaded[name] = item
            if index % 20 == 0:
                print(json.dumps({'public_registered_files_verified': index, 'total': len(names)}), flush=True)
    receipt = {'status': 'PASS_REMOTE_REGISTERED_STORED_STATE_COMPARISON_GO',
               'input_scope': 'ALL_FOUR_RETAINED_INCOMING_CAPSULES_AT_FIXED_ANCHOR',
               'public_repository': 'https://github.com/maldonado-research/HDblast',
               'repository_checkpoint_path': PREFIX, 'freeze_commit': a.commit,
               'registration_sha256': a.registration_sha256,
               'all_registered_public_blobs_independently_downloaded': True,
               'new_target_evaluations_before_public_verification': zeros,
               'remote_registered_files': downloaded,
               'remote_files_verified': len(downloaded),
               'public_visibility': {'repository_private': False, 'anonymous_immutable_README_bytes_verified': True},
               'external_peer_review': False,
               'review_scope': 'Root acceptance of independent internal source, mathematical and security reviews'}
    raw = (json.dumps(receipt, sort_keys=True, indent=2) + '\n').encode()
    output = Path(a.output)
    require(root not in output.resolve().parents, 'External receipt required')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('xb') as f:
        f.write(raw)
    print(json.dumps({'status': receipt['status'], 'files': len(downloaded),
                      'public_go_sha256': hashlib.sha256(raw).hexdigest()}), flush=True)

if __name__ == '__main__':
    main()
