#!/usr/bin/env python3
"""Read-only proof that the local Git public-freeze object binds all frozen bytes."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess

from active_integrity import digest, require, verify_frozen


def verify(root, repository, registration_sha256, freeze_commit):
    inputs = verify_frozen(root, registration_sha256, freeze_commit)
    repository = Path(repository).resolve()
    relative = Path(root).resolve().relative_to(repository).as_posix()
    def git(*args):
        return subprocess.check_output(['git', '-C', str(repository), *args])
    require(git('cat-file', '-t', freeze_commit).strip() == b'commit', 'Freeze object is not a commit')
    names = {'FULL_REGISTRATION.json': registration_sha256, **inputs['frozen_files']}
    for name, pin in names.items():
        require(digest(git('show', freeze_commit + ':' + relative + '/' + name)) == pin,
                'Public freeze Git bytes differ: ' + name)
    return {'status': 'PASS_LOCAL_FREEZE_GIT_BYTES', 'public_freeze_commit': freeze_commit,
            'registration_sha256': registration_sha256, 'git_files_verified': len(names),
            'science_evaluated': False, 'network_operations': 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint-root', type=Path, required=True)
    parser.add_argument('--repository-root', type=Path, required=True)
    parser.add_argument('--registration-sha256', required=True)
    parser.add_argument('--freeze-commit', required=True)
    args = parser.parse_args()
    import json
    print(json.dumps(verify(args.checkpoint_root, args.repository_root,
                            args.registration_sha256, args.freeze_commit), sort_keys=True))


if __name__ == '__main__':
    main()
