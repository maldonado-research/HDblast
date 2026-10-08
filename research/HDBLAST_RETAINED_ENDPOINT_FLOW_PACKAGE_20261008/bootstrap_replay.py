#!/usr/bin/env python3
"""Trusted-source companion bootstrap: capture external-pinned control bytes.

Use this code from the independently reviewed release/workflow, not an unknown
unverified extracted copy. It creates no GO, installs nothing and decodes no
array. The captured helper then authenticates the complete source bundle.
"""
import sys
if not sys.flags.isolated or not sys.flags.dont_write_bytecode:
    raise SystemExit('Trusted bootstrap requires Python -I -B')
sys.dont_write_bytecode = True
import argparse
import hashlib
import os
from pathlib import Path
import re
import stat
import types

parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
parser.add_argument('--package', required=True)
parser.add_argument('--payload-manifest-sha256', required=True)
parser.add_argument('--helper-sha256', required=True)
parser.add_argument('--replay-pins-sha256', required=True)
parser.add_argument('helper_args', nargs=argparse.REMAINDER)
args = parser.parse_args()
if args.helper_args and args.helper_args[0] != '--':
    raise SystemExit('After -- supply fixed helper operation and original repository root')
operation = argparse.ArgumentParser(prog='authenticated replay operation', allow_abbrev=False)
mode = operation.add_mutually_exclusive_group(required=True)
mode.add_argument('--verify-only', action='store_true')
mode.add_argument('--output-root')
scope = operation.add_mutually_exclusive_group(required=True)
scope.add_argument('--repository-root')
scope.add_argument('--manufactured-control', action='store_true')
remaining = args.helper_args[1:] if args.helper_args else []
allowed = {'--verify-only', '--output-root', '--repository-root', '--manufactured-control'}
counts = {}
for item in remaining:
    if item.startswith('--'):
        name = item.split('=', 1)[0]
        if name not in allowed:
            raise SystemExit('Bootstrap trust roots and unknown helper options cannot be overridden: ' + name)
        counts[name] = counts.get(name, 0) + 1
        if counts[name] != 1:
            raise SystemExit('Repeated helper operation option forbidden: ' + name)
selected = operation.parse_args(remaining)
root = Path(os.path.abspath(args.package))
if root.resolve() != root or not root.is_dir() or any(p.is_symlink() for p in (root, *root.parents)):
    raise SystemExit('Real package root without symlink ancestors required')
controls = {'replay_checkpoint.py': args.helper_sha256,
            'PAYLOAD_MANIFEST.json': args.payload_manifest_sha256,
            'REPLAY_PINS.json': args.replay_pins_sha256}
captured = {}
root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
identity = lambda st: (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns, st.st_nlink)
try:
    root_before = os.fstat(root_fd)
    for name, expected in controls.items():
        if re.fullmatch('[0-9a-f]{64}', expected) is None:
            raise SystemExit('Independent final external SHA256 required: ' + name)
        fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC, dir_fd=root_fd)
        try:
            before = os.fstat(fd)
            if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1 or before.st_size > 20 * 1024 * 1024:
                raise SystemExit('Bounded regular single-link control file required: ' + name)
            data = bytearray()
            while chunk := os.read(fd, 65536):
                data.extend(chunk)
                if len(data) > 20 * 1024 * 1024:
                    raise SystemExit('Control grew beyond bound: ' + name)
            after = os.fstat(fd)
            named = os.stat(name, dir_fd=root_fd, follow_symlinks=False)
            if identity(before) != identity(after) or identity(after) != identity(named):
                raise SystemExit('Control changed while hashing: ' + name)
            if len(data) != before.st_size or hashlib.sha256(data).hexdigest() != expected:
                raise SystemExit('External control SHA256 mismatch: ' + name)
            captured[name] = bytes(data)
        finally:
            os.close(fd)
    named_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
    try:
        named = os.fstat(named_fd)
        if (root_before.st_dev, root_before.st_ino) != (named.st_dev, named.st_ino):
            raise SystemExit('Package root replaced')
    finally:
        os.close(named_fd)
finally:
    os.close(root_fd)
entry = root / 'replay_checkpoint.py'
module = types.ModuleType('hdblast_authenticated_replay_helper')
module.__file__ = str(entry)
sys.modules[module.__name__] = module
# Execute the already captured authenticated source; never reload helper disk.
exec(compile(captured['replay_checkpoint.py'], str(entry), 'exec'), module.__dict__)
argv = ['--package', str(root), '--payload-manifest-sha256', args.payload_manifest_sha256,
        '--helper-sha256', args.helper_sha256, '--replay-pins-sha256', args.replay_pins_sha256,
        *(['--verify-only'] if selected.verify_only else ['--output-root', selected.output_root]),
        *(['--manufactured-control'] if selected.manufactured_control else ['--repository-root', selected.repository_root])]
raise SystemExit(module.main(argv))
