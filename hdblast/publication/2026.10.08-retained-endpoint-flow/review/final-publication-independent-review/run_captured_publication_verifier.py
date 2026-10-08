"""Run captured, externally pinned read-only witness verifier under -I -B."""
import sys
if not sys.flags.isolated or not sys.dont_write_bytecode:
    raise SystemExit('isolated no-bytecode Python required')
from pathlib import Path
import hashlib
import os
import stat

source = Path(__file__).resolve().parent / 'verifier-closure/verify_publication.py'
fd = os.open(source, os.O_RDONLY | os.O_NOFOLLOW)
try:
    before = os.fstat(fd)
    if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1 or before.st_size != 37342:
        raise SystemExit('invalid captured verifier source')
    with os.fdopen(os.dup(fd), 'rb') as stream:
        raw = stream.read(37343)
    after = os.fstat(fd)
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        raise SystemExit('verifier source changed')
finally:
    os.close(fd)
if hashlib.sha256(raw).hexdigest() != 'a6d7e0db5f8a68cb893c37e24fde2aa2e37d12dad96ec13ad3c9983570b9e9bf':
    raise SystemExit('external verifier SHA256 mismatch')
sys.argv[0] = str(source)
exec(compile(raw, str(source), 'exec'), {'__name__': '__main__', '__file__': str(source)})
