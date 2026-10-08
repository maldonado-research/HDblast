"""Execute the exact independently reviewed bootstrap bytes, captured once."""
import sys
if not sys.flags.isolated or not sys.flags.dont_write_bytecode:
    raise SystemExit('Independent invocation requires Python -I -B')
import hashlib
import os
from pathlib import Path
import stat
expected='438a0478fa185f46ca06372036caa7aee1f9e9659a2778b9d1a9e98651122f0a'
path=Path(__file__).absolute().with_name('trusted_bootstrap_replay.py')
fd=os.open(path,os.O_RDONLY|os.O_NOFOLLOW|os.O_CLOEXEC)
try:
    before=os.fstat(fd)
    if not stat.S_ISREG(before.st_mode) or before.st_nlink!=1 or before.st_size!=5229:
        raise SystemExit('Reviewed bootstrap file identity/type differs')
    raw=os.read(fd,65536);after=os.fstat(fd)
    if (before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=(after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns):
        raise SystemExit('Bootstrap changed during capture')
finally:os.close(fd)
if len(raw)!=5229 or hashlib.sha256(raw).hexdigest()!=expected:
    raise SystemExit('Independent immutable bootstrap hash mismatch')
sys.argv[0]=str(path)
exec(compile(raw,str(path),'exec'),{'__name__':'__main__','__file__':str(path)})
