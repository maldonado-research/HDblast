"""Harmless regression for directories hidden from Path.rglob in source closure."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
BASE=Path(__file__).resolve().parent
source=Path(sys.argv[1]);registration=Path(sys.argv[2]);label=sys.argv[3]
out=BASE/label;out.mkdir();candidate=out/'candidate';shutil.copytree(source,candidate)
for copied in [candidate,*candidate.rglob('*')]:copied.chmod(0o755 if copied.is_dir() else 0o644)
package=candidate/'source'/'later_source';package.mkdir()
marker='HIDDEN_UNREGISTERED_PACKAGE_EXECUTED'
(package/'__init__.py').write_text('print('+repr(marker)+')\n')
package.chmod(0o111)
program='''import sys
from pathlib import Path
import hashlib
import types
root=Path(sys.argv[1]);reg=Path(sys.argv[2])
module=types.ModuleType('review_guard')
exec(compile((root/'execution/registration_guard.py').read_bytes(),'captured guard','exec'),module.__dict__)
module.authenticate(root,reg,hashlib.sha256(reg.read_bytes()).hexdigest(),True)
print('FULL_CLOSURE_AUTHENTICATION_ACCEPTED')
sys.path.insert(0,str(root/'source'))
import later_source
'''
try:
 p=subprocess.run([sys.executable,'-I','-B','-c',program,str(candidate),str(registration)],capture_output=True,text=True)
 log=p.stdout+p.stderr;(out/'probe.log').write_text(log)
 result={'status':'FAIL_UNREADABLE_DIRECTORY_HIDES_IMPORTABLE_UNREGISTERED_PACKAGE' if marker in log else 'PASS_UNREADABLE_DIRECTORY_REJECTED',
  'authentication_accepted':'FULL_CLOSURE_AUTHENTICATION_ACCEPTED' in log,'hidden_marker_executed':marker in log,
  'package_directory_probe_mode':'0111','exit_code':p.returncode,'real_decodes':0,'physical_source_calls':0,
  'registered_guard_sha256':hashlib.sha256((source/'execution/registration_guard.py').read_bytes()).hexdigest(),
  'registration_sha256':hashlib.sha256(registration.read_bytes()).hexdigest()}
 (out/'RESULT.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
finally:
 package.chmod(0o555)
