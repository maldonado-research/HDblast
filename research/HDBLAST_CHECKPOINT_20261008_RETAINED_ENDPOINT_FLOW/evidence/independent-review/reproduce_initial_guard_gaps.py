"""Negative-control reproduction; pure metadata, no worker or source import."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
BASE = Path(__file__).resolve().parent
ROOT = BASE.parent / 'implementation'
def digest(raw): return hashlib.sha256(raw).hexdigest()
def write(path, value):
    raw = (json.dumps(value, sort_keys=True, separators=(',', ':'))+'\n').encode()
    with path.open('xb') as out: out.write(raw)
    return digest(raw)

snapshot = BASE/'initial_registration_guard.py.snapshot'
raw = (ROOT/'execution/registration_guard.py').read_bytes()
with snapshot.open('xb') as out: out.write(raw)
namespace = {'__name__': 'reviewed_initial_guard'}
exec(compile(raw, str(snapshot), 'exec'), namespace)
registration = json.loads((BASE.parent/'manufactured-registration-002.json').read_bytes())
results = []

# The all-.py closure stays byte-identical; only consumed JSON is omitted.
del registration['files']['provenance_decoder/INPUT_SPEC.json']
path = BASE/'omitted_input_spec_registration.json'
pin = write(path, registration)
try:
    namespace['authenticate'](ROOT, path, pin, True)
except Exception as error:
    results.append({'case':'omitted_consumed_input_spec','rejected':True,'error':str(error)})
else:
    results.append({'case':'omitted_consumed_input_spec','rejected':False})

# A self-consistent but irrelevant tree must not claim the endpoint scope.
fake = BASE/'fabricated_subset_root'
fake.mkdir()
payload = b'"""No endpoint implementation exists in this fabricated tree."""\n'
(fake/'irrelevant.py').write_bytes(payload)
registration['files'] = {'irrelevant.py':{'bytes':len(payload),'sha256':digest(payload)}}
path = BASE/'fabricated_subset_registration.json'
pin = write(path, registration)
try:
    namespace['authenticate'](fake, path, pin, True)
except Exception as error:
    results.append({'case':'irrelevant_subset_source_tree','rejected':True,'error':str(error)})
else:
    results.append({'case':'irrelevant_subset_source_tree','rejected':False})
result = {'status':'INITIAL_GUARD_NEGATIVE_CONTROLS', 'source_guard_sha256':digest(raw),
          'controls':results,'retained_decodes':0,'physical_source_calls':0}
write(BASE/'INITIAL_GUARD_GAPS.json',result)
print(json.dumps(result,sort_keys=True))
