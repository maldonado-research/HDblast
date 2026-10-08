"""Build local manufactured-only registration; supplies no physical GO."""
from pathlib import Path
import hashlib
import json
import sys
import types

sys.dont_write_bytecode=True
root=Path(__file__).resolve().parents[1]
guard_path=root/'execution/registration_guard.py'
guard=types.ModuleType('manufactured_registration_guard');guard.__file__=str(guard_path)
exec(compile(guard_path.read_bytes(),str(guard_path),'exec'),guard.__dict__)
files={}
for relative in sorted(guard.source_tree(root)):
    raw=guard.safe_file(root,relative).read_bytes()
    files[relative]={'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
guard.source_tree(root,files)
contract={'arithmetic_bits':1024,'momentum_degree':2048,'source_degree':24,'source_coefficient_bits':512,
    'source_cells':64,'export_bits':96,'export_gate':'1/1000000000000000000','anchor':'-9/2','endpoints':['-4','-7/2'],
    'selected_decode_members':['k.npy','momentum_weights.npy','observation_eta.npy','u_1.npy','w_1.npy',
                               'u_2.npy','w_2.npy','u_3.npy','w_3.npy'],
    'wall_seconds':900,'rss_kib':524288,'output_bytes':134217728}
data={'schema_version':1,'scope':'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE','files':files,'contract':contract}
raw=(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n').encode()
with Path(sys.argv[1]).open('xb') as out:out.write(raw)
digest=hashlib.sha256(raw).hexdigest()
guard.authenticate(root,Path(sys.argv[1]),digest,True)
print(digest)
