import hashlib,json,zipfile
from pathlib import Path
ROOT=Path('/workspace/hdblast-research-work/publish-HDblast-20261014')
CAP=ROOT/'research/HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011'
EXT=ROOT/'hdblast/continuation/2026.10.11-scalar-response-monotonicity/SCIENTIFIC_GITHUB_RELEASE.json'
PUB=ROOT/'hdblast/publication/2026.10.11-scalar-response-monotonicity/proof/PUBLIC_INVENTORY.json'
READ=ROOT/'hdblast/publication/2026.10.11-scalar-response-monotonicity/proof/IMMUTABLE_PUBLIC_SOURCE_READBACK.json'
def sha(b):return hashlib.sha256(b).hexdigest()
def authenticate():
 manifest=json.loads((CAP/'MANIFEST.json').read_text()); external=json.loads(EXT.read_text()); read=json.loads(READ.read_text()); inv=json.loads(PUB.read_text())
 mh=sha((CAP/'MANIFEST.json').read_bytes())
 if mh!='917480172ae98b2bf39f29d37fe0ad3cf279ad0d8508c98efd469efdba26c0ec' or mh!=external['capsule_manifest_sha256']:raise ValueError('external manifest pin mismatch')
 expected={'MANIFEST.json'}|{r['path'] for r in manifest['files']}
 actual={str(p.relative_to(CAP)) for p in CAP.rglob('*') if p.is_file()}
 if actual!=expected:raise ValueError(('inventory mismatch',actual^expected))
 rows=[]
 for r in manifest['files']:
  p=CAP/r['path']; b=p.read_bytes()
  if p.is_symlink() or len(b)!=r['bytes'] or sha(b)!=r['sha256']: raise ValueError(('payload mismatch',r['path']))
  rows.append(r)
 for p,h in manifest['frozen_pins'].items():
  if sha((CAP/p).read_bytes())!=h:raise ValueError(('frozen mismatch',p))
 readrows=next(v for v in read.values() if isinstance(v,list) and v and isinstance(v[0],dict) and 'path' in v[0])
 for r in readrows:
  b=(CAP/r['path']).read_bytes()
  if len(b)!=r['bytes'] or sha(b)!=r['sha256']:raise ValueError(('immutable source mismatch',r['path']))
 archive=next(r for r in inv['files'] if r['filename']=='HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011.zip')
 archive_path=Path(__file__).with_name('HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011.zip')
 archive_replayed=False
 if archive_path.exists():
  ab=archive_path.read_bytes()
  if len(ab)!=archive['bytes'] or sha(ab)!=archive['sha256'] or hashlib.md5(ab).hexdigest()!=archive['md5']:raise ValueError('external archive bytes mismatch')
  with zipfile.ZipFile(archive_path) as z:
   names=z.namelist();prefix=CAP.name+'/'
   if len(names)!=len(set(names)) or {n[len(prefix):] for n in names if n.startswith(prefix)}!=expected or len(names)!=len(expected):raise ValueError('archive exact inventory mismatch')
   for n in names:
    if not n.startswith(prefix) or z.read(n)!=(CAP/n[len(prefix):]).read_bytes():raise ValueError(('archive payload mismatch',n))
  archive_replayed=True
 return {'status':'PASS_COMPLETE_SELECTED_OLDER_BRANCH_CLOSURE','payload_count':len(rows),'payload_bytes':sum(r['bytes'] for r in rows),'manifest_sha256':mh,'external_manifest_pin':{'path':str(EXT),'sha256':sha(EXT.read_bytes()),'source_commit':external['source_commit']},'independent_immutable_source_readback':{'path':str(READ),'sha256':sha(READ.read_bytes()),'checked_files':len(readrows)},'external_archive_pin':archive,'archive_bytes_replayed':archive_replayed,'archive_pin_source':{'path':str(PUB),'sha256':sha(PUB.read_bytes())},'selected_capsule':str(CAP),'all_payloads':rows,'scope':'Retained complete local capsule authenticated against independently retained external manifest and immutable-public-source pins. Archive bytes are checked when present, including all44 exact capsule members. No numerical module imports occurred.'}
if __name__=='__main__':
 out=authenticate(); Path(__file__).with_name('BRANCH_DEPENDENCY_AUTHENTICATION.json').write_text(json.dumps(out,indent=2)+'\n'); print(out['status'],out['payload_count'])
