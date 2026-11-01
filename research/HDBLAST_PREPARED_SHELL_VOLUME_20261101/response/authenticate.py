"""Standard-library-only response input guards; no scientific imports."""
import hashlib,json
from pathlib import Path,PurePosixPath

class Failure(ValueError):
 def __init__(self,witness):super().__init__(witness);self.witness=witness

def need(ok,witness):
 if not ok:raise Failure(witness)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def regular(p):return p.is_file() and not p.is_symlink() and all(not x.is_symlink()for x in p.absolute().parents)

def relative(s):
 p=PurePosixPath(s)
 need(isinstance(s,str) and bool(s) and p.as_posix()==s and not p.is_absolute()and '..'not in p.parts and s!='.','CANONICAL_RELATIVE_PATH')
 return Path(*p.parts)

def closed(root,pin):
 need(root.is_dir()and not root.is_symlink(),'OLD_CAPSULE_ROOT')
 manifest=root/'MANIFEST.json';need(regular(manifest)and sha(manifest)==pin,'OLD_MANIFEST_PIN')
 rows=json.loads(manifest.read_text())['files']
 if isinstance(rows,dict):rows=[{'path':k,**v}for k,v in rows.items()]
 seen={'MANIFEST.json'}
 for row in rows:
  rel=relative(row['path']);need(str(rel)not in seen,'OLD_DUPLICATE_MEMBER');seen.add(str(rel));p=root/rel
  need(regular(p)and p.stat().st_size==row['bytes']and sha(p)==row['sha256'],'OLD_MEMBER_BYTES')
 actual=set()
 for p in root.rglob('*'):
  need(not p.is_symlink(),'OLD_SYMLINK')
  if p.is_file():actual.add(str(p.relative_to(root)))
  else:need(p.is_dir(),'OLD_SPECIAL_FILE')
 need(actual==seen,'OLD_CLOSED_INVENTORY')
 return {'manifest_sha256':pin,'members_including_manifest':len(seen)}

def authenticate(repo,cap):
 pinsfile=cap/'response/SOURCE_PINS.json';need(regular(pinsfile),'RESPONSE_PINS_REGULAR')
 pins=json.loads(pinsfile.read_text());old={name:closed(repo/'research'/relative(name),pin)for name,pin in sorted(pins['source_capsules'].items())}
 selected={}
 for name,row in sorted(pins['new_files'].items()):
  p=cap/relative(name);need(regular(p)and p.stat().st_size==row['bytes']and sha(p)==row['sha256'],'NEW_SOURCE_BYTES')
  selected[name]=row
 retained=json.loads((repo/'research/HDBLAST_FORCED_RESPONSE_BOUNDS_20261023/golden/BOUNDS.json').read_text())['exact']
 need(retained['F_H_upper']=='17387513156573'and retained['full_G_H_upper']=='3851400','SOURCE_CAP_BINDING')
 premises=json.loads((repo/'research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json').read_text())
 need(premises['coefficient_norm_upper']=={'A':'1/10000','B':'1/10000','D':'1/1000000','M':'1/2000'}and premises['mean_gram']==[1,2,1,3,1],'OPERATOR_CAP_BINDING')
 return {'source_pins_sha256':sha(pinsfile),'old_closures':old,'new_files':selected}

def output_guard(output,sources):
 p=output.absolute();need(all(not q.is_symlink()for q in [p,*p.parents]),'OUTPUT_SYMLINK')
 p=p.resolve();need(not p.exists(),'FRESH_OUTPUT_DIRECTORY')
 need(all(p!=r.resolve()and r.resolve()not in p.parents for r in sources),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
