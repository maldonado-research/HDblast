"""Copy public immutable input blobs and authenticate original retained pins."""
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path(__file__).resolve().parent.parent
OLD=Path('/workspace/hdblast-research-work/quadratic-shell-forcing-20261015')
GIT=OLD/'publish-evidence'
COMMIT='97fe9d2828b93b237b571e221a2422fd4366ef32'
def sha(raw): return hashlib.sha256(raw).hexdigest()
def blob(path): return subprocess.check_output(['git','-C',str(GIT),'show',COMMIT+':'+path])
prior=json.loads((OLD/'inputs/INPUT_PINS.json').read_text())
entries=[]
for entry in prior['files']:
    raw=blob(entry['original_public_source']['path'])
    assert len(raw)==entry['bytes'] and sha(raw)==entry['sha256']
    assert raw==(OLD/entry['path']).read_bytes()
    path=ROOT/entry['path']; path.parent.mkdir(parents=True,exist_ok=True); path.write_bytes(raw)
    entries.append(dict(entry,verified_at_commit=COMMIT))
extra=[
 ('BACKGROUND_DERIVATIVE_BOUNDS.json','HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015','archive/BACKGROUND_DERIVATIVE_BOUNDS.json'),
 ('INERTIAL_WAVEFORM_MANIFEST.json','HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015','MANIFEST.json'),
 ('PHYSICAL_SHELL_RESPONSE.md','HDBLAST_PHYSICAL_SHELL_RESPONSE_20261015','PHYSICAL_SHELL_RESPONSE.md'),
 ('PHYSICAL_SHELL_RESPONSE_MANIFEST.json','HDBLAST_PHYSICAL_SHELL_RESPONSE_20261015','MANIFEST.json'),
 ('FULL_CONTACT_MANIFEST.json','HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014','MANIFEST.json'),
 ('GH_BOUNDARY_SOURCE_AND_CONSTRAINT_COMPLETION.md','HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014','boundary/GH_BOUNDARY_SOURCE_AND_CONSTRAINT_COMPLETION.md'),
 ('FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md','HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014','FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md'),
 ('COMPLETE_L4_SCALAR_BULK.md','HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014','astra-bulk/COMPLETE_L4_SCALAR_BULK.md'),
 ('FULL_L4_SCALAR_SOURCE_AND_OPERATOR.md','HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014','sol-bulk/FULL_L4_SCALAR_SOURCE_AND_OPERATOR.md')]
for name,cap,rel in extra:
    source='research/'+cap+'/'+rel; raw=blob(source); dst='inputs/reference/'+name
    path=ROOT/dst;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
    entries.append({'path':dst,'bytes':len(raw),'sha256':sha(raw),'original_public_source':{'repository':'maldonado-research/HDblast','commit':COMMIT,'capsule':cap,'path':source,'url':'https://github.com/maldonado-research/HDblast/blob/'+COMMIT+'/'+source},'verified_at_commit':COMMIT})
for manifest_name,cap in [('PHYSICAL_SHELL_RESPONSE_MANIFEST.json','HDBLAST_PHYSICAL_SHELL_RESPONSE_20261015'),('FULL_CONTACT_MANIFEST.json','HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014'),('INERTIAL_WAVEFORM_MANIFEST.json','HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015')]:
    m=json.loads((ROOT/'inputs/reference'/manifest_name).read_text()); by={x['path']:x for x in m['files']}
    for e in entries:
        src=e['original_public_source']
        if src['capsule']==cap and src['path'].rsplit('/',1)[-1]!='MANIFEST.json':
            rel=src['path'].split('research/'+cap+'/',1)[1]
            assert e['sha256']==by[rel]['sha256']
record={'schema':'hdblast-coordinate-shell-sources-input-pins-v1','status':'PASS','immutable_public_commit':COMMIT,'scope':'Exact public Git blobs; previous original pins preserved; no prior scientific program imported or executed. Physical reconstruction and contact notes checked against their retained manifests.','external_inertial_manifest_sha256':prior['external_inertial_manifest_sha256'],'external_full_contact_manifest_sha256':prior['external_full_contact_manifest_sha256'],'files':entries}
(ROOT/'inputs/INPUT_PINS.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'status':'PASS','input_files':len(entries),'immutable_commit':COMMIT}))
