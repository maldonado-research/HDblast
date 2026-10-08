#!/usr/bin/env python3
"""Copy only frozen additive publication preparations; no real-witness capture."""
from pathlib import Path
import hashlib
import json
import stat

W3 = Path('/workspace/hdblast-research-work/continuation-trajectory-20261008')
HERE = W3/'publication-delivery-candidate'
PACKET = HERE/'hdblast/publication/2026.10.08-retained-endpoint-flow'
SOURCE = W3/'release-assets/publication-verifier'
SEALED = W3/'release-assets/sealed-candidate-001'
MANIFEST_SHA = 'f0f094a5e68f2ed02b7d5810cb928ed178783327eadaaf4fac167976fbb34ca1'
VERIFIER_SHA = 'a6d7e0db5f8a68cb893c37e24fde2aa2e37d12dad96ec13ad3c9983570b9e9bf'

def sha(raw): return hashlib.sha256(raw).hexdigest()
def read(path):
    require = path.is_file() and all(not p.is_symlink() for p in (path,*path.parents))
    if not require or not stat.S_ISREG(path.stat().st_mode) or path.stat().st_nlink != 1: raise ValueError('Real frozen source required')
    return path.read_bytes()
def copy(raw,name):
    destination = PACKET/name
    if destination.exists():
        if read(destination) != raw: raise ValueError('Existing preparation differs: '+name)
    else:
        destination.parent.mkdir(parents=True,exist_ok=True); destination.write_bytes(raw)
    return {'bytes':len(raw),'sha256':sha(raw)}

raw = read(SOURCE/'MANIFEST.json')
if sha(raw) != MANIFEST_SHA: raise ValueError('Independent accepted verifier closure differs')
manifest = json.loads(raw); leaves = {}
for name,pin in manifest['files'].items():
    raw = read(SOURCE/name)
    if len(raw) != pin['bytes'] or sha(raw) != pin['sha256']: raise ValueError('Verifier source closure leaf differs:'+name)
    leaves[name] = copy(raw,name)
if leaves['verify_publication.py']['sha256'] != VERIFIER_SHA: raise ValueError('Accepted verifier source differs')
leaves['VERIFIER_SOURCE_CLOSURE_MANIFEST.json'] = copy(read(SOURCE/'MANIFEST.json'),'VERIFIER_SOURCE_CLOSURE_MANIFEST.json')
leaves['FINAL_VERIFIER_PREPARATION_RECEIPT.json'] = copy(read(SOURCE/'FINAL_PREPARATION_RECEIPT.json'),'FINAL_VERIFIER_PREPARATION_RECEIPT.json')

for name in ('HDBLAST_LATER_STATE_AND_BULK_RESPONSE_PUBLICATION_SUPPLEMENT_20261008.zip',
             'HDBLAST_LATER_STATE_AND_BULK_RESPONSE_OVERVIEW_20261008.md',
             'ASSEMBLY_AND_SEAL_RECEIPT.json','CONTAINED_LEAF_MANIFEST.json',
             'EXTERNAL_FULL_REPLAY_REFERENCE.json','OMITTED_LEAVES_AND_INPUT_DEPENDENCIES.json'):
    leaves['assets/'+name] = copy(read(SEALED/name),'assets/'+name)
for source,name in (
    (SEALED/'SEALED_UPLOAD_INVENTORY.json','proof/inventory/FINAL_UPLOAD_INVENTORY.json'),
    (SEALED/'SEALED_METADATA_MODERN.json','proof/inventory/FINAL_METADATA_MODERN.json'),
    (W3/'release-assets/FINAL_RELEASE_INPUTS.json','preparation/FINAL_RELEASE_INPUTS.json'),
    (W3/'root/SCIENTIFIC_RELEASE_ACCEPTANCE.json','review/SCIENTIFIC_RELEASE_ACCEPTANCE.json'),
    (W3/'root/INDEPENDENT_RELEASE_ASSET_REVIEW.json','review/INDEPENDENT_RELEASE_ASSET_REVIEW.json'),
    (W3/'root/INDEPENDENT_31_FILE_VERIFIER_REVIEW.json','review/INDEPENDENT_31_FILE_VERIFIER_REVIEW.json'),
    (W3/'root/INDEPENDENT_VERIFIER_CONTROLS_NORMAL.json','review/INDEPENDENT_VERIFIER_CONTROLS_NORMAL.json'),
    (W3/'root/INDEPENDENT_VERIFIER_CONTROLS_OPTIMIZED.json','review/INDEPENDENT_VERIFIER_CONTROLS_OPTIMIZED.json'),
    (W3/'publication-planning/new_edition_controller_v2.py','controller/new_edition_controller_v2.py'),
    (W3/'publication-planning/CONTROLLER_V2_ADAPTATION.diff','controller/CONTROLLER_V2_ADAPTATION.diff'),
    (W3/'publication-planning/CONTROLLER_V2_FINAL_DELIVERY.json','controller/CONTROLLER_V2_FINAL_DELIVERY.json'),
    (W3/'publication-independent-review/FINAL_V2_ACCEPTANCE.json','controller/FINAL_V2_ACCEPTANCE.json')):
    leaves[name] = copy(read(source),name)
state = {'schema_version':1,'status':'STATIC_PREPARATION_COPIED_REAL_PUBLICATION_WITNESS_PENDING',
         'files':leaves,'verifier_sha256':VERIFIER_SHA,'verifier_closure_manifest_sha256':MANIFEST_SHA,
         'live_journal_copied':False,'new_record_or_DOI_announced':False,'repository_edits':0,
         'remote_writes':0,'network_calls':0,'array_decodes':0,'physical_source_callbacks':0}
(HERE/'STATIC_PREPARATION_MANIFEST.json').write_text(json.dumps(state,sort_keys=True,indent=2)+'\n')
print(json.dumps({'status':state['status'],'copied_files':len(leaves),'copied_bytes':sum(p['bytes'] for p in leaves.values())}))
