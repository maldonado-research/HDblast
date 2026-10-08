#!/usr/bin/env python3
"""Finalize a public-safe packet only after genuine publication and review."""
from pathlib import Path
import hashlib
import json

W3=Path('/workspace/hdblast-research-work/continuation-trajectory-20261008')
HERE=W3/'publication-delivery-candidate'
PACKET=HERE/'hdblast/publication/2026.10.08-retained-endpoint-flow'
REVIEW=W3/'final-publication-independent-review'
VERIFIER='a6d7e0db5f8a68cb893c37e24fde2aa2e37d12dad96ec13ad3c9983570b9e9bf'
WITNESS='534a07bcd2d7e6c5ce1d21642724cccbb51fcdd93f47ff55ba5f526217bae398'
INVENTORY='2c6d0a3d4c4e227c3a8aa5f6e48dd3ff8494ac4e2f1b9045c603cf35617ced3c'
METADATA='43feec5ca6bd9ab8c5c6a2cd925f6f42ff13418bd90b334de7ab604367dc385c'
FINAL_REVIEW='024a6d9d3f757f2b700bff80fcd823c347e7557f43898624afcafd4cdbaa621a'
ORIGINAL_REVIEW_MANIFEST='bf3e3245c682dbd247200d9d8f958bccbba5d2783031bd815c50161e4bc19195'
PRIVATE_BODY='journal-audit/final-observation-001/owner_page_1.json'

def sha(raw):return hashlib.sha256(raw).hexdigest()
def enc(value):return (json.dumps(value,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+'\n').encode()
def put(name,raw):
    path=PACKET/name;path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        if path.read_bytes()!=raw:raise ValueError('Existing packet leaf differs:'+name)
    else:path.write_bytes(raw)
def require(ok,why):
    if not ok:raise ValueError(why)

original_raw=(REVIEW/'FINAL_PUBLICATION_REVIEW_MANIFEST.json').read_bytes()
require(sha(original_raw)==ORIGINAL_REVIEW_MANIFEST,'Final independent review manifest differs')
original=json.loads(original_raw);public_files={};omitted={}
for name,pin in original['files'].items():
    source=REVIEW/name;require(source.is_file() and not source.is_symlink(),'Real independent review source required')
    raw=source.read_bytes();require(len(raw)==pin['bytes'] and sha(raw)==pin['sha256'],'Final review leaf differs:'+name)
    if name==PRIVATE_BODY:
        omitted[name]=dict(pin,reason='PRIVATE_AUTHENTICATED_OWNER_LIST_BODY_RETAINED_LOCALLY_NOT_PUBLISHED')
        require(not (PACKET/'review/final-publication-independent-review'/name).exists(),'Private owner body entered packet')
        continue
    require(not any(word in name.lower() for word in ('owner_page','account_response','access_token','__pycache__')),'Unexpected private response/cache leaf')
    put('review/final-publication-independent-review/'+name,raw);public_files[name]=pin
require(len(public_files)==197 and set(omitted)=={PRIVATE_BODY},'Exact public-safe197 review selection required')
put('review/final-publication-independent-review/FINAL_PUBLICATION_REVIEW_MANIFEST.json',original_raw)
portable={'schema_version':1,'scope':'EXACT197PUBLIC_SAFE_SELECTED_REVIEW_LEAVES_EXCLUDING_THIS_AND_PROVENANCE_OMISSION_MANIFESTS','files':public_files,'file_count':197,'total_bytes':sum(p['bytes'] for p in public_files.values()),'original198_manifest_sha256':ORIGINAL_REVIEW_MANIFEST,'original198_closure_standalone_complete':False,'private_original_body_omitted':True}
put('review/final-publication-independent-review/PORTABLE_PUBLIC_REVIEW_MANIFEST.json',enc(portable))
omission={'schema_version':1,'original_review_manifest_sha256':ORIGINAL_REVIEW_MANIFEST,'original_selected_files':198,'public_selected_files':197,'omitted':omitted,'limitation':'The unchanged original198manifest binds the privately retained full review. Its completeness cannot be authenticated from this public197leaf selection alone. The separate portable manifest binds all197 present selected leaves. The direct31file publication witness closure and offline verifier remain complete and unchanged.'}
put('review/final-publication-independent-review/PRIVATE_BODY_OMISSION.json',enc(omission))
review_raw=(REVIEW/'FINAL_PUBLICATION_REVIEW_RECEIPT.json').read_bytes();require(sha(review_raw)==FINAL_REVIEW,'Final review receipt differs');review=json.loads(review_raw)
require(review['status']=='PASS_INDEPENDENT_FINAL_PUBLICATION_WITNESS_AND_LIVE_PRESERVATION_REVIEW' and review['record_id']=='23244754' and review['file_count']==31 and review['total_bytes']==453646943,'Genuine complete publication review required')
inventory=json.loads((PACKET/'proof/inventory/FINAL_UPLOAD_INVENTORY.json').read_bytes())
external=json.loads((PACKET/'assets/EXTERNAL_FULL_REPLAY_REFERENCE.json').read_bytes())
facts={'schema_version':1,'status':'PASS_PUBLISHED_31_FILE_EDITION_AND_INDEPENDENT_SAVED_WITNESS','record_id':'23244754','record_url':'https://zenodo.org/records/23244754','version_doi':'10.5281/zenodo.23244754','concept_id':'17088132','concept_doi':'10.5281/zenodo.17088132','semantic_version':'2026.10.08-later-state-and-bulk-response','scientific_checkpoint_date_pacific':'2026-10-08','full_public_content_verified_utc':review['publication_times_utc']['public_31_stream_verification_complete'],'publish_response_observed_utc':review['publication_times_utc']['publish_response'],'independent_final_review_sealed_utc':review['sealed_utc'],'file_count':31,'total_bytes':453646943,'parent_record_id':'23228395','inherited_file_count':29,'inherited_bytes':451598148,'new_files':[{k:x[k] for k in ('filename','bytes','md5','sha256')} for x in inventory['new_files']],'inventory_sha256':INVENTORY,'metadata_sha256':METADATA,'witness_manifest_sha256':WITNESS,'verifier_sha256':VERIFIER,'offline_verification_receipt_sha256':'bc4f5872a518f13c3e53a585d9d9cac6c15154bf9bf1aaee44444d712fdff010','independent_final_review_sha256':FINAL_REVIEW,'original_independent_review_manifest_sha256':ORIGINAL_REVIEW_MANIFEST,'portable_public_review_manifest_sha256':sha(enc(portable)),'exact_journal_events':197,'complete_draft_rounds':2,'complete_draft_streams':62,'complete_public_rounds':1,'complete_public_streams':31,'successful_recorded_writes':9,'publish_intents':1,'unknown_failed_truncated_write_outcomes':0,'normal_optimized_saved_witness_receipts_identical':True,'preserved_records':['23228395','23225288','23114217','22347452','23111008'],'preserved_concept_families':{'17088132':['23228395','23225288','23114217','22347452'],'22922927':['23111008']},'scientific_release_status':'PASS_COMPLETE_LATER_STATE_AND_BULK_RESPONSE_REPLAY','later_endpoint_nodes':49152,'later_endpoint_comparisons':98304,'finite_weighted_cases':24,'signed_finite_intervals_excluding_zero':48,'full_replay_reference':external,'compact_supplement_standalone_replay':False,'large_node_exports_omitted':2,'external_original_dependencies':11,'external_original_dependency_bytes':213114777,'source_freeze_commit':'3fc078553d1c9b10f11099870195c7a317951a71','full_science_commit':'52126b98d3e4b24680017fa2537b8449bd75d94f','scientific_limits':review['scientific_limits'],'automatic_GitHub_to_Zenodo_archiving':'HISTORICAL_OFF_UNTOUCHED_NOT_FRESHLY_INSPECTED','publication_limits':['Dated client readbacks; no server-side atomic publication compare-and-swap claim.','No cryptographic hash chain or account-wide publisher exclusion claim.','The independent reviewer fetched public metadata/file lists, not another full set of attachment streams.','Internal AI-assisted reviews are not external peer review or proof-assistant formalization.']}
put('PUBLICATION_FACTS.json',enc(facts))
readme=f'''# Retained endpoint flow and causal bulk response: publication proof

Published [record23244754](https://zenodo.org/records/23244754), DOI **10.5281/zenodo.23244754**, retains the established concept DOI **[10.5281/zenodo.17088132](https://doi.org/10.5281/zenodo.17088132)**. Complete public readback at **2026-10-08T17:52:21.331322+00:00** verifies **31 files /453,646,943 bytes**, every filename, integer size, MD5 and complete-content SHA256, the full approved metadata and latest owner/version/concept identity. Independent final saved-witness and live preservation review passes. This guide was created after genuine publication and completed checks.

All **29 files /451,598,148 bytes** from record23228395 are retained. The additions are the **2,040,745-byte compact scientific supplement** and **8,050-byte original overview**. Scientific date: **8 October2026, America/Los_Angeles**. Semantic version: **2026.10.08-later-state-and-bulk-response**. The publish response was observed at2026-10-08T17:49:19.454259+00:00; actual observations use UTC timestamps.

## Scientific result and limits

The new comparison encloses retained later-state errors relative to exact forced evolution beginning at the same exact represented incoming state at eta=-9/2. All **49,152 capsule-node occurrences**, **98,304 endpoint comparisons** and **24 finite weighted cases** are enclosed. Maximum normalized Cartesian complex L1 errors satisfy **U&lt;2.520e-16** and **W&lt;3.158e-15**. All48 signed finite R/P intervals exclude zero. Exact rational [results](../../../research/HDBLAST_CHECKPOINT_20261008_RETAINED_ENDPOINT_FLOW/RESULTS.md) and the [24-row table](../../../research/HDBLAST_CHECKPOINT_20261008_RETAINED_ENDPOINT_FLOW/EXACT_ENDPOINT_METRICS_24.csv) are authoritative. The complete target export radii are U&lt;7.574e-29 and W&lt;1.263e-28; the1e-18 gate concerns target uncertainty, rather than historical endpoint error. Earlier incoming Bunch-Davies mismatch is separate.

Saved and fresh portable normal/optimized executions reproduce all nine stable scientific files exactly. Independent source-free readback verifies serialized states, representation, enclosures and finite aggregates; it does not independently decode the original archives or reconstruct another flow. Source/flow truth also requires authenticated implementation and reviewed analytic proof. The unsaved continuous historical solver path cannot be recovered from retained endpoints. Full continuous pressure/contact/time, momentum quadrature and ultraviolet certification remain **UNRESOLVED**, historical metric calibration **FAIL**, higher-dimensional Big Bang origin **NOT_ESTABLISHED**, external novelty **NOT_ASSESSED**.

The separate [stable half-space scalar model](../../../research/HDBLAST_CHECKPOINT_20261008_CAUSAL_BULK_RESPONSE/README.md) derives a causal response/noise testbed from a declared flat4+1-dimensional action. An explicit four-dimensional continuum of free fields reproduces its Gaussian response/noise laws, so those observables alone do not identify extra-dimensional geometry. The [incident extension](../../../research/HDBLAST_CHECKPOINT_20261008_CAUSAL_BULK_RESPONSE/incident/README.md) proves unit reflection and vanishing late local boundary energy for the specified smooth pure-continuum packets with both bound-mode projections zero. This static linear single-port mechanism does not permanently capture/reheat the boundary from those packets; the conclusion depends on that action/state. Gravity, cosmological backreaction, nonlinear energy transfer and observational likelihood remain open. KMS equilibrium is a state premise, and raw fixed-k noise needs ultraviolet filtering.

The bounded literature update records12 queries and8 selected primary-source readings, including verified2026 papers. Only original synthesis, compact provenance and original figures are delivered. Raw papers and extracted third-party article bodies are excluded. Internal AI-assisted reviews are not external peer review, proof-assistant formalization or observational confirmation.

## Compact supplement and complete replay

The compact ZIP contains **166 safe leaves**, each read back against its leaf manifest. It retains source/proof, compact actual results and source/error receipts, reviews, both bulk/incident models, bounded literature notes and original figures. It omits two large NODE_CERTIFICATES.jsonl.gz exports and11 external original input dependencies. It is **not a standalone numerical replay package**.

The [complete103-file source-replay package at immutable Git commit52126b98d3e4b24680017fa2537b8449bd75d94f]({external['package_url']}) contains the full node exports and captured-source bootstrap. All103 package leaves were independently streamed from that immutable public commit before release. The11 original dependencies,213,114,777 bytes, remain separately pinned at older commit13a30ef46f5c90c1b01ae83b609db65fd0f8a709. The older commit supplies input bytes; the new commit contains the new science. Authenticate external pins before executing the full package bootstrap and follow its README.

The [external replay reference](assets/EXTERNAL_FULL_REPLAY_REFERENCE.json) and [omission manifest](assets/OMITTED_LEAVES_AND_INPUT_DEPENDENCIES.json) preserve exact input/output identities. Payload manifest:976e92a7a14ecb13ed04cbb9f7ed15fff3a285bd40175a28929dff53677d0752. Helper:5de25511aea5579b7e0ebe94208ee1584d0869ca248fa8939cdbb29ed25ff6fd. Replay pins:1407dc4bcf0a312d1bc23a29e93e3c1149071f4ce2705ccbbeef430e82446d68.

## Publication witnesses and preservation

Inventory SHA256: `{INVENTORY}`. Full metadata SHA256: `{METADATA}`. The [exact witness manifest](proof/WITNESS_MANIFEST.json) SHA256 is **`{WITNESS}`**. Accepted verifier SHA256 is **`{VERIFIER}`**. Proof has31 exact witness leaves plus that manifest,11 roles and5 preserved record groups. Generate verification outputs outside proof/ and manufactured-seeds/.

The terminal journal independently passes **197 events**, **9 successful writes**, **two complete31-file draft stream rounds**, **one publish intent** and **one complete31-file public stream round**. All93 streams check the full453,646,943-byte inventory once per round. No unknown, failed or truncated write occurred in this edition. Normal/optimized saved-witness receipts match exactly. The independently rerun manufactured suite passes **7 positive /175 negative controls per mode**; the separate journal checker passes17 controls per mode.

The final reviewer independently fetched public metadata/file lists before and after publication for23228395,23225288,23114217,22347452 and companion23111008. Full editable metadata, complete file listings and immutable file/version IDs match pinned baselines. The companion remains in concept22922927. The new record's public identity, metadata and file list also match an independent live GET. The reviewer authenticated root-owned31 public attachment streams from saved evidence and did not duplicate them. Complete Notes history, the historical literal Notes **FAIL**, earlier numerical/packaging failures and rejected controls remain preserved.

The original final review198-leaf manifest remains unchanged as dated provenance. A private authenticated owner-list body is omitted from this public packet. The separate [197-leaf portable review manifest](review/final-publication-independent-review/PORTABLE_PUBLIC_REVIEW_MANIFEST.json) binds all present selected review leaves; the [hashed omission record](review/final-publication-independent-review/PRIVATE_BODY_OMISSION.json) identifies the privately retained missing leaf. This public selection cannot fully authenticate the original198-leaf closure without that private original. The direct31-leaf publication proof and offline CI closure remain complete.

These dated client checks do not establish server-side atomic publication compare-and-swap, a cryptographic journal hash chain or account-wide publisher exclusion. Automatic GitHub-to-Zenodo archiving remains a historical **OFF** observation, untouched; no fresh account configuration inspection is claimed. Archived controller/state files are provenance. Do not execute an archived publication controller or reuse a completed publisher state. Startup and saved-witness checks perform no publication action.

## Reproduce the saved-witness check

From this publication directory, authenticate captured source before executing it. Use Python3.12.14 and a fresh external output filename:

```bash
python -I -S -B - <<'PYVERIFY'
from pathlib import Path
import hashlib, sys
p = Path('verify_publication.py').resolve()
raw = p.read_bytes()
if hashlib.sha256(raw).hexdigest() != '{VERIFIER}':
    raise SystemExit('Accepted verifier source differs')
sys.argv = [str(p), '--root', 'proof', '--witness-manifest', 'proof/WITNESS_MANIFEST.json',
    '--witness-manifest-sha256', '{WITNESS}', '--expected-record-id', '23244754',
    '--inventory-sha256', '{INVENTORY}', '--metadata-sha256', '{METADATA}',
    '--output', '/tmp/HDBLAST_FRESH_31_FILE_PUBLICATION_RECEIPT.json']
exec(compile(raw, str(p), 'exec'), {{'__name__': '__main__', '__file__': str(p)}})
PYVERIFY
```

Require PASS_OFFLINE_31_FILE_PUBLICATION_WITNESSES, **31 files /453,646,943 bytes**, with zero network calls, numerical imports, source callbacks and array decodes. Adding -O produces the same receipt. These checks authenticate dated saved observations; current live state requires a separate readback. [PUBLICATION_FACTS.json](PUBLICATION_FACTS.json) provides verified structured facts, and [UPDATE_GUIDE.md](UPDATE_GUIDE.md) gives downstream update instructions.
'''
for a,b in {'record232':'record 232','/453,646':'/ 453,646','/451,598':'/ 451,598','October2026':'October 2026','at2026':'at 2026','All48':'All 48','the1e-18':'the 1e-18','flat4+1':'flat 4+1','records12':'records 12','and8':'and 8','verified2026':'verified 2026','and11':'and 11','complete103-file':'complete 103-file','commit521':'commit 521','All103':'All 103','The11':'The 11','dependencies,213':'dependencies, 213','commit13a':'commit 13a','complete31-file':'complete 31-file','All93':'All 93','full453,646':'full 453,646','/175':'/ 175','passes17':'passes 17','for232':'for 232','companion231':'companion 231','concept229':'concept 229','root-owned31':'root-owned 31','review198-leaf':'review 198-leaf','original198-leaf':'original 198-leaf','direct31-leaf':'direct 31-leaf','has31 exact':'has 31 exact','manifest,11':'manifest, 11','and5':'and 5','Python3.12':'Python 3.12'}.items():readme=readme.replace(a,b)
put('README.md',readme.encode())
guide='''# Verified downstream update contract

Use PUBLICATION_FACTS.json and the frozen proof/inventory files: record 23244754, version DOI 10.5281/zenodo.23244754, concept DOI 10.5281/zenodo.17088132, semantic version 2026.10.08-later-state-and-bulk-response, **31 files / 453,646,943 bytes**. Complete public readback finished at17:52:21 UTC on8 October2026.

Update current GitHub, website, status and registry fields while preserving all dated29-file and older publication evidence. Link the new immutable full science package52126b98d3e4b24680017fa2537b8449bd75d94f and its external replay pins. Original dependencies remain separately pinned to13a30ef46f5c90c1b01ae83b609db65fd0f8a709. The compact supplement is not standalone replay.

Run the exact new31-file saved-witness check in normal and optimized Python, plus authenticated7-positive/175-negative manufactured controls, before announcing matching deployment. Preserve old29-file checks as historical checks. Startup and CI must not run a publisher or alter archiving configuration. OFF is a historical account observation, rather than a fresh inspection.

State the science accurately:49,152 retained nodes,98,304 endpoint comparisons,24 finite cases, nine scientific artifacts reproduced and conditional bulk/replica/incident results. Preserve full continuous pressure/contact/time/momentum/UV UNRESOLVED, metric calibration FAIL, higher-dimensional Big Bang origin NOT_ESTABLISHED, novelty NOT_ASSESSED and historical literal Notes FAIL. Internal independent AI-assisted review is not external peer review.

Preserve unrelated concurrent projects, including shared-site TMD content. Keep the private archive's existing interactive upload tool. Keep proof/ and manufactured-seeds/ exact; put generated verification receipts elsewhere. The unchanged198-leaf original final-review manifest includes one privately retained owner response. Use the separate197-leaf portable review manifest and hashed omission record for the public selection; do not claim the public selection fully contains the original198-leaf closure. Do not reuse the archived completed publisher.
'''
for a,b in {'at17':'at 17','on8':'on 8','October2026':'October 2026','dated29':'dated 29','package521':'package 521','to13a':'to 13a','new31':'new 31','authenticated7':'authenticated 7','old29':'old 29','accurately:49':'accurately: 49',',98,304':', 98,304',',24 finite':', 24 finite','unchanged198':'unchanged 198','separate197':'separate 197','original198':'original 198'}.items():guide=guide.replace(a,b)
put('UPDATE_GUIDE.md',guide.encode())
put('preparation/finalize_packet.py',Path(__file__).read_bytes())
files={f.relative_to(PACKET).as_posix():{'bytes':f.stat().st_size,'sha256':sha(f.read_bytes())} for f in sorted(PACKET.rglob('*')) if f.is_file()}
require('PUBLICATION_DELIVERY_MANIFEST.json' not in files,'Packet already sealed')
packet_manifest={'schema_version':1,'status':'PASS_COMPLETED_PUBLICATION_DELIVERY_PACKET_SEALED','scope':'Exact safe additive public/private delivery packet excluding this externally pinned manifest; original198review closure provenance differs from explicit197public selection','record_id':'23244754','file_count':31,'total_bytes':453646943,'packet_files':len(files),'packet_bytes':sum(pin['bytes'] for pin in files.values()),'witness_manifest_sha256':WITNESS,'verifier_sha256':VERIFIER,'inventory_sha256':INVENTORY,'metadata_sha256':METADATA,'independent_final_review_sha256':FINAL_REVIEW,'public_private_response_bodies':0,'files':files}
put('PUBLICATION_DELIVERY_MANIFEST.json',enc(packet_manifest));digest=sha(enc(packet_manifest))
(HERE/'FINAL_PACKET_SEAL.json').write_bytes(enc({'status':'PASS_FINAL_ADDITIVE_PUBLICATION_PACKET_SEALED','packet_path':str(PACKET),'packet_manifest_sha256':digest,'packet_files_excluding_manifest':len(files),'packet_bytes_excluding_manifest':packet_manifest['packet_bytes'],'record_id':'23244754','file_count':31,'total_bytes':453646943,'no_further_packet_edits_planned':True,'repository_mutations':0,'remote_writes':0,'additional_attachment_streams':0,'array_decodes':0,'source_callbacks':0}))
print(json.dumps({'status':packet_manifest['status'],'manifest_sha256':digest,'packet_files':len(files),'packet_bytes':packet_manifest['packet_bytes'],'record_id':'23244754'}))
