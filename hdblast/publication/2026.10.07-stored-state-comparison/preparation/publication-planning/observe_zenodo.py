#!/usr/bin/env python3
"""GET-only baseline observation. Never sends tokens via URLs or logs."""
from pathlib import Path
from urllib import request,error,parse
from datetime import datetime,timezone
import hashlib,importlib.util,json,os,ssl
OLD=Path('/workspace/hdblast-research-work/continuation-next-20261007/release-planning')
OUT=Path(__file__).resolve().parent/'read-only-observation'
SRC=OLD/'new_edition_controller.py'
PIN='7dade987e77ff166e5849db06cd47d019f6ab16d67bd18a22c2120564c480110'
assert hashlib.sha256(SRC.read_bytes()).hexdigest()==PIN
spec=importlib.util.spec_from_file_location('historical_controller',SRC);m=importlib.util.module_from_spec(spec)
import sys
sys.modules[spec.name]=m;spec.loader.exec_module(m)
token=os.environ.get('ZENODO_ACCESS_TOKEN','')
t=m.Transport(token)
receipts=[]
def get(path,label,public=False,accept=m.VENDOR):
 url=m.ORIGIN+path
 r=t.request('GET',url,headers={'Accept':accept,'_public':public})
 if r.status in (301,302,303,307,308):
  target=next((v for k,v in r.headers.items() if k.lower()=='location'),None)
  m.safe_url(target)
  if not target.startswith(m.ORIGIN+'/api/records/') or not target.rsplit('/',1)[-1].isdigit():raise m.Stop('LATEST_REDIRECT_INVALID')
  receipts.append({'label':label+'_redirect','status':r.status,'response_sha256':m.sha(r.body)})
  r=t.request('GET',target,headers={'Accept':accept,'_public':public})
 if r.status!=200:raise m.Stop('GET_'+label+'_HTTP_'+str(r.status))
 if token.encode() in r.body:raise m.Stop('TOKEN_IN_RESPONSE_REFUSE_SAVE')
 obj=r.json()
 out=OUT/(label+'.raw.json');out.write_bytes(r.body);out.chmod(0o600)
 receipts.append({'label':label,'status':r.status,'bytes':len(r.body),'sha256':m.sha(r.body),'authenticated':not public,'etag':next((v for k,v in r.headers.items() if k.lower()=='etag'),None)})
 return obj
p=get('/api/records/23225288','CURRENT_PRIOR')
latest=get('/api/records/23225288/versions/latest','CURRENT_LATEST',public=True)
files=get('/api/records/23225288/files','CURRENT_FILES',public=True,accept='application/json')
legacy=get('/api/deposit/depositions/23225288','OWNED_CURRENT_PRIOR',accept='application/json')
owner_pages=[];family=[];seen=set()
for page in range(1,21):
 body=get('/api/deposit/depositions?'+parse.urlencode({'page':page,'size':100,'sort':'mostrecent'}),'OWNER_PAGE_'+str(page),accept='application/json')
 items=body if isinstance(body,list) else body.get('hits',{}).get('hits')
 m.require(isinstance(items,list) and len(items)<=100,'OWNER_PAGE_SHAPE')
 ids=[str(x.get('id')) for x in items];m.require(all(i.isdigit() for i in ids) and len(set(ids))==len(ids) and not set(ids)&seen,'OWNER_PAGE_IDENTITY');seen.update(ids)
 family += [x for x in items if str(x.get('conceptrecid',x.get('parent',{}).get('id','')))=='17088132']
 owner_pages.append({'page':page,'count':len(items)})
 if len(items)<100:break
else:raise m.Stop('OWNER_LIST_NOT_EXHAUSTED')
fmap=m.file_map(files)
previous_file_bytes=(OLD/'new-edition-execution/NEW_PUBLIC_FILES.json').read_bytes()
previous_files=m.file_map(json.loads(previous_file_bytes))
m.require(set(fmap)==set(previous_files),'PRIOR_FILE_MEMBERSHIP_CHANGED')
for name,f in fmap.items():
 previous=previous_files[name]
 m.require(f.get('file_id',f.get('id'))==previous.get('file_id',previous.get('id')) and f.get('version_id')==previous.get('version_id'),'CURRENT_IMMUTABLE_FILE_IDENTITY_CHANGED')
i=json.loads((OLD/'FINAL_UPLOAD_INVENTORY.json').read_text());pins=i['inherited_files']+i['new_files']
m.assert_pins(fmap,pins)
m.require(len(fmap)==27 and sum(f.get('size',f.get('filesize')) for f in fmap.values())==448387915,'CURRENT_COUNT_BYTES_CHANGED')
m.require(p.get('id')=='23225288' and p.get('parent',{}).get('id')=='17088132' and p.get('is_published') is True,'CURRENT_IDENTITY_CHANGED')
m.metadata_equal(p,json.loads((OLD/'FINAL_METADATA_MODERN.json').read_text()))
m.require(latest.get('id')=='23225288' and latest.get('is_published') is True,'LATEST_CHANGED_FROM_EXPECTED')
owner=legacy.get('owner');owner=owner.get('id') if isinstance(owner,dict) else owner
m.require(str(owner)=='1386319' and str(legacy.get('conceptrecid'))=='17088132' and legacy.get('submitted') is True,'CURRENT_OWNER_IDENTITY_CHANGED')
rows=[]
for pin in pins:
 f=fmap[pin['filename']]
 rows.append({'filename':pin['filename'],'status':'PASS_PRIOR_COMPLETE_BYTES_WITH_FRESH_IMMUTABLE_IDENTITY','observed':{k:pin[k] for k in ['bytes','md5','sha256']},'remote_file_id':f.get('file_id',f.get('id')),'remote_version_id':f.get('version_id'),'sha256_basis':'PINNED_PRIOR_COMPLETE_PUBLIC_27_FILE_STREAM_RECEIPT_NOT_REPEATED_IN_THIS_GET_ONLY_OBSERVATION'})
 m.require(rows[-1]['remote_file_id'] and rows[-1]['remote_version_id'],'CURRENT_FILE_IDS_MISSING')
prior_receipt=json.loads((OLD/'new-edition-execution/VERIFIED_PUBLIC.json').read_text())
m.require(prior_receipt['published'] is True and prior_receipt['draft_id']=='23225288' and prior_receipt['file_count']==27 and prior_receipt['total_bytes']==448387915 and all(x['basis']=='FRESH_COMPLETE_CONTENT_STREAM' for x in prior_receipt['content']),'PRIOR_COMPLETE_RECEIPT_DIFFERS')
m.require({x['filename']:x['sha256'] for x in prior_receipt['content']}=={x['filename']:x['sha256'] for x in pins},'PRIOR_COMPLETE_PIN_DIFFERS')
owned=[x for x in family if str(x.get('id'))=='23225288'];m.require(len(owned)==1 and owned[0].get('submitted') is True,'OWNER_CURRENT_NOT_LISTED')
for x in family:
 owner=x.get('owner');owner=owner.get('id') if isinstance(owner,dict) else owner
 m.require(str(owner)=='1386319','FAMILY_OWNER_DIFFERS')
drafts=[str(x.get('id')) for x in family if x.get('submitted') is False or x.get('is_draft') is True]
summary={'utc':datetime.now(timezone.utc).isoformat(),'status':'PASS_GET_ONLY_CURRENT_27_FILE_BASELINE','current_record':'23225288','current_latest_record':latest['id'],'concept_id':'17088132','owner':1386319,'file_count':27,'total_bytes':448387915,'metadata_equal_to_completed_release':True,'same_family_owned_unpublished_draft_ids':drafts,'owner_pages':owner_pages,'no_remote_mutations':True,'sha256_content_basis':'Prior completed root publication full streams; this observation repeats metadata/listing/identity GETs only','source_controller_sha256':PIN,'source_inventory_sha256':m.sha((OLD/'FINAL_UPLOAD_INVENTORY.json').read_bytes()),'source_metadata_sha256':m.sha((OLD/'FINAL_METADATA_MODERN.json').read_bytes()),'source_complete_public_files_sha256':m.sha(previous_file_bytes),'source_complete_public_receipt_sha256':m.sha((OLD/'new-edition-execution/VERIFIED_PUBLIC.json').read_bytes()),'requests':receipts}
OUT.mkdir(parents=True,exist_ok=True)
m.atomic_json(OUT/'READ_ONLY_OBSERVATION.json',summary)
m.atomic_json(OUT/'PRIOR27_CONTENT_RECEIPT.json',{'status':'PRIOR_FULL_STREAM_ATTESTATION_WITH_FRESH_GET_IDENTITY','rows':rows,'basis':summary})
baseline={'current_record':23225288,'concept_id':17088132,'concept_doi':'10.5281/zenodo.17088132','current_version_doi':'10.5281/zenodo.23225288','current_semantic_version':p['metadata']['version'],'expected_files':27,'expected_bytes':448387915,'files':[{k:x[k] for k in ['filename','bytes','md5','sha256']} for x in pins]}
m.atomic_json(OUT/'CURRENT_PUBLISHED_INVENTORY.json',baseline)
print(json.dumps({k:summary[k] for k in ['status','utc','current_record','current_latest_record','concept_id','owner','file_count','total_bytes','same_family_owned_unpublished_draft_ids','no_remote_mutations']}))
