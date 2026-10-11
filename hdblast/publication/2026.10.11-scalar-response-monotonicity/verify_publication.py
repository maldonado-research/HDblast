#!/usr/bin/env python3
"""Verify pinned saved public witnesses offline; never contacts a server."""
import argparse,datetime,hashlib,json,re
from pathlib import Path

def require(test,message):
    if not test:raise ValueError(message)
def sha(data):return hashlib.sha256(data).hexdigest()
def read(path):return json.loads(path.read_text())
def exact(a,b):return json.dumps(a,sort_keys=True,separators=(',',':'),allow_nan=False)==json.dumps(b,sort_keys=True,separators=(',',':'),allow_nan=False)
def filemap(obj):
    rows=obj['entries'] if 'entries' in obj else obj['files']['entries']
    if isinstance(rows,dict):
        require(all(type(key) is str and type(value) is dict and value.get('key')==key for key,value in rows.items()),'FILE_MAPPING_KEY_DIFFERS')
        rows=list(rows.values())
    require(isinstance(rows,list),'INVALID_FILE_LIST')
    names=[r['key'] for r in rows];require(all(type(n) is str for n in names) and len(names)==len(set(names)),'DUPLICATE_FILE_NAMES')
    return dict(zip(names,rows))
def vocabulary(obj):
    require(type(obj) is dict and set(obj)<={'id','title','description','icon','props'} and type(obj.get('id')) is str,'UNRECOGNIZED_VOCABULARY_SHAPE')
    return {'id':obj['id']}
def normal_metadata(metadata):
    value=json.loads(json.dumps(metadata))
    for item in value['related_identifiers']:
        for key in ('relation_type','resource_type'):
            if key in item:item[key]=vocabulary(item[key])
    return value
def compare_metadata(record,expected):
    require(exact(normal_metadata(record['metadata']),normal_metadata(expected['metadata'])),'COMPLETE_METADATA_DIFFERS')
    require(exact(record['custom_fields'],expected['custom_fields']),'COMPLETE_CUSTOM_FIELDS_DIFFERS')
    access=dict(record['access']);status=access.pop('status',None)
    require(status in (None,'open') and exact(access,expected['access']),'ACCESS_DIFFERS')
def main():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--manifest-sha256',required=True);p.add_argument('--output');a=p.parse_args()
    root=Path(a.root).resolve();raw=(root/'WITNESS_MANIFEST.json').read_bytes()
    require(re.fullmatch('[a-f0-9]{64}',a.manifest_sha256) and sha(raw)==a.manifest_sha256,'WITNESS_MANIFEST_PIN_DIFFERS')
    manifest=json.loads(raw);expectedpaths=set(manifest['files'])
    present={str(x.relative_to(root)) for x in root.rglob('*') if x.is_file()}
    require(present==expectedpaths|{'WITNESS_MANIFEST.json'},'WITNESS_FILE_MEMBERSHIP_DIFFERS')
    for name,pin in manifest['files'].items():
        path=root/name
        require(not path.is_symlink() and path.resolve().is_relative_to(root),'UNSAFE_WITNESS_PATH')
        data=path.read_bytes();require(type(pin['bytes']) is int and len(data)==pin['bytes'] and sha(data)==pin['sha256'],'WITNESS_BYTES_DIFFER_'+name)
    record=read(root/'NEW_PUBLIC_RECORD.json');listing=read(root/'NEW_PUBLIC_FILES.json');inventory=read(root/'PUBLIC_INVENTORY.json');verified=read(root/'VERIFIED_PUBLIC.json');streams=read(root/'PUBLIC_STREAM_OBSERVATIONS.json');writes=read(root/'WRITE_SUMMARY.json')
    rid=inventory['record_id'];require(type(rid) is str and rid.isdigit() and rid!='23244754','NEW_RECORD_ID_INVALID')
    require(record['id']==rid and record['parent']['id']=='17088132' and record['is_published'] is True and record['is_draft'] is False,'RECORD_IDENTITY_STATE_DIFFERS')
    require(record['pids']['doi']['identifier']=='10.5281/zenodo.'+rid and record['parent']['pids']['doi']['identifier']=='10.5281/zenodo.17088132','DOI_FAMILY_DIFFERS')
    require(record['versions']['is_latest'] is True and str(record['parent']['access']['owned_by']['user'])=='1386319','LATEST_OR_OWNER_DIFFERS')
    compare_metadata(record,read(root/'FINAL_METADATA_MODERN.json'))
    preceding=streams['preceding_public_file_list_event']
    require(streams['record_id']==rid and preceding['kind']=='GET' and preceding['url']=='https://zenodo.org/api/records/'+rid+'/files' and type(preceding['status']) is int and preceding['status']==200 and preceding['authenticated_request'] is False,'PUBLIC_STREAM_ROUND_SCOPE_DIFFERS')
    require(streams['source_journal_sha256']==writes['source_journal_sha256'] and re.fullmatch('[0-9a-f]{64}',streams['source_journal_sha256']),'PUBLIC_JOURNAL_BINDING_DIFFERS')
    index=preceding['source_journal_index'];require(type(index) is int and index>=0,'PUBLIC_LIST_EVENT_INDEX_INVALID')
    previous_time=datetime.datetime.fromisoformat(preceding['utc'])
    require(previous_time.tzinfo is not None,'PUBLIC_EVENT_TIMEZONE_REQUIRED')
    for offset,event in enumerate(streams['files'],1):
        require(event['kind']=='FULL_CONTENT_STREAM_VERIFIED' and type(event['source_journal_index']) is int and event['source_journal_index']==index+offset,'PUBLIC_STREAM_EVENT_ORDER_DIFFERS')
        observed_time=datetime.datetime.fromisoformat(event['utc']);require(observed_time.tzinfo is not None and observed_time>=previous_time,'PUBLIC_STREAM_TIME_ORDER_DIFFERS');previous_time=observed_time
    files=filemap(listing);pins={r['filename']:r for r in inventory['files']};observed={r['filename']:r for r in streams['files']};vrows={r['filename']:r for r in verified['content']}
    require(len(inventory['files'])==len(pins)==len(streams['files'])==len(observed)==len(verified['content'])==len(vrows)==33,'FILE_COUNT_OR_DUPLICATE_DIFFERS')
    require(set(files)==set(pins)==set(observed)==set(vrows),'FILE_MEMBERSHIP_DIFFERS')
    for name,pin in pins.items():
        e=files[name];obs=observed[name];v=vrows[name]
        require(type(pin['bytes']) is int and type(e['size']) is int and e['size']==pin['bytes'] and e['checksum']=='md5:'+pin['md5'] and e['status']=='completed','FILE_PIN_DIFFERS_'+name)
        require(exact(obs['observed'],{k:pin[k] for k in ('bytes','md5','sha256')}) and type(obs['utc']) is str and bool(obs['utc']),'COMPLETE_STREAM_OBSERVATION_DIFFERS_'+name)
        require(v['sha256']==pin['sha256'] and v['basis']=='FRESH_COMPLETE_CONTENT_STREAM' and e['file_id']==v['remote_file_id'] and e['version_id']==v['remote_version_id'],'VERIFIED_IDENTITY_DIFFERS_'+name)
    require(inventory['total_bytes']==sum(x['bytes'] for x in pins.values())==454129865 and inventory['file_count']==33,'TOTAL_BYTES_DIFFERS')
    require(verified['published'] is True and verified['draft_id']==rid and verified['parent_id']=='17088132' and verified['file_count']==33 and verified['total_bytes']==454129865,'VERIFICATION_RECEIPT_IDENTITY_DIFFERS')
    require(verified['metadata_sha256']==sha((root/'FINAL_METADATA_MODERN.json').read_bytes()),'APPROVED_METADATA_PIN_DIFFERS')
    require(verified['inventory_sha256']==inventory['derived_from_frozen_private_inventory_sha256'],'APPROVED_INVENTORY_PROVENANCE_DIFFERS')
    require(writes['publish_intent_count']==1 and writes['publish_response_count']<=1 and writes['draft_id']==rid,'PUBLICATION_INTENT_COUNTS_DIFFER')
    preserved=read(root/'PRESERVED_RECORDS.json')
    expected_groups={'23244754':'17088132','23228395':'17088132','23225288':'17088132','23114217':'17088132','22347452':'17088132','23111008':'22922927'}
    require(len(preserved['records'])==6 and {g['record_id']:g['concept_id'] for g in preserved['records']}==expected_groups,'EXACT_SIX_PRESERVED_RECORDS_REQUIRED')
    for group in preserved['records']:
        for key in ('before_record','after_record','before_files','after_files'):
            relative=group[key];require(type(relative) is str and relative.startswith('preserved/') and not Path(relative).is_absolute() and '..' not in Path(relative).parts and relative in expectedpaths,'UNSAFE_PRESERVATION_WITNESS_PATH')
        before=read(root/group['before_record']);after=read(root/group['after_record']);bf=read(root/group['before_files']);af=read(root/group['after_files'])
        require(all(before[k]==after[k] for k in ('metadata','custom_fields','access')) and filemap(bf)==filemap(af),'PRESERVED_RECORD_CHANGED_'+group['record_id'])
        require(all(str(r['id'])==group['record_id'] and str(r['parent']['id'])==group['concept_id'] and r['is_published'] is True for r in (before,after)),'PRESERVED_IDENTITY_CHANGED')
    result={'status':'PASS_OFFLINE_SAVED_33_FILE_PUBLICATION_WITNESSES','record_id':rid,'concept_id':'17088132','file_count':33,'total_bytes':454129865,'witness_files':len(expectedpaths),'preserved_records':len(preserved['records']),'manifest_sha256':a.manifest_sha256,'network_calls':0,'scope':'Authenticates these pinned saved client observations; does not independently re-download attachment bytes, prove mathematics, or establish current server state.'}
    if a.output:
        output=Path(a.output).resolve();require(not output.is_relative_to(root),'OUTPUT_MUST_BE_OUTSIDE_WITNESS_ROOT');output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
