#!/usr/bin/env python3
"""Read-only serialized metadata review; stdlib only, no verifier imports."""
import hashlib
import html.parser
import json
from pathlib import Path
import re

ROOT=Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
PACKET=ROOT/'publication-supplement-priority/genuine-publication-packet-001'
OUT=ROOT/'publication-supplement-priority/genuine-publication-packet-001-metadata-review'
MANIFEST_SHA='d751185656bb6b11757770482c9c2df8bc58af70b190227ce4b765f8d324d88b'
METADATA_SHA='70d4eb5751e344323305a10187b89b6c1357744aa37b067aeee6da3aeb77b640'

def require(ok,reason):
    if not ok:raise RuntimeError(reason)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()

class Markup(html.parser.HTMLParser):
    VOID={'area','base','br','col','embed','hr','img','input','link','meta','param','source','track','wbr'}
    def __init__(self):
        super().__init__(convert_charrefs=True);self.items=[];self.depth=0
    def handle_starttag(self,tag,attrs):
        self.items.append(('start',tag,tuple(sorted(attrs))))
        if tag not in self.VOID:self.depth+=1
    def handle_endtag(self,tag):self.items.append(('end',tag));self.depth=max(0,self.depth-1)
    def parse_endtag(self,index):
        end=self.rawdata.find('>',index)
        if end>=0:require(re.fullmatch(r'</[A-Za-z][A-Za-z0-9:-]*\s*>',self.rawdata[index:end+1]) is not None,'malformed closing tag')
        return super().parse_endtag(index)
    def handle_startendtag(self,tag,attrs):self.items.append(('empty',tag,tuple(sorted(attrs))))
    def handle_data(self,data):
        if self.depth==0 and data and all(c in '\t\n\f\r ' for c in data):self.items.append(('edge_space',data))
        elif self.items and self.items[-1][0]=='text':self.items[-1]=('text',self.items[-1][1]+data)
        else:self.items.append(('text',data))
    def handle_comment(self,data):self.items.append(('comment',data))
    def handle_decl(self,data):self.items.append(('declaration',data))
    def handle_pi(self,data):self.items.append(('processing_instruction',data))
    def unknown_decl(self,data):self.items.append(('unknown_declaration',data))
    def events(self):
        result=[]
        for i,event in enumerate(self.items):
            if event[0]=='edge_space':
                if i==0 or i+1==len(self.items):continue
                event=('text',event[1])
            if result and result[-1][0]==event[0]=='text':result[-1]=('text',result[-1][1]+event[1])
            else:result.append(event)
        return result

VOCAB={('metadata','resource_type'),('metadata','rights'),('metadata','languages'),('metadata','creators','role'),('metadata','contributors','role'),('metadata','related_identifiers','relation_type'),('metadata','related_identifiers','resource_type'),('metadata','additional_descriptions','type'),('metadata','additional_descriptions','lang'),('custom_fields','code:programmingLanguage'),('custom_fields','code:developmentStatus')}
def fold(actual,expected,path=()):
    require(type(actual) is type(expected),'JSON type differs '+repr(path))
    if type(actual) is dict:
        words=tuple(x for x in path if not isinstance(x,int))
        drop={'title','description','icon','props'} if 'id' in actual and 'id' in expected and words in VOCAB else set()
        return {key:fold(value,expected.get(key),path+(key,)) for key,value in actual.items() if key not in drop}
    if type(actual) is list:return [fold(value,expected[i] if i<len(expected) else None,path+(i,)) for i,value in enumerate(actual)]
    if type(actual) is str and path and path[-1]=='description':
        parser=Markup();parser.feed(actual);parser.close();return parser.events()
    return actual
def same(actual,expected,path):return canonical(fold(actual,expected,path))==canonical(fold(expected,expected,path))
def editable(actual,expected):
    for key in ['metadata','custom_fields']:
        require(key in actual and key in expected and same(actual[key],expected[key],(key,)),key+' differs')
    require(type(actual['access']) is dict and set(actual['access'])<=set(expected['access'])|{'status'},'unexpected access fields')
    require(set(expected['access'])<=set(actual['access']),'access keys missing')
    require(same({key:actual['access'][key] for key in expected['access']},expected['access'],('access',)),'access differs')
def filemap(doc):
    raw=doc.get('files',doc)
    if type(raw) is dict:raw=raw.get('entries',raw)
    rows=list(raw.values()) if type(raw) is dict else raw
    require(type(rows) is list,'file list shape')
    result={}
    for row in rows:
        name=row.get('key',row.get('filename'));require(type(name) is str and name not in result,'duplicate/missing filename');result[name]=row
    return result
def fid(row):return row.get('file_id',row.get('id'))
def size(row):return row.get('size',row.get('filesize'))
def md5(row):return row.get('checksum','').removeprefix('md5:')

manifest_raw=(PACKET/'proof/WITNESS_MANIFEST.json').read_bytes();require(sha(manifest_raw)==MANIFEST_SHA,'manifest pin differs')
manifest=json.loads(manifest_raw)
for relative,pin in manifest['files'].items():
    raw=(PACKET/'proof'/relative).read_bytes();require(len(raw)==pin['bytes'] and sha(raw)==pin['sha256'],'witness pin differs '+relative)
require(manifest['historical_original_literal_Notes_guard']=='FAIL_PRESERVED','historical guard changed')
approved_raw=(PACKET/'proof/inventory/FINAL_METADATA_MODERN.json').read_bytes();require(sha(approved_raw)==METADATA_SHA,'approved70 metadata differs')
require(approved_raw==(ROOT/'publication-supplement-priority/FINAL_METADATA_MODERN.json').read_bytes(),'approved metadata origin bytes differ')
approved=json.loads(approved_raw)
public=json.loads((PACKET/'proof/publication/NEW_PUBLIC_RECORD.json').read_bytes())
latest=json.loads((PACKET/'proof/publication/CURRENT_LATEST.json').read_bytes())
editable(public,approved);editable(latest,approved)
prior=json.loads((PACKET/'proof/preserved/23225288-baseline-record.json').read_bytes())
protected=sorted(set(prior['metadata'])-{'title','description','version','additional_descriptions','related_identifiers'})
for key in protected:require(same(approved['metadata'][key],prior['metadata'][key],('metadata',key)),'prior protected field differs '+key)
require(same(approved['custom_fields'],prior['custom_fields'],('custom_fields',)),'prior custom_fields changed')
require(same(approved['access'],{key:prior['access'][key] for key in approved['access']},('access',)),'prior editable access changed')
old_rel=prior['metadata']['related_identifiers'];new_rel=approved['metadata']['related_identifiers']
require(same(new_rel[:len(old_rel)],old_rel,('metadata','related_identifiers')),'prior related identifiers changed/reordered')
old_notes=prior['metadata']['additional_descriptions'][0]['description'];new_notes=public['metadata']['additional_descriptions'][0]['description']
require(old_notes in new_notes,'prior Notes0 literal history missing')
require(same(public['metadata']['additional_descriptions'][1:],prior['metadata']['additional_descriptions'][1:],('metadata','additional_descriptions')),'historical Other description changed')
origins=json.loads((PACKET/'verification/independent-review/DATED_BASELINE_PROVENANCE_final.json').read_bytes())
origin_rows=[]
for row in origins['rows']:
    raw=Path(row['original_path']).read_bytes();require(len(raw)==row['bytes'] and sha(raw)==row['sha256'],'dated origin differs')
    name=manifest['preserved_records'][row['record_id']][row['role']]
    require(raw==(PACKET/'proof'/name).read_bytes(),'packet dated baseline differs from origin')
    origin_rows.append({'record_id':row['record_id'],'role':row['role'],'original_path':row['original_path'],'bytes':len(raw),'sha256':sha(raw)})
preserved=[]
for record_id,group in manifest['preserved_records'].items():
    objs={role:json.loads((PACKET/'proof'/name).read_bytes()) for role,name in group.items()}
    baseline,current=objs['baseline'],objs['record']
    if all(key in baseline for key in ['metadata','custom_fields','access']):editable(current,{key:baseline[key] for key in ['metadata','custom_fields','access']})
    else:require(same(current['metadata'],baseline['metadata'],('metadata',)),'legacy preserved metadata differs')
    old,embedded,dated,fresh=map(filemap,[baseline,current,objs['baseline_files'],objs['files']])
    require(set(old)==set(embedded)==set(dated)==set(fresh),'preserved file membership differs')
    for name,entry in old.items():
        for row in [embedded[name],dated[name],fresh[name]]:require(type(size(row)) is int and size(row)==size(entry) and md5(row)==md5(entry),'preserved size/md5 differs')
        require(fid(fresh[name])==fid(dated[name]) and fresh[name]['version_id']==dated[name]['version_id'],'dated immutable pair differs')
        require(fid(embedded[name])==fid(entry)==fid(fresh[name]),'embedded file ID differs')
        if 'version_id' in entry:require(embedded[name]['version_id']==entry['version_id'],'baseline embedded version differs')
        if 'version_id' in embedded[name]:require(embedded[name]['version_id']==fresh[name]['version_id'],'current embedded version differs')
    preserved.append({'record_id':record_id,'file_count':len(old),'editable_metadata_preserved':True,'file_membership_size_md5_and_immutable_ids_preserved':True,'record_baseline_sha256':sha((PACKET/'proof'/group['baseline']).read_bytes()),'record_current_sha256':sha((PACKET/'proof'/group['record']).read_bytes()),'dated_files_sha256':sha((PACKET/'proof'/group['baseline_files']).read_bytes()),'current_files_sha256':sha((PACKET/'proof'/group['files']).read_bytes())})
history=[]
original_verification=ROOT/'publication-final/verification'
for path in sorted((PACKET/'verification').rglob('*')):
    if not path.is_file():continue
    relative=path.relative_to(PACKET/'verification');origin=original_verification/relative
    require(origin.is_file() and path.read_bytes()==origin.read_bytes(),'verification history copy differs '+str(relative))
    history.append({'path':str(relative),'sha256':sha(path.read_bytes())})
result={'status':'PASS_GENUINE_SERIALIZED_METADATA_AND_PRESERVATION_REVIEW','scope':'Metadata and four preserved records/files only; no content streaming or scientific replay','witness_manifest_sha256':MANIFEST_SHA,'approved_metadata_sha256':METADATA_SHA,'public_record_sha256':sha((PACKET/'proof/publication/NEW_PUBLIC_RECORD.json').read_bytes()),'accepted_verifier_source_sha256':sha((PACKET/'verification/verify_publication.py').read_bytes()),'public_and_latest_match_approved_metadata_under_strict_semantics':True,'protected_prior_metadata_fields':protected,'custom_fields_access_preserved':True,'prior_related_identifier_count':len(old_rel),'additional_related_identifiers':new_rel[len(old_rel):],'old_Notes0_literal':{'chars':len(old_notes),'sha256':sha(old_notes.encode()),'contained_in_public_new_Notes0':True},'historical_original_literal_Notes_guard':'FAIL_PRESERVED','historical_Other_exact_text_sha256':sha(prior['metadata']['additional_descriptions'][1]['description'].encode()),'preserved_records':preserved,'dated_origin_rows':origin_rows,'full_verification_history_unchanged_file_count':len(history),'full_verification_history_unchanged_files':history,'network_calls':0,'controller_or_numerical_imports':0,'source_asset_repository_packet_edits':0,'blockers':[]}
(OUT/'METADATA_PRESERVATION_ACCEPTANCE.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':result['status'],'preserved_record_count':len(preserved),'dated_origins':len(origin_rows),'history_files':len(history),'blockers':result['blockers']}))
