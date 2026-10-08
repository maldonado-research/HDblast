#!/usr/bin/env python3
"""Staged controller for one additive HDBLAST version; default stage is GET-only.

Prospective 27 inherited + 2 added edition. No mutation has been run.
The completed 23225288 publisher and its state remain untouched.
All writable inputs must be frozen by
explicit SHA256 arguments. Token values are read only from the environment.
"""
from __future__ import annotations
import argparse
import fcntl
import hashlib
import html.parser
import json
import os
from pathlib import Path
import re
import ssl
import sys
import time
from datetime import datetime, timezone
from urllib import error, parse, request
from dataclasses import dataclass

HERE = Path(__file__).resolve().parent
ORIGIN = "https://zenodo.org"
VENDOR = "application/vnd.inveniordm.v1+json"
PRIOR_ID = "23225288"
INHERITED_COUNT = 27
ADDITION_COUNT = 2
FINAL_COUNT = INHERITED_COUNT + ADDITION_COUNT
PRIOR_SEMANTIC_VERSION = "2026.10.07-bd-prehistory-target"
SCIENTIFIC_REPLAY_STATUS = "PASS_COMPLETE_STORED_STATE_COMPARISON_REPLAY"
PARENT_ID = "17088132"
EXPECTED_OWNER = 1386319
PRIOR_RECEIPT = HERE/'read-only-observation/PRIOR27_CONTENT_RECEIPT.json'
PRIOR_RECEIPT_SHA = '2664ec6a867d806e6146847fcf677d2e5c20061a400262acb8563d3ca552529c'
MAX_JSON = 8 * 1024 * 1024
MAX_ASSET = 256 * 1024 * 1024
BASELINE_SHA = '94ec06030315e99c6465c784c3dcae741ed028dd89e276335b656482e11c63b4'
PUBLISHED_BASELINE_SHA = '261d81302396c1683f5c0e1185b071189cf83bbb5fde5fb4a098ca64e8fe47ca'

class Stop(RuntimeError):
    pass

class EvidenceLock:
    """Hold a nonblocking exclusive process lock for one entire CLI stage."""
    def __init__(self,evidence):
        self.directory=Path(evidence);self.fd=None
    def __enter__(self):
        self.directory.mkdir(mode=0o700,parents=True,exist_ok=True)
        self.fd=os.open(self.directory/'CONTROLLER.lock',os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW,0o600)
        try:
            fcntl.flock(self.fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            os.close(self.fd);self.fd=None
            raise Stop('EVIDENCE_DIRECTORY_ALREADY_LOCKED') from None
        return self
    def __exit__(self,*args):
        if self.fd is not None:
            fcntl.flock(self.fd,fcntl.LOCK_UN);os.close(self.fd);self.fd=None

def require(condition, code):
    if not condition:
        raise Stop(code)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def canonical_bytes(value):
    return json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf-8')

def revision_if_match(raw_etag,record):
    require(isinstance(raw_etag,str),'ETAG_NOT_CANONICAL_QUOTED_DECIMAL')
    match=re.fullmatch(r'"(0|[1-9][0-9]*)"',raw_etag)
    require(match is not None,'ETAG_NOT_CANONICAL_QUOTED_DECIMAL')
    revision=record.get('revision_id')
    require(type(revision) is int and revision>=0 and str(revision)==match.group(1),
            'ETAG_REVISION_ID_DIFFERS')
    return match.group(1)

def bounded_error_json(body):
    def reject_constant(value):
        raise ValueError('NONSTANDARD_JSON_CONSTANT')
    value=json.loads(body,parse_constant=reject_constant)
    stack=[(value,0)];count=0
    while stack:
        item,depth=stack.pop();count+=1
        if depth>64 or count>100000:
            raise ValueError('ERROR_JSON_STRUCTURE_LIMIT')
        if isinstance(item,dict):
            stack.extend((x,depth+1) for x in item.values())
        elif isinstance(item,list):
            stack.extend((x,depth+1) for x in item)
    canonical_bytes(value)
    return value

def utc():
    return datetime.now(timezone.utc).isoformat()

def safe_url(url, query=False):
    require(isinstance(url, str) and not any(ord(c) < 32 for c in url), 'UNSAFE_URL')
    p = parse.urlsplit(url)
    require(p.scheme == 'https' and p.hostname == 'zenodo.org' and p.port is None
            and p.username is None and p.password is None and not p.fragment
            and '\\' not in url and p.path.startswith('/api/'), 'UNSAFE_DESTINATION')
    require(not p.query or query, 'UNEXPECTED_QUERY')
    if p.query:
        pairs = parse.parse_qsl(p.query, keep_blank_values=True)
        require(all(k in ('page', 'size', 'sort', 'q', 'expand') for k, v in pairs), 'UNSAFE_QUERY')
    return url

def redact(obj, token=''):
    if isinstance(obj, dict):
        return {k: ('[REDACTED]' if k.lower() in ('authorization', 'access_token', 'token', 'secret')
                    else redact(v, token)) for k, v in obj.items()}
    if isinstance(obj, list):
        return [redact(x, token) for x in obj]
    if isinstance(obj, str):
        value = obj.replace(token, '[REDACTED]') if token else obj
        def clean_url(match):
            url=match.group(0)
            try:
                p = parse.urlsplit(url)
            except ValueError:
                return '[REDACTED_MALFORMED_URL]'
            if p.username or p.password:
                return '[REDACTED_CREDENTIAL_URL]'
            if p.query:
                return parse.urlunsplit((p.scheme, p.netloc, p.path, '[REDACTED_QUERY]', ''))
            return url
        return re.sub(r"https?://[^\s<>\"']+",clean_url,value)
    return obj

def atomic_json(path, obj, token=''):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    with tmp.open('w', encoding='utf-8') as f:
        json.dump(redact(obj, token), f, ensure_ascii=False, indent=2)
        f.write('\n'); f.flush(); os.fsync(f.fileno())
    tmp.chmod(0o600); os.replace(tmp, path)

def frozen_json(path, pin):
    data = Path(path).read_bytes()
    require(re.fullmatch('[0-9a-f]{64}', pin or '') and sha(data) == pin, 'FROZEN_INPUT_SHA256_DIFFERS')
    return json.loads(data)

class HtmlEvents(html.parser.HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True); self.events=[]; self.depth=0
    def handle_starttag(self, tag, attrs):
        self.events.append(('start', tag, tuple(sorted(attrs)))); self.depth += 1
    def handle_endtag(self, tag):
        self.events.append(('end', tag)); self.depth = max(0, self.depth-1)
    def handle_startendtag(self, tag, attrs):
        self.events.append(('empty', tag, tuple(sorted(attrs))))
    def handle_data(self, data):
        if self.depth == 0 and not data.strip():
            return
        if self.events and self.events[-1][0] == 'text':
            self.events[-1] = ('text', self.events[-1][1] + data)
        else:
            self.events.append(('text', data))
    def handle_comment(self, data):
        self.events.append(('comment', data))

def vocabulary_path(path):
    words=tuple(x for x in path if not isinstance(x,int))
    return words in (
        ('metadata','resource_type'),('metadata','rights'),('metadata','languages'),
        ('metadata','creators','role'),('metadata','contributors','role'),
        ('metadata','related_identifiers','relation_type'),
        ('metadata','related_identifiers','resource_type'),
        ('metadata','additional_descriptions','type'),
        ('metadata','additional_descriptions','lang'),
        ('custom_fields','code:programmingLanguage'),
        ('custom_fields','code:developmentStatus'),
    )

def normalized(obj, expected=None, path=()):
    if isinstance(obj, dict):
        drop = {'title','description','icon','props'} if (isinstance(expected,dict)
                and 'id' in expected and 'id' in obj and vocabulary_path(path)) else set()
        return {k: normalized(v,expected.get(k) if isinstance(expected,dict) else None,path+(k,))
                for k,v in obj.items() if k not in drop}
    if isinstance(obj, list):
        return [normalized(v,expected[i] if isinstance(expected,list) and i<len(expected) else None,path+(i,))
                for i,v in enumerate(obj)]
    if isinstance(obj,str) and isinstance(expected,str) and path and path[-1]=='description':
        p=HtmlEvents(); p.feed(obj); p.close(); return p.events
    return obj

def value_equal(actual,expected,path):
    return canonical_bytes(normalized(actual,expected,path))==canonical_bytes(normalized(expected,expected,path))

def metadata_equal(actual, expected):
    require(set(expected) == {'metadata', 'custom_fields', 'access'}, 'METADATA_BODY_FIELDS_DIFFER')
    for key in ('metadata', 'custom_fields'):
        require(value_equal(actual.get(key,{}),expected[key],(key,)), 'COMPLETE_' + key.upper() + '_READBACK_DIFFERS')
    access=actual.get('access', {})
    require(isinstance(access,dict) and set(access)<=set(expected['access'])|{'status'},'ACCESS_UNEXPECTED_FIELDS')
    projected={k: access.get(k) for k in expected['access']}
    require(value_equal(projected,expected['access'],('access',)), 'ACCESS_READBACK_DIFFERS')
    return True

def file_map(document):
    if isinstance(document.get('files'), dict):
        entries=document['files'].get('entries', {})
    else:
        entries=document.get('entries', document.get('files', []))
    if isinstance(entries, dict):
        entries=list(entries.values())
    require(isinstance(entries, list), 'INVALID_FILE_LIST')
    answer={}
    for f in entries:
        name=f.get('key', f.get('filename'))
        require(isinstance(name, str) and name not in answer, 'DUPLICATE_OR_MISSING_REMOTE_FILENAME')
        answer[name]=f
    return answer

def immutable_file_identity(entry):
    file_id=entry.get('file_id',entry.get('id'));version_id=entry.get('version_id')
    require(isinstance(file_id,str) and bool(file_id) and isinstance(version_id,str) and bool(version_id),
            'COMPLETE_REMOTE_FILE_AND_VERSION_IDS_REQUIRED')
    return file_id,version_id

def assert_pins(files, pins, exact=True, completed=True):
    expected={p['filename'] for p in pins}
    require(set(files) == expected if exact else expected <= set(files), 'FILE_MEMBERSHIP_DIFFERS')
    for p in pins:
        f=files[p['filename']]
        size=f.get('size',f.get('filesize'))
        require(type(size) is int and type(p['bytes']) is int and size==p['bytes'], 'FILE_SIZE_DIFFERS')
        checksum=f.get('checksum', '')
        require(checksum in (p['md5'], 'md5:'+p['md5']), 'FILE_MD5_DIFFERS')
        if completed and 'status' in f:
            require(f['status'] == 'completed', 'FILE_NOT_COMPLETED')
    return True

@dataclass
class Response:
    status: int
    headers: dict
    body: bytes
    body_truncated: bool = False
    body_read_failed: bool = False
    def json(self):
        require(len(self.body) <= MAX_JSON, 'JSON_RESPONSE_TOO_LARGE')
        return json.loads(self.body)

class NoRedirect(request.HTTPRedirectHandler):
    def redirect_request(self, *args, **kwargs):
        return None

class Transport:
    def __init__(self, token):
        require(token and '\n' not in token and '\r' not in token, 'TOKEN_BINDING_REQUIRED')
        self.token=token
        self.opener=request.build_opener(NoRedirect(), request.HTTPSHandler(context=ssl.create_default_context()))
    def request(self, method, url, data=None, headers=None):
        safe_url(url, query=method=='GET')
        supplied=headers or {}
        h={'User-Agent':'HDBLAST-additive-edition-controller/20261007'}
        if not supplied.get('_public'):
            h['Authorization']='Bearer '+self.token
        h.update({k:v for k,v in supplied.items() if k!='_public'})
        if data is not None:
            h['Content-Length']=str(len(data))
        req=request.Request(url, data=data, method=method, headers=h)
        try:
            with self.opener.open(req, timeout=60) as r:
                require(r.geturl()==url, 'FINAL_URL_DIFFERS')
                body=r.read(MAX_JSON+1)
                require(len(body)<=MAX_JSON, 'JSON_RESPONSE_TOO_LARGE')
                return Response(r.status, dict(r.headers.items()), body)
        except error.HTTPError as e:
            try:
                body=e.read(MAX_JSON+1)
                return Response(e.code,dict(e.headers.items()),body,body_truncated=len(body)>MAX_JSON)
            except Exception:
                return Response(e.code,dict(e.headers.items()),b'',body_read_failed=True)
        except Stop:
            raise
        except Exception:
            raise Stop('TRANSPORT_OUTCOME_UNKNOWN') from None
    def stream(self, url, headers=None):
        safe_url(url)
        supplied=headers or {}
        h={'Accept-Encoding':'identity'}
        if not supplied.get('_public'):
            h['Authorization']='Bearer '+self.token
        h.update({k:v for k,v in supplied.items() if k!='_public'})
        try:
            with self.opener.open(request.Request(url, method='GET', headers=h), timeout=60) as r:
                require(r.status==200 and r.geturl()==url, 'CONTENT_GET_NOT_FULL_200')
                require(r.headers.get('Content-Range') is None and r.headers.get('Content-Encoding','identity')=='identity', 'PARTIAL_OR_ENCODED_CONTENT')
                while chunk:=r.read(1024*1024):
                    yield chunk
        except Stop:
            raise
        except Exception:
            raise Stop('CONTENT_STREAM_FAILED') from None

class Controller:
    def __init__(self, evidence, transport, token='', inventory=None, metadata=None,
                 inventory_sha=None, metadata_sha=None, baseline=None, prior_receipt=None,published_baseline=None):
        self.dir=Path(evidence); self.dir.mkdir(mode=0o700, parents=True, exist_ok=True)
        self.transport=transport; self.token=token
        self.baseline=baseline or frozen_json(HERE/'read-only-observation/CURRENT_PUBLISHED_INVENTORY.json', BASELINE_SHA)
        self.published_baseline=published_baseline or frozen_json(HERE/'read-only-observation/CURRENT_PRIOR.raw.json',PUBLISHED_BASELINE_SHA)
        self.old=self.baseline['files']
        require(len(self.old)==INHERITED_COUNT and type(self.baseline.get('expected_bytes')) is int
                and self.baseline['expected_bytes']>0
                and sum(p['bytes'] for p in self.old)==self.baseline['expected_bytes']
                and self.baseline['concept_id']==int(PARENT_ID) and self.baseline['current_record']==int(PRIOR_ID), 'BASELINE_IDENTITY_DIFFERS')
        self.prior=prior_receipt or frozen_json(PRIOR_RECEIPT, PRIOR_RECEIPT_SHA)
        self.prior_rows={r['filename']:r for r in self.prior['rows']}
        require(set(self.prior_rows)=={p['filename'] for p in self.old}, 'PRIOR_SHA_RECEIPT_MEMBERSHIP_DIFFERS')
        for p in self.old:
            r=self.prior_rows[p['filename']]
            require(r['status']=='PASS_PRIOR_COMPLETE_BYTES_WITH_FRESH_IMMUTABLE_IDENTITY'
                    and canonical_bytes(r['observed'])==canonical_bytes({k:p[k] for k in ('bytes','md5','sha256')})
                    and r.get('remote_file_id') and r.get('remote_version_id'), 'PRIOR_SHA_RECEIPT_DIFFERS')
        self.inventory=inventory; self.metadata=metadata; self.inventory_sha=inventory_sha; self.metadata_sha=metadata_sha
        self.state_path=self.dir/'STATE.json'
        binding={'prior_record_id':PRIOR_ID,'parent_id':PARENT_ID,'controller_sha256':sha(Path(__file__).read_bytes()),
                 'baseline_sha256':sha(canonical_bytes(self.baseline)),
                 'published_baseline_sha256':sha(canonical_bytes(self.published_baseline)),
                 'prior_receipt_sha256':sha(canonical_bytes(self.prior))}
        if self.state_path.exists():
            self.state=json.loads(self.state_path.read_text())
            require(self.state.get('continuation_binding')==binding,'JOURNAL_BASELINE_OR_CONTROLLER_BINDING_DIFFERS')
        else:
            require(not (self.dir/'JOURNAL.jsonl').exists(),'ORPHAN_JOURNAL_REQUIRES_MANUAL_REVIEW')
            self.state={'pending':None,'draft_id':None,'published':False,'continuation_binding':binding}
        self.add=[]
        if inventory is not None:
            self.validate_inputs()
        for key,value in (('inventory_sha256',inventory_sha),('metadata_sha256',metadata_sha)):
            if value:
                require(not self.state.get(key) or self.state[key]==value, 'JOURNAL_INPUT_PINS_DIFFER')
                self.state[key]=value
        self.save()
    def validate_inputs(self):
        i=self.inventory
        require(i.get('sealed_upload_manifest') is True and i.get('concept_doi')=='10.5281/zenodo.17088132', 'UNSEALED_OR_WRONG_FAMILY_INVENTORY')
        old=i.get('inherited_files',[]); self.add=i.get('new_files',[])
        require(len(old)==INHERITED_COUNT and len(self.add)==ADDITION_COUNT, 'INVENTORY_COUNTS_DIFFER')
        require(i.get('scientific_replay_status')==SCIENTIFIC_REPLAY_STATUS,'COMPLETED_SCIENTIFIC_REPLAY_REQUIRED')
        require(canonical_bytes([{k:p[k] for k in ('filename','bytes','md5','sha256')} for p in old])
                ==canonical_bytes([{k:p[k] for k in ('filename','bytes','md5','sha256')} for p in self.old]), 'INHERITED_PINS_DIFFER')
        names=[p['filename'] for p in old+self.add]
        require(len(names)==len(set(names))==FINAL_COUNT, 'INVENTORY_NAMES_NOT_UNIQUE')
        for p in self.add:
            require(Path(p['filename']).name==p['filename'] and '/' not in p['filename'] and '\\' not in p['filename'], 'UNSAFE_ASSET_FILENAME')
            require(type(p['bytes']) is int and 0<p['bytes']<=MAX_ASSET
                    and re.fullmatch('[0-9a-f]{32}',p.get('md5','')) and re.fullmatch('[0-9a-f]{64}',p.get('sha256','')), 'INVALID_ASSET_PIN')
            b=Path(p['local_path']).read_bytes()
            require(len(b)==p['bytes'] and hashlib.md5(b).hexdigest()==p['md5'] and sha(b)==p['sha256'], 'LOCAL_ASSET_PIN_DIFFERS')
        total=sum(p['bytes'] for p in old+self.add)
        require(type(i.get('total_bytes_final')) is int and i['total_bytes_final']==total
                and type(i.get('expected_total_files_final')) is int
                and i['expected_total_files_final']==FINAL_COUNT, 'FINAL_INVENTORY_TOTAL_DIFFERS')
        require(self.metadata is not None and set(self.metadata)=={'metadata','custom_fields','access'}, 'FULL_METADATA_REQUIRED')
        require(self.metadata['metadata'].get('version') not in (None,'',PRIOR_SEMANTIC_VERSION), 'NEW_SEMANTIC_VERSION_REQUIRED')
        require(not any(k in self.metadata for k in ('id','parent','pids','files','versions','links')), 'SERVER_ASSIGNED_METADATA_FIELDS')
        # Preserve all inherited metadata except explicitly reviewed update fields.
        prior=self.published_baseline
        for key in ('custom_fields','access'):
            expected=prior.get(key,{})
            if key=='access':expected={k:v for k,v in expected.items() if k!='status'}
            require(value_equal(self.metadata[key],expected,(key,)),'INHERITED_'+key.upper()+'_DIFFERS')
        old_meta=prior['metadata'];new_meta=self.metadata['metadata']
        require(set(new_meta)==set(old_meta),'INHERITED_METADATA_FIELD_SET_DIFFERS')
        editable={'title','description','additional_descriptions','publication_date','version','related_identifiers'}
        for key in set(old_meta)-editable:
            require(value_equal(new_meta[key],old_meta[key],('metadata',key)),'INHERITED_METADATA_'+key.upper()+'_DIFFERS')
        before=old_meta.get('additional_descriptions',[]);after=new_meta.get('additional_descriptions',[])
        require(isinstance(before,list) and isinstance(after,list) and len(after)==len(before)>0,'NOTES_SHAPE_DIFFERS')
        require(value_equal(after[1:],before[1:],('metadata','additional_descriptions')),'HISTORICAL_NOTES_DIFFERS')
        require(set(after[0])==set(before[0]) and all(value_equal(after[0][k],before[0][k],('metadata','additional_descriptions',0,k)) for k in set(before[0])-{'description'}),'CURRENT_NOTES_TYPE_OR_LANGUAGE_DIFFERS')
        previous=old_meta.get('related_identifiers',[]);related=new_meta.get('related_identifiers',[])
        require(isinstance(previous,list) and isinstance(related,list) and len(related)>=len(previous),'INHERITED_REFERENCES_REMOVED')
        require(value_equal(related[:len(previous)],previous,('metadata','related_identifiers')),'INHERITED_REFERENCES_DIFFERS')
        allowed_append={'identifier':'10.5281/zenodo.'+PRIOR_ID,'relation_type':{'id':'references'},'resource_type':{'id':'software'},'scheme':'doi'}
        require(len(related)==len(previous) or (len(related)==len(previous)+1 and value_equal(related[-1],allowed_append,('metadata','related_identifiers',len(previous)))),'UNREVIEWED_REFERENCE_ADDITION')

    def save(self):
        atomic_json(self.state_path,self.state,self.token)
    def event(self, kind, **values):
        obj=redact({'utc':utc(),'kind':kind,**values},self.token)
        with (self.dir/'JOURNAL.jsonl').open('a',encoding='utf-8') as f:
            f.write(json.dumps(obj,ensure_ascii=False)+'\n');f.flush();os.fsync(f.fileno())
    def capture(self, label, obj):
        atomic_json(self.dir/(label+'.json'),obj,self.token)
    def get(self,url,label,accept=VENDOR,public=False):
        r=self.transport.request('GET',safe_url(url,query=True),headers={'Accept':accept,'_public':public})
        self.event('GET',url=url,status=r.status,response_sha256=sha(r.body),authenticated_request=not public)
        require(r.status==200, 'GET_'+label+'_FAILED')
        o=r.json(); self.capture(label,o)
        return o,r.headers
    def write(self,kind,url,data=None,headers=None,**context):
        require(self.inventory is not None and self.metadata is not None, 'FROZEN_INPUTS_REQUIRED_FOR_WRITE')
        require(self.state.get('pending') is None, 'AMBIGUOUS_WRITE_REQUIRES_READ_ONLY_RECONCILIATION')
        safe_url(url)
        body=data if isinstance(data,bytes) else (json.dumps(data,ensure_ascii=False).encode() if data is not None else b'')
        intent={'kind':kind,'url':url,'body_sha256':sha(body),'body_bytes':len(body),'context':context,'utc':utc()}
        self.state['pending']=intent;self.save()
        self.event('WRITE_INTENT',operation=kind,**{k:v for k,v in intent.items() if k!='kind'})
        accept='application/json' if kind in ('initialize_file','put_content','commit_file') else VENDOR
        h={'Accept':accept,'Content-Type':'application/json'};h.update(headers or {})
        try:
            r=self.transport.request('PUT' if kind in ('put_content','save_metadata') else 'POST',url,data=body,headers=h)
        except Exception:
            self.event('WRITE_OUTCOME_UNKNOWN',operation=kind);raise Stop('WRITE_OUTCOME_UNKNOWN_RECONCILE_ONLY') from None
        self.event('WRITE_RESPONSE',operation=kind,status=r.status,response_sha256=sha(r.body),
                   body_truncated=r.body_truncated,body_read_failed=r.body_read_failed)
        if r.status>=400:
            diagnostic={'operation':kind,'http_status':r.status,'captured_body_bytes':len(r.body),
                        'captured_body_sha256':sha(r.body),'body_truncated':r.body_truncated,
                        'body_read_failed':r.body_read_failed}
            if r.body_read_failed:
                diagnostic['format']='BODY_READ_FAILED_NO_RAW_BODY_SAVED'
            elif r.body_truncated or len(r.body)>MAX_JSON:
                diagnostic['format']='OVERSIZED_BODY_NO_RAW_BODY_SAVED'
            else:
                try:
                    diagnostic['server_error']=bounded_error_json(r.body)
                    diagnostic['format']='JSON'
                except Exception:
                    diagnostic['format']='NON_JSON_OR_MALFORMED_NO_RAW_BODY_SAVED'
            label='error_'+kind+'_'+str(time.time_ns())
            self.capture(label,diagnostic)
            self.event('ERROR_DIAGNOSTIC_CAPTURED',operation=kind,http_status=r.status,
                       path=str(self.dir/(label+'.json')),format=diagnostic['format'])
        if 400<=r.status<500:
            self.state['pending']=None;self.save();raise Stop('WRITE_REJECTED_'+str(r.status))
        require(r.status in (200,201,202,204), 'WRITE_OUTCOME_UNKNOWN_RECONCILE_ONLY')
        o=r.json() if r.body else {}
        self.capture('response_'+kind+'_'+str(time.time_ns()),o)
        self.state['pending']=None;self.save()
        return o,r.headers
    def owner_list(self,max_pages=20,page_size=100):
        family=[];seen=set(); exhausted=False
        for page in range(1,max_pages+1):
            url=ORIGIN+'/api/deposit/depositions?'+parse.urlencode({'page':page,'size':page_size,'sort':'mostrecent'})
            r=self.transport.request('GET',url,headers={'Accept':'application/json'})
            self.event('OWNER_PAGE_GET',page=page,status=r.status,response_sha256=sha(r.body))
            require(r.status==200,'OWNER_LIST_FAILED')
            body=r.json();items=body if isinstance(body,list) else body.get('hits',{}).get('hits')
            require(isinstance(items,list) and len(items)<=page_size,'OWNER_LIST_SHAPE_DIFFERS')
            ids=[str(x.get('id')) for x in items]
            require(all(x.isdigit() for x in ids),'OWNER_LIST_IDS_INVALID')
            require(len(set(ids))==len(ids) and not (set(ids)&seen),'OWNER_PAGINATION_REPEATED_IDS');seen.update(ids)
            for x in items:
                parent=str(x.get('conceptrecid',x.get('parent',{}).get('id','')))
                if parent==PARENT_ID:
                    owner=x.get('owner')
                    owner=owner.get('id') if isinstance(owner,dict) else owner
                    require(str(owner)==str(EXPECTED_OWNER),'OWNED_FAMILY_OWNER_DIFFERS')
                    require(type(x.get('submitted')) is bool or type(x.get('is_draft')) is bool,'OWNED_FAMILY_STATE_UNKNOWN')
                    family.append(x)
            self.capture('owner_page_'+str(page),{'page':page,'returned_count':len(items),'family_entries':[x for x in items if str(x.get('conceptrecid',x.get('parent',{}).get('id','')))==PARENT_ID]})
            if len(items)<page_size:exhausted=True;break
        require(exhausted,'OWNER_PAGINATION_NOT_EXHAUSTED')
        candidates=[x for x in family if x.get('submitted') is False or x.get('is_draft') is True]
        require(len(candidates)<=1,'MULTIPLE_OWNED_FAMILY_DRAFTS')
        prior_owned=any(str(x.get('id'))==PRIOR_ID and x.get('submitted') is True for x in family)
        journaled_draft=self.state.get('draft_id')
        journal_published=self.state.get('published')
        matching=[x for x in family if str(x.get('id'))==journaled_draft]
        latest_only_match=(isinstance(journaled_draft,str) and journaled_draft.isdigit()
                           and journaled_draft!=PRIOR_ID and type(journal_published) is bool
                           and len(matching)==1 and matching[0].get('submitted') is journal_published
                           and ('is_draft' not in matching[0]
                                or matching[0]['is_draft'] is (not journal_published)))
        require(prior_owned or latest_only_match,'PUBLISHED_PRIOR_NOT_IN_OWNED_LIST')
        if not prior_owned:
            self.event('OWNED_LATEST_ONLY_JOURNALED_RECORD_MATCH',record_id=journaled_draft,
                       published=journal_published,original_published_record_omitted_from_owner_listing=True)
        return candidates[0] if candidates else None
    def current(self):
        o,h=self.get(ORIGIN+'/api/records/'+PRIOR_ID,'CURRENT_PRIOR',public=True)
        require(o.get('id')==PRIOR_ID and o.get('parent',{}).get('id')==PARENT_ID
                and o.get('is_published') is True,'CURRENT_PRIOR_IDENTITY_DIFFERS')
        require(str(o.get('parent',{}).get('access',{}).get('owned_by',{}).get('user'))==str(EXPECTED_OWNER),'CURRENT_PRIOR_OWNER_DIFFERS')
        metadata_equal(o,{k:self.published_baseline[k] for k in ('metadata','custom_fields','access')})
        latest_url=o['links']['latest']
        require(latest_url==ORIGIN+'/api/records/'+PRIOR_ID+'/versions/latest','LATEST_LINK_DIFFERS')
        r=self.transport.request('GET',latest_url,headers={'Accept':VENDOR,'_public':True})
        self.event('LATEST_LINK_GET',url=latest_url,status=r.status,response_sha256=sha(r.body),authenticated_request=False)
        if r.status in (301,302,303,307,308):
            target=next((v for k,v in r.headers.items() if k.lower()=='location'),None)
            safe_url(target)
            require(re.fullmatch(r'https://zenodo\.org/api/records/\d+',target),'LATEST_REDIRECT_DESTINATION_DIFFERS')
            self.event('EXPLICIT_PUBLIC_LATEST_RESOLUTION',source=latest_url,target=target)
            latest,lh=self.get(target,'CURRENT_LATEST',public=True)
        else:
            require(r.status==200,'GET_CURRENT_LATEST_FAILED')
            latest,lh=r.json(),r.headers;self.capture('CURRENT_LATEST',latest)
        if latest.get('id')==self.state.get('draft_id') and latest.get('is_published') is True and self.metadata:
            metadata_equal(latest,self.metadata);assert_pins(file_map(latest),self.old+self.add)
            self.state['published']=True;self.save()
        if self.state.get('published'):
            require(latest['id']==self.state['draft_id'],'LATEST_PUBLISHED_DIFFERENT')
        else:
            pending_publish=(self.state.get('pending') or {}).get('kind')=='publish'
            require(latest['id']==PRIOR_ID or (pending_publish and latest['id']==self.state.get('draft_id')),'LATEST_CHANGED_FROM_FROZEN_PRIOR')
        require(o['links']['files']==ORIGIN+'/api/records/'+PRIOR_ID+'/files','PRIOR_FILES_LINK_DIFFERS')
        files,fh=self.get(o['links']['files'],'CURRENT_PRIOR_FILES',accept='application/json',public=True)
        original_files=file_map(files);assert_pins(original_files,self.old)
        for p in self.old:
            f=original_files[p['filename']];r=self.prior_rows[p['filename']]
            require(f.get('file_id',f.get('id'))==r['remote_file_id']
                    and f.get('version_id')==r['remote_version_id'],'ORIGINAL_PUBLISHED_IMMUTABLE_FILE_IDENTITY_CHANGED')
        self.current_prior_metadata=o['metadata']
        return o
    def draft(self,id_):
        require(str(id_).isdigit() and str(id_)!=PRIOR_ID,'PRIOR_RECORD_CANNOT_BE_NEW_EDITION_DRAFT')
        o,h=self.get(ORIGIN+'/api/records/'+str(id_)+'/draft','DRAFT_'+str(id_))
        require(o.get('id')==str(id_) and o.get('parent',{}).get('id')==PARENT_ID
                and o.get('is_draft') is True and o.get('is_published') is False,'DRAFT_FAMILY_OR_STATE_DIFFERS')
        require(str(o.get('parent',{}).get('access',{}).get('owned_by',{}).get('user'))==str(EXPECTED_OWNER),'DRAFT_OWNER_DIFFERS')
        require(o['links']['self']==ORIGIN+'/api/records/'+str(id_)+'/draft','DRAFT_SELF_LINK_DIFFERS')
        require(o['links']['files']==o['links']['self']+'/files','DRAFT_FILES_LINK_DIFFERS')
        if self.metadata:
            require(o.get('metadata',{}).get('title') in (self.current_prior_metadata['title'],self.metadata['metadata']['title']),'OWNED_DRAFT_DIFFERENT_TITLE')
            for key in ('creators','rights','resource_type','publisher','copyright','languages'):
                require(value_equal(o.get('metadata',{}).get(key),self.metadata['metadata'].get(key),('metadata',key)),'DRAFT_CORE_METADATA_DIFFERS')
            require(o.get('metadata',{}).get('version') in (None,'',PRIOR_SEMANTIC_VERSION,self.metadata['metadata']['version']),'OWNED_DRAFT_DIFFERENT_SEMANTIC_VERSION')
        files,fh=self.get(o['links']['files'],'DRAFT_FILES_'+str(id_),accept='application/json')
        return o,h,file_map(files)
    def reconcile(self):
        current=self.current(); candidate=self.owner_list()
        if candidate:
            id_=str(candidate['id']);o,h,files=self.draft(id_)
            allowed={p['filename'] for p in self.old+self.add}
            require(set(files)<=allowed,'OWNED_DRAFT_HAS_UNEXPECTED_FILES')
            assert_pins(files,self.old,exact=False)
            require(not self.state.get('draft_id') or self.state['draft_id']==id_,'OWNED_DRAFT_ID_DIFFERS_FROM_JOURNAL')
            self.state['draft_id']=id_;self.save()
        else:
            o=h=files=None
            require(not self.state.get('draft_id') or self.state.get('published') or (self.state.get('pending') or {}).get('kind')=='publish','JOURNALED_DRAFT_NOT_IN_OWNED_LIST')
        pending=self.state.get('pending')
        if pending:
            kind=pending['kind'];ctx=pending.get('context',{});recovered=False
            if kind=='create_version' and candidate:recovered=True
            elif kind=='save_metadata' and o:
                metadata_equal(o,self.metadata);recovered=True
            elif kind in ('initialize_file','commit_file','put_content') and files and ctx.get('filename') in files:
                f=files[ctx['filename']];p=next(x for x in self.add if x['filename']==ctx['filename'])
                if kind=='initialize_file':recovered=True
                elif kind=='commit_file':assert_pins({p['filename']:f},[p]);recovered=True
                else:
                    self.stream_check(f['links']['content'],p)
                    self.state.setdefault('content_verified_files',[]).append(p['filename']);recovered=True
            elif kind=='publish' and self.state.get('draft_id'):
                p,ph=self.get(ORIGIN+'/api/records/'+self.state['draft_id'],'RECONCILE_PUBLISH',public=True)
                require(p.get('is_published') is True and p.get('parent',{}).get('id')==PARENT_ID,'PUBLISH_OUTCOME_STILL_UNKNOWN')
                metadata_equal(p,self.metadata);assert_pins(file_map(p),self.old+self.add)
                self.state['published']=True;recovered=True
            require(recovered,'AMBIGUOUS_WRITE_NOT_RECONCILED')
            self.event('WRITE_RECONCILED_FROM_GETS',operation=kind);self.state['pending']=None;self.save()
        return {'current_prior':PRIOR_ID,'draft_id':self.state.get('draft_id'),'published':self.state.get('published'),'owner_listing_exhausted':True}
    def create(self,retry_create_406=False):
        result=self.reconcile()
        if result['draft_id']:return result
        if retry_create_406:
            require(self.state.get('pending') is None and not result['draft_id']
                    and self.state.get('published') is False,'CREATE_406_RETRY_REQUIRES_NO_PENDING_DRAFT_OR_PUBLICATION')
            require(self.state.get('create_post_attempted') is True,'CREATE_406_RETRY_REQUIRES_ORIGINAL_ATTEMPT')
            require(not self.state.get('create_406_retry_attempted'),'CREATE_406_RETRY_ALREADY_ATTEMPTED')
            journal=[json.loads(line) for line in (self.dir/'JOURNAL.jsonl').read_text().splitlines()]
            outcomes=[x for x in journal if x.get('operation')=='create_version'
                      and x.get('kind') in ('WRITE_INTENT','WRITE_RESPONSE','WRITE_OUTCOME_UNKNOWN')]
            require(len(outcomes)==2 and outcomes[0].get('kind')=='WRITE_INTENT'
                    and outcomes[-1].get('kind')=='WRITE_RESPONSE'
                    and type(outcomes[-1].get('status')) is int and outcomes[-1]['status']==406,
                    'CREATE_406_RETRY_REQUIRES_DEFINITIVE_LAST_406')
            require(outcomes[0].get('url')==ORIGIN+'/api/deposit/depositions/'+PRIOR_ID+'/actions/newversion'
                    and type(outcomes[0].get('body_bytes')) is int and outcomes[0]['body_bytes']==0
                    and outcomes[0].get('body_sha256')==sha(b''),'CREATE_406_RETRY_ORIGINAL_INTENT_DIFFERS')
        else:
            require(not self.state.get('create_post_attempted'),'NEWVERSION_POST_ALREADY_ATTEMPTED_RECONCILE_ONLY')
        legacy,lh=self.get(ORIGIN+'/api/deposit/depositions/'+PRIOR_ID,'LEGACY_PRIOR',accept='application/json')
        require(str(legacy.get('conceptrecid'))==PARENT_ID and legacy.get('submitted') is True,'LEGACY_PRIOR_DIFFERS')
        url=legacy['links']['newversion']
        require(url==ORIGIN+'/api/deposit/depositions/'+PRIOR_ID+'/actions/newversion','NEWVERSION_ACTION_DIFFERS')
        if retry_create_406:
            self.state['create_406_retry_attempted']=True;self.save()
            self.event('REVIEWED_406_CREATE_RETRY_INTENT',original_attempt_preserved=True,
                       reason='DEFINITIVE_HTTP_406_LEGACY_ACCEPT_CORRECTION')
        else:
            self.state['create_post_attempted']=True;self.save()
        response,rh=self.write('create_version',url,headers={'Accept':'application/json'})
        id_=str(response.get('id',''))
        if id_==PRIOR_ID:
            link=safe_url(response.get('links',{}).get('latest_draft',''))
            match=re.fullmatch(r'https://zenodo\.org/api/deposit/depositions/(\d+)',link)
            require(match and match.group(1)!=PRIOR_ID,'NEW_DRAFT_LINK_NOT_IDENTIFIED')
            id_=match.group(1)
        require(id_.isdigit() and id_!=PRIOR_ID,'NEW_DRAFT_ID_NOT_IDENTIFIED')
        o,h,files=self.draft(id_);assert_pins(files,self.old)
        self.state['draft_id']=id_;self.save()
        return self.reconcile()
    def prepare(self):
        self.reconcile();id_=self.state.get('draft_id');require(id_,'CREATE_OR_IDENTIFY_NEW_DRAFT_FIRST')
        o,h,files=self.draft(id_)
        for p in self.add:
            name=p['filename']
            if name in files and files[name].get('status')=='completed':
                assert_pins({name:files[name]},[p]);continue
            if name not in files:
                self.write('initialize_file',o['links']['files'],[{'key':name}],filename=name)
                o,h,files=self.draft(id_)
            f=files[name]
            prefix=ORIGIN+'/api/records/'+id_+'/draft/files/'+parse.quote(name,safe='')
            require(f['links']['content']==prefix+'/content' and f['links']['commit']==prefix+'/commit','UPLOAD_LINKS_DIFFER')
            payload=Path(p['local_path']).read_bytes()
            require(len(payload)==p['bytes'] and sha(payload)==p['sha256'] and hashlib.md5(payload).hexdigest()==p['md5'],'ASSET_CHANGED_BEFORE_UPLOAD')
            if name not in self.state.get('content_verified_files',[]):
                self.write('put_content',f['links']['content'],payload,headers={'Content-Type':'application/octet-stream'},filename=name)
            else:
                self.stream_check(f['links']['content'],p)
            self.write('commit_file',f['links']['commit'],filename=name)
            o,h,files=self.draft(id_);assert_pins({name:files[name]},[p])
        assert_pins(files,self.old+self.add)
        o,h,files=self.draft(id_)
        etag=next((v for k,v in h.items() if k.lower()=='etag'),None)
        require(etag,'DRAFT_ETAG_REQUIRED_FOR_METADATA_WRITE')
        try:
            metadata_equal(o,self.metadata)
        except Stop:
            self.write('save_metadata',o['links']['self'],self.metadata,
                       headers={'If-Match':revision_if_match(etag,o)})
        return self.verify()
    def stream_check(self,url,pin,public=False):
        safe_url(url);count=0;md5=hashlib.md5();digest=hashlib.sha256()
        for chunk in self.transport.stream(url,headers={'_public':public}):
            count+=len(chunk);require(count<=pin['bytes'],'CONTENT_EXCEEDS_FROZEN_SIZE');md5.update(chunk);digest.update(chunk)
        observed={'bytes':count,'md5':md5.hexdigest(),'sha256':digest.hexdigest()}
        require(canonical_bytes(observed)==canonical_bytes({k:pin[k] for k in ('bytes','md5','sha256')}),'CONTENT_COMPLETE_MD5_SHA256_DIFFERS')
        self.event('FULL_CONTENT_STREAM_VERIFIED',filename=pin['filename'],observed=observed)
        return observed
    def verify(self,published=None):
        if published is None:
            published=bool(self.state.get('published'))
        require(self.inventory is not None and self.state.get('draft_id'),'FINAL_INPUTS_AND_NEW_ID_REQUIRED')
        require(self.state.get('pending') is None,'PENDING_WRITE_BLOCKS_VERIFY')
        self.current();id_=self.state['draft_id']
        if published:
            o,h=self.get(ORIGIN+'/api/records/'+id_,'NEW_PUBLIC_RECORD',public=True)
            require(o.get('is_published') is True and o.get('parent',{}).get('id')==PARENT_ID,'PUBLIC_NEW_RECORD_IDENTITY_DIFFERS')
            require(o['links']['files']==ORIGIN+'/api/records/'+id_+'/files','NEW_PUBLIC_FILES_LINK_DIFFERS')
            listing,lh=self.get(o['links']['files'],'NEW_PUBLIC_FILES',accept='application/json',public=True);files=file_map(listing)
        else:
            o,h,files=self.draft(id_)
        metadata_equal(o,self.metadata);assert_pins(files,self.old+self.add)
        content=[]
        for p in self.old+self.add:
            f=files[p['filename']]
            prefix=ORIGIN+'/api/records/'+id_+('/files/' if published else '/draft/files/')+parse.quote(p['filename'],safe='')+'/content'
            require(f['links']['content']==prefix,'CONTENT_LINK_DIFFERS_FROM_VERIFIED_RECORD_AND_FILENAME')
            file_id,version_id=immutable_file_identity(f)
            self.stream_check(f['links']['content'],p,public=published)
            content.append({'filename':p['filename'],'sha256':p['sha256'],'basis':'FRESH_COMPLETE_CONTENT_STREAM',
                            'remote_file_id':file_id,'remote_version_id':version_id})
        require(len(content)==FINAL_COUNT and all(x['basis']=='FRESH_COMPLETE_CONTENT_STREAM' for x in content),'ALL_EDITION_FILES_REQUIRE_FRESH_STREAMS')
        receipt={'status':'PASS_COMPLETE_'+str(FINAL_COUNT)+'_FILE_EDITION_READBACK','utc':utc(),'draft_id':id_,'parent_id':PARENT_ID,'inventory_sha256':self.inventory_sha,'metadata_sha256':self.metadata_sha,'published':published,'file_count':FINAL_COUNT,'total_bytes':self.inventory['total_bytes_final'],'content':content,'etag':next((v for k,v in h.items() if k.lower()=='etag'),None)}
        self.capture('VERIFIED_PUBLIC' if published else 'VERIFIED_DRAFT',receipt)
        if not published:
            self.state['verified_receipt_sha256']=sha((self.dir/'VERIFIED_DRAFT.json').read_bytes());self.save()
        return receipt
    def publish(self):
        require(not self.state.get('published'),'EDITION_ALREADY_PUBLISHED_NO_REPEAT_POST')
        require(self.state.get('verified_receipt_sha256'),'VERIFY_STAGE_REQUIRED_BEFORE_PUBLISH')
        require(sha((self.dir/'VERIFIED_DRAFT.json').read_bytes())==self.state['verified_receipt_sha256'],'VERIFICATION_RECEIPT_CHANGED')
        self.reconcile();require(not self.state.get('published'),'EDITION_ALREADY_PUBLISHED_NO_REPEAT_POST')
        require(not self.state.get('publish_post_attempted'),'PUBLISH_POST_ALREADY_ATTEMPTED_RECONCILE_ONLY')
        verified=self.verify();o,h,files=self.draft(self.state['draft_id'])
        metadata_equal(o,self.metadata);assert_pins(files,self.old+self.add)
        verified_identities={r['filename']:[r['remote_file_id'],r['remote_version_id']] for r in verified['content']}
        fresh_identities={p['filename']:list(immutable_file_identity(files[p['filename']])) for p in self.old+self.add}
        require(canonical_bytes(fresh_identities)==canonical_bytes(verified_identities),
                'DRAFT_CONTENT_IDENTITIES_CHANGED_SINCE_COMPLETE_VERIFICATION')
        url=o['links']['publish'];require(url==o['links']['self']+'/actions/publish','PUBLISH_LINK_DIFFERS')
        etag=next((v for k,v in h.items() if k.lower()=='etag'),None);require(etag,'DRAFT_ETAG_REQUIRED_FOR_PUBLISH')
        require(etag==verified['etag'],'DRAFT_CHANGED_SINCE_COMPLETE_VERIFICATION')
        revision_header=revision_if_match(etag,o)
        self.state['publish_post_attempted']=True;self.save()
        self.write('publish',url,headers={'If-Match':revision_header})
        self.state['published']=True;self.save()
        for attempt in range(20):
            r=self.transport.request('GET',ORIGIN+'/api/records/'+self.state['draft_id'],headers={'Accept':VENDOR,'_public':True})
            self.event('PUBLICATION_GET_POLL',attempt=attempt+1,status=r.status,response_sha256=sha(r.body))
            if r.status==200 and r.json().get('is_published') is True:
                break
            require(r.status in (200,404,503),'PUBLICATION_READBACK_HTTP_ERROR')
            if attempt<19:
                time.sleep(1)
        else:
            raise Stop('PUBLICATION_READBACK_PENDING_RETRY_VERIFY_ONLY')
        result=self.verify(published=True)
        self.current()
        return result

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=['reconcile','create','prepare','verify','publish'],default='reconcile')
    p.add_argument('--evidence-dir',required=True)
    p.add_argument('--inventory');p.add_argument('--inventory-sha256')
    p.add_argument('--metadata');p.add_argument('--metadata-sha256')
    p.add_argument('--mutate',action='store_true');p.add_argument('--publish',action='store_true')
    p.add_argument('--retry-create-406',action='store_true')
    a=p.parse_args(argv)
    require(a.stage not in ('create','prepare','publish') or a.mutate,'MUTATION_STAGE_REQUIRES_EXPLICIT_MUTATE_FLAG')
    require(a.stage!='publish' or a.publish,'PUBLISH_STAGE_REQUIRES_EXPLICIT_PUBLISH_FLAG')
    require(not a.publish or a.stage=='publish','PUBLISH_FLAG_ONLY_VALID_FOR_PUBLISH_STAGE')
    require(not a.retry_create_406 or (a.stage=='create' and a.mutate),'CREATE_406_RETRY_FLAG_ONLY_FOR_EXPLICIT_CREATE_MUTATION')
    require(bool(a.inventory)==bool(a.metadata),'BOTH_FROZEN_INPUTS_REQUIRED')
    inventory=frozen_json(a.inventory,a.inventory_sha256) if a.inventory else None
    metadata=frozen_json(a.metadata,a.metadata_sha256) if a.metadata else None
    require(a.stage=='reconcile' or inventory is not None,'FROZEN_FINAL_INPUTS_REQUIRED')
    token=os.environ.get('ZENODO_ACCESS_TOKEN','')
    with EvidenceLock(a.evidence_dir):
        c=Controller(a.evidence_dir,Transport(token),token,inventory,metadata,a.inventory_sha256,a.metadata_sha256)
        result=c.create(retry_create_406=a.retry_create_406) if a.stage=='create' else getattr(c,a.stage)()
        print(json.dumps(redact({'stage':a.stage,'status':'COMPLETE','record_id':c.state.get('draft_id'),'summary':result},token),ensure_ascii=False))
    return 0

if __name__=='__main__':
    try:
        raise SystemExit(main())
    except Stop as e:
        print(json.dumps({'status':'STOPPED','reason':str(e)}));raise SystemExit(2)
    except Exception:
        print(json.dumps({'status':'STOPPED','reason':'UNEXPECTED_FAILURE_DETAILS_NOT_OUTPUT'}));raise SystemExit(3)
