#!/usr/bin/env python3
"""Bounded mechanical replay/watch candidate; no publishing or model requests."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone, timedelta
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from scientific_projection import POLICY as SCIENTIFIC_HASH_POLICY, project as scientific_projection

REPOSITORY = 'maldonado-research/HDblast'
CHECKPOINT = 'research/HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER'
CLASSIFICATIONS = {'LEDGER_ERROR_DEMONSTRATED', 'NO_GATE_SCALE_ATTRIBUTION', 'CONSISTENCY_FAILURE'}
FIELDS = {'schema_version','state','repository','checkpoint_relative_path','science_commit',
          'freeze_commit','registration_sha256','requirements_sha256','replay_driver_sha256',
          'expected_classification','expected_old_metric_status','expected_scientific_sha256',
          'scientific_hash_policy'}
MAX_CAMPAIGN_DAYS = 30
MAX_RESPONSE_BYTES = 4*1024*1024
QUERY_SIZE = 20
QUERIES = {
    'adiabatic_scalar': 'find date > 2023 and (title adiabatic or title renormalization) and (title scalar or title semiclassical or title cosmology)',
    'ward_response': 'find date > 2023 and (title "linear response" or title Ward or title conservation) and (title semiclassical or title scalar or title cosmology)',
}

def require(ok, message):
    if not ok:
        raise ValueError(message)

def stamp():
    return datetime.now(timezone.utc).isoformat()

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def load_json(path):
    def pairs(items):
        obj={}
        for key,value in items:
            require(key not in obj, 'Duplicate JSON key')
            obj[key]=value
        return obj
    def constant(value):
        raise ValueError('Nonfinite JSON literal')
    return json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=constant)

def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')

def hex_string(value,length):
    return isinstance(value,str) and re.fullmatch('[0-9a-f]{%d}'%length,value) is not None

def safe_relative(value):
    require(isinstance(value,str) and value and '\\' not in value and ':' not in value,
            'Unsafe relative path')
    parts=value.split('/')
    require(not PurePosixPath(value).is_absolute() and all(p not in ('','.','..') for p in parts),
            'Unsafe relative path')
    return parts

def safe_file(root,name):
    root=Path(root).absolute()
    require(root.is_dir() and not root.is_symlink() and all(not p.is_symlink() for p in root.parents),
            'Missing or symlink source root')
    path=root
    for part in safe_relative(name):
        path=path/part
        require(not path.is_symlink(), 'Symlink source path')
    require(path.is_file() and path.resolve().is_relative_to(root.resolve()),'Missing or escaped source file')
    return path

def fresh_output(path):
    path=Path(path).absolute()
    require(not path.exists() and not path.is_symlink() and all(not p.is_symlink() for p in path.parents),
            'Output must be fresh and have no symlink ancestors')
    path.mkdir(parents=True)
    return path

def validate_target(target):
    require(isinstance(target,dict) and set(target)==FIELDS,'Unexpected target schema fields')
    require(type(target['schema_version']) is int and target['schema_version']==1,'Unsupported target schema')
    require(target['scientific_hash_policy']==SCIENTIFIC_HASH_POLICY,'Unexpected scientific comparison policy')
    require(target['repository']==REPOSITORY and target['checkpoint_relative_path']==CHECKPOINT,
            'Wrong repository or checkpoint')
    if target['state']=='PENDING_COMPLETED_ACTIVE_CHECKPOINT':
        require(all(target[k] is None for k in ('science_commit','freeze_commit','registration_sha256',
             'requirements_sha256','replay_driver_sha256','expected_classification')),
             'Pending target must not contain guessed pins')
        require(target['expected_scientific_sha256']=={} and target['expected_old_metric_status']=='FAIL',
                'Pending target has unexplained results')
        return False
    require(target['state']=='COMPLETED_REVIEWED_CHECKPOINT','Unrecognized target state')
    for name in ('science_commit','freeze_commit'):
        require(hex_string(target[name],40),'Expected immutable full Git SHA: '+name)
    for name in ('registration_sha256','requirements_sha256','replay_driver_sha256'):
        require(hex_string(target[name],64),'Expected SHA256: '+name)
    require(target['expected_classification'] in CLASSIFICATIONS and target['expected_old_metric_status']=='FAIL',
            'Unknown scientific classification or changed historical status')
    expected=target['expected_scientific_sha256']
    require(isinstance(expected,dict) and {'primary/diagnostic.json','independent/diagnostic.json'}<=set(expected),
            'Both completed scientific diagnostics must be pinned')
    for name,value in expected.items():
        safe_relative(name)
        require(name.split('/')[0] in {'primary','independent'} and name.endswith(('.json','.jsonl')),
                'Expected file is outside scientific route outputs')
        require(hex_string(value,64),'Expected scientific SHA256')
    return True

def parse_utc(value):
    require(isinstance(value,str) and value.endswith('Z'),'UTC campaign timestamps must end in Z')
    parsed=datetime.fromisoformat(value[:-1]+'+00:00')
    require(parsed.utcoffset()==timedelta(0),'Expected UTC timestamp')
    return parsed

def campaign_ok(start,end,now):
    start,end=parse_utc(start),parse_utc(end)
    require(timedelta(0)<end-start<=timedelta(days=MAX_CAMPAIGN_DAYS),'Campaign must last at most 30 days')
    require(start<=now<end,'Campaign is not currently active')

def gate(target,context,now):
    require(context['repository']==REPOSITORY,'Wrong workflow repository')
    require(context['branch']==context['default_branch'],'Monitor may run only on the default branch')
    require(context['attempt']=='1','Automatic unchanged reruns are disabled')
    require(context['event'] in {'schedule','workflow_dispatch'},'Unsupported active monitor event')
    if context['event']=='schedule':
        require(context['enabled']=='true','Explicit follow-up enable flag is absent')
    campaign_ok(context['start'],context['end'],now)
    ready=validate_target(target)
    if not ready:
        return {'status':'INACTIVE_PENDING_COMPLETED_TARGET','ready':False,'do_watch':False,'do_replay':False}
    mode=context['mode']
    require(mode in {'watch','replay','all'},'Unknown manual mode')
    scheduled=context['event']=='schedule'
    return {'status':'ACTIVE_BOUNDED_MECHANICAL_MONITOR','ready':True,
            'do_watch':scheduled or mode in {'watch','all'},
            'do_replay':(scheduled and now.weekday()==0) or (not scheduled and mode in {'replay','all'}),
            'science_commit':target['science_commit'],'checkpoint_relative_path':CHECKPOINT}

def normalize_records(data):
    require(isinstance(data,dict) and isinstance(data.get('hits'),dict),'Missing INSPIRE hits schema')
    hits=data['hits'].get('hits');total=data['hits'].get('total')
    require(type(total) is int and total>=0 and isinstance(hits,list) and len(hits)<=QUERY_SIZE,
            'Unexpected hit count or first-page limit exceeded')
    records=[];seen=set()
    for hit in hits:
        require(isinstance(hit,dict) and isinstance(hit.get('metadata'),dict),'Malformed metadata hit')
        identifier=str(hit.get('id',''))
        require(identifier.isdigit() and identifier not in seen,'Missing or repeated INSPIRE ID')
        seen.add(identifier);m=hit['metadata'];titles=m.get('titles')
        require(isinstance(titles,list) and titles and isinstance(titles[0],dict) and isinstance(titles[0].get('title'),str),
                'Missing primary title')
        authors=m.get('authors',[]);arxiv=m.get('arxiv_eprints',[]);dois=m.get('dois',[])
        require(all(isinstance(x,list) for x in (authors,arxiv,dois)),'Malformed bibliographic list')
        for collection,key in ((authors,'full_name'),(arxiv,'value'),(dois,'value')):
            require(all(isinstance(x,dict) and isinstance(x.get(key),str) for x in collection),
                    'Malformed bibliographic item')
        date=m.get('preprint_date')
        require(date is None or isinstance(date,str),'Malformed preprint date')
        records.append(dict(id=identifier,title=titles[0]['title'],authors=[x['full_name'] for x in authors],
                            arxiv=sorted(x['value'] for x in arxiv),doi=sorted(set(x['value'] for x in dois)),preprint_date=date))
    return sorted(records,key=lambda r:int(r['id']))

def bibliographic_diff(previous,current):
    old={r['id']:r for r in previous};new={r['id']:r for r in current}
    require(len(old)==len(previous) and len(new)==len(current),'Duplicate baseline/current IDs')
    return dict(added=[new[k] for k in sorted(new.keys()-old.keys())],
                absent_from_current_first_page=[old[k] for k in sorted(old.keys()-new.keys())],
                modified=[{'previous':old[k],'current':new[k]} for k in sorted(old.keys()&new.keys()) if old[k]!=new[k]])

class SameOriginRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self,req,fp,code,msg,headers,newurl):
        # Reject even same-origin redirects to keep exactly one HTTP attempt
        # per query; no response can send the monitor to another destination.
        raise ValueError('Metadata redirect rejected')

def public_get(url):
    parsed=urllib.parse.urlsplit(url)
    require(parsed.scheme=='https' and parsed.hostname=='inspirehep.net' and parsed.path=='/api/literature'
            and parsed.port in (None,443) and not parsed.username and not parsed.password,'Unapproved public endpoint')
    opener=urllib.request.build_opener(SameOriginRedirect())
    req=urllib.request.Request(url,headers={'User-Agent':'HDBLAST-bounded-followup-monitor/1.0','Accept':'application/json'})
    with opener.open(req,timeout=35) as response:
        raw=response.read(MAX_RESPONSE_BYTES+1)
        require(len(raw)<=MAX_RESPONSE_BYTES,'Response exceeds 4 MiB limit')
        require(response.status==200 and 'application/json' in response.headers.get('Content-Type',''),
                'Unexpected response status/type')
        return raw,dict(status=response.status,final_url=response.url,bytes=len(raw),
                        sha256=hashlib.sha256(raw).hexdigest(),response_date=response.headers.get('Date'))

def watch(baseline,output,getter=public_get):
    require(baseline.get('schema_version')==1 and set(baseline.get('queries',{}))==set(QUERIES),
            'Unexpected reviewed baseline schema')
    state={'status':'RUNNING','started_utc':stamp(),'scope':'Two fixed first-page public metadata GETs, at most 40 records; no full-text, novelty or scientific-truth assessment.','queries':{}}
    failed=False
    for label,query in QUERIES.items():
        url='https://inspirehep.net/api/literature?'+urllib.parse.urlencode({'q':query,'size':QUERY_SIZE,'sort':'mostrecent'})
        item={'query':query,'url':url,'requested_at_utc':stamp(),'method':'GET','attempts':1}
        try:
            require(baseline['queries'][label]['query']==query,'Baseline query differs')
            raw,receipt=getter(url);item.update(receipt)
            require(len(raw)<=MAX_RESPONSE_BYTES,'Getter exceeded response bound')
            data=json.loads(raw);current=normalize_records(data)
            item.update(total_hits=data['hits']['total'],returned_hits=len(current),records=current,
                        diff=bibliographic_diff(baseline['queries'][label]['records'],current))
        except Exception as error:
            failed=True;item.update(status='FAILED_READ',error_type=type(error).__name__,http_status=getattr(error,'code',None))
            # Do not save raw error bodies, headers or exception URLs.
        item['completed_utc']=stamp();state['queries'][label]=item
        write_json(output/'WATCH.json',state)
    state.update(status='FAILED_METADATA_READ' if failed else 'COMPLETED_METADATA_READ',completed_utc=stamp(),requests_attempted=2)
    write_json(output/'WATCH.json',state)
    return state

def verify_completed_replay(target,output):
    require(validate_target(target),'Cannot compare an uncomputed/pending target')
    replay=load_json(safe_file(output,'REPLAY.json'))
    require(replay.get('status')=='COMPLETED_REPLAY' and replay.get('plan_only') is False,'Replay is not complete physical reproduction')
    require(replay.get('classification')==target['expected_classification'] and replay.get('old_metric_status')=='FAIL',
            'Classification differs from completed reference or history changed')
    require(replay.get('physical_routes_completed')==2 and replay.get('successful_commands')==4,'Incomplete route/validator execution')
    require(replay.get('source_verification_before')=='PASS' and replay.get('source_verification_after')=='PASS','Source authentication failed')
    checks={}
    for name,digest in target['expected_scientific_sha256'].items():
        path=safe_file(Path(output)/'fresh',name)
        raw_digest=sha(path)
        if name in {'primary/diagnostic.json','independent/diagnostic.json'}:
            measured,observations=scientific_projection(load_json(path),name.split('/')[0])
            checks[name]={'comparison_policy':SCIENTIFIC_HASH_POLICY,'raw_sha256':raw_digest,
                          'retained_variable_observations':observations,
                          'expected_sha256':digest,'actual_sha256':measured,'equal':digest==measured}
        else:
            measured=raw_digest
            checks[name]={'comparison_policy':'RAW_SCIENTIFIC_FILE_SHA256',
                          'expected_sha256':digest,'actual_sha256':measured,'equal':digest==measured}
    require(all(c['equal'] for c in checks.values()),'Retained scientific output bytes differ from completed reference')
    return dict(status='VERIFIED_IDENTICAL_COMPLETED_REPLAY',classification=replay['classification'],
                old_metric_status='FAIL',scientific_checks=checks,
                limitation='Resources are observed separately; unchanged scientific hashes do not establish the full hypothesis.')

def verify_checkout(target,source):
    require(validate_target(target),'Cannot replay an uncomputed/pending target')
    current=subprocess.run(['git','rev-parse','HEAD'],cwd=source,check=True,capture_output=True,text=True).stdout.strip()
    require(current==target['science_commit'],'Checkout differs from immutable completed science commit')
    subprocess.run(['git','merge-base','--is-ancestor',target['freeze_commit'],'HEAD'],cwd=source,check=True)
    root=Path(source)/CHECKPOINT
    require(sha(safe_file(root,'requirements-replay.txt'))==target['requirements_sha256'],'Replay requirements changed')
    driver=safe_file(root,'code/replay_active_source.py')
    require(sha(driver)==target['replay_driver_sha256'],'Replay driver changed')
    require(sha(safe_file(root,'FULL_REGISTRATION.json'))==target['registration_sha256'],'Full registration changed')
    registration=load_json(safe_file(root,'FULL_REGISTRATION.json'));files=registration.get('files')
    require(registration.get('schema_version')==1 and isinstance(files,dict) and 0<len(files)<=1000,
            'Unexpected full registration schema or size')
    require({'EXPERIMENT.json','requirements-replay.txt','code/replay_active_source.py',
             'code/verify_public_freeze.py'}<=set(files),'Required replay controls are not frozen')
    for name,digest in files.items():
        require(hex_string(digest,64) and sha(safe_file(root,name))==digest,'Frozen source bytes changed')
    permitted_data={'.json','.jsonl','.log','.txt','.md','.csv','.tsv','.npy','.npz','.png','.svg','.pdf'}
    for path in root.rglob('*'):
        require(not path.is_symlink(),'Checkpoint contains a symlink')
        if not path.is_file():continue
        name=path.relative_to(root).as_posix()
        if name in files or name in {'FULL_REGISTRATION.json','FREEZE_RECEIPT.json','MANIFEST.json'}:continue
        require(name.split('/')[0] in {'outputs','reports','figures','evidence'} and path.suffix in permitted_data,
                'Unregistered source or unsupported postfreeze file')
    return root,driver

def replay(target,source,output):
    root,driver=verify_checkout(target,source)
    args=['--checkpoint-root',str(root),'--registration-sha256',target['registration_sha256'],
          '--freeze-commit',target['freeze_commit']]
    with (output/'reproduction.stdout.log').open('xb') as stdout,(output/'reproduction.stderr.log').open('xb') as stderr:
        verifier=safe_file(root,'code/verify_public_freeze.py')
        subprocess.run([sys.executable,str(verifier),*args,'--repository-root',str(source)],check=True,stdout=stdout,stderr=stderr,timeout=300)
        # The frozen driver enforces the entire-route budgets and normal Python.
        subprocess.run([sys.executable,str(driver),*args,'--output-dir',str(output/'complete-replay')],check=True,
                       stdout=stdout,stderr=stderr,timeout=3900)
    result=verify_completed_replay(target,output/'complete-replay')
    write_json(output/'REPRODUCTION_COMPARISON.json',result)
    return result

def write_summary(path,state):
    lines=['### HDBLAST bounded mechanical follow-up','',
           'Status: `'+state['status']+'`.','',
           'This job checks a completed checkpoint or bibliographic metadata. It does not conduct autonomous AI research or publish scientific claims.']
    for label,item in state.get('queries',{}).items():
        diff=item.get('diff',{})
        lines.append('- '+label+': '+str(item.get('status'))+'; added '+str(len(diff.get('added',[])))+', modified '+str(len(diff.get('modified',[])))+', absent from current first page '+str(len(diff.get('absent_from_current_first_page',[])))+'.')
    if state.get('classification') in CLASSIFICATIONS:
        lines.extend(['','Reproduced classification: `'+state['classification']+'`; historical metric status remains `FAIL`.'])
    lines.extend(['','Inspect the artifact for complete receipts/diffs. Missing first-page records are not retractions; metadata changes are not evidence for HDBLAST.',''])
    with Path(path).open('a') as f:f.write('\n'.join(lines))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['gate','watch','verify-source','replay'])
    p.add_argument('--target',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--baseline',type=Path);p.add_argument('--baseline-sha256');p.add_argument('--source',type=Path)
    p.add_argument('--github-output',type=Path);p.add_argument('--summary',type=Path);args=p.parse_args()
    output=fresh_output(args.output);state={'status':'FAILED_MONITOR','started_utc':stamp()};code=1
    try:
        target=load_json(args.target)
        if args.mode!='gate':
            # Recheck when each actual job starts; a long runner queue must not
            # extend the reviewed campaign beyond its declared end.
            campaign_ok(os.environ.get('MONITOR_START_UTC',''),os.environ.get('MONITOR_UNTIL_UTC',''),
                        datetime.now(timezone.utc))
        if args.mode=='gate':
            context={k:os.environ.get(n,'') for k,n in dict(repository='GITHUB_REPOSITORY',branch='GITHUB_REF_NAME',default_branch='MONITOR_DEFAULT_BRANCH',attempt='GITHUB_RUN_ATTEMPT',event='GITHUB_EVENT_NAME',enabled='MONITOR_ENABLED',start='MONITOR_START_UTC',end='MONITOR_UNTIL_UTC',mode='MONITOR_MODE').items()}
            state=gate(target,context,datetime.now(timezone.utc));code=0
            if args.github_output:
                values={k:str(state.get(k,False)).lower() for k in ('ready','do_watch','do_replay')}
                values.update(science_commit=state.get('science_commit',''),checkpoint_relative_path=state.get('checkpoint_relative_path',''))
                with args.github_output.open('a') as f:
                    for name,value in values.items():f.write(name+'='+value+'\n')
        elif args.mode=='watch':
            require(validate_target(target),'Cannot activate monitor before completed target review')
            require(args.baseline and hex_string(args.baseline_sha256,64),'Reviewed baseline pin required')
            require(sha(args.baseline)==args.baseline_sha256,'Reviewed baseline bytes changed')
            state=watch(load_json(args.baseline),output);code=0 if state['status']=='COMPLETED_METADATA_READ' else 1
        elif args.mode=='verify-source':
            require(args.source is not None,'Science source required')
            root,_=verify_checkout(target,args.source.resolve())
            state={'status':'VERIFIED_COMPLETED_SOURCE_BYTES','science_commit':target['science_commit'],
                   'registration_sha256':target['registration_sha256'],'checkpoint_relative_path':CHECKPOINT,
                   'registered_files':len(load_json(root/'FULL_REGISTRATION.json')['files']),
                   'physical_science_evaluated':False};code=0
        else:
            require(args.source is not None,'Science source required')
            state=replay(target,args.source.resolve(),output);code=0
    except Exception as error:
        state.update(status='FAILED_MONITOR',error_type=type(error).__name__)
        # Known local guard errors are safe fixed messages; transport bodies/URLs are not emitted.
        if isinstance(error,ValueError):state['guard_message']=str(error)
    state['completed_utc']=stamp();write_json(output/'MONITOR.json',state)
    if args.summary:write_summary(args.summary,state)
    print(json.dumps({'status':state['status']}));return code

if __name__=='__main__':raise SystemExit(main())
