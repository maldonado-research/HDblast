#!/usr/bin/env python3
"""Replay an externally authenticated, explicit positive/rejection roster."""
import argparse
from decimal import Decimal
import hashlib
import json
from pathlib import Path,PurePosixPath
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parent


def require(ok,message):
    if not ok:raise ValueError(message)


def digest(raw):return hashlib.sha256(raw).hexdigest()


def canonical(name):
    require(isinstance(name,str),'Member path must be a string')
    p=PurePosixPath(name)
    require(bool(name) and name!='.' and '\\' not in name and not p.is_absolute()
            and '..' not in p.parts and p.as_posix()==name,'Noncanonical member')
    return p


def complete_json(raw):
    def pairs(items):
        result={}
        for key,value in items:
            require(key not in result,'Duplicate JSON key')
            result[key]=value
        return result
    def reject(value):raise ValueError('Nonfinite JSON constant: '+value)
    return json.loads(raw,object_pairs_hook=pairs,parse_float=Decimal,parse_constant=reject)


def identical(a,b):
    if type(a) is not type(b):return False
    if isinstance(a,dict):return a.keys()==b.keys() and all(identical(a[k],b[k]) for k in a)
    if isinstance(a,list):return len(a)==len(b) and all(identical(x,y) for x,y in zip(a,b))
    return a==b


def verify(expected):
    require(isinstance(expected,str) and re.fullmatch('[0-9a-f]{64}',expected)
            and expected!='0'*64,'Concrete external manifest pin required')
    raw=(ROOT/'MANIFEST.json').read_bytes()
    require(digest(raw)==expected,'External manifest hash differs')
    manifest=complete_json(raw);listed=set()
    for item in manifest['files']:
        name=item['path'];canonical(name)
        require(name not in listed and name!='MANIFEST.json','Duplicate/reserved member')
        source=ROOT/name
        require(source.is_file() and not any(p.is_symlink() for p in (source,*source.parents)),'Nonregular source member')
        body=source.read_bytes()
        require(type(item['bytes']) is int and len(body)==item['bytes']
                and digest(body)==item['sha256'],'Source differs: '+name)
        listed.add(name)
    walked=list(ROOT.rglob('*'))
    require(all(not p.is_symlink() and (p.is_file() or p.is_dir()) for p in walked),'Special source inventory entry')
    actual={p.relative_to(ROOT).as_posix() for p in walked if p.is_file()}
    require(actual==listed|{'MANIFEST.json'},'Source inventory differs')
    require({'REPLAY_ROSTER.json','reproduce.py'}<=listed,'Runner and roster must be authenticated')
    return listed


def register(roster,listed):
    require(isinstance(roster.get('scope'),str) and bool(roster['scope'].strip()),'Explicit roster scope required')
    labels=set();expanded={'positive':[],'rejection':[]}
    for kind,key in (('positive','jobs'),('rejection','rejections')):
        rows=roster.get(key,[])
        require(isinstance(rows,list),'Roster rows must be a list')
        for job in rows:
            require(isinstance(job.get('name'),str) and re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*',job['name']), 'Safe job name required')
            canonical(job['script']);require(job['script'] in listed,'Job script outside authenticated source')
            arguments=job.get('arguments',[])
            require(isinstance(arguments,list) and all(isinstance(v,str) for v in arguments),'Job arguments must be strings')
            require(not any('{output}' in v and v!='{output}' for v in arguments),'Output token must be a whole argument')
            if '{output}' in arguments:
                require(arguments.count('{output}')==1,'Exactly one output token required')
            else:
                require(not any(v=='--output' or v.startswith('--output=') for v in arguments),'Runner owns the output argument')
            modes=job.get('modes',['normal','optimized'])
            require(isinstance(modes,list) and bool(modes) and len(set(modes))==len(modes)
                    and all(v in ('normal','optimized') for v in modes),'Valid unique job modes required')
            timeout=job.get('timeout_seconds',240)
            require(type(timeout) is int and timeout>0,'Positive job timeout required')
            if kind=='rejection':
                require(isinstance(job.get('expected_reason'),str) and bool(job['expected_reason']),'Specific rejection reason required')
            for mode in modes:
                label=job['name']+'.'+mode
                require(label not in labels,'Duplicate expanded job label');labels.add(label)
                golden=None
                if kind=='positive':
                    golden=job['goldens'][mode];canonical(golden)
                    require(golden in listed,'Golden outside authenticated source')
                    complete_json((ROOT/golden).read_bytes())
                expanded[kind].append({'label':label,'mode':mode,'job':job,'golden':golden,'timeout':timeout})
    require(bool(expanded['positive']),'Nonempty positive replay roster required')
    return expanded


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest-sha256',required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args()
    listed=verify(args.manifest_sha256)
    roster_raw=(ROOT/'REPLAY_ROSTER.json').read_bytes();roster=complete_json(roster_raw)
    expanded=register(roster,listed)
    out=args.output_dir.resolve()
    require(out!=ROOT and ROOT not in out.parents and not out.exists(),'New replay output outside immutable source required')
    out.mkdir(parents=True,exist_ok=False);records=[];rejections=[]
    for kind in ('positive','rejection'):
        for row in expanded[kind]:
            label=row['label'];job=row['job'];report=out/(label+'.json')
            flags=['-O'] if row['mode']=='optimized' else []
            arguments=job.get('arguments',[])
            output_arguments=([str(report) if v=='{output}' else v for v in arguments]
                              if '{output}' in arguments else [*arguments,'--output',str(report)])
            command=[sys.executable,'-I','-S',*flags,'-B',str(ROOT/job['script']),*output_arguments]
            try:
                run=subprocess.run(command,cwd=ROOT,capture_output=True,timeout=row['timeout'])
            except subprocess.TimeoutExpired as error:
                (out/(label+'.stdout')).write_bytes(error.stdout or b'')
                (out/(label+'.stderr')).write_bytes(error.stderr or b'')
                raise ValueError('Checker timed out: '+label) from error
            (out/(label+'.stdout')).write_bytes(run.stdout);(out/(label+'.stderr')).write_bytes(run.stderr)
            if kind=='positive':
                require(run.returncode==0,'Checker failed: '+label)
                require(report.is_file() and not report.is_symlink(),'Checker produced no regular report: '+label)
                actual_raw=report.read_bytes();golden_raw=(ROOT/row['golden']).read_bytes()
                require(identical(complete_json(actual_raw),complete_json(golden_raw)),'Complete report differs: '+label)
                records.append({'name':label,'returncode':run.returncode,'report_sha256':digest(actual_raw),
                                'golden':row['golden'],'golden_sha256':digest(golden_raw),
                                'complete_typed_JSON_matches':True,'raw_bytes_match':actual_raw==golden_raw})
            else:
                require(run.returncode!=0,'Rejection unexpectedly exited zero: '+label)
                require(job['expected_reason'].encode() in run.stderr,'Rejection reason differs: '+label)
                require(not report.exists() and not report.is_symlink(),'Rejection left output: '+label)
                rejections.append({'name':label,'returncode':run.returncode,'expected_reason':job['expected_reason'],
                                   'reason_matched':True,'no_output':True})
    verify(args.manifest_sha256)
    result={'status':'PASS_COMPLETE_REPORT_AGREEMENT','manifest_sha256':args.manifest_sha256,
            'roster_sha256':digest(roster_raw),'processes':len(records),'rejection_processes':len(rejections),
            'records':records,'rejections':rejections,'source_unchanged':True,'scope':roster['scope']}
    (out/'REPRODUCTION_RECORD.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'processes':len(records),'rejection_processes':len(rejections)},indent=2))


if __name__=='__main__':
    try:main()
    except (ValueError,KeyError,OSError) as error:sys.exit('REJECT: '+str(error))
