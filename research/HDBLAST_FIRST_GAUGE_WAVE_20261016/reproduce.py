#!/usr/bin/env python3
"""Replay the explicit exact-check roster from an authenticated source capsule."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def verify(expected):
    raw = (ROOT/'MANIFEST.json').read_bytes()
    require(digest(raw) == expected, 'External manifest hash differs')
    inventory = json.loads(raw)['files']
    listed = set()
    for entry in inventory:
        name = entry['path']
        p = Path(name)
        require(not p.is_absolute() and '..' not in p.parts and name == p.as_posix(), 'Noncanonical member')
        require(name not in listed and name != 'MANIFEST.json', 'Duplicate/reserved member')
        listed.add(name)
        source = ROOT/p
        require(not any(q.is_symlink() for q in [source,*source.parents]), 'Symlink source')
        body = source.read_bytes()
        require(len(body) == entry['bytes'] and digest(body) == entry['sha256'], 'Source differs: '+name)
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()}
    require(actual == listed|{'MANIFEST.json'}, 'Source inventory differs')
    return digest(raw)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest-sha256', required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    args = p.parse_args()
    manifest = verify(args.manifest_sha256)
    out = args.output_dir.resolve()
    require(out != ROOT and ROOT not in out.parents, 'Write replay outside immutable source')
    out.mkdir(parents=True, exist_ok=False)
    roster = json.loads((ROOT/'REPLAY_ROSTER.json').read_text())
    records = []
    for job in roster['jobs']:
        for mode in job.get('modes',['normal','optimized']):
            require(mode in ('normal','optimized'), 'Unknown execution mode')
            flags = ['-O'] if mode == 'optimized' else []
            label = job['name']+'.'+mode
            report = out/(label+'.json')
            command = [sys.executable,'-I',*flags,'-B',str(ROOT/job['script']),*job.get('arguments',[]),'--output',str(report)]
            run = subprocess.run(command,cwd=ROOT,capture_output=True,timeout=job.get('timeout_seconds',240))
            (out/(label+'.stdout')).write_bytes(run.stdout)
            (out/(label+'.stderr')).write_bytes(run.stderr)
            require(run.returncode == 0, 'Checker failed: '+label)
            actual = json.loads(report.read_text())
            expected = json.loads((ROOT/job['goldens'][mode]).read_text())
            require(actual == expected, 'Complete report differs: '+label)
            records.append({'name':label,'returncode':run.returncode,'report_sha256':digest(report.read_bytes()),'golden':job['goldens'][mode],'complete_report_matches':True})
    verify(manifest)
    result = {'status':'PASS_COMPLETE_REPORT_AGREEMENT','manifest_sha256':manifest,'processes':len(records),'records':records,'source_unchanged':True,'scope':'Finite exact checks and recorded controls; conditional first gauge-wave construction and coordinate derivative estimates, no full second-order solution or external peer review'}
    (out/'REPRODUCTION_RECORD.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'processes':len(records)},indent=2))

if __name__ == '__main__':
    main()
