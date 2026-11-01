"""Standard-library authentication and non-overwriting portable output policy."""
from pathlib import Path
import argparse,hashlib,json,sys
CAP='HDBLAST_PREPARED_SHELL_VOLUME_20261101'
HERE=Path(__file__).resolve().parent
def fail(stage,reason):
    print(json.dumps({'status':'REJECT','stage':stage,'reason':reason},sort_keys=True))
    raise SystemExit(2)
class Parser(argparse.ArgumentParser):
    def error(self,message):fail('ARGUMENT','INVALID')
def arguments(mutants=()):
    p=Parser()
    p.add_argument('--repo-root',type=Path)
    p.add_argument('--initial-report',type=Path)
    p.add_argument('--conditioned-report',type=Path)
    p.add_argument('--response-report',type=Path)
    p.add_argument('--output',default='-')
    p.add_argument('--mutant',choices=list(mutants))
    a=p.parse_args();a.repo_root=(a.repo_root or (HERE.parents[2] if len(HERE.parents)>2 else Path.cwd())).resolve()
    cap=a.repo_root/'research'/CAP
    a.initial_report=(a.initial_report or cap/'initial/references/RAW_BOXES.json').resolve()
    a.conditioned_report=(a.conditioned_report or cap/'initial/references/CONDITIONED_RAW_BOXES.json').resolve()
    a.response_report=(a.response_report or cap/'response/references/FINAL_NORMAL.json').resolve()
    if a.output!='-':
        a.output=Path(a.output).resolve()
        if a.output.exists():fail('OUTPUT','EXISTS')
        for root in (a.repo_root,HERE):
            if a.output.is_relative_to(root):fail('OUTPUT','PROTECTED_TREE')
        for source in (a.initial_report,a.conditioned_report,a.response_report):
            if a.output==source:fail('OUTPUT','SOURCE_COLLISION')
        if not a.output.parent.is_dir():fail('OUTPUT','PARENT_MISSING')
    return a
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def authenticate(a):
    try:
        pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
        for row in pins['repository_files']:
            p=a.repo_root/row['path']
            if not p.is_file() or p.is_symlink() or p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:
                fail('AUTH_SOURCE',row['id'])
        for row in pins['reports']:
            p=getattr(a,row['argument'])
            if not p.is_file() or p.is_symlink() or p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:
                fail('AUTH_REPORT',row['id'])
        local=json.loads((HERE/'LOCAL_PINS.json').read_text())
        for row in local['files']:
            p=HERE/row['path']
            if not p.is_file() or p.is_symlink() or p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:
                fail('AUTH_LOCAL',row['path'])
    except (OSError,KeyError,ValueError) as error:fail('AUTH_SCHEMA',type(error).__name__)
    a._authenticated=True
def finish(a,result):
    authenticate(a)
    data=json.dumps(result,sort_keys=True,indent=2)+'\n'
    if a.output=='-':sys.stdout.write(data)
    else:
        with a.output.open('x') as f:f.write(data)
