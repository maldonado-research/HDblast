#!/usr/bin/env python3
"""Build the editable four-page working note and record exact selected hashes."""
import datetime,hashlib,json,os,re,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
source=ROOT/'INPUT_PRECISION.tex';initial=ROOT/'INITIAL_RESULT.tex'
if 'Draft placeholder' in initial.read_text():
 raise SystemExit('Initial result is not finalized; refusing final build receipt.')
env=os.environ.copy();env['SOURCE_DATE_EPOCH']=str(int(datetime.datetime(2026,10,23,tzinfo=datetime.timezone.utc).timestamp()));env['FORCE_SOURCE_DATE']='1';env['TZ']='UTC'
commands=[]
for n in range(2):
 cmd=['pdflatex','-interaction=nonstopmode','-halt-on-error','INPUT_PRECISION.tex']
 p=subprocess.run(cmd,cwd=ROOT,env=env,capture_output=True,text=True)
 (ROOT/f'build-pass{n+1}.txt').write_text(p.stdout+p.stderr)
 commands.append({'argv':cmd,'returncode':p.returncode})
 if p.returncode:raise SystemExit('pdflatex failed')
log=(ROOT/'INPUT_PRECISION.log').read_text();bad=[x for x in log.splitlines() if re.search(r'Overfull|undefined|^!|Rerun to get',x)]
if bad:raise SystemExit('Unresolved final TeX diagnostics: '+repr(bad))
info=subprocess.run(['pdfinfo','INPUT_PRECISION.pdf'],cwd=ROOT,capture_output=True,text=True,check=True).stdout
pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M).group(1))
if pages not in (3,4):raise SystemExit('Expected3–4 pages, got '+str(pages))
def record(name):
 b=(ROOT/name).read_bytes();return {'path':name,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
ver=subprocess.run(['pdflatex','--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0]
receipt={'status':'PASS_EDITABLE_TEX_PDF_BUILD','title':'A sharper regular-shell enclosure and correlated initial-data errors','author':'Ricardo Maldonado','affiliation':'Independent Researcher, Henderson, Nevada, United States','commands':commands,'compiler':ver,'source_date_epoch':env['SOURCE_DATE_EPOCH'],'pages':pages,'overfull_boxes':0,'undefined_references':0,'final_log_error_lines':bad,'files':[record(s) for s in ['INPUT_PRECISION.tex','INITIAL_RESULT.tex','INPUT_PRECISION.pdf','INPUT_PRECISION.log','build_note.py']],'scope':'New conditional working note. Frozen older manuscripts are not modified. No journal submission.'}
(ROOT/'BUILD_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
