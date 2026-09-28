#!/usr/bin/env python3
"""Which sources in theory_sources.py were already in the HDBLAST programme's own records?

Scans (read-only) Markdown/text files of: D-Blast 3 folders 146-152 (Chats 9-14 and controls), the
DBlast 4 vector-store handoffs, and the 22 Sept 2026 checkpoint package. For each source it searches for
its arXiv identifiers (from URL and venue) and for its title (case-insensitive). Output: prior_mentions.json.
A miss means "not found in the scanned files", not "unknown to the author".
"""
import importlib.util, json, os, re, time
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ts', os.path.join(HERE, 'theory_sources.py'))
ts = importlib.util.module_from_spec(spec); spec.loader.exec_module(ts)
DB3 = "/home/user/unified-theory-maldonado/new-files/D-Blast 3"
ROOTS = [os.path.join(DB3, f"untitled folder {n}") for n in range(146, 153)] + [
    os.path.join(DB3, "VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store"),
    "/home/user/unified-theory-maldonado/new-files/latest-work/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922"]
t0 = time.time()
texts = {}
for r in ROOTS:
    for dp, _, fs in os.walk(r):
        for f in fs:
            if f.endswith(('.md', '.txt')):
                p = os.path.join(dp, f)
                try:
                    if os.path.getsize(p) < 200_000_000:
                        texts[p] = open(p, errors='ignore').read().lower()
                except OSError:
                    pass
def ids_of(s):
    out = set()
    blob = s['url'] + ' ' + s.get('venue', '')
    out.update(re.findall(r'\b\d{4}\.\d{4,5}\b', blob))
    out.update(re.findall(r'\b(?:hep-th|hep-ph|gr-qc|astro-ph)/\d{7}\b', blob))
    return sorted(out)
res = {}
for s in ts.SOURCES:
    keys = ids_of(s)
    title = re.sub(r'\s+', ' ', s['title'].split(';')[0].split(' (')[0]).lower().strip()
    hits = []
    for p, txt in texts.items():
        if any(k.lower() in txt for k in keys) or (len(title) > 18 and title in txt):
            hits.append(os.path.relpath(p, "/home/user/unified-theory-maldonado/new-files"))
    res[s['id']] = dict(identifiers=keys, title_key=title, n_files=len(hits), files=sorted(hits)[:6])
summary = dict(files_scanned=len(texts), seconds=round(time.time() - t0, 1),
               sources_with_prior_mention=sum(1 for v in res.values() if v['n_files'] > 0),
               sources_total=len(res), roots=ROOTS, results=res)
json.dump(summary, open(os.path.join(HERE, 'prior_mentions.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in summary.items() if k != 'results'}, indent=1))
print(' '.join(f"{k}:{v['n_files']}" for k, v in res.items()))
