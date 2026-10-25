#!/usr/bin/env python3
"""Post-publication source/ZIP and live-page readback on public owned URLs."""
import argparse
import concurrent.futures
import hashlib
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

def need(ok,message):
    if not ok:
        raise ValueError(message)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def fetch(url):
    request=urllib.request.Request(url,headers={"User-Agent":"HDBLAST-publication-check/1.0",
                                               "Cache-Control":"no-cache"})
    with urllib.request.urlopen(request,timeout=20) as response:
        need(response.status==200,"HTTP200")
        return response.read()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root",type=Path,required=True)
    p.add_argument("--pins",type=Path,required=True)
    p.add_argument("--pins-sha256",required=True)
    p.add_argument("--research-sha",required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    need(not a.output.exists(),"fresh publication readback")
    need(len(a.research_sha)==40 and all(c in "0123456789abcdef" for c in a.research_sha),
         "exact source commit")
    root=a.repo_root.resolve()
    pinraw=a.pins.read_bytes()
    need(sha(pinraw)==a.pins_sha256,"external pins")
    pins=json.loads(pinraw)
    evidence=a.pins.parent
    additions_raw=(evidence/"SITE_ADDITIONS.json").read_bytes()
    need(sha(additions_raw)==pins["site_additions_sha256"],"site insertion pins")
    additions=json.loads(additions_raw)["files"]
    capsule=root/pins["capsule"]["path"]
    manifest_raw=(capsule/"MANIFEST.json").read_bytes()
    need(sha(manifest_raw)==pins["capsule"]["manifest_sha256"],"capsule manifest")
    members=json.loads(manifest_raw)["files"]
    jobs=[]
    for entry in members+[{"path":"MANIFEST.json","sha256":sha(manifest_raw),
                           "bytes":len(manifest_raw)}]:
        path=pins["capsule"]["path"]+"/"+entry["path"]
        url="https://raw.githubusercontent.com/maldonado-research/HDblast/"+a.research_sha+"/"+urllib.parse.quote(path,safe="/")
        jobs.append((path,url,entry["sha256"],entry["bytes"]))
    archive=pins["archive"]
    jobs.append((archive["path"],"https://raw.githubusercontent.com/maldonado-research/HDblast/"+a.research_sha+"/"+archive["path"],archive["sha256"],archive["bytes"]))
    for row in additions:
        repository=row["repository"]
        need(repository in ("maldonado-research/HDblast",
                            "maldonado-research/maldonado-research.github.io"),"owned repository")
        ref=a.research_sha if repository=="maldonado-research/HDblast" else "main"
        url="https://raw.githubusercontent.com/"+repository+"/"+ref+"/"+row["path"]
        jobs.append((repository+":"+row["path"],url,row["new_sha256"],None))
    raw_reports=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        futures=[pool.submit(fetch,j[1]) for j in jobs]
        for (path,url,expected,size),future in zip(jobs,futures):
            raw=future.result()
            need(sha(raw)==expected and (size is None or len(raw)==size),
                 "published source bytes:"+path)
            raw_reports.append({"path":path,"url":url,"sha256":sha(raw),"bytes":len(raw)})
    live_urls={
        "site-en":"https://maldonado-research.github.io/HDblast/",
        "site-es":"https://maldonado-research.github.io/HDblast/es/",
        "shared-site":"https://maldonado-research.github.io/projects/hdblast/",
    }
    live_reports=[]
    attempts=[]
    for attempt in range(1,7):
        found=[]
        errors=[]
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            work=[(row,pool.submit(fetch,live_urls[row["label"]]+"?checkpoint="+a.research_sha[:12]))
                  for row in additions if row["label"] in live_urls]
            for row,future in work:
                try:
                    raw=future.result()
                    need(sha(raw)==row["new_sha256"],"whole live HTML hash")
                    need(raw.count(row["block"].strip().encode())==1,"one complete live insertion")
                    for marker,count in row["prior_marker_counts"].items():
                        need(raw.count(marker.encode())==count,"prior live marker")
                    found.append({"url":live_urls[row["label"]],"sha256":sha(raw),
                                  "bytes":len(raw),"whole_html_equals_reviewed_source":True,
                                  "prior_marker_types":len(row["prior_marker_counts"])})
                except Exception as exc:
                    errors.append({"page":row["label"],"error":str(exc)})
        attempts.append({"attempt":attempt,"errors":errors})
        if not errors:
            live_reports=found
            break
        if attempt<6:
            time.sleep(10)
    need(len(live_reports)==3,"all live pages must match after bounded deployment wait")
    result={"status":"PASS_PUBLISHED_SOURCE_ARCHIVE_AND_LIVE_PAGES",
            "research_source_commit":a.research_sha,"raw_file_count":len(raw_reports),
            "raw_files":raw_reports,"live_pages":live_reports,"bounded_live_attempts":attempts,
            "scope":"Selected public source files including exact ZIP bytes, all four edited source files, three complete live HTML pages and retained markers. Canonical merge/PR identities are verified separately."}
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("HDB_LIVE_READBACK:"+json.dumps(result,sort_keys=True,separators=(",",":")),flush=True)

if __name__=="__main__":
    main()
