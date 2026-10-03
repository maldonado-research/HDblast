"""Bounded public metadata retrieval; no project numerical inputs are loaded."""
import concurrent.futures
import datetime
import hashlib
import json
import pathlib
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = pathlib.Path(__file__).resolve().parent
NOW = lambda: datetime.datetime.now(datetime.timezone.utc).isoformat()
QUERIES = [
    ("ball_integration", 'ti:"integration" AND all:"ball arithmetic"', 10),
    ("recent_verified_oscillatory", 'submittedDate:[202401010000 TO 202610022359] AND all:"oscillatory" AND (all:"validated" OR all:"rigorous" OR all:"certified")', 15),
    ("recent_taylor_integration", 'submittedDate:[202401010000 TO 202610022359] AND all:"Taylor models" AND all:"integration"', 10),
    ("recent_brane_reheating", 'submittedDate:[202501010000 TO 202610022359] AND all:"brane" AND (all:"reheating" OR all:"particle production")', 10),
    ("recent_ward_scalar", 'submittedDate:[202605010000 TO 202610022359] AND all:"Ward" AND all:"scalar" AND (all:"gravity" OR all:"stress")', 10),
]

def fetch(label, url, kind="atom"):
    entry = {"label": label, "url": url, "requested_at_utc": NOW()}
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "HDBLAST bounded primary-literature provenance check", "Accept": "application/atom+xml, application/json"})
        with urllib.request.urlopen(req, timeout=50) as r:
            body = r.read()
            entry.update(status=r.status, final_url=r.geturl(), content_type=r.headers.get("Content-Type"), response_date=r.headers.get("Date"))
        suffix = ".xml" if kind == "atom" else ".json"
        name = f"raw/{label}{suffix}"
        (ROOT / name).write_bytes(body)
        entry.update(retrieved_at_utc=NOW(), bytes=len(body), sha256=hashlib.sha256(body).hexdigest(), saved_path=name)
        if kind == "atom":
            ns = {"a": "http://www.w3.org/2005/Atom", "o": "http://a9.com/-/spec/opensearch/1.1/", "x": "http://arxiv.org/schemas/atom"}
            xml = ET.fromstring(body)
            entry["total_hits"] = xml.findtext("o:totalResults", namespaces=ns)
            records = []
            for e in xml.findall("a:entry", ns):
                records.append({"id": e.findtext("a:id", namespaces=ns), "title": " ".join(e.findtext("a:title", default="", namespaces=ns).split()), "authors": [x.findtext("a:name", namespaces=ns) for x in e.findall("a:author", ns)], "published": e.findtext("a:published", namespaces=ns), "updated": e.findtext("a:updated", namespaces=ns), "doi": e.findtext("x:doi", namespaces=ns), "journal_ref": e.findtext("x:journal_ref", namespaces=ns), "summary": " ".join(e.findtext("a:summary", default="", namespaces=ns).split()), "links": [x.attrib for x in e.findall("a:link", ns)]})
            entry["records"] = records
        else:
            entry["record"] = json.loads(body)
    except Exception as e:
        entry.update(error=f"{type(e).__name__}: {e}", retrieved_at_utc=NOW())
    return entry

def main():
    urls = []
    for label, query, count in QUERIES:
        url = "https://export.arxiv.org/api/query?" + urllib.parse.urlencode({"search_query": query, "start": 0, "max_results": count, "sortBy": "submittedDate", "sortOrder": "descending"})
        urls.append((label, url, "atom"))
    urls.extend([
        ("selected_method_ids", "https://export.arxiv.org/api/query?" + urllib.parse.urlencode({"id_list": "1610.01503,1802.07946,2404.11448,2503.08169,2607.13580", "max_results": 10}), "atom"),
        ("flint_latest_release", "https://api.github.com/repos/flintlib/flint/releases/latest", "json"),
        ("flint_main_ref", "https://api.github.com/repos/flintlib/flint/git/ref/heads/main", "json"),
    ])
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        entries = list(pool.map(lambda args: fetch(*args), urls))
    ledger = {"schema_version": 1, "scope": "Five bounded first-page arXiv metadata searches (maximum 55 hits including duplicates), one selected-ID check, and two FLINT public release/ref checks; selected primary documents are a separate reading stage. No entire-web coverage claim.", "client_cutoff_date": "2026-10-02 America/Los_Angeles", "completed_at_utc": NOW(), "queries": entries}
    (ROOT / "SEARCH_LEDGER.json").write_text(json.dumps(ledger, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps([{ "label": e["label"], "status": e.get("status"), "error": e.get("error"), "total_hits": e.get("total_hits"), "titles": [{"id": r["id"], "title": r["title"]} for r in e.get("records", [])], "tag_name": e.get("record", {}).get("tag_name"), "sha": e.get("record", {}).get("object", {}).get("sha") } for e in entries], indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
