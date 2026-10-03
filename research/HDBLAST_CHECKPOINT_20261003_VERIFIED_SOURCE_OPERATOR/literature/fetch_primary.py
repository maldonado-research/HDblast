"""Download selected public primary text and pin the successful response bytes."""
import concurrent.futures
import datetime
import hashlib
import json
import pathlib
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent
SHA = "8d5454b96761fafe4d5a9da76a369a602f500f49"
SOURCES = [
    ("johansson_arb_1611_02831v1.pdf", "https://arxiv.org/pdf/1611.02831v1", "primary preprint"),
    ("johansson_integration_1802_07942v1.pdf", "https://arxiv.org/pdf/1802.07942v1", "primary preprint"),
    ("borinsky_collider_2609_10673v1.pdf", "https://arxiv.org/pdf/2609.10673v1", "primary preprint"),
    ("mottola_ward_2607_18180v1.pdf", "https://arxiv.org/pdf/2607.18180v1", "primary preprint"),
    ("chattopadhyay_brane_2609_37346v1.pdf", "https://arxiv.org/pdf/2609.37346v1", "primary preprint"),
    ("flint_v3_6_0_acb_calc.rst", f"https://raw.githubusercontent.com/flintlib/flint/{SHA}/doc/source/acb_calc.rst", "upstream API documentation at released commit"),
    ("flint_v3_6_0_arb.rst", f"https://raw.githubusercontent.com/flintlib/flint/{SHA}/doc/source/arb.rst", "upstream API documentation at released commit"),
    ("flint_v3_6_0_acb.rst", f"https://raw.githubusercontent.com/flintlib/flint/{SHA}/doc/source/acb.rst", "upstream API documentation at released commit"),
    ("flint_v3_6_0_acb_calc_integrate.c", f"https://raw.githubusercontent.com/flintlib/flint/{SHA}/src/acb_calc/integrate.c", "upstream integration implementation at released commit"),
]
def fetch(args):
    name, url, kind = args
    e = {"filename": name, "url": url, "kind": kind, "requested_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "HDBLAST selected-primary-literature provenance check"})
        with urllib.request.urlopen(req, timeout=55) as r:
            body = r.read()
            e.update(status=r.status, final_url=r.geturl(), content_type=r.headers.get("Content-Type"), response_date=r.headers.get("Date"))
        (ROOT / "raw" / name).write_bytes(body)
        e.update(bytes=len(body), sha256=hashlib.sha256(body).hexdigest(), saved_path=f"raw/{name}")
        if name.endswith(".pdf") and not body.startswith(b"%PDF-"):
            raise ValueError("response is not a PDF")
    except Exception as err:
        e["error"] = f"{type(err).__name__}: {err}"
    e["retrieved_at_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    return e
def main():
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(fetch, SOURCES))
    doc = {"schema_version": 1, "release_tag": "v3.6.0", "release_commit": SHA, "scope": "Five selected preprints and four pinned upstream FLINT documentation/implementation files. Retrieval is not full-paper or code verification; reviewed passages are specified separately.", "retrievals": rows}
    (ROOT / "PRIMARY_ACCESS_PROVENANCE.json").write_text(json.dumps(doc, indent=2) + "\n")
    print(json.dumps([{k: r.get(k) for k in ("filename", "status", "bytes", "sha256", "error")} for r in rows], indent=2))
if __name__ == "__main__":
    main()
