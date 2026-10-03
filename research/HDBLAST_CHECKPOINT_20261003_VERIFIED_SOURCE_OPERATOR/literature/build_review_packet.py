"""Build literature metadata from retrieved public records; no numerical experiment."""
import datetime
import hashlib
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent
PRIOR = pathlib.Path("/workspace/HDblast/research/HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER/literature/BIBLIOGRAPHY.json")

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write(name, content):
    (ROOT / name).write_text(json.dumps(content, indent=2, ensure_ascii=False) + "\n")

def main():
    search = json.loads((ROOT / "SEARCH_LEDGER.json").read_text())
    follow = json.loads((ROOT / "FOLLOWUP_METADATA.json").read_text())
    targeted = json.loads((ROOT / "TARGETED_FOLLOWUP_SEARCH.json").read_text())
    access = json.loads((ROOT / "PRIMARY_ACCESS_PROVENANCE.json").read_text())
    candidates = {}
    for q in search["queries"] + follow["records"] + targeted["records"]:
        for rec in q.get("records", []):
            candidates[rec["id"].rsplit("/", 1)[-1]] = rec
    doi_metadata = json.loads((ROOT / "DOI_METADATA.json").read_text())
    confirmed_dois = {}
    for q in doi_metadata["records"]:
        for rec in q.get("record", {}).get("message", {}).get("items", []):
            if rec.get("author") != [{"given": "Fredrik", "family": "Johansson", "sequence": "first", "affiliation": [], "role": [{"vocabulary": "crossref", "role": "author"}]}]:
                continue
            title = rec.get("title", [""])[0].lower()
            if title == "numerical integration in arbitrary-precision ball arithmetic":
                confirmed_dois["1802.07942v1"] = rec["DOI"]
            elif title == "arb: efficient arbitrary-precision midpoint-radius interval arithmetic":
                confirmed_dois["1611.02831v1"] = rec["DOI"]
    expected_dois = {"1802.07942v1": "10.1007/978-3-319-96418-8_30", "1611.02831v1": "10.1109/tc.2017.2690633"}
    if confirmed_dois != expected_dois:
        raise ValueError("title/author-confirmed DOI set differs from reviewed citation set")
    reviewed = {
        "1611.02831v1": ("johansson_arb_1611_02831v1.pdf", "Sec.2.2, printed pp.2–3", "Ball arithmetic and actual-radius precision control; not a validated trajectory."),
        "1802.07942v1": ("johansson_integration_1802_07942v1.pdf", "Secs.2, 2.1, 3.2, 3.3, printed pp.2–3 and 6–7", "Callback and complex-domain premises, tolerance limits, separate tails."),
        "2607.18180v1": ("mottola_ward_2607_18180v1.pdf", "Selected introduction; Sec.II, Eqs.(2.14)–(2.21), printed pp.7–8", "Full metric-variation contact terms; flat-space vacuum scope."),
        "2609.10673v1": ("borinsky_collider_2609_10673v1.pdf", "Secs.3.1–3.2, Eqs.(43)–(62), printed pp.12–14", "Analytic elimination of oscillatory scale; numerical validation remains empirical."),
        "2609.37346v1": ("chattopadhyay_brane_2609_37346v1.pdf", "Selected introduction; Sec.8.1, Eqs.(8.4)–(8.22), printed pp.23–25", "Conditional trace coupling and branching/abundance/thermalization premises in a different thick-brane model."),
    }
    old = json.loads(PRIOR.read_text())
    old_ids = {r.get("arxiv_version", r.get("arxiv")) for r in old["records"]}
    rows = []
    access_by_name = {r["filename"]: r for r in access["retrievals"]}
    for key, (filename, passages, relevance) in reviewed.items():
        rec = candidates[key]
        a = access_by_name[filename]
        rows.append({"arxiv_version": key, "title": rec["title"], "authors": rec["authors"], "first_submitted_at_utc": rec["published"], "version_date_at_utc": rec["updated"], "doi": confirmed_dois.get(key, rec.get("doi")), "primary_url": f"https://arxiv.org/abs/{key}", "primary_pdf_url": a["url"], "retrieved_at_utc": a["retrieved_at_utc"], "pdf_bytes": a["bytes"], "pdf_sha256": a["sha256"], "primary_read_status": "SELECTED_PASSAGES_READ_NOT_WHOLE_PAPER_VERIFICATION", "inspected_passages": passages, "supported_relevance_and_limit": relevance, "present_in_inspected_prior_bibliography": key in old_ids})
    inherited = []
    for key in ("2607.13580v1", "2503.08169v3", "2404.11448v2"):
        rec = candidates[key]
        inherited.append({"arxiv_version": key, "title": rec["title"], "authors": rec["authors"], "doi": rec.get("doi"), "primary_url": f"https://arxiv.org/abs/{key}", "reading_this_round": "CURRENT_METADATA_CHECK_ONLY_PRIOR_SELECTED_PRIMARY_REVIEW_RETAINED", "present_in_inspected_prior_bibliography": key in old_ids})
    docs = []
    for filename in ("flint_v3_6_0_arb.rst", "flint_v3_6_0_acb_calc.rst"):
        a = access_by_name[filename]
        docs.append({"title": filename, "release_tag": access["release_tag"], "release_commit": access["release_commit"], "url": a["url"], "sha256": a["sha256"], "bytes": a["bytes"], "retrieved_at_utc": a["retrieved_at_utc"], "read_status": "SELECTED_CONTRACT_AND_INTEGRATION_SECTIONS_READ_NOT_WHOLE_LIBRARY_VERIFICATION"})
    write("BIBLIOGRAPHY.json", {"schema_version": 1, "client_cutoff_date": "2026-10-02 America/Los_Angeles", "inspected_prior_bibliography": str(PRIOR), "inspected_prior_bibliography_sha256": digest(PRIOR), "selected_primary_readings": rows, "selected_upstream_documentation_readings": docs, "inherited_metadata_confirmations": inherited, "retrieved_without_selected_reading": ["raw/flint_v3_6_0_acb.rst", "raw/flint_v3_6_0_acb_calc_integrate.c"], "identifier_correction_disclosure": {"recalled_arxiv_ids_excluded": ["1610.01503", "1802.07946"], "recalled_doi_excluded_http404": "10.1145/3095149", "policy": "Only retrieved title/author-confirmed IDs/DOIs appear in selected bibliography."}})
    query_rows = search["queries"][:5] + targeted["records"]
    write("REVIEW_LIMITS.json", {"schema_version": 1, "status": "BOUNDED_SELECTED_PRIMARY_REVIEW_NO_SCIENTIFIC_EVALUATION", "first_page_searches": 7, "maximum_search_results_including_duplicates": 75, "returned_search_results_including_duplicates": sum(len(r.get("records", [])) for r in query_rows), "unique_search_record_ids": len({r["id"] for q in query_rows for r in q.get("records", [])}), "selected_preprint_readings": 5, "selected_documentation_readings": 2, "arxiv_api_access_success": all(q.get("status") == 200 for q in query_rows), "whole_web_coverage": False, "whole_paper_verification": False, "external_peer_review_of_this_packet": False, "external_novelty_assessment": "NOT_ASSESSED", "new_fundamental_math": "NOT_CLAIMED", "source_callback_execution": False, "project_array_decode": False, "scientific_experiment": False, "remote_mutation": False, "paper_pdf_text_extraction_limit": "pdftotext -layout emitted font warnings; selected passages locate support but every equation was not independently transcribed.", "library_qualification_limit": "Release-pinned C docs do not establish the runtime version/guarantees of a separately installed Python binding.", "public_payload_recommendation": "Review plus bounded metadata/provenance; external raw third-party document cache is not a repository payload."})
    package_files = ["PRIMARY_LITERATURE_REVIEW.md", "BIBLIOGRAPHY.json", "REVIEW_LIMITS.json", "SELECTED_PASSAGES.json", "SEARCH_LEDGER.json", "TARGETED_FOLLOWUP_SEARCH.json", "FOLLOWUP_METADATA.json", "DOI_METADATA.json", "PRIMARY_ACCESS_PROVENANCE.json", "TEXT_EXTRACTION_PROVENANCE.json", "fetch_metadata.py", "fetch_primary.py", "build_selected_passages.py", "build_review_packet.py"]
    for name in package_files:
        if not (ROOT / name).is_file():
            raise ValueError(f"missing staged public-provenance candidate: {name}")
    write("MANIFEST.json", {"schema_version": 1, "scope": "Literature/provenance candidate only; third-party raw cache excluded", "files": {name: {"bytes": (ROOT/name).stat().st_size, "sha256": digest(ROOT/name)} for name in package_files}})
    print(json.dumps({"status": "PASS_LOCAL_SELECTED_PRIMARY_LITERATURE_PACKET", "path": str(ROOT), "review_sha256": digest(ROOT / "PRIMARY_LITERATURE_REVIEW.md"), "bibliography_sha256": digest(ROOT / "BIBLIOGRAPHY.json"), "manifest_sha256": digest(ROOT / "MANIFEST.json"), "selected_papers": len(rows), "selected_docs": len(docs), "inherited_metadata_only": len(inherited), "new_physical_evaluations": 0}, indent=2))

if __name__ == "__main__":
    main()
