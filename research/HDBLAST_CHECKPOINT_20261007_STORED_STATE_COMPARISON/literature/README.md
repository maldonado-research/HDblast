# Literature refresh working evidence

Read [LITERATURE_REFRESH_AND_MODEL_TEST.md](LITERATURE_REFRESH_AND_MODEL_TEST.md) for the physical implications and limits. `PRIMARY_READING_EVIDENCE.json` records eight primary PDF selected-passage readings; it does not assert complete proof or pipeline verification. `COVERAGE_AND_LIMITS.json` records the bounded search size and distinguishes six newly screened items from two existing snippet citations upgraded to primary text.

Current runtime HTTPS responses for arXiv/export, INSPIRE and Crossref are HTTP200. Every network operation was a verified-TLS GET. GitHub authentication was used only in HTTPS headers to api.github.com; no credential values or headers were saved. No callable scholarly browser/search connector appeared in the tool catalog. No remote mutation or physical callback/quantum-array evaluation occurred.

`ALL_NETWORK_RECEIPTS.json` contains 46 read-only network receipts. `EVIDENCE_VALIDATION.json` confirms every recorded request returned HTTP 200 and the retained response/PDF/text/decoded-document hashes match. These are provenance checks, not scientific replication. `TMD_SINCE_PREVIOUS_PIN.json` confirms two commits and 45 changed files since the old TMD pin; the six-document reading is a selected methods audit, not full inspection of all 45 files.

The source caches are private working references. They should not be copied into a public checkpoint wholesale. Original summaries, typed provenance and relevant source links are the candidate portable outputs.

The helper scripts reproduce read-only retrieval or provenance construction; rerunning them overwrites this working evidence and may observe later source heads. They are not physical solvers. `build_primary_evidence.py` performs byte/page bindings after the selected source files have been retrieved. No additional Python dependencies were installed: the initial optional `requests` import was unavailable, and retrieval was implemented with the standard library.
