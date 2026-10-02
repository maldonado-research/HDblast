# Bounded search and limits

This 2026-10-02 UTC search used thirteen single-page INSPIRE public API requests. Exact URLs, retrieval timestamps, response hashes, limits and returned counts are in LITERATURE_SEARCH_METADATA.json. No recursive pagination or all-owner recursive source crawl was performed. Returned titles/metadata selected primary papers, after which the particular primary passages listed in LITERATURE_REVIEW.md were actually read.

| Request focus | API total | Returned | Limit | Interpretation |
|---|---:|---:|---:|---|
| Recent de Sitter titles with scalar/vacuum/effective words | 119 | 20 | 20 | Recent screen; includes unrelated black-hole/AdS records |
| Recent de Sitter vacuum/effective/stress titles | 6 | 6 | 30 | Narrow title screen; limited recall |
| Recent adiabatic renormalization/regularization titles | 6 | 6 | 20 | Mostly different spins/dispersion; existing PSAR recovered |
| Recent scalar/de Sitter/renormalization broad terms | 127 | 25 | 25 | Low precision; selected 2026 effective-potential/nonlocal/matching leads |
| Broad scalar stress/effective potential since 2010 | 968 | 25 | 25 | Low precision; unrelated titles excluded |
| Older de Sitter scalar/stress/adiabatic titles | 390 | 35 | 35 | One recent-first page only; far from complete historical coverage |
| de Sitter tadpole/propagator/renormalization titles | 78 | 35 | 35 | One recent-first page; many different fields/geometries |
| de Sitter renormalized/zeta/effective-Lagrangian titles | 14 | 14 | 30 | Selected Chodos–Kaiser; excluded Kitamoto–Kitazawa after primary scope reading |
| Exact Ferreiro–Torrenti arXiv lookup | 1 | 1 | 30 | Bibliographic confirmation |
| Bunch + Davies before 1980 | 4 | 4 | 20 | Classic metadata lookup |
| Dowker + Critchley | 9 | 9 | 20 | Classic metadata lookup |
| Loops in de Sitter space title lookup | 2 | 2 | 5 | 2403.13145 retained as cited pointer, primary body not read |
| Unfielded de Sitter/digamma query | 1,598,264 | 25 | 25 | Query failure: returned a near-database-wide set; discarded |

The broad full-record and phrase searches match many irrelevant uses of de Sitter, including anti-de Sitter, black holes and other spins. API totals are not counts of relevant papers. The unfielded `find de Sitter and digamma` request failed to constrain its intended terms; it was discarded without a scientific inference. No result dated beyond 2026-10-02 was adopted. No search result means that all relevant 2022–2026 literature was considered.

Bibliography following from Bonanno–Cacciatori–Moschella identified García-Consuegra–Rajantie 2511.23076v2, whose May 2026 revision was then downloaded and primary-checked. Known-method lookup identified Markkanen–Rajantie 1607.00334v3. Each of the seven adopted new primary PDFs was fetched once through its unversioned URL, then retrieved again through its explicit version URL; all seven SHA-256 hashes matched. Primary checks covered selected sections and equation pages, not an asserted end-to-end referee review.

The earlier eight-paper SMOOTH_FRW review and metadata were inspected before this search. Seven are retained as prior evidence without unchanged re-review; the de Sitter section of PSAR is the sole new passage check within that prior eight-paper set. Ferreiro–Torrenti 2212.01078 was previously inherited-only and is upgraded to primary-passage checked. None of the recent interacting/IR methods was implemented by this literature task.

Publisher PDF requests for Bunch–Davies and Dowker–Critchley failed with HTTP/tunnel 403. INSPIRE metadata was available; their primary formulas are not claimed to have been checked. Other cited articles in the newly read bibliographies remain pointers unless explicitly listed as primary checked. Raw PDFs/text, rendered article pages and raw API bodies remain outside the public payload. The curated payload contains original assessment, bibliographic metadata, request receipts/hashes, and public-owner context only. No private handoff or unavailable local Mac folder was used as public evidence.
