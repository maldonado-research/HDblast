# Bounded search log for the causal variance calibration

Searches and retrievals were performed on 2026-10-02 UTC, also 2026-10-02
in America/Los_Angeles. Machine timestamps below use UTC explicitly.
Companion JSON preserves exact query URLs, result counts, response hashes,
selected bibliographic records, and evidence status. Raw responses and
article bodies remain outside the public payload.

The available tool catalog contained no dedicated web search or browser
tool. Public verified-HTTPS requests to INSPIRE and arXiv were made from
the execution environment. No authentication or private-body retrieval
was used.

Nine INSPIRE requests were bounded to one page each:

1. `find eprint 0712.2282` — 1 result, limit 15; primary paper selected.
2. `find title semiclassical and date > 2025 and (keyword "linear response" or keyword "de Sitter")` — 1 result, limit 15; inherited holographic-stability source, metadata only this round.
3. `find title "nonlocal corrections" and title "de Sitter"` — 1 result, limit 15; 2601.22644 selected and reread.
4. `find eprint 2212.01078` — 1 result, limit 15; inherited source reread.
5. `find eprint 0907.0823` — 1 result, limit 15; inherited source reread, additional appendix passages checked.
6. `find date > 2024 and title "de Sitter" and (title response or title memory or title nonlocal)` — 10 results, limit 20; 2601.22644 retained; remaining nearby titles metadata-screened only.
7. `find date > 2024 and (title "Wick square" or title "composite operator") and (keyword curved or keyword cosmology)` — 0 results, limit 20; no completeness implication.
8. `find title "de Sitter" and (keyword "linear response" or keyword "initial conditions")` — 2 results, limit 20; neither added as equation evidence. The narrow keyword field did not retrieve the already known Anderson source, illustrating the query's incompleteness.
9. `find eprint gr-qc/0103074` — 1 result, limit 1; primary paper selected.

Requests 2, 6, 7, and 8 were topical screens. The other five were targeted
bibliographic retrievals. All returned successfully; no exhaustive search,
new-source discovery, or publication-priority inference follows.

Unversioned arXiv abstract pages were used only to verify submission
histories. The downloaded bodies were pinned to 0712.2282v1,
0907.0823v2, gr-qc/0103074v2, 2212.01078v2, and 2601.22644v1.
Versions and hashes are recorded separately from publication metadata.
INSPIRE supplied authoritative bibliographic cross-checks; a DOI does not
mean the publisher body was read.

The literature review identifies the precise passages read and two pages
visually inspected. Five PDF retrievals do not imply five whole-paper
validations. Prior checkpoint reviews were read to distinguish inherited
evidence from this round's passage checks. No new GitHub repository search
or private-body audit was conducted.
