# Bounded search provenance

Six one-page INSPIRE queries were run on 2026-10-02 UTC. Their exact URL, query, UTC retrieval time, response SHA-256, returned titles and count are in [LITERATURE_SEARCH_METADATA.json](LITERATURE_SEARCH_METADATA.json). Five selected primary reads and hashes are in [LITERATURE_METADATA.json](LITERATURE_METADATA.json).

| Query key | Returned | Reported total | Limit |
| --- | ---: | ---: | ---: |
| recent_semiclassical | 29 | 29 | 30 |
| recent_brane | 19 | 19 | 30 |
| recent_species | 10 | 10 | 30 |
| anderson_mottola | 6 | 6 | 30 |
| brane_quantum | 30 | 34 | 30 |
| starobinsky | 8 | 8 | 30 |

The old brane query returned 30 of 34 records; the page was not completed by pagination. Brane New World was retrieved directly as a known primary lead. Recent-query results include some older preprints because INSPIRE date expressions can match publication dates. None of the totals is a relevant-paper count.

The five version-pinned PDF retrievals succeeded and matched the corresponding unversioned PDF hashes. No additional body was read for Starobinsky 1980, or for recent titles screened only in metadata/abstracts. This bound deliberately favors close primary methods over an expansive novelty claim.

No browsing tool was available. Python urllib over verified HTTPS accessed public INSPIRE/arXiv endpoints. All article bodies, extracted text, raw source responses and page renders remain outside the public payload. No private chat, credentials or third-party full text are included.
