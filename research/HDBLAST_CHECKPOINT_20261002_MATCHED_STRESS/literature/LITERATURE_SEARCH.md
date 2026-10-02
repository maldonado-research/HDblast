# Bounded stress-response literature update

Date: 2026-10-02 in America/Los_Angeles. All machine retrieval and file
verification timestamps are UTC and appear in the companion JSON. The
available tool catalog had no dedicated web search/browser tool. Requests
used the public INSPIRE API and arXiv abstract pages over verified HTTPS.

Six topical queries were each limited to one page of at most 25 records:

| Query | Total hits | Records screened |
| --- | ---: | ---: |
| `find date > 2024 and (title "stress tensor" or title "stress-energy") and (keyword "de Sitter" or keyword scalar or keyword renormalization)` | 6 | 6 |
| `find date > 2024 and (title "linear response" or title "Ward identities") and (keyword scalar or keyword semiclassical or keyword cosmology)` | 1 | 1 |
| `find date > 2024 and title "adiabatic" and (title renormalization or title regularization)` | 2 | 2 |
| `find date > 2024 and (title "time-dependent mass" or title "time dependent mass") and (keyword quantum or keyword cosmology)` | 0 | 0 |
| `find date > 2024 and title "linear response"` | 54 | 25 |
| `find date > 2024 and (title renormalization or title regularization) and title scalar` | 21 | 21 |

The broader fifth query recovered the known Lai–Ota source, which the narrow
keyword query missed. This illustrates a limitation of keyword filtering.
It also found Wang–Cui–Wu, *Covariant linear response theory for a photon gas
in curved spacetime*, arXiv:2609.02615. Its abstract specifies Boltzmann
relaxation-time transport and null-particle phase-space geometry. It was
excluded as a different observable and approximation, not cited as a
vacuum stress or renormalization result. No body was read for that lead.

Other nearby records concerned modified dispersion relations, spinors,
time-dependent compact dimensions, AdS energy conditions, gravitational
tensor stress, and primordial spectra. Those screened records were not
promoted to equation evidence. The compact-dimension and primordial-spectrum
papers were already inherited pointers. The date filter also returned
several pre-2025 preprints by their later publication dates. Neither these
hits nor the zero-result query establishes completeness or novelty.

Five additional targeted INSPIRE requests used `find eprint` with size=1
for 2311.08986, gr-qc/9908037, 1607.00334, 2606.16296, and 0907.0823.
Their unversioned arXiv abstract pages were retrieved to check version
histories. The actual primary bodies read were the previously downloaded
2311.08986v2, gr-qc/9908037v1, 1607.00334v3, 2606.16296v1, and
0907.0823v2 PDFs. Their recomputed hashes matched prior metadata. No new
primary PDF was downloaded in this round.

All eleven INSPIRE requests and five abstract-page requests succeeded.
Exact URLs, timestamps, response hashes, hit counts, and selected metadata
are preserved in `LITERATURE_SEARCH_METADATA.json`; body/version/hash and
passage provenance is in `LITERATURE_METADATA.json`. PDFs, text extraction,
rendered pages, and raw response bodies remain private. There was no new
owner-repository audit, external peer review, or physical numerical test.
