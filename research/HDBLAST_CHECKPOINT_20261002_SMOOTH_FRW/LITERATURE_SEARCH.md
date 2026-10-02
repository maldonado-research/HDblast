# Bounded recent methods search

The search covered submissions from 2022-01-01 through 2026-10-02 UTC and combined recent metadata screening with older authoritative primary-method checks. Four arXiv Atom API queries returned bounded metadata sets. Article metadata/abstracts were used to select leads; the in-body claims in [LITERATURE_REVIEW.md](LITERATURE_REVIEW.md) come from the eight downloaded primary papers explicitly identified there. Version/hash metadata is in [LITERATURE_METADATA.json](LITERATURE_METADATA.json); exact queries and response hashes are in [LITERATURE_SEARCH_METADATA.json](LITERATURE_SEARCH_METADATA.json).

| Query focus | API total | Returned and screened |
|---|---:|---:|
| adiabatic + renormalization, recent date range | 123 | 20 |
| time-dependent mass + renormalization | 1 | 1 |
| stress + particle production | 5 | 5 |
| gr-qc/hep-th, adiabatic/renormalization titles, scalar/cosmological | 199 | 30 |

The first and fourth sets contain many unrelated uses of adiabatic/RG terminology; their totals are not counts of relevant cosmological papers. A driven Gross–Neveu coupling item was rejected as a different fermionic/integrability model. Other potentially relevant titles about inhomogeneous GW production, modified dispersion, static-shell vacuum and de Sitter infrared matching remain metadata-only leads; no in-body result from them was adopted.

Two newly screened sources were selected for primary passage checks: Ballesteros–Gambín Egea–Riccardi [2607.06042v1](https://arxiv.org/abs/2607.06042v1), for analytic UV/IR separation and finite local terms in later perturbation methods, and Gao–Anderson–Link [2308.11040v1](https://arxiv.org/abs/2308.11040v1), for state-compatible treatment of higher derivatives in future backreaction. Neither method is implemented in the present checkpoint; their different observables/couplings are described in the review.

Initial source requests returned HTTP/tunnel 403. After networking changed, primary PDF/HTML and metadata requests succeeded. A later complex category/variable-mass query timed out; a coherence query and a narrower exact-title adiabatic query returned HTTP 429. The search stopped without repeated rate-limit retries. These failures are retained in the metadata and limit coverage; they are not scientific null results.

This targeted search cannot establish that every relevant 2022–2026 paper was considered, that no competing method exists, or that the checkpoint has physical or mathematical priority. Seven earlier pointers remain explicitly marked inherited-only in the review. Full article bodies and raw response files are not included; source URLs, response hashes and query outcomes document what was actually checked.
