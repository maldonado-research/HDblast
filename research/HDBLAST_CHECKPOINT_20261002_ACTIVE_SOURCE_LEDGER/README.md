# Source-active conservation diagnostic

This checkpoint tests the fixed saved interval `[-4.5,-3.5]` for both compact
metric sources, both inherited resolutions and `K=64,128,256`: twelve cases.
The saved inputs and earlier active-source failures were already published.
The new reference-flow, contact and direct-integral attribution diagnostic is
frozen prospectively before evaluation, with that prior knowledge explicit.
The completed source-free diagnostic established a Simpson-ledger error only
after the source ended. It cannot establish conservation during forcing.

The new protocol keeps separately defined density and pressure operators,
the full metric/subtraction contacts, the saved momentum rules and the source
work. It compares a phase-aware Green-function ledger with an independently
integrated direct stress ledger. The saved final density never defines either
integral. Signed flow, momentum/contact, operator and reconstruction terms are
reported separately, with fixed midpoint and endpoint audits.

All earlier metric experiments retain **FAIL**, including the 29 endpoint and
30 ledger-refinement failures. The source-free result and its scope are
unchanged. This diagnostic supplies no evidence of an extra dimension, blast,
particle yield, thermal history or a hot Big Bang.

Read `PROTOCOL.md`, `REPRODUCTION.md`, the formal proofs in `theory/`, and the
bounded primary-literature review in `literature/`. The complete prospective
registration must be published and authenticated before physical evaluation.
Subsequent execution evidence and outcomes belong in `reports/` and `outputs/`;
they do not change the frozen protocol, algorithms, inputs or gates.

Native long-double source/phase calculations and numerical quadrature controls
are empirical checks. The 80/100-digit reductions and contact evaluations have
the explicitly narrower scopes stated in the implementation notes. They do
not certify the full trajectory or its inherited initial state.

Publication follows the existing HDBLAST concept DOI
`10.5281/zenodo.17088132`. Automatic GitHub archiving is OFF; the existing
main-family draft and the separate static companion are preserved. A reserved
draft DOI is not a published version. Current publication state is recorded in
`publication/CURRENT_PUBLICATION_STATE.json`.
