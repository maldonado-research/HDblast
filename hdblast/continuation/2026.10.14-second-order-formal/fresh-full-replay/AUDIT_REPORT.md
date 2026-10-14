# Fresh second-order localized checkpoint audit

**Result: pass.** All ten hardened default-mode predicates pass on freshly
generated outputs. The original eight scientific scripts completed with exit
code 0, and the independent angular symbolic verifier and fresh `f₀`
comparison also passed. This is a combined audited replay, rather than a claim
that the revised runner executed the complete set in one uninterrupted process.

Completed 14 October 2026, 12:06:35 PM America/Los_Angeles
(`2026-10-14T19:06:35Z`). Python 3.12.14, SymPy 1.14.0, mpmath 1.3.0.

## Provenance and execution

The uploaded original capsule remains unchanged. All 75 original manifest
entries matched both byte count and SHA-256. Execution took place in a copy.
Every uploaded pickle was removed from that copy; `sources_cc.pkl` was
regenerated from inspected text code in 212.557 seconds. All five first-order
constraint residuals vanished symbolically during regeneration.

All old JSON outputs for the original eight-step replay were removed before
execution. Their timings and exit codes are preserved in
`REPRODUCTION_RECEIPT.json`:

| Script | Seconds | Exit |
|---|---:|---:|
| `new/kids_tensor_only.py` | 194.5 | 0 |
| `new/kids_symbolic_l01.py` | 208.5 | 0 |
| `new/l0_full_sector.py` | 63.9 | 0 |
| `new/wgg_canonical.py` | 30.7 | 0 |
| `new/exact_random_point.py` | 505.4 | 0 |
| `new/fast_verify.py` | 467.7 | 0 |
| `independent/indep_check.py` | 209.6 | 0 |
| `independent/indep_l0.py` | 2867.4 | 0 |

The enclosing original replay took 4548.737 seconds. Its orchestration stopped
after recording successful completion of those eight scripts and before
`compare_f0.py`; the precise stopping cause is not established. Finalization
was resumed honestly, preserving the completed costly runs. The prior `f₀`
comparison receipt was removed and the comparison rerun with exit 0 in
12.334 seconds. All ten fresh outputs were then validated with the revised
runner predicates. Detailed execution and resumed-finalization records are in
`INDEPENDENT_REPLAY_RECEIPT.json`; the compact final record is
`AUDIT_COMPLETION_RECEIPT.json`.

## Scientific verification

The independent numerical code checked ten points on two radii and two
profile cases. Maximum absolute background Hamiltonian residual was
`2.5714e-39`; maximum first-order residual was `2.8928e-39`; maximum relative
momentum residual was `5.3155e-41`. The angular Hamiltonian remainder was
constant to the receipt's precision on each sphere. The four fresh `f₀`
comparisons had maximum relative difference `1.4235e-20`.

The independent radial quadrature gave weighted-identity relative difference
`3.2212e-25`, maximum relative radial-gauge residual `1.6976e-36`, and
scalar-only weighted-integral relative residual `7.5294e-30`.

Three seeded exact-rational samples gave all 15 Hamiltonian/momentum residuals
exactly zero, and all nine corresponding wrong-formula detections passed.
The symbolic momentum verifier detected its doubled-scalar-momentum control.
No probabilistic error bound is claimed for the finite samples.

A new angular-polynomial reduction proves `P_H[h_T] = S_H − f₀` symbolically
for the regenerated source, with no sampling. It uses the independently
reviewed conformal identity for tangential trace-free tensors and reduces
polynomials modulo `|x|² = s` before radial coefficient substitution. Normal
and optimized execution passed (15.553 and 9.789 seconds). A wrong-sign `C₁`
coefficient caused exit 1 in both modes while preserving the reconstruction
relations, demonstrating that the result does not depend on Python assertions.
The source itself still comes from the original constraint Hessian generator;
its independent derivation is checked numerically by `indep_check.py`.

## Corrections and limits

Prepared text corrections fix the obsolete `basis/results_snapshot` path,
remove unsupported negligible-probability language, and distinguish the
Hamiltonian full-Hessian argument in Delay Lemma 8.1 from the momentum
divergence estimate recorded in Delay §9.3. Existing mathematical formulas
and the eight replayed scientific scripts remain unchanged; the existing
exact-sample script's proposed change is its docstring and result method label.
The new angular verifier strengthens the Hamiltonian verification.

The revised runner removes stale receipts, always regenerates source caches,
validates all outputs, rejects non-finite numbers, requires fresh successful
full-verifier receipts, and includes the `f₀` comparison. Its `--core` mode
omits only the expensive independent radial quadrature and labels that exclusion.
Five targeted corrupted-output controls and five lightweight orchestration
scenarios passed. The original slow `--full` verifier was not rerun and is not
counted as completed; earlier stalled runs remain part of the historical record.

The verified mathematical scope is localized second-order initial constraint
data with vanishing initial junction sources near the shell. Shell-contact
evolution, higher-order solvability, nonlinear solution existence, physical
energy positivity, and cosmological interpretation are not established by this
audit. No GitHub mutation or external submission was performed by the auditor.
