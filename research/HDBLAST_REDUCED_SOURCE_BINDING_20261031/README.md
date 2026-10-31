# Reduced gauge sources, constraint propagation and actual third jets

31 October 2026. Exact implementation checks for the retained admissible
quadrupole, with an additive derivative adapter. Existing analytic
reconstruction/subsidiary proofs, original source maps, numerical candidates,
error certificates and frozen manuscripts are preserved.

The new calculation verifies the full-metric gauge lowering, gauge-gradient,
PC1 and reduced-source roles across their eight covector and 20 tensor modes.
It checks all eight general reduced Bianchi rows with first residuals and
PC1 retained before specialization. Separate arithmetic/frame implementations
derive the covector wave directly and agree on every exported modal coefficient,
including the radial scalar-field p² curvature contribution.

The new actual-field adapter supplies ordinary physical y/t profile jets
through third order: 70 metric and 10 scalar entries. It preserves the earlier
lower two-jets and all 24 existing time-containing third intervals, adding
eight pure radial yyy entries. Third-order intervals are coarse uniform
magnitude bounds under the inherited actual higher-jet certificates. They
are not smaller approximation-error bounds.

Read the [combined conditional constraint-propagation chain](theory/CONSTRAINT_PROPAGATION_CHAIN.md)
for the equations, premises and remaining limits.

## Reproduce

Use a complete checkout and Python 3.12.14. All distributed checks use the
standard library. From the repository root, choose a fresh output directory:

```sh
python -I -S -B research/HDBLAST_REDUCED_SOURCE_BINDING_20261031/checks/run_checkpoint.py \
  --repo-root . --output-dir /tmp/hdblast-reduced-source-fresh
```

The runner authenticates this closed capsule before and after use and
reproduces six complete canonical reports in normal and optimized Python.
It rejects 150 actual mathematical corruptions and 12 additional source,
import or readback corruptions, for 162 corrupted subprocesses total.
Mathematical controls require the intended nonzero witness/mismatch and no
successful output. The source controls use temporary copies and check that
changed dependencies are rejected before use or on readback.

| Canonical report | Scope | Mathematical controls per mode |
| --- | --- | ---: |
| [REDUCED](reference/REDUCED.json) | Source roles and general eight-row reduced Bianchi identity | 20 |
| [WAVE](reference/WAVE.json) | Independent rough-wave contraction, angular reconstruction and modal export | 13 |
| [INDEPENDENT_ROLES](reference/INDEPENDENT_ROLES.json) | Own arithmetic/frame source-role and wave derivation | 20 |
| [WITNESS](reference/WITNESS.json) | Coherent pointwise-premise counterexamples | 2 |
| [ACTUAL_THIRD](reference/ACTUAL_THIRD.json) | Same-actual-field third-jet adapter and prior interval joins | 10 |
| [INDEPENDENT_THIRD](reference/INDEPENDENT_THIRD.json) | Nonlinear coordinate substitution and independent raw composition | 10 |

The final row's 10 mathematical cases and six additional source cases per
mode run through the separate real-disk corruption driver. Counts describe
the execution roster, not independent theorem counts. All six complete
canonical report bytes must agree between modes and with the references.
Commands, timings and failure witnesses are generated separately in the fresh
output tree; no machine-specific runtime receipt is a canonical reference.

Sparse checkouts need this capsule, the complete unchanged
`HDBLAST_NATURAL_FIELD_BINDING_20261030` and
`HDBLAST_GEOMETRIC_BULK_SOURCES_20261016` roots, plus the selected files listed
in `actual_first_third_jets.py::SELECTED` from
`HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030` and
`HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030`, under `research/`.
Those two older roots have no manifest; the new code pins its selected
static inputs explicitly. It does not claim to authenticate their entire
capsules or to replay their higher-jet certificate proofs. The executed
geometric numeric/import closure is authenticated before importing it.

## Derivations and independent internal review

- [Source roles and reduced Bianchi identity](theory/REDUCED_SOURCE_BINDING.md)
- [Covector wave and exported modal coefficients](theory/COVECTOR_WAVE_BINDING_NOTE.md)
- [Actual third-jet adapter and input certificates](theory/ACTUAL_FIRST_THIRD_JET_ADAPTER.md)
- [Independent source-role and wave review](review/INDEPENDENT_REDUCED_REVIEW.md)
- [Analytic sign and prolongation review](review/ANALYTIC_SIGN_REVIEW.md)
- [Independent actual third-jet review](review/INDEPENDENT_THIRD_JET_REVIEW.md)
- [Bounded primary-literature comparison](literature/PRIMARY_METHODS_COMPARISON.md)

The source producer shares audited arithmetic/sphere primitives with the
wave producer. The Sol source/wave reviewer uses separate arithmetic and
frame calculus, importing no producer scientific code; it compares producer
wave output as authenticated data. The third-jet reviewer uses direct
nonaffine-coordinate substitution and separate product/inverse arithmetic.
These are internal AI-assisted reviews, not external specialist assessments.

## Scientific scope

Generic identities are checked within the registered seven-profile metric
quadrupole and scalar profile, with coherent ordinary derivatives. Full
angular reconstructions are retained. The eight covector Bianchi rows
annihilate TT tensors and cannot by themselves validate all natural source
components; the all-20 role joins and prior geometric checks have separate
responsibilities.

The homogeneous subsidiary statement retains the first natural equations
as neighbourhood identities, their needed derivatives, and the separate
harmonic conditions/prolongations. Matching QC values at a point is
insufficient inside a wave operator. The coherent witnesses illustrate this
logical distinction; they are not newly constructed physical perturbations
of the registered actual solution.

Actual background/master/gauge solutions, higher-jet certificates, physical
initial/shell/corner compatibility and subsidiary uniqueness remain the
identified analytic inputs. The numerical adapter encloses derivatives of
those same functions, with physical-side one-sided shell traces. Arbitrary
selections inside its interval boxes need not form a coherent jet. No
selected interval or unknown Taylor remainder is differentiated.

No improved response candidate, smaller error, exact finite-amplitude
nonlinear family, physical force, observational fit or higher-dimensional
Big Bang origin is established. The existing quantitative error limitations
remain. The next useful target is a distinct, quantitatively better response
candidate with explicit absolute and local error comparisons; this checkpoint
does not promise that improvement.

No new Zenodo version or journal submission accompanies this checkpoint.
Public supplements contain research material only; raw third-party paper
retrievals remain private.
