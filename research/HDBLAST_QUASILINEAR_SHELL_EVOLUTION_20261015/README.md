# Exact nonlinear shell evolution at fixed finite regularity

**Conditional mathematical result:** the selected aligned quadrupole has a
small exact Einstein–scalar family on both copies of the controlled
positive-cone region for a prescribed finite shell-interaction interval,
under the stated background, support and exact-initial-data hypotheses.
The physical shell laws and constraints hold. Evolution has fixed
`H^14 × H^13` regularity, with `C²` parameter dependence in the lower
`X_10` topology and hence a lower piecewise spacetime `C⁴` norm.

Read the [main theorem and proof](EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md)
for the exact domain, gauge, smallness and regularity qualifications.
This supplies existence across the reference interaction window; an actual
first-arrival time and nonzero gauge-invariant shell response remain to be
calculated. No global stability, common spatial `C∞` parameter family,
external novelty, external peer review, observational agreement or
higher-dimensional Big Bang origin is established.

## Proof and review

- [Variable-coefficient construction](astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md): positive energy, smooth frozen solver, tame estimates, exact initial jets, whole-interval Picard construction and finite Sobolev completion.
- [Independent derivation](sol/INDEPENDENT_QUASILINEAR_MIXED_PROOF.md) and [linear construction audit](sol/linear_audit/LINEAR_VARIABLE_WAVE_AUDIT.md).
- [Initial harmonic completion and subsidiary constraints](gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md), with the [finite-regularity parameter completion](gauge-parameter/FINITE_REGULARITY_VARIATION_COMPLETION.md).
- [Consolidated theorem review](referee/CONSOLIDATED_THEOREM_ACCEPTANCE.md) and [analytic cross-review](referee/FINAL_ANALYTIC_PROOF_ACCEPTANCE.md). These are independent internal agent reviews, not external peer review.
- [Primary literature review](literature/QUASILINEAR_CONTACT_PRIMARY_LITERATURE.md), with its source ledger. Full retrieved papers and private source archives are not redistributed.

The preserved [exact initial-data theorem](inputs/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md)
and [exact boundary-flow calculation](inputs/NONLINEAR_BOUNDARY_STRAIGHTENING.md)
are retained inputs. Their historical scope statements describe their
original checkpoints; the current claim is stated in the new main theorem.
`INPUT_SNAPSHOT_PROVENANCE.json` records their unchanged bytes. The earlier
specific second Taylor coefficient need not be the second derivative of
this exact initial family; its initial second correction can differ.

## Reproduce the finite checks

The required runtime is Python 3.12.14 and `sympy==1.14.0`. The repository's
`.github/workflows/hdblast_quasilinear_shell_evolution.yml` and the separate
publication receipt distribute external SHA-256 pins for this manifest,
ZIP, runner and verifier. Authenticate those bytes before executing source
code. `MANIFEST.json` then authenticates every regular capsule file, including
all proofs and review records; the matching ZIP contains exactly that tree.
The manifest cannot authenticate its own external trust anchor.

After that external authentication, from this directory run:

```sh
python -I -S -B verify_manifest.py --expected '<full external manifest SHA-256>'
python -I -B -m pip install --only-binary=:all: -r requirements.txt
python -I -B reproduce.py --manifest-sha256 '<full external manifest SHA-256>' --output-dir /tmp/hdblast-quasilinear-fresh
```

Use a new or empty output directory outside the capsule. The runner copies
each checker to a separate fresh directory and executes it in normal and
optimized Python. It compares every scientific output with the frozen
reference; only the gauge checker's intentional optimization-mode field
is excluded from scientific equality.

Expected verification:

| Check category | In each execution mode | Both modes |
| --- | ---: | ---: |
| Successful independent checker processes | 3 | 6 |
| Finite positive check groups | 19 | 38 |
| Deliberately corrupted formula processes rejected | 7 | 14 |
| Additional embedded gauge wrong-formula controls rejected | 6 | 12 |

The categories are distinct; the embedded controls are not separate
processes. The checks cover exact variable-energy identities, flux
cancellation and missing cross terms, finite derivative-count bookkeeping,
nonzero-shift harmonic completion fixtures, field-dependent-tension Codazzi
cancellation, the subsidiary Robin factor, and the trace-reversed initial
Hamiltonian projection. Their success does **not** computer-certify the PDE
existence proof. `REPRODUCTION_RECEIPT.json` is generated outside this frozen
source tree and records fresh outputs, logs, executable hashes and outcomes.

## Scope and next calculation

The proof constructs the whole prescribed interval from the original open
background shell collar. It does not claim a general finite-Sobolev restart
theorem at arbitrary later nonzero corners. The positive energy is an
auxiliary analytic norm, not a physical gravitational energy or stability
criterion. Constants and the small parameter interval depend on the fixed
background and region.

The next target is a gauge-invariant scalar or geometric response on the
physical shell, with a nonvanishing argument and controlled errors. Its
first and second variations must use the exact family's actual initial
completion. GitHub carries this checkpoint; no new Zenodo version or APS
submission accompanies it.
