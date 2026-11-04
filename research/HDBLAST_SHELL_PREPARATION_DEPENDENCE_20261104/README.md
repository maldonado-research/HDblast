# Preparation dependence of intrinsic shell volume

4 November 2026. Under the retained analytic assumptions, the same first-order
Einstein–scalar packet admits compact second-order preparations with opposite
signs of the proper-time volume and expansion changes. The initial constraints,
initial shell collar, first tangent, first quadratic forcing and second TT pair
are preserved. This is a concrete comparison of two preparations, not a new
general principle of perturbation theory.

At the physical shell and `tau_star=rho_b*10^-6*(5/4)`, the second parameter
derivatives have these safe outward bounds:

| Prepared family | Relative volume `v2=V2/V0` | Expansion `Theta2` |
|---|---|---|
| Earlier preparation, preserved | `[-5.4520e16,-2.6521e16]` | `[-9.0208e21,-4.4458e21]` |
| New preparation, `lambda=-10^21` | `[1.0343e17,1.8442e17]` | `[8.3831e21,1.8496e22]` |

For sufficiently small nonzero perturbations of either sign, both exact families
still have positive volume and positive expansion. One is below the background
and the other above it at that same anchored proper clock. The parameter range
is qualitative, uncomputed and family dependent. These are model-normalized
second derivatives, not finite-amplitude field values or SI predictions.

## What is specified and checked

The original compact mean profile `f` is replaced by `f_lambda=f+lambda*h`, where
`h=rho^-2*chi_y`, `chi=ell*rho_b^3*Q(r)` and Q is the same fixed flat packet
bump. Its constraint moment vanishes exactly. An explicit compact metric
correction preserves the complete second constraint equation. The change is
physical initial-data freedom; it is distinct from the gauge conditioning used
to calculate it.

The difference evolves through the full coupled five-mean system with the
retained physical and gauge boundary rows. The first-field quadratic sources,
first proper-clock and measure corrections, and first-only volume remainder
cancel exactly when comparing these two second responses. They remain in the
absolute earlier calculation. No coupling is set to zero.

Two separately checked sufficient methods certify the signs:

- Whole-cell rational enclosures rebuild the new profile over 241 continuous
  cells and one target, then bound the complete homogeneous resolvent and the
  proper-time Green histories. Independent Fraction code reconstructs the
  profiles, norms, histories, slopes and finite-parameter witness.
- An elementary analytic certificate uses broad authenticated background caps,
  proved continuous bump-derivative bounds and exact rational inequalities.
  It obtains the same simultaneous sign reversal without quadrature.

The target, direction, coverage and sign gates were fixed before the detailed
evaluation. The finite `lambda` is chosen by a declared exact algebraic rule
from the interval endpoints; it changes no profile, target or error method.
Neither method makes new background or evolved-field oracle calls. Reused
background enclosures and analytic PDE inputs remain shared premises.

Read the [theorem and proof](DERIVATION.md), [elementary certificate](theory/PROOF.md),
[source comparison](REFERENCES.md), and [next targets](NEXT_STEP.md). The working
paper is provided as [PDF](manuscript/PREPARATION_DEPENDENCE.pdf) and
[editable TeX](manuscript/PREPARATION_DEPENDENCE.tex) for author and specialist
review. Compilation is documented separately; PDF recompilation is not a runtime
dependency of the mathematical replay.

## Reproduce

Use Python 3.12.14 from a checkout containing the pinned dependencies:

```sh
python -I -S -B research/HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-preparation-dependence-fresh
```

Choose a fresh output directory outside the checkout. The root replay
authenticates complete new and inherited inventories before scientific imports
and afterward. It reproduces new canonical reports in normal and optimized
Python and separates rejected subprocesses, successful controls and in-process
formula checks. Authenticating the earlier source closures does not replay
their proofs or prove their analytic hypotheses.

## Scientific meaning

The earlier sign was a result for its declared preparation. It is not determined
by the first packet and constraints alone. A physical model must select a second
preparation, or establish a prediction robust within a physically defined class.
One direction suffices for the counterexample; no independent two-observable
control or rank-two result is claimed.

The result is conditional on the retained background, localized-constraint
inverse, exact-family development, parameter regularity, mixed linear uniqueness
and anchored-observable premises. Internal AI-assisted mathematical review,
exact arithmetic and CI are distinct from external specialist review, a
proof-assistant formalization and a demonstration of originality.

The higher-dimensional origin of the Big Bang remains unverified. Matter and
radiation production, a physical scale, nonlinear remainder bounds and comparison
with observations remain open. Earlier capsules and manuscripts are preserved.
No Zenodo edition or journal submission accompanies this checkpoint.
