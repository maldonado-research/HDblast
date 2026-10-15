# Explicit compact second initial constraint jet

15 October 2026. This continuation constructs all physical second initial
Cauchy fields for the registered aligned Einstein–scalar quadrupole so that
the Hamiltonian and four momentum constraints hold through second order.
The correction stays in the original packet annulus. It follows the
[interior derivative and bulk-source checkpoint](https://github.com/maldonado-research/HDblast/tree/a85ba66bdde23ece81648d1826849829d43e0876/research/HDBLAST_INTERIOR_MASTER_JETS_20261015).

Two results are established under the retained actual-background assumptions:

- Zero second initial data fail: at the packet maximum, the quadratic
  Hamiltonian source is strictly below `-1e19` in the stated model units.
- Explicit local formulas cancel all five constraint sources. Allowing the
  second scalar and normal scalar-velocity fields to vary makes this possible
  without a radial boundary-value solve or an enlarged support annulus.

For scalar angular eigenvalues `lambda=0,8,24`, the Hamiltonian correction is

```
f_l = -rho^2 S_l / [2(lambda_l-3)],
phi2_l = -3 partial_y f_l / (2p).
```

Its three denominators are `-6,10,42`; the absent dipole eigenvalue would be
the exceptional value. The pure-trace extrinsic-curvature correction and
normal scalar-velocity correction cancel both angular and radial momentum
sources. All operations are local differentiation and multiplication by
smooth background factors with `p` bounded away from zero. The resulting
profiles are smooth, flat at their support endpoints, and vanish on an open
initial shell collar on both identical copies.

Read the [complete derivation and conventions](EXPLICIT_COMPACT_SECOND_INITIAL_JET.md).
The path uses physical Cauchy fields `(gamma,K,phi,Pi)`, with K an extrinsic
curvature tensor and Pi a unit-normal scalar derivative. They are not ADM
momentum densities. The derivation explicitly converts this path to spacetime
coordinate velocities, including the lapse cross terms. Omitting that
conversion would leave an incorrect momentum constraint.

The [independent implementation and review](review/INDEPENDENT_INITIAL_JET_REVIEW.md)
and fresh source-hashed review accompany exact standard-library checkers.
The producer has 48 generic/sector checks and rejects seven specified errors
in each Python mode; the independent implementation checks 75 groups and the
actual nonzero witness. These finite activities support the displayed analytic
identities and conditional bounds. They establish neither external peer review
nor novelty.

## Reproduce

Use Python 3.12.14 with no third-party dependencies or network. From a fresh
capsule extraction, choose a new external output directory:

```sh
python -I -B reproduce.py --manifest-sha256 ACTUAL_EXTERNAL_MANIFEST_HASH --output-dir /tmp/hdblast-initial-jet-replay-new
```

The adjacent publication record and CI workflow externally pin the manifest
and matching ZIP. The runner authenticates all bytes, compares complete normal
and optimized reports, executes all fourteen intended-error subprocesses,
and checks that the source remains unchanged. `REPLAY_ROSTER.json` records
the exact jobs. Both implementations authenticate the actual background report.
The included LICENSE retains the repository's code/text licensing.

## Remaining evolution problem

These formulas solve the differentiated initial constraints for this chosen
second jet. They do not solve the constraints exactly at finite epsilon or
construct an exact nonlinear initial-data family. Quantitative norms of the
second initial fields, the first relative-harmonic gauge bridge and transformed
second data, the complete common-gauge boundary forcing, and a validated
second-order evolution remain open. The earlier frozen checkpoints and their
conditional constructions retain their stated scope. No invariant force,
nonlinear remainder, cosmological origin, Zenodo version or journal submission
is established by this GitHub checkpoint.
