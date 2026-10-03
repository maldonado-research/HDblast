# Conditional certificate for the registered source and Duhamel operators

Completed **3 October 2026, Pacific time**. Both independently implemented
registered calculations pass the complete numerical enclosure gate. This
establishes a bounded integration prerequisite in the fixed-geometry model.
The full twelve-case pressure/contact conservation certificate remains open.

Each route returns **1,170 moment rows and 130 source-primitive rows**. Every corresponding
real and imaginary enclosure overlaps across the two implementations. The
largest whole-interval hull radius is exactly

`1459446311 / 10633823966279326983230456482242756608`,

which is **less than 1.373×10^-28**, below the unchanged registered gate `1e-26`.
The radius is the sum of the outward real and imaginary component half-widths.
It includes source approximation, phase/defect, coefficient/arithmetic,
propagation and endpoint serialization error. It is an absolute error in model
units, not a relative error or an observational sensitivity.

## Exact target and coverage

On `a=-9/2`, `b=-7/2`, define `z=eta+4`,

`B=exp(-z^2/(1-z^2))`, `h=B` or `h=zB`, `L=-1/eta`,

`g=4 L^2 h - 2 L h' - h''`.

The registered source operators are

\[
 M_0=\int_a^b g(s)\,ds,\qquad
 M_{\exp}(k)=\int_a^b e^{2ik(b-s)}g(s)\,ds,
\]
\[
 M_u(k)=\int_a^b\Phi_k(b-s)g(s)\,ds,\qquad
 \Phi_k(t)=\frac{e^{2ikt}-1}{2ik},\quad\Phi_0(t)=t.
\]

The same moments are enclosed on each of 64 equal panels. The separate real
primitive `J_Lg=integral Lg` is also enclosed on every panel and the whole
interval. This primitive is not the action-derived pressure/contact work.
The exact momentum probes are `0,2^-40,2^-12,1/4,1,16,64,128,256`.

The analytic source and local phase bounds hold throughout the 64 panels and
real `k` in `[0,256]`. **The complete computed numerical widths apply to the
nine registered probes**, not every momentum in that continuous range. This
interior interval excludes the bump's nonanalytic support cutoffs. Neither
incoming quantum-state arrays nor archived physical samples are inputs.

## Two enclosure routes

The primary method uses 256-bit Arb balls, degree 24 source jets and degree 96 
entire local kernels. Its global transport uses validated exponential and
regularized small-phase drift evaluation. The independent method uses its own
exact rational source recurrences, a degree 200 alternating scalar-exponential
enclosure, 512-bit chosen coefficients and a degree 96 polynomial ODE defect.
It bounds the actual residual against the true frequency and carries every
prefix error, point displacement and upward radius rounding. It imports no
primary numerical helper. Zero momentum is treated through exact real-target
identities, without unstable division by a tiny momentum.

Closed complex disks of radius `1/8` give common bounds `|g|<64`, `|Lg|<32`.
The degree 24 uniform tails are `1/(15*2^90)` and half that amount respectively.
The explicit phase tails and every scale/transport factor are retained in the
positive error evidence. The primary polynomial-only arithmetic baseline is
an inspection aid, not an additive partition of all rounding error. Small
signed residuals and cross-method agreement are not substituted for bounds.

| Route | Authoritative wall time | Peak RSS | Registered source constructions |
|---|---:|---:|---:|
| Primary Arb | 2.504 seconds | 56,556 KiB |128 |
| Independent rational defect | 30.947 seconds | 46,372 KiB |128 |

Both runs complete within the fixed 120-second/128-MiB budget, as nonroot
processes with one numeric thread. Every source attempt is recorded before
the callback. Frozen sources are unchanged before and after the calculation.

## Provenance and reproduction

The [public prospective freeze](https://github.com/maldonado-research/HDblast/commit/ead803c9ba575877fbda000cf2452f3ea35e272b)
preceded both real-source runs. All 155 registered files, including the
registration itself, were independently downloaded from public Git blobs and
matched by path, byte length, Git hash and SHA256 before execution.

- [Immutable protocol](PROTOCOL.md) and [full registration](FULL_REGISTRATION.json)
- [Public readback receipt](FREEZE_RECEIPT.json)
- [Exact result and resource review](RESULT_REVIEW.json)
- [Primary complete output](results/primary/OUTPUT.json) and [independent output](results/independent/OUTPUT.json)
- [Mathematical source proof](primary/cauchy-proof/CAUCHY_REMAINDER_PROOF.md)
- [Independent method and proof review](independent/bound-review/SOURCE_OPERATOR_FORMULA_REVIEW.md)
- [Standalone package and fresh reproduction](../HDBLAST_VERIFIED_SOURCE_OPERATOR_PACKAGE_20261003/RESULTS_READ_FIRST.md)
- [Primary literature review and bounded-search provenance](literature/PRIMARY_LITERATURE_REVIEW.md)

The preparation preserves disclosed series-cap, schema, radius-label,
serialization and guard defects and their repairs before real-source evaluation.
Earlier sampling/metric outcomes were already known; this is not described as
a blind confirmation. Internal independent-agent checks are not external peer
review or proof-assistant formalization. The certificate is conditional on the
archived analytic identities/proofs, numerical-library contracts and custodian
readback chronology. External mathematical novelty is **NOT_ASSESSED**.

## What this enables, and what remains open

For zero incoming forcing state, `u'=w`, `w'=2ik w-epsilon*g` implies
`w(b)=-epsilon*Mexp` and `u(b)=-epsilon*Mu`. These enclosures can supply a
controlled forcing contribution to a later mode/stress calculation. Initial
state, momentum integration, renormalization, action pressure/contact terms
and geometry errors still require their own analyses. The original 29 endpoint
and 30 refinement metric failures remain **FAIL**; they were not rerun or repaired.

An analytic sensitivity corollary follows from `|g|<=64`: on unit interval
length, `|dMexp/dk|<=64` and `|dMu/dk|<=64/3`, by differentiating the entire
kernels and integrating absolute bounds. These are standard consequences of
the model bound, not a new physical law or a sparse-probe continuum certificate.

The next registered study should propagate the bounded source into the full
action-derived finite-band stress/work ledger and account for its independent
state and momentum errors before another unchanged metric calibration. The
recent trace-coupling literature also makes physical energy transfer and
thermalization explicit requirements. Nonlinear persistence, a hot radiation
era, agreement with cosmological data and a higher-dimensional cause of the
Big Bang remain **NOT_ESTABLISHED**.
