# Exact finite-band stress kernels and a complete conservation-error budget

Ricardo Maldonado · 5 October 2026 UTC · mathematical methods checkpoint

The unchanged linear metric-response model now has exact finite-band kernels
for the independently defined energy density and pressure, and an explicit
error budget separating the incoming state, differential-equation residual,
source approximation, contact terms and quadrature. Applied to the inherited
degree-24 source Taylor remainder, these kernels give uniform bounds on the
**analytic truncation contribution alone** throughout the declared time interval.
The full twelve-case pressure/contact certificate remains **UNRESOLVED**.

The conservation identity itself was established in an earlier checkpoint.
This continuation integrates the forcing-error kernels in momentum exactly,
transports state and residual errors, and supplies a constructive obstruction
to certification from sparse momentum probes. It does not claim a new physical
law, externally assessed mathematical novelty, or a higher-dimensional cause
of the Big Bang.

## Fixed scope and conditional premises

The model retains conformal time η ∈ [−9/2,−7/2], L=−1/η, L′=L²,
minimal coupling, physical mass squared and fixed subtraction reference x=r=2,
and the two registered smooth metric sources. With normalized amplitudes
U=u/ε and W=w/ε, the exact flow is U′=W, W′=2ikW−g, with
g=4L²h−2Lh′−h″. Density and pressure keep the inherited a₀⁴/ε scaling.
All results use the same declared momentum measure k²dk/(2Pi²); the rational
upper bounds need only Pi≥3 and cover the inherited represented Pi constant.

The proof is conditional on the archived source definitions, Cauchy remainder
proof, prescribed geometry, and complete action/subtraction contact inventory.
Pinned copies of seven inherited analytic documents are supplied in
[`premises/`](premises/); their original repository paths and hashes are in
[`theory/INPUT_PINS.json`](theory/INPUT_PINS.json). The independent derivation
checks the arbitrary-jet W0/W2/W4 identities, but the programs are exact algebra
checks, not proof-assistant formalizations of all analytic premises.

## Continuous-band source-truncation contribution

The inherited exact Taylor approximation has 64 equal panels, degree 24,
complex Cauchy radius 1/8 and common bound 64. Its absolute remainder integrates
to at most α*=1/(390·2⁹⁰); the envelope's first time moment is at most α*/2.
The finite-band kernels integrate this forcing error without momentum
quadrature. The following outward upper bounds hold for every time in the
interval, with identical incoming data and unchanged exact contact terms.

| Cutoff K | Scaled density contribution | Scaled pressure contribution | Integrated source-mismatch work |
|---|---:|---:|---:|
| 64 | <3.4117173×10⁻²⁹ | <3.3854823×10⁻²⁷ | <3.3666285×10⁻²⁹ |
| 128 | <1.3556692×10⁻²⁸ | <2.6948292×10⁻²⁶ | <1.3466514×10⁻²⁸ |
| 256 | <5.4046411×10⁻²⁸ | <2.1504587×10⁻²⁵ | <5.3866056×10⁻²⁸ |

Exact rational values, rather than the displayed decimals, determine these
bounds: [`EXACT_CHECKS_FINAL_NORMAL.json`](EXACT_CHECKS_FINAL_NORMAL.json).
They concern the continuous interval k∈[0,K], not just the nine earlier probes.
They are analytic truncation bounds, not measured total errors or acceptance
criteria. Rounded Taylor coefficients, phase/arithmetic uncertainty, incoming
state, contact evaluation, time integration, geometry and the ultraviolet
tail above K still need separate bounds. This result does not extend the
earlier nine-probe arithmetic-width certificate to every momentum.

## What the error identity retains

For a reconstructed real forcing f and approximate state, define
rᵤ=Û′−Ŵ and r𝓌=Ŵ′−2ikŴ+f. The direct density obeys

\[
\widehat R'-\widehat F
=\frac{(2k^2+3L^2)\operatorname{Re}r_u-k\operatorname{Im}r_w
       -L\operatorname{Re}r_w}{2k}
 +\frac{L(f-g)}{2k}+\zeta_C,
\qquad
\widehat F=L(\widehat R-3\widehat P)-3h'\widehat B_0.
\]

Here ζC is the independently evaluated contact defect. Pressure comes from
the action-derived operator; the identity is a check on it. Dropping the
canonical normalization drift, the source mismatch or the direct contact
terms can produce an incomplete conservation test. Changing from a continuum
measure to a discrete quadrature also produces a known measure-mismatch term.
Piecewise forcing is allowed with continuous Duhamel states; any discontinuous
contact reconstruction needs an explicit interface-jump contribution.

The proof gives exact retarded kernels and Duhamel transport for the incoming
state and residuals. It identifies where to add time quadrature and endpoint
serialization errors, without changing the inherited 2×10⁻⁸ full-integral gate.
See [`theory/FINITE_BAND_ERROR_TRANSPORT.md`](theory/FINITE_BAND_ERROR_TRANSPORT.md)
and the independent [`INDEPENDENT_WARD_DERIVATION.md`](independent/INDEPENDENT_WARD_DERIVATION.md).

## Why the initial state needs a separate bound

A smooth real A(k), supported strictly between k=2 and k=4, can be invisible
at all nine registered probes. The homogeneous perturbation
δU=A(k)exp[2ik(η−a)], δW=2ikδU preserves the flow, first-order Wronskian and
Ward identity. Nevertheless its initial integrated density changes by

\[
\delta R_K(a)=\frac{3L(a)^2}{4\mathrm{Pi}^2}
             \int_2^4 kA(k)\,dk>0.
\]

A normalized Bogoliubov completion exists with the basis phase included in
β=εA(k)exp(−2ika). Thus conservation, normalization and a finite set of probes
do not certify the incoming continuum state. The next numerical certificate
must include an independently justified state enclosure or spectral regularity
assumption. This is a methodological counterexample within the fixed model.

## Verification and literature

Three separately implemented routes pass in ordinary and optimized Python:

- 41 standard-library exact polynomial/rational checks; both output files
  reproduce byte for byte.
- 55 independent SymPy identities/rational checks and 12 rejected mutations.
- 42 independent action/subtraction identities and 10 rejected omission controls.

No physical source callbacks, archived quantum-array decoding or trajectory
runs occur in these checks. The final independent publication review reports
no blocking mathematical finding for this limited scope. It is internal
AI-assisted verification, not external peer review.

The bounded [`literature review`](literature/PRIMARY_LITERATURE_REVIEW.md)
reads selected primary passages from six papers and records provenance,
limitations and a recent real-time ultraviolet/contact calculation. Related
work supports care with contact terms, state preparation and renormalization;
it supplies no missing project result or evidence for the cosmological hypothesis.

## Reproduction and next step

[`REPRODUCE.md`](REPRODUCE.md) documents the portable, manifest-checked replay.
The [`standalone package`](../HDBLAST_FINITE_BAND_WARD_PACKAGE_20261005/RESULTS_READ_FIRST.md)
contains the proof, exact verifiers, results and premises. Scientific receipts
are compared after excluding only explicitly named runtime fields, and frozen
payload bytes must remain unchanged.

Next, enclose the actual incoming state, coefficient/arithmetic uncertainty,
full contact inventory and time/momentum integrals for the unchanged twelve
cases. Only then rerun the unchanged metric calibration. Earlier endpoint and
refinement results remain **FAIL**. Coupled evolution, heating, cosmological
agreement and the higher-dimensional Big Bang hypothesis remain **NOT_ESTABLISHED**;
external novelty remains **NOT_ASSESSED**.

This dated result is published on GitHub. It is not a new Zenodo version;
[`current publication status`](../../hdblast/CURRENT_STATUS.md) tracks the
existing same-family draft and its metadata-recovery requirement separately.
