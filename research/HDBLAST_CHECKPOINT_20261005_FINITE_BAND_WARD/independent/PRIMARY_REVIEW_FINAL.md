# Independent final review of the finite-band analytic checkpoint

Assessment: **PASS for scientific publication of the stated math-only
checkpoint**. The approved scientific scope is an exact finite-band
forcing-error kernel, complete off-shell action/Ward budget, and explicit
all-time analytic source-truncation component bounds. The full twelve-case
numerical pressure/contact certificate remains **UNRESOLVED**; external
mathematical novelty remains **NOT_ASSESSED**.

The independent action/bilinear derivation and its first exact symbolic
receipt preceded reading the new primary implementation or derivation.
Full approximate-state residual and closed-band coefficients were sent to
the root agent before the primary coefficients were read. The independent
SymPy route passes 42 exact identities and 10 omission controls in both normal
and optimized Python, with matching stable receipts.

The refreshed standard-library verifier `../verify_finite_band.py` passes 41
checks. `../EXACT_CHECKS_FINAL_NORMAL.json` and its optimized counterpart
are byte-identical, SHA256
`71fe37795f2c5c7a9dad2219f415604bfea588b9c9a6434d9abf2464178f7353`.
The current checker and theory hashes are bound in
`FINAL_PUBLICATION_REVIEW_RECEIPT.json`.

The full raw and epsilon-normalized state residuals, signed forcing mismatch,
missing canonical drift, independent density/pressure kernel extraction,
finite-K moment corrections, Pi convention and rational K=64/128/256 source
bounds agree with the independent derivation. The primary's joint-complex
state-error coefficient correctly sharpens the triangle component bound by
Cauchy--Schwarz. Its integrations and zero-lag limits are valid.

The reflected Cauchy-envelope first moment alpha/2 is valid on the full
domain and every prefix. The finite 41 checks verify the reflected centers
and exact rational constants; the documented calculus proof supplies the
infinite-series integration and reflection step. Similarly, entireness of
the stress kernels follows from their compact momentum integral definitions,
while the finite series tests verify the implemented removable limits and
coefficients. No proof-assistant formalization is claimed.

The interface qualification is necessary and correct. A bounded piecewise
Taylor forcing produces continuous absolutely continuous Duhamel amplitudes
and mode stresses; the differential identity holds almost everywhere and
the integrated identity has no state jump. An arbitrary approximation with
discontinuous density or contact data instead requires the explicit signed
interface jumps, or their absolute allowances in a bound. An almost-everywhere
contact residual alone cannot bound distributional interface impulses.

The added homogeneous counterexample is stronger than a real constant
amplitude: `U=A(k) exp(2ik(t-a)), W=2ikU` satisfies both homogeneous equations,
the first-order Wronskian and zero Ward defect. For A nonnegative, smooth and
supported in `(2,4)`, it is invisible at all nine inherited probes but has
strictly positive nonzero integrated initial density at every registered
cutoff. Thus conservation, Wronskian and sparse probes do not certify the
incoming state or continuum stress. In the specified global plane-wave basis
its exact normalized completion uses
`beta=epsilon A(k) exp(-2ika)`, `alpha=sqrt(1+|beta|^2)`.

No blocking algebra, scope or gate defect remains. The density and pressure
numbers are declared **analytic truncation-only** bounds for exact degree 24
Taylor coefficients with identical incoming states and unchanged contacts.
They exclude coefficient/phase arithmetic, state/geometry/contact errors,
time/momentum quadrature and cutoff tails. The pressure bound is not promoted
to the distinct `1e-26` source-moment gate. No real source callback, physical
array, physical trajectory or remote operation was used by this review.

Publication of these conditional analytic claims supplies a reviewable
methods checkpoint. It does not establish particle heating, thermalization,
nonlinear Einstein evolution, a hot radiation era, agreement with data, or
a higher-dimensional cause of the Big Bang. All inherited metric failures
remain unchanged.
