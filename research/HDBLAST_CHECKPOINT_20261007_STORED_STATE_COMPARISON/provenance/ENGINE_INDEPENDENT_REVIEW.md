# Independent generic target-engine review

The generic entire-moment engine is mathematically valid conditional on its
caller supplying the authenticated prescribed-source polynomial family and
certified uniform errors. This review grants no real-source or retained-array
GO. The engine contains no I/O and executes no source callback.

With a=-9/2 and y=2(a-s)=-9-2s, the computed moments are
M_n=integral p(s)y^n ds. The degree-832 coefficients are exactly
`W_n=-i^n*M_n/n!` and `U_n=-i^n*M_(n+1)/(2*(n+1)!)`. The normalized affine
substitution uses the fixed parent center and halfwidth. All 113 ordinary
source coefficients survive; successive affine products and polynomial
antiderivatives give outward Arb512 enclosures of the moments. No decimal or
binary64 conversion enters this calculation.

The exact polynomial L1 majorant controls all moments because 0<=y<=1. The
geometric tail ratios are bounded by 256/834 for W and 256/835 for U, giving
the stated factorial-tail formulas. These are errors of the finite polynomial
kernel representation. Source uniform errors add E to W and E/2 to U, with
E=sum 2*half*uniform_error. They must also remain distinct from the omitted
flat-cap error.

The cap formulas agree with the previously proved integration-by-parts
identities in the frozen prior independent method. Their terms include both
endpoint derivatives and the interior 6L^2 term. They apply throughout the
real band and retain the cap as a positive analytic error. The cap is neither
evaluated nor assumed zero by this review.

Every component endpoint is taken from an outward Arb bound, converted to an
exact rational, expanded by its full analytic error, then floored/ceiled to
the 2^-96 grid. The complete exported L1 radius is tested afterward against
1e-18 for every U and W. At k=0, real forcing proves the exact target
imaginary components vanish; the engine checks its arithmetic imaginary
component is exact zero before this narrowing. It does not project any
stored quantity.

`INDEPENDENT_ENGINE_CONTROL_REEXECUTION.json` records a separate reexecution
of the nineteen fabricated controls and nine rejecting mutations under the
repository's research environment. They pass. The controls check constant
finite-source exact solutions at zero, tiny and high k; retention of degree
112; nonzero stored normalization; kA sign and factor; and fail-closed
geometry/degree/precision/error/gate requirements. The author's 49,152-node
fabricated resource receipt is separate evidence, rather than an independently
rerun benchmark or a physical result.

The essential remaining integration condition is the actual guarded caller.
It must authenticate the frozen prior independent constructor and all its
dependencies, the prior BUDGET and source proof, and the new complete public
registration before any of the twenty-two source constructions. For each
ordered source/parent it must compare the rebuilt canonical 113-coefficient
vector SHA256 to the prior recorded SHA256, verify its exact geometry, and
use the safe certified uniform error
`max(recorded_uniform_error, recorded_analytic_tail+recorded_coefficient_error)`.
A digest alone supplies no coefficients and proves no source approximation.
The real caller must reject `source=None`, require the correct registered
source label, and retain all four capsule identities and all 49,152 nodes.

The generic library deliberately does not enforce this external provenance
contract itself. Until the caller and public freeze are reviewed, its passing
controls support only the conditional generic engine. The full twelve-case
pressure certificate remains **UNRESOLVED**, metric calibration remains
**FAIL_UNCHANGED**, and the proposed higher-dimensional cause remains
**NOT_ESTABLISHED**.
