# Independent S4 source review: prospective dimensionless math check

Registered before numerical source values are computed, 2026-10-02 UTC.
Repository baseline: `0205cc651bfb614c32e39dfe833d93d957264229`.
This staged work is a mathematical benchmark, not a choice of physical HDBLAST
parameters, a field evolution, or a publication claim.

The independent implementation integrates the proper-time expression directly,
using the S4 eigenvalues l(l+3) and multiplicities (l+1)(l+2)(2l+3)/6.
It does not import the producer's Hurwitz-zeta determinant implementation.
Its small-proper-time heat coefficients come from the half-integer
Euler–Maclaurin expansion (Bernoulli polynomials); its remaining integral uses
explicit spectral sums. Q and rho use separately differentiated integrands.

Grid: fixed r=1, x in {1,2}, z=H² in {1/2,1}. All masses and curvatures
are strictly positive; these are dimensionless test inputs only.

Primary settings: 40 decimal digits, UV join r*s0=0.01, heat expansion
through bracket power 10, upper proper time r*L=60.
Refinement: 55 digits, r*s0=0.005, heat expansion through bracket power 12,
r*L=80. Quadrature uses mpmath Gauss–Legendre on geometric subintervals.
Spectral truncation uses a conservative Gaussian polynomial tail envelope.
IR tails use explicit incomplete-gamma/exponential-integral bounds.

Registered gates for W/r², rho/r², Q/r: refinement and cross-method
difference at most 1e-9 + 1e-7*abs(value). Trace residual uses the same gate.
Any failed gate remains in evidence. Negative controls flip the current sign,
flip the curvature response in metric variation, and change the 29/15 heat
coefficient. The executable must reject all three. No physical parameter
search or unregistered adjustment is permitted.

Analytic spectral and IR tail bounds do not certify total numerical error.
The UV expansion is asymptotic, and its omitted remainder is tested by
refinement but has no proved enclosing constant here. mpmath quadrature and
arithmetic are not interval-certified. Cross-method agreement and precision
refinement therefore constitute numerical evidence, not a rigorous enclosure.
