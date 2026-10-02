# Independent review of the logarithmic-memory follow-up

Prepared 2 October 2026, before any follow-up physical evaluations. This is an
analytic implementation review by a separate project agent, not external peer
review. The original registered failure and its partial output remain in
`metric-original-complete-replay-001`: SciPy's weighted-log quadrature raised
its roundoff warning at the first active positive_B observation eta=-4.5.
This review does not reinterpret that failed experiment as a pass.

The proposal changes the arithmetic and endpoint treatment of the continuum
logarithmic-memory integral in a separately registered follow-up. It changes
no metric profile, amplitude, initial state, mass, finite reference, observation,
cutoff, physical stress definition, numerical acceptance gate, or independent
mode producer. Any additional numerical change needs explicit registration.

For g=4L^2h-2Lh'-h'', the required memory terms are

    F_n(eta) = integral[-5,min(eta,-3)] g^(n)(t) ln(eta-t) dt
               + c(eta) g^(n-1)(eta), n=1,2,3,
    c(eta) = ln(sqrt(2)/(-eta))+EulerGamma+1.

All g jets vanish at the source-free lower boundary. The last term follows
from the integral of g^(n); it must use the current derivative of g, not h.
At or before eta=-5 the result is exactly zero, without taking a logarithm.

For -5<eta<=-3, let ell=eta+5>0 and t=eta-ell*y. Exact endpoint subtraction
and the further substitution y=z^2 give

    I_n = ell*g^(n)(eta)*(ln ell-1)
          + 2ell integral[0,1] z
            [g^(n)(eta-ell*z^2)-g^(n)(eta)]
            [ln ell+2ln z] dz.

The endpoint integrand has a continuous value zero at z=0 and behaves as
O(z^3 ln z). Define that endpoint value before evaluating ln z. The analytic
constant term is indispensable; simply assigning zero to a divergent original
endpoint integrand does not implement this subtraction. The transformed
integrand's first and second endpoint derivatives also have limit zero.
At eta=-3 all source derivative jets are flat at the upper endpoint, so both
this expression and the post-source expression below have regular limits.

For eta>=-3 one can instead use the fixed source coordinate t=-5+2z:

    I_n = 2 integral[0,1] g^(n)(-5+2z) ln(eta+5-2z) dz.

At eta=-3, define the endpoint value at z=1 as zero before taking its logarithm;
the C-infinity flat source makes this the exact continuous extension. For
eta>-3 its logarithm has no endpoint singularity.

The proposed numerical settings are fixed panels [0,1/2,1], mpmath tanh-sinh
at 50 and 70 decimal digits, and maximum degree 10. They require a new public
freeze before any pulse, response or quadrature is evaluated. The two
precisions and endpoint transformation are implementation settings, not
permission to loosen the unchanged response or error gates.

True arbitrary precision requires evaluating the source's integer polynomial
coefficients, Horner recurrence, denominator, exponential, metric derivative
recurrence, eta, L, pi, Euler constant, scale logarithm, endpoint terms and q
contacts inside each precision context. Calling the existing double source
function at mpmath nodes would retain its double precision errors. Construct
observation eta from its exact registered decimal representation. Keep the
support endpoints and derivative coefficients exact. Cast only the final
reported float values if the existing response schema requires floats.

The exact forcing recurrence is

    g_n = 4 sum[j=0..n] binomial(n,j) (n-j+1)!
            L^(n-j+2) h_j
          - 2 sum[j=0..n] binomial(n,j) (n-j)!
            L^(n-j+1) h_(j+1) - h_(n+2), n=0..3.

In particular

    g_3=96L^5h_0+60L^4h_1+12L^3h_2
        -2L^2h_3-2Lh_4-h_5.

The same derivative polynomials and compact profiles are used; high precision
is not a different source. The convolution differentiation identity is
F_n'=F_(n+1)+L g_(n-1), which gives the unchanged canonical-memory jets

    q=-F_1/(8pi^2),
    q'=-(F_2+L g_0)/(8pi^2),
    q''=-(F_3+2L g_1+L^2g_0)/(8pi^2).

This derivative identity follows by writing the integral in the retarded
lag coordinate and differentiating its smooth compact source under the
integrable logarithm. It avoids a naive boundary term at ln(0).

Retain each tanh-sinh error estimate and the actual 50/70-digit difference in
the run evidence. Compute their difference within the higher-precision
context. Preserve enough decimal digits with an explicit serialization width;
plain string conversion outside the context may print only default precision.
Use these empirical estimates in the existing error propagation. Include
final float conversion and local-coefficient arithmetic in the stated scope
of numerical uncertainty; do not describe precision agreement as a certified
interval enclosure. Keep the frozen wall budget active during callbacks and
fail rather than silently exceed the quadrature degree, precision, or budget.

The companion pure proof verifies 34 exact identities and five symbolic
mutation witnesses, normally and with Python optimization. It checks both
endpoint transformations using arbitrary formal polynomial coefficients,
not physical profile values; the inherited integer derivative recurrence;
g_0 through g_3; and the unchanged memory derivative contacts. Its input
copies and source hashes are recorded. No physical source sampling, mode,
response value, or physical quadrature has been evaluated in this review.

A second independent review executes only the actual pure forcing recurrence
with formal symbolic jets and compares the actual AST integrand and contact
expressions with the analytic identities. Twenty source-structure/algebra
checks pass normally and with optimization for candidate module SHA256
`3f56906cad28258f988d7002ca66b8cbf4d249ea6121774591fc701393d2d8ae`.
No high-precision physical entry point is invoked. This is code preparation
evidence, not a numerical result or an empirical convergence measurement.

The final helper adds finite-value guards, separate raw quadrature diagnostics,
and upward serialization of the complete empirical qjet allowance into the
legacy memory-error field. The serializer uses the unchanged validator's exact
binary64 division and explicitly rejects failure to meet its allowance
postcondition. These allowances are not certified errors. The true raw
mp.quad estimates remain separately recorded. Exact AST comparisons confirm
that source, forcing, logarithmic-history and complete qjet expressions are
unchanged from the first reviewed helper. The earlier actual source and review
captures are preserved in `development/initial-helper-9817-review/`.
The final helper/producer/validator bindings and review scope are recorded in
`FINAL_SOURCE_REVIEW_BINDING.json`.
