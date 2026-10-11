# Independent review: uniform scalar-branch improvement

**Outcome: PASS within the stated mathematical scope.**

Reviewed proof: `JACOBI_ENDPOINT_MONOTONICITY.md`.

SHA-256: `f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8`.

Review date: 11 October 2026. This is an internal review performed by a separate
AI agent. It is not an external referee report, machine-checked formal proof,
novelty determination, or certification of cosmological evidence.

I read the full frozen proof and independently reconstructed its key estimates
and identities. The earlier analytic derivations are recorded in
`ANALYTIC_REVIEW_NOTES.md`. The independently implemented exact arithmetic
checker uses SymPy polynomial algebra and rational arithmetic; it neither
imports nor executes the certificate producer.

| Review item | Finding |
|---|---|
| General derivative scale | Green derivative normalization and all metric/scalar product factors are consistent. |
| Automatic strict bulk convexity | Follows from the contraction gate, `M3 B b <= D/2 < m^2`; no extra theorem hypothesis is needed. |
| Regular-center expansions | Leading Jacobi coefficients and error orders suffice under the stated differentiability assumptions. |
| Central variation signs | The damped cooperative system preserves `h>=1`, `h'>0`, and the specified metric sign. |
| First Jacobi identity | `j'+4aj=-2 rho^-2 v` follows from the differentiated constraint. Its integrating-factor center term vanishes. |
| Second Jacobi identity | `D'+4aD=4p(rho^-2 h+p j)` is correct. Although `D(0)` need not vanish, `rho^4 D` does. |
| Central-to-endpoint parameter | `h(T)>0` gives nondegeneracy and the correct normalization of every endpoint derivative. |
| Moving endpoint | The chain-rule formula holds along the full scalar-junction locus. Finite-endpoint `C^1` dependence precedes its use. |
| Metric residual sign | `R_T<-H^2` for nonconstant scalars and equality for the exact scalar-free branch. Only the existing positive scalar-junction derivative is needed; no shell-convexity hypothesis is added. |
| Shell uniqueness | The inherited endpoint signs plus strict decrease prove uniqueness within the stated tube and `T>=1/k`. No global uniqueness is asserted. |
| Riccati comparison | The ratio starts at zero and has a strict upper first-crossing barrier for every finite positive radius. |
| Sharper nonlinear endpoint constant | Automatic convexity gives a nonnegative metric source and `R<=1`, justifying `C_H=lambda/[k(w-3k)]`. |
| Tight scalar-root bound | `L=|f'(0)|/w` follows from the existing scalar-junction derivative bound. The zero-slope case is handled exactly. |
| Remainder refinement | Polynomial terms are bounded exactly, including the improved `min(L,|alpha|+K_eta delta_bar)` factor. No missing truncation term was found. |
| Registered corollary | Expansion coefficients, exact-decimal model specification, readable upper bounds, and separate comparator are consistent. |

The proof therefore removes the old absolute slope gate, improves the endpoint
estimate, and supports the registered exact-decimal model throughout
`0 < delta <= 1/500`, subject to the separately checked sufficient inequalities.
The documented rational computations give `K_eta < 3.196`, `K_H < 0.254`, and
a strict curvature-comparison margin greater than `0.000423` in the stated
normalization.

The producer also checks explicit polynomial bulk and shell convexity. These
are redundant sufficient diagnostics for the evaluated model; the structural
theorem does not require shell convexity. A conservative endpoint-factor variant
is a comparison of estimates, not a claim about nonconvex bulk potentials.

The comparator changes the full boundary theory by making the detuning
function scalar independent, while retaining `2W` in the tension. It is not
a constant-scalar solution of the original scalar-dependent detuning theory.

Remaining scope limits are substantive: the result does not certify the error
of individual historical floating-point trajectories, prove perturbative
stability, complete an energy-transfer mechanism, establish a distinctive
observational signature, or demonstrate a higher-dimensional origin of the
Big Bang. Original author-approved submission files were not changed by this
review.
