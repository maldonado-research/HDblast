# Infrared corollary for the registered analytic prehistory target

The two completed registered routes enclose the same zero-anchored analytic
BD prehistory target at `a=-9/2`. Their saved, complete `k=0` rectangles have
exactly zero imaginary components and imply the following outward coarse
rational bounds. Each displayed interval contains both route intervals;
no midpoint averaging is used.

| Source | Exact rational enclosure, written as terminating decimals, for W(a,0) |
|---|---|
| `positive_B` | `[1.550744117255, 1.550744117256]` |
| `signed_uB` | `[-0.053498425349, -0.053498425348]` |

These are post-result deductions from the saved JSON, conditional on the
registered numerical proofs and runtime premises. The public source freeze
was `a253f50158051881f18bbdec81c12dbbb4f49eba`, with full registration SHA256
`1d55b6310631fcf600a8165ec1d326b9880a0de346b531c9d9659336778f765a`.
The independent review verifies the saved execution and entry receipts,
every output byte pin, and the unchanged frozen files. It reproduces all
36 route-rectangle intersections, all 72 component intersections, all 666
independent error-trace rows, and the continuous source/coefficient
certificate using JSON and exact rational arithmetic only.

For the audited real forcing,

\[
 W(a,k)=-\int_{-5}^{a}e^{2ik(a-s)}g(s)\,ds,
 \qquad
 U(a,k)=-\int_{-5}^{a}\Phi_k(a-s)g(s)\,ds.
\]

Compact integrable forcing makes these functions entire in complex `k`.
Here `Phi_0(tau)=tau`; both zero-frequency values are removable limits.
In particular `W(a,k)=W(a,0)+O(k)`. Because each certified `W(a,0)` is
nonzero, the exact normalized amplitude `A(k)=W(a,k)/(2ik)` has a simple
pole:

\[
 \lim_{k\to0^+}kA(k)=\frac{W(a,0)}{2i}
                       =-\frac{i}{2}W(a,0)\ne0.
\]

The residue is negative imaginary for `positive_B` and positive imaginary
for `signed_uB`. The full rational residue enclosures are recorded in
`ACTUAL_RESULTS_INDEPENDENT_REVIEW_NORMAL.json`. A uniform bounded
`|A|` hypothesis would therefore be false for these exact targets. The
weighted finite-band theorem allows a `1/k` envelope and remains applicable;
this pole is compatible with its finite linear stress integrals.

For each `k>0`, decompose the first-order variation pair at the matching time
in the free basis `v0=e^(-ika)/sqrt(2k)` and its conjugate. With the audited
definitions `delta_v=epsilon*v0*U` and
`delta_v_prime=epsilon*v0*(W-ikU)`, elementary two-by-two inversion gives

\[
 \delta\beta=\epsilon e^{-2ika}A,
 \qquad
 \delta\alpha=\epsilon(U-A).
\]

These are matching-time linear coefficients; the forcing is still active
at `a`, so no source-free future evolution is being assumed for the whole
physical variation. The represented epsilon is the nonzero rational
`3777893186295716171/37778931862957161709568`. Therefore

\[
 \lim_{k\to0^+}k|\delta\beta|
 =\lim_{k\to0^+}k|\delta\alpha|
 =\frac{|\epsilon|\,|W(a,0)|}{2}>0.
\]

By continuity, `|W(a,k)| >= |W(a,0)|/2` for sufficiently small positive
`k`, which already proves `|delta_beta| >= |epsilon W(a,0)|/(4k)` there.
Thus at fixed nonzero epsilon these first-order coefficients are not
uniformly small deformations of the unperturbed pair `(alpha,beta)=(1,0)`
as `k` approaches zero. Their opposite poles cancel in
`delta_alpha+exp(2ika)*delta_beta=epsilon*U`, so this statement does not
assert a pole in the relative field variation `delta_v/v0`. It also does
not, by itself, estimate a higher-order epsilon remainder or prove that a
linearized stress approximation fails.

The `k=0` numerical entries are analytic limits of normalized `U,W`;
they do not supply a finite physical plane-wave mode at exactly zero
momentum, since its `1/sqrt(2k)` normalization is singular. The identity
`Re(U)-Im(W)/(2k)=0` is a first-order canonical condition. It does not
upgrade the linear pair into an exactly normalized finite-epsilon state.
No particle occupation, heating, nonlinear dynamics, UV tail, or
cosmological origin follows from the residue calculation.

The source/coefficient certificate concerns the exact Duhamel solution of
the single code-identified real polynomial family for all momenta through
`K=256`. Numerical phase/ODE evaluation is separately enclosed at the nine
registered probes. Neither result compares the original stored binary80
incoming arrays to their exact target. That error remains `NOT_ENCLOSED`,
the full twelve-case certificate remains `UNRESOLVED`, metric calibration
remains `FAIL`, and the proposed Big Bang cause remains `NOT_ESTABLISHED`.
External novelty is `NOT_ASSESSED`.
