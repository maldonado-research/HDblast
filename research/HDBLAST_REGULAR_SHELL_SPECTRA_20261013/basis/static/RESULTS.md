# Result, comparison, and remaining work

## What has been proved

In the regular two-sided Einstein–scalar shell construction, let `T` be the
interior endpoint, `beta(T)` the scalar boundary value selected by the scalar
junction, `a = rho'/rho`, and `sigma = 2W + delta f` the tension. Within the
constructed small-field family,

\[
\frac{d}{dT}\left(a(T;\beta(T))-\frac{\sigma(\beta(T))}{6}\right)
\leq -\rho(T;\beta(T))^{-2}<0.
\]

The inequality is strict on the left for a nonconstant scalar. The previously
used contraction condition already supplies strict bulk convexity on the
realized tube. Coupled scalar/metric Jacobi fields then have controlled signs.
Their exact identities retain a cancellation lost by the former separate
absolute-value estimates. Together with the inherited endpoint signs, this
proves uniqueness of the constructed junction for `T >= 1/k`. It is not a
global uniqueness theorem outside that tube.

A Riccati comparison also improves the endpoint constant to
`C_H = (w-4k)/(k(w-3k))`. The full hypotheses, eleven sufficient construction
inequalities, differentiability assumptions, and proof are in the
[technical note source](theory/JACOBI_ENDPOINT_MONOTONICITY.md).
Separate shell convexity is not a hypothesis of the structural argument.

## Registered example and explicit error

Use the dimensionless conventions of the proof:

\[
k=\frac19,\qquad w=2,\qquad W(\eta)=\frac13+\eta^2+\frac{\eta^3}{3},
\qquad f(\eta)=1+c+c\eta,
\]
\[
c=\frac{2}{1.0357712571566784}-\frac43
=\frac{96714321425692}{161839258930731},
\]

where the displayed decimal is interpreted exactly, not as a binary
floating-point input. For every `0 < delta <= 1/500`, the existence and
uniqueness conditions hold with the new norm and tube, and

\[
\left|\beta+\frac{9c}{64}\delta\right|<3.196\,\delta^2,
\]
\[
\left|H^2-\frac{1+c}{27}\delta
-\left[\frac{(1+c)^2}{36}-\frac{c^2}{384}\right]\delta^2\right|
<0.254\,\delta^3.
\]

At `delta = 0.002`, the displayed outward bounds are respectively
`1.2784e-5` and `2.032e-9`. These are analytic expansion errors, not estimates
of an ODE solver's roundoff or truncation error.

For the **separate boundary theory** with `f_metric = 1+c`, while retaining
`2W` in the tension, the scalar-free comparison branch obeys

\[
H_{\rm metric}^2=\frac{1+c}{27}\delta+\frac{(1+c)^2}{36}\delta^2.
\]

Using the exact polynomial remainder constant `K_H = 0.2533589973...`, the
uniform comparison is

\[
H^2-H_{\rm metric}^2<-0.000423\,\delta^2<0.
\]

The strict `0.000423` margin uses the exact certificate; it must not be
recomputed from the looser displayed upper bound `0.254`. This comparator
changes the full detuning function. It is not a constant-scalar solution of
the original scalar-dependent boundary theory.

| Quantity | 8 October sufficient estimate | 11 October sufficient estimate |
|---|---:|---:|
| Upper detuning endpoint | `1.0915542913e-7` | `0.002` |
| Endpoint remainder constant `C_H` | `1442/3` | `42/5` |
| Additional absolute shell-slope inequality | Required by that proof | Removed by structural sign identity |

The old result is preserved. The new proof supplies stronger estimates; it
does not show that the earlier solution ceases to exist above its old bound.
The new interval reaches all three previously sampled detunings
(`0.0005`, `0.001`, `0.002`). Distinct binary64 coefficient conventions are
checked as explicitly separate rational models, without claiming certified
numerical trajectories.

## Validation and limits

The exact producer uses `fractions.Fraction`. A separately implemented SymPy
checker reconstructs 56 rational constants and seven polynomial expressions
for each of eight variants, checks the sufficient inequalities, readable
rounding and suppression margin, and rejects six corrupted certificates.
An additional verifier checks 36 exact algebraic identities, including center
boundary terms, both scalar signs, the moving endpoint, and the Riccati barrier.
Normal and optimized Python executions are compared. Hashes identify exactly
which proof and source files were checked; hashes do not establish mathematical
truth by themselves.

The existing quadratic response coefficient and general curved-wall
Hamilton–Jacobi framework are prior work. The present project contribution is
the stated conditional argument and its controlled domain. The accompanying
literature report does not establish global originality or APS eligibility.

The next substantive target is a **validated branch-wide perturbative stability
bound**. The [27 September gauge-invariant checkpoint](https://github.com/maldonado-research/HDblast/tree/30c2ef623a1b2a322ebf245402433eca57ce24f9/research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant)
already supplies coupled equations, 55 symbolic checks, and numerical spectral
diagnostics at six detunings on the corresponding near-`phi=1` branch. The shift
`phi=1+eta` maps its potential and detuning to the present conventions. Those
floating-point diagnostics are not a rigorous bound over this full interval;
their binary64 parameter convention and the singular-looking `Psi/phi'`
variable as `delta` approaches zero need explicit treatment. Earlier stability
criteria for other branches/tensions cannot be transferred automatically.
A static solution
and a sign-definite curvature response do not prove stable propagation,
energy capture, a hot Big Bang, or an observational preference for additional
dimensions. The project's broader continuum certificate remains unresolved;
the recorded metric calibration remains failed. Neither status is repaired
by this theorem.
