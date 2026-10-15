# Final independent acceptance of the conditional first-variation waveform

15 October 2026 UTC. Internal independent AI-assisted analytic, arithmetic,
artifact and visual review. This is not external peer review, novelty
certification, or a formal proof of the retained PDE/exact-family hypotheses.
The companion JSON pins the exact reviewed bytes.

**Decision: accept the stated conditional first-variation waveform enclosure,
actual global-peak exclusion windows, and sharpened leading-response interval.
No substantive blocker remains.**

The accepted state is the same exact-decimal near-`phi=+1` branch at
`delta=1/1000`, aligned bump width `ell=1e-10`, simultaneous nonzero
displacement and velocity, and observation coordinate time `3ell` used by
the preceding quantitative checkpoint. The accepted statements are

\[
H(x)=\ell D(\ell x),\quad
\widehat H(x)=Q'(x)+\ell e^{-\ell x}Q(x),\quad
|H-\widehat H|<1.0000004\,10^{-7}\quad(0\le x\le3),
\]

\[
H=0\quad(0\le x\le1),\qquad
|H|<4.752\,10^{-14}\quad(2\le x\le3),
\qquad 6.06\,10^{15}<\mathcal C_{3\ell}<7.44\,10^{15}.
\]

The full physical scalar first variation includes the angular factor
`Y=n1^2-n2^2`. Its rescaled harmonic amplitude is
`ell S(ell x)=G cosh(ell x)^(-3/2)H(x)`, with **negative**
`G=3/(rho_b^3 p_b)`. The leading proper-time coefficient is
`C=(K/ell) integral H^2 dx`, `K=3*pi^2/(rho_b^2*p_b^2)>0`.

## Analytic and normalization review

The pointwise theorem is supplied by incoming/reflected characteristics and
the positive shell boundary kernel. It is not inferred from the preceding
`L2` estimate. The large derivative input is integrated by parts against
the damped kernel; its two nonnegative terms each lie in `[0,1/mu]`, so
their difference is bounded by `1/mu`. The positive two-variable majorant,
its determinant and inverse, and the true output `D=J-v-gamma*q_b` were
independently derived/checked. An independent error-field derivation in
this directory obtains the same sharp determinant formula and uniform
bound without importing either primary implementation.

The tail estimate retains the boundary memory after the free pulse has
ended. Its derivative-input memory is at most `ell/m^2`, its velocity
input at most `ell/m`, and the remaining potential/boundary terms are
explicit. Combining exact prefront zero, one support interval of length
one and one tail interval of length one justifies the support-aware norm
error. All characteristic arguments remain in the actual coefficient
collar; no artificial outer-boundary condition is imposed.

The retained background enclosure applies to the same actual analytic
branch. The signed gain and norm prefactor use the exact shell constraint
and scalar junction, preserving their correlation through the shell
scalar parameter. Their derivative signs, endpoint extrema and polynomial
identities were checked. The time-conjugation range
`1-t^2<=cosh(t)^(-3/2)<=1` is safe here. The norm Jacobian is `1/ell`;
the proper-time weight cancels exactly as in the retained response proof.

The nonlinear correction annulus and evolution interval are explicitly
distinguished from the linear energy domain. The displayed strict margins
and qualitative nearby-speed bound are consistent with the retained
parameter/exact-family hypotheses. They do not produce an effective
nonlinear parameter interval.

## Continuous exported cells and arithmetic

The production 96-bit interval code was reviewed: rational floor/ceil
rounding, signed products, squares crossing zero, exponential range
reduction, Taylor remainder, repeated squaring and the large-argument tail
are outward. Continuous reciprocal ranges and the monotonicity of
`v^2 exp(4-v)` give full real-cell bounds, including flat support edges.
The independent reciprocal positive-Taylor checks are a distinct
exponential implementation, rather than a call to the production code.

An independent checker reconstructed all **4098** exported CSV cells,
including exact-zero `[0,1]`, all 4096 support cells, and the nonzero tail
envelope `[2,3]`. It verified coverage, integer endpoint order, global-error
widening, signed physical gain and clock multiplication, the complete
quadrature sum, an independent 128-bit norm recheck, and both response
coefficient endpoints. The reported interval lies strictly inside the
published outward interval above. Source hashes and gain-input hashes
match the exported receipt.

The referee's primary checker replay passed 44 groups per normal/optimized
isolated Python mode, with identical parsed results. The independent
characteristic audit passed 20 groups per mode and the independent
CSV/norm audit passed 23 per mode. These are finite analytic-bookkeeping
and arithmetic checks, not a machine-checked PDE proof.

## Peak claims and presentation

The comparator's peak and zero geometry, interval Taylor/square-root
arithmetic, root sign brackets and monotonicity argument were reviewed.
Both endpoint height gaps exceed `2*(1.000001e-7)+ell`. That value gap,
together with the uniform error, proves that every global maximizer of
**H** lies in `(1.3033,1.3034)` and every global minimizer in
`(1.6966,1.6967)`. The maximum/minimum height intervals transfer by the
supremum norm. A read-only replay matched the saved peak certificate in
both modes and rejected all four built-in wrong-formula controls.

This does not transfer the comparator's precise critical-point locations,
number of zeros, or uniqueness of local/global extrema to `H`, and does
not identify which lobe dominates the absolute maximum. The primary note
correctly preserves these limits and labels the peak windows as those of
`H`, whose physical scalar orientation is reversed by the negative gain.

The rendered PNG was inspected. It correctly shows the physical sign
orientation, conditional first-variation label, model-coordinate axis,
midpoint-gain illustration, full continuous-cell band and unresolved small
tail band. Floating-point plotting is a display step; the CSV and rational
JSON contain the certificate data.

## Scope retained

This is a conditional resolved waveform approximation for the linear first
variation of one very narrow normalized state. It is not an exact nonlinear
waveform, finite-epsilon prediction, effective nonlinear Taylor remainder,
physical-unit time calibration, gravitational energy measurement, global
stability result, or cosmological-origin test. High-precision comparator
geometry and large normalized coefficients do not establish any of those
claims. The accepted scientific next target is a wider practical-scale
state with validated coupled reconstruction, followed by an effective
nonlinear remainder estimate.
