# Endpoint-inertia and Volterra methods: bounded primary-source comparison

15 October 2026 UTC. The target is the retained halfline equation
\(w_{tt}-w_{zz}+W(t,z)w=0\), with
\(-w_z+\mu w_{tt}+\beta(t)w_b=0\), \(\mu>0\), and physical output
\(D=w_z-\gamma w_b\). The proposed state has \(\ell=10^{-6}\) and
\(T=3\ell\) in model coordinates.

This packet records three discovery attempts: two successful paper-index
searches and one web search rejected by Firecrawl with HTTP 402 before it ran.
Four targeted primary passage reads succeeded. The failed web search was not
retried. This is a bounded review, not a whole-web search or a novelty claim.
The previously reviewed dynamic-boundary and Arb papers were not reread.

## 1. Causal boundary memory

[*Green's Function For Linear Differential Operators In One Variable*,
arXiv:1301.1004](https://arxiv.org/abs/1301.1004), §2 equation (11), constructs
causal Green functions through iterated integrals. Corollary 2, equations
(19)–(20), gives constant-coefficient kernels through roots and exponential
convolutions. The returned passages support the causal-kernel method. They do
not supply a verified numerical remainder for the present endpoint-mass PDE.

An elementary comparison construction useful here can be obtained directly,
without importing an external existence theorem. Let \(f=q_\ell\),
\(g=e^{-t}f\), and \(j_0=f'+g\). For a free comparison problem with
\(W=0\), \(\beta=0\), and the **same positive mass** \(\mu\), put

\[
v_I(t)=\frac1\mu\int_0^t e^{-(t-s)/\mu}j_0(s)\,ds,\qquad
b_I(t)=\int_0^t[1-e^{-(t-s)/\mu}]j_0(s)\,ds.
\]

These solve \(\mu v_I'+v_I=j_0\), \(b_I'=v_I\), with zero initial
boundary data. If \(u_D\) is the free Dirichlet solution for the given
initial data, then
\(w_I(t,z)=u_D(t,z)+\mathbf1_{t\ge z}b_I(t-z)\) is the corresponding
free endpoint-inertia comparison, and
\(D_I=j_0-v_I-\gamma b_I\). Smooth compatible corner data justify the
boundary lift. This is a direct model calculation, not a claim that the
physical \(W,\beta\) vanish: their differences must remain in the error
equations.

Integrating only the **known** derivative input gives

\[
v_I(t)=\frac{f(t)}\mu+
\frac1\mu\int_0^t e^{-(t-s)/\mu}
\left[g(s)-\frac{f(s)}\mu\right]ds,
\]
\[
b_I(t)=\frac1\mu\int_0^t e^{-(t-s)/\mu}f(s)ds+
\int_0^t[1-e^{-(t-s)/\mu}]g(s)ds.
\]

For \(t\ge2\ell\), the point term \(f(t)\) is zero and only integrals
over the initial support remain. The memory is generally nonzero; a zero
Dirichlet tail does not imply a zero endpoint-inertia tail. These exact
identities expose cancellation before absolute values are taken and avoid
estimating the convolution merely by \(\|f'\|_\infty=O(\ell^{-1})\).
They require their own sign, normalization and arithmetic checks in HDB.

## 2. Iterated paths and controlled truncation

[*A path integral method for solution of the wave equation with
continuously-varying coefficients*, arXiv:1901.04158](https://arxiv.org/abs/1901.04158),
§4.1 Lemmas 3–4 and Theorem 3, proves convergence of characteristic-path
integral series for its spatially varying acoustic medium. Time-ordered path
volumes give factorial denominators; the retrieved text also controls
derivative series and gives finite truncation estimates. This supports
keeping the causal integration geometry when estimating an iterated series.
The bound can be conservative, as the paper itself observes.

Its sharper alternating remainder in Theorem 4, equations (47)–(48), assumes
monotone impedance, restricted impedance contrast and particular step data.
Those conditions are not established for the HDB potential/boundary problem.
We therefore import neither that alternating-tail bound nor its numerical
constants. A Born/Volterra expansion for HDB must use its own coupled bulk
and boundary operator, including the output trace. A small bulk operator norm
does not alone bound a derivative observation.

The odd-reflection identity is particularly useful for the Dirichlet source
component. For a positive-side source \(F\), with odd extension \(F_o\),

\[
\frac12\int_0^t\!\int_{z-(t-s)}^{z+(t-s)}F_o(s,\xi)\,d\xi\,ds
=\frac12\int_0^t\!\int_{|z-(t-s)|}^{z+t-s}F(s,\xi)\,d\xi\,ds.
\]

The remaining spatial width is \(2\min\{z,t-s\}\), rather than the
full unreduced cone width. Preserve this cancellation, and the actual pulse
support where justified, before bounding the source. The physical solution
also has a boundary-memory lift, so this identity does not assert a vanishing
post-pulse field or tail. The actual HDB remainder proof must account for
both pieces and every retained coefficient.

## 3. What nonlinear Einstein boundary theory does and does not provide

[Kreiss, Reula, Sarbach and Winicour, *Well-posed initial-boundary value
problem for the harmonic Einstein equations using energy estimates*,
arXiv:0707.4188](https://arxiv.org/abs/0707.4188), derives bulk and boundary
energy estimates for a constraint-preserving harmonic boundary hierarchy.
Sections 2–3 treat the frozen-coefficient form; the introduction explicitly
points to Appendix 1 for the local quasilinear extension. The displayed scalar
prototype has a first-order Sommerfeld-type boundary law with a strict
coefficient inequality. This is not a dynamic Israel junction condition.

[*The initial boundary value problem for the Einstein equations with totally
geodesic timelike boundary*, arXiv:2006.01498](https://arxiv.org/abs/2006.01498),
uses the special condition that the boundary second fundamental form vanishes.
The abstract identifies loss of normal-derivative control and an anisotropic
Sobolev framework. Retrieved Proposition 4.8 distinguishes \(H^s\), \(s\ge3\),
for the Cauchy case from its \(B^s\), \(s\ge7\), boundary regularity; these
are different spaces, not interchangeable derivative counts. The retrieved
proof closes the constraint/curvature energy using boundary-flux cancellation.
Its equations are the Einstein **vacuum** equations. Its totally-geodesic
boundary cannot be substituted for a scalar-bearing Israel shell.

Neither read supplies an effective second parameter derivative or nonlinear
Taylor-remainder bound for this Einstein–scalar shell model. That gap remains
a separate scientific target. For example, the general Banach-space estimate

\[
\|U_\varepsilon-U_0-\varepsilon U'_0\|_X
\le\tfrac12\varepsilon^2
\sup_{|s|\le|\varepsilon|}\|U''_s\|_X
\]

becomes effective only after the norm, a valid parameter interval and the
right-hand bound are quantitatively controlled for the actual equations,
constraints and boundary conditions. Scalar matter, the shell embedding,
intrinsic clock, measure and physical pullback all enter the corresponding
observable estimate. Qualitative smooth dependence or a verified linear
waveform does not provide those numerical constants.

## 4. Recent discovery with limited evidence

The second search returned [*Optimal error bounds on the exponential
integrator for dispersive equations with highly concentrated potential*,
arXiv:2602.22068](https://arxiv.org/abs/2602.22068). Only its index abstract was
read. It reports cancellation-sensitive iterated-Duhamel analysis for a
one-dimensional **dispersive** model and an exponential integrator. Its body,
authors and applicability to a hyperbolic endpoint-mass boundary were not
verified. It is a methodological lead, not an imported error bound. No new
quantum, string-theory or observational connection follows from this search.

## Provenance and limits

[INERTIAL_WAVEFORM_SOURCE_LEDGER.json](INERTIAL_WAVEFORM_SOURCE_LEDGER.json)
records locators, evidence levels and transfer limits.
[BOUNDED_SEARCH_RECEIPT.json](BOUNDED_SEARCH_RECEIPT.json) records exact tool
arguments, result counts, the failed request and separate internal result-file
hashes. The raw internal responses contain copyrighted passages and must be
excluded from publication. No downloaded code was executed, and no external
messages or Git writes were made.

The endpoint-mass-specific web discovery remains incomplete because that
request did not run. No source identified in this bounded review is represented
as a complete certificate for the proposed HDB first variation or its nonlinear
completion. The HDB proofs and computations require independent verification.
