# Why the finite-band surrogate cannot be extended unchanged

This analysis concerns the interpolation and state recipes, not a singularity established in the underlying HDBLAST solution. The previous calculation was explicitly finite-band and remains within that scope.

Let L=ln a, H=Ldot and x=G²(phi-phi_star)². In conformal time,
\[
u''+[k^2+U]u=0,\quad u=a\chi,\quad
U=a^2[x-\dot H-2H^2].
\]
For a frequency step, continuity of u and u' gives
\[
\beta=\frac{\omega_+-\omega_-}{2\sqrt{\omega_+\omega_-}},
\quad
n=\frac{(\Delta U)^2}{4\omega_+\omega_-(\omega_++\omega_-)^2}.
\]
At large k, beta=DeltaU/(4k²) up to a phase. For several distinct interior jumps,
\[
\beta_k\sim\frac1{4k^2}\sum_j\Delta U_j e^{-2ik\eta_j}.
\]
The diagonal excitation-energy contribution is
\[
\rho_{\rm excitation}(K)
=\frac{\sum_j(\Delta U_j)^2}{32\pi^2a_f^4}\ln K+O(1)
\]
in the leading particle-energy term with an admissible incoming UV state. Distinct-knot cross terms integrate as differences of cosine integrals and cannot cancel this logarithmic coefficient. Even two opposite jumps with zero total jump retain a positive sum of squares.

C1 Hermite interpolation of L has DeltaU=-a² DeltaHdot. A C2 cubic reconstruction makes U continuous but generally gives
\[
\Delta U'_\eta=-a^3\Delta\ddot H,\qquad
\beta\sim\Delta U'_\eta/(8ik^3).
\]
Its particle energy is integrable, but its phase-sensitive pressure and fourth-order local curvature variations still require care at the knots. C4 or C6 splines improve finite regularity; they do not prove a smooth physical solution or an all-orders Hadamard state.

For a degree d spline with C^(d-1) continuity, the first generic potential derivative jump is order n=d-2, with A=-a^d DeltaL^(d). Its reflection amplitude falls as k^-(n+2). The corresponding leading diagonal omitted energy is
\[
\rho_{\rm tail}(K)\sim
\frac{\sum_j A_j^2}{2^{2n+5}\pi^2a_f^4(2n)K^{2n}},\qquad n>0.
\]
These are asymptotic coefficients, not certified finite-K error bounds. A steeper power does not determine where the asymptotic regime begins; large high derivatives matter.

## The independent initial-state issue

The prior amplitude-WKB recipe maps to
\[
u=(2w)^{-1/2},\quad
u'=(-iw-w'/2w)u,\quad w^2=k^2+a^2x.
\]
It omits the curvature term in the canonical frequency. Relative to curvature-aware UV data, the magnitude of its leading beta coefficient is |(a''/a)_0|/(4k²). If this curvature is nonzero, the unchanged infinite-band extension carries its own logarithmic excitation-energy coefficient
\[
C_{\rm initial}=\frac{[(a''/a)_0]^2}{32\pi^2a_f^4}.
\]
This is distinct from the interior jumps; do not count an initial boundary both ways.

The physical-Hamiltonian data map instead to u'=aH u-iwu. They generically have beta of order i aH/(2k) relative to an admissible curved UV reference, giving a quadratic particle-energy tail. A final instantaneous particle basis can show this tail even for an admissible state. Its spectrum alone is not an invariant diagnosis of particle production or renormalized energy.

Local vacuum counterterms cannot subtract the later history/state-dependent logarithmic excitation term while retaining a common fixed local action. Regularity, state preparation and covariant matching must be addressed separately.

## Curvature sensitivity

For a unit action integral sqrt(-g) R² with T=-2 deltaS/delta g, the regular FRW pieces are
\[
\rho_{R^2}=6R\dot H-12H\dot R,\qquad
p_{R^2}=2R\dot H+4\ddot R+8H\dot R.
\]
Here R=6(Hdot+2H²). C1 curvature jumps generate knot distributions; C2 curvature-derivative jumps still generate a distributional pressure term. The plotted regular pieces omit those distributions and do not choose a quantum matching coefficient.

## Smooth-width and UV limits

For the separate analytic toy step U=1+1.5[1+tanh(eta/epsilon)], the exact smooth spectrum approaches the sudden spectrum for epsilon*k small, but is exponentially suppressed for epsilon*k large. At large k its ratio to the sudden spectrum behaves as
\[
\frac{n_{\rm smooth}}{n_{\rm sudden}}
\sim \left[\frac{\pi\epsilon k}{\sinh(\pi\epsilon k)}\right]^2.
\]
Thus the UV integral is finite at each positive smoothing width, while taking epsilon toward zero is nonuniform:
\[
E(\epsilon)=\frac{(\Delta U)^2}{32\pi^2}\ln(1/\epsilon)+O(1).
\]
The post hoc exact-spectrum quadrature tests this asymptotic behavior for DeltaU=3. It does not determine a physical smoothing width for HDBLAST or a shell radiation budget. This is a known scattering mechanism applied to the project's numerical-readiness question, not a claim of new fundamental mathematics.


## Phase-independent archive cross-term check

Integration by parts in the defining Ci integral gives |Ci(z)| <= 2/z for z>0. Consequently the absolute interference integral over [K_L,K_U] is bounded by

    2 (1/K_L + 1/K_U) sum_{i<j} |d_i d_j| / |eta_i-eta_j|,

where d_i=Delta U_i. For the declared fine-archive band this is 6.20724e-4 of the diagonal logarithm, already below the original 1e-3 target without resolving the oscillatory phases. This analytic inequality is evaluated on floating archived coefficients and times; it is not an interval-certified physical error bound, a validation of the effective theory at arbitrarily high momentum, or a bound on subleading exact-mode corrections.
