# Independent prospective review

This file registers an independent reviewer implementation and acceptance logic.
No field evolution or new numerical research result was executed while preparing
it. The implementation imports no primary numerical or subtraction source.

## Physical modes

Use r=4, T=1, A in {0,0.2}, and k in {0,0.5,2,8,24,48}. Start the exact static
incoming vacuum at eta=0. Solve the four-real-component equation with SciPy Radau,
its analytic Jacobian, rtol=2e-12, atol=2e-14, and max_step=min(1/100,0.1/w_future).
Evaluate eta={0,.125,.25,.375,.5,.625,.75,.875,1,1.25}. Evolve to eta=1 and use
the exact constant-frequency propagator thereafter. No Wronskian rescaling is
permitted. Compare primary and independent raw u,u',rho,p,Q at identical nodes.

The independently implemented compact step uses g=-1/s+1/(1-s), B=logistic(g)
and analytic g',g''. Computing B(1-B) naively fails near the future endpoint once
B rounds to 1; use exp(-abs(g))/(1+exp(-abs(g)))^2 for that product. Endpoint
branches set exact static values. This supplies the actual background curvature,
not finite differences or a spline fit.

The independent solver clips only min(s,1-s)<1e-6, where its two derivative
magnitudes admit conservative bounds |B_s|<2*10^12 exp(-999998) and
|B_ss|<8*10^24 exp(-999998). Both are below the smallest float64 magnitude, as
is the profile departure from the exact static branch. This is a stated floating
point representation bound, not an approximation to a resolved transition.

Predeclared comparisons: |Delta u| and |Delta u'| <= 2e-8; observable differences
<= 2e-7+2e-8|primary bare observable|. Raw internal Wronskian error <1e-8.
Occupation comparisons require absolute errors as well as relative errors; do not
reject harmless relative disagreement when |beta| is below the solver's absolute
error floor. Parent integration/protocol may impose stricter gates.

## Exact identities and negative controls

With z=u'-Hcal u and q=k^2+a^2x, z'=-q u-Hcal z. Hence the bare mode identity is
rho'+3Hcal(rho+p)-x'Q/2=0 for rho=(|z|^2+q|u|^2)/(2a^4),
p=(|z|^2-(k^2/3+a^2x)|u|^2)/(2a^4), Q=|u|^2/a^2.
The implementation differentiates the unreduced quadratic expression using the
physical equation. Its algebraic residual has a 1e-10 mixed absolute/scale gate;
this is an algebra and normalization check, not an independent conservation
accuracy estimate for the field evolution.

Dropping the current term must give a nonzero residual at some registered
interior node of every selected case, exceeding 1e-5. Dropping pressure must give
a residual exceeding 1e-5 for A=.2; this control is inapplicable at A=0 because
Hcal=0. The corresponding primary counterterm control should drop p4 or Q2,
which directly checks the subtraction pairing.

Incoming static vacuum observables minus analytic static counterterms must vanish
to 1e-12(1+k). In the exactly static future, with phase origin eta=T, verify

rho_ren = w |beta|^2/a^4,
p_ren = [k^2|beta|^2/3-(2k^2/3+a^2r)Re(alpha beta* exp(-2iw(eta-T)))]/(a^4 w),
Q_ren = [|beta|^2+Re(alpha beta* exp(-2iw(eta-T)))]/(a^2 w).

Compare these to direct bare-minus-static subtraction at both future nodes with
2e-10(1+|bare observable|) gate. Store the diagonal-only expressions separately
so coherent pressure/current cannot be silently replaced by occupations.

## Cutoff and cancellation review

A fixed COMOVING cutoff preserves an exact modewise exchange identity upon
integration. Its residual is therefore not evidence that the cutoff integral
has converged, and not a proof of local covariance. A changing physical cutoff
introduces a moving-boundary term that must be included in the ledger.

On smooth generic interior profiles, after stress order4/Q order2 subtraction,
remaining radial integrands are generically O(k^-3), with O(K^-2) cutoff errors.
The exact static future cancels these local remainder terms, but pressure and Q
retain phase-sensitive oscillatory interference. Momentum quadrature must resolve
the oscillation phase, independently of the small occupation spectrum.

Test the registered K=24,48,96 ladder and conditional K=192; refine momentum
quadrature and ODE tolerances at each relevant cutoff. Require shell increments
to shrink and to be stable under those independent refinements before claiming
a tail estimate. Do not assume 1/K^2 asymptotics before the shell behavior supports
it. If higher K worsens solver/roundoff sensitivity, report a numerical ceiling.

Longdouble arithmetic in the final subtraction reduces subtraction roundoff but
cannot restore double-precision ODE information. A normalization error epsilon
can produce integrated contamination scaling as epsilon K^4. Preserve raw
Wronskians and solve again more accurately; do not renormalize modes post hoc.
An absolute error floor must accompany any relative residual near zero crossings.

The full subtraction needs derivatives through a'''' and x'' for pressure. The
energy density has lower highest derivatives, so the exchange identity can use
the same derivative budget if differentiated after the complete cancellations.
High derivatives must come from analytic expressions or controlled jets. Endpoints must be checked
using stable exact-static branches, not by extrapolating near-singular formulas.

## Finite covariant action review

Use signature -+++ and R=6(Hdot+2H^2), T_munu=-2 delta S/delta g^munu. For the
declared local action density -V(x)+F(x)R+alpha R^2, independently varied pieces are

rho = V-6F H^2-6H Fdot+alpha(36Hdot^2-216H^2Hdot-72H Hddot),
p = -V+2F(2Hdot+3H^2)+2Fddot+4H Fdot
    +alpha(108Hdot^2+216H^2Hdot+144H Hddot+24H'''),
Q = 2V_x-2F_x R.

All derivatives in this paragraph are cosmic-time derivatives. These formulas
are regular at H=0 and should be checked directly there; deriving pressure by
division by H is insufficient. A Weyl-squared term has vanishing first variation
on FRW; the Euler term is topological and total derivatives require boundary
conventions. FRW data therefore cannot measure or prove the general-spacetime
Weyl/Euler coefficients. A covariant heat-kernel bridge can supply that extension.

PV matching must name regulator weights, verify cancellation of sum c_j,
sum c_j r_j, sum c_j r_j^2, distinguish physical masses x+j Lambda^2 from reference
masses r+j Lambda^2, and take the combined momentum integral before removing the
regulators. The proposed alpha_PV=-sum c_j ln(r_j/mu^2)/(2304 pi^2) has the correct
sign for this Lorentzian action convention. Heavy-field decoupling should remain
a qualified asymptotic statement requiring bounded smooth derivatives and the
specified incoming vacuum, unless a uniform remainder bound is demonstrated.

No prescription removes the physical freedom to add finite local gravitational
or x-dependent action terms. Declaring their coefficients relative to a proved
base scheme is a legitimate benchmark; it is not a renormalization-choice-free
source. Passing this control supplies neither a self-consistent HDBLAST solution
nor backreaction, decay, thermalization, or observational evidence.
