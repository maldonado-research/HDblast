# Independent shell-gradient geometry and exact-family join

Read-only source repository HEAD
`b51908d7e62e7ec80cef82438edd00923ff94723`; merged commit
`c2293ee83945365f24547bd83dbf828478e99ab5`.

## What is inherited and what this check adds

The observable Q and its pointwise/mean second derivative are **already in**
`research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/shell-invariant/SHELL_INTRINSIC_SCALAR.md`,
lines 9–68, as I. Its old report also supplies two pairs of mean zeros, in
[5/4,11/8] and [13/8,7/4]. No novelty claim for that formula, standard gauge
machinery, scalar-gradient interpretation or original numerical mean zeros
is justified. The useful new joins are (i) conditional exact-family and
completion-independent identification, (ii) a specified intrinsically
anchored shell proper clock/measure, and (iii) the angular measure formula
independently derived in `ANGULAR_FRACTION.md`.

The exact-family source has common domain `0<=x<=3/2,0<=r<=6`, not through
x=7/4. Thus only the first old bracket lies inside that source-bound exact
family. A fresh analytic domain argument would be needed to transfer its
second bracket to the selected exact family.

## Completion-independent quadratic derivative

Let E_e embed the same physical shell; set h_e=E_e^*g_e (Lorentzian),
Phi_e=E_e^*phi_e. Assume e -> (h_e,Phi_e) is C2 in a topology controlling
h_e and dPhi_e, and Phi_0 is constant on the **whole shell**. Write
Phi_e=Phi_0+e eta+e² omega/2+o(e²) in C1. Then

    Q_e=h_e^{-1}(dPhi_e,dPhi_e),
    Q_0=Q'_0=0,
    Q''_0=2h_0^{-1}(deta,deta).

Only the zeroth metric and first physical-shell scalar pullback occur.
This statement holds for every sufficiently regular kinematical family
with the same h_0 and eta, and hence conditionally for every exact solution
family having those data. It requires no matching of the second scalar,
second metric, nonlinear constraint correction or TT initial pair. This
contrasts with the source-bound TT completion counterexample; that
counterexample does not obstruct this observable.

For any regular event/clock pullback A_e with a common zeroth map A_0,

    (A_e^*Q_e)''_0=A_0^*Q''_0.

Indeed Q_e=e²Q''_0/2+o(e²), so first/second variations of A_e have no
second-order contribution. Likewise for a regular finite measure mu_e
with common mu_0 and positive total mass,

    [ integral Q_e dmu_e / integral dmu_e ]''_0
       = integral Q''_0 dmu_0 / integral dmu_0.

The background event, foliation, observer and measure are part of the
observable. Their first and second variations disappear; replacing their
zeroth data is not an invariance statement.

For a regular future unit shell vector u_e, define

    K_e=(u_e Phi_e)^2,
    S_e=(h_e^{-1}+u_e tensor u_e)(dPhi_e,dPhi_e),
    Q_e=S_e-K_e.

Then K''_0=2(u_0 eta)^2 and S''_0=Q''_0+K''_0. The same event/measure
argument applies to K and S. These formulas describe scalar-gradient
causal character and its observer split, not a shell stress tensor,
matter energy, reheating or cosmological origin mechanism.

## Selected background and exact moments

Use h_0=rho_b²[-dt²+cosh²(t) kappa], eta=psi(t,0)Y,
Y=n0²-n1², the anchored background observer u_0=rho_b^{-1}partial_t,
and round probability measure dvol_kappa/(2pi²). Then

    Q''_0=(2/rho_b²)[-psi_t²Y²+psi²|DY|²/cosh²(t)],
    K''_0=2psi_t²Y²/rho_b²,
    S''_0=2psi²|DY|²/[rho_b²cosh²(t)].

Here psi_t is ordinary unscaled t derivative; x=t/ell implies
psi_t=psi_x/ell. Probability moments are <Y²>=1/6, <|DY|²>=4/3, so

    <Q''_0>=[-psi_t²+8psi²/cosh²(t)]/(3rho_b²),
    <K''_0>=psi_t²/(3rho_b²),
    <S''_0>=8psi²/[3rho_b²cosh²(t)].

These are second derivatives; the coefficient of e² is one half.

## Physical shell and anchor join

The exact-family source, sections 3–4, uses a common chart and **fixed
physical shell z=0**, identifies the full actual first solution by the
linearized bulk/boundary system and complete initial one-jet, and retains
f|shell=0. In a moving-shell description eta=chi_b+p_b zeta; under a plus
pullback chi_b -> chi_b+p_b f and zeta -> zeta-f. The terms cancel.
Consequently eta=psi(t,0)Y in that fixed-shell representative. A physical
shell displacement added without transforming the scalar would change eta.
Tangential relabelings cannot alter eta because Phi_0 is constant.
The ambient norm g^{-1}(dphi,dphi) is a different observable: its background
normal derivative p_b is nonzero, and the quadratic cancellation above
does not apply to it.

For the shell proper-time prescription, use the initial **physical shell
section** t=0 and intrinsic h_e-normal geodesics, as specified in the old
`HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/observable/SHELL_PROPER_TIME_OBSERVABLE.md`.
These are shell geodesics; they need not be ambient spacetime geodesics.
The background map is explicitly F_0(tau,Omega)=(tau/rho_b,Omega), with
metric -d tau²+rho_b²cosh²(tau/rho_b)kappa. On any compact history with
0<=tau/rho_b<=t_max<3ell/2 it is an embedding with strictly positive
spatial Gram matrix. C2 parameter dependence in spacetime C4, retained by
the exact-family source, gives geodesic and label dependence and stability
of the embedding on such a compact domain for a common sufficiently small
parameter interval. At t=0 use one-sided forward dependence. Compactness,
background injectivity, and local inverse-map stability give injectivity
after decreasing the parameter interval. This is qualitative; it supplies
no effective parameter radius.

The selected exact initial correction support is ell/2<z<5ell/2, and the
first seed support is ell<=z<=2ell. All exact canonical initial fields
therefore equal the background on an open shell collar. The prescribed
lapse N_e=rho(1+e uY) and zero shift also equal the background there.
The initial shell section and its normal are consequently fixed physical
anchor data. Under a gauge change transport the shell and anchor jointly.
The old volume observable depends on second metric/clock data because
its nonzero background and first differential permit them; Q,K,S all
vanish to first order, which is why their present join is simpler.

For either fixed initial round-label averaging or normalized physical
slice volume, mu_0 is round at a fixed tau. The second derivatives above
are identical. This does not identify the finite-e means or their higher
coefficients.

## Necessary qualifications and exact controls

A different zeroth observer changes K and S: in a flat shell with deta=dt,
u_0=(1,0,0,0) gives K''=2, whereas the unit boost (5/3,4/3,0,0) gives
K''=50/9. A different zeroth angular measure changes the moments: normalized
weight 1+Y² gives <Y²>_weighted=29/140 rather than 1/6. Thus “any clock or
measure” must mean variations with the specified common zeroth data.

There is no common finite parameter radius based only on the first
pullback: phi_e=e t+e²c x on a fixed flat shell has Q=-e²+e⁴c² and
Q''_0=-2 for every c. For any proposed nonzero radius, choosing c large
enough reverses Q there. For each fixed C2 family a strictly signed
quadratic coefficient does imply its sign sufficiently near e=0, without
an effective or family-uniform radius. No all-angle negative sign follows:
if both psi and psi_t are nonzero, explicit positive and negative angles
exist, and Z=0 always gives Q''=0.

`check_geometry.py` independently checks generic symmetric first/second
metric jets, scalar second jets, observer variations, event and normalized
measure variations, exact S3 moments and these necessary-hypothesis
counterexamples. It imports no source or producer code and evaluates no
actual first oracle. The separate numeric agent is responsible for actual
psi and ordinary t-derivative enclosures.

## Strict angular fraction and qualitative finite-family consequence

Set R=<S''_0>/<K''_0>=8psi²/[cosh²(t)psi_t²] when psi_t!=0.
The angular threshold in the independent derivation is r=R/(R+2).
If a frozen actual-jet certificate supplies R<=10^-12, then exactly

    r < 1/(2*10^12) < [3/(4*10^6)]²,
    1-F(r) <= (4/pi)sqrt(r) < (4/3)sqrt(r) < 10^-6.

The r=0 case has F=1 and obeys the conclusion directly. Thus more than
99.9999% of round background angles have a strictly timelike quadratic
scalar gradient. This is a corollary of a frozen numerical target; the
geometric derivation itself does not certify the actual R bound.

There is also a qualitative consequence for each retained exact C2 family.
The strict fraction gap and inner regularity yield a compact angular set
A contained in {Q''_0<0} with mu_0(A)>1-10^-6. Continuity gives a strict
negative coefficient margin on A. Uniform C2 Taylor convergence of the
proper-clock pullback gives Q_e/e² -> Q''_0/2 uniformly on A. Hence Q_e<0
on A for all sufficiently small nonzero e. The normalized physical slice
volume has a positive density converging uniformly to its round background
density, so mu_e(A)>1-10^-6 after another parameter shrinkage. This proves
an actual proper-volume fraction greater than 99.9999% for sufficiently
small nonzero e, with **no effective or family-uniform parameter radius**.
The same argument works for the fixed initial-label measure.

For uniformity on a compact time bracket, require the continuous R bound,
continuous coefficients, and |psi_t| bounded away from zero throughout.
For example choose one rational threshold

    r_c = 17/(32*10^12),
    1/(2*10^12) < r_c < 9/(16*10^12).

The compact angular set {Z cos²(2phi)>=r_c}, understood by continuity as
{Y²>=r_c Z, Z>=r_c}, lies strictly inside every negative coefficient set,
and its round mass exceeds 1-10^-6. The amplitude lower bound and threshold
gap give a common negative coefficient margin on the whole time product.
Uniform Taylor and measure convergence then give one qualitative parameter
interval for the complete bracket. Without the continuous hypotheses only
the certified individual times inherit this assertion. The shell-clock
history must still fit the retained x<3/2 family domain.
