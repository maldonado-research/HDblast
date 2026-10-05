# Finite-band stress kernels and a composable conservation-error theorem

Analytic continuation of the verified source/operator checkpoint, prepared
5 October 2026 UTC. This document evaluates no source callback, saved array,
physical trajectory or numerical source integral. The results below are exact
linear-model identities and conditional error bounds. They do not complete the
twelve-case numerical pressure/contact certificate. The method is established
linear ODE, action-variation and integral-inequality mathematics; external
novelty is not assessed.

## 1. Fixed model, normalization and inherited premises

Use the unchanged interior interval `a=-9/2`, `b=-7/2`, background
`L=-1/eta`, `L'=L^2`, `H=1`, physical mass squared and fixed positive
subtraction reference `x=r=2`, and minimal coupling. The metric direction is
`a_metric=a0(1+epsilon h)` with `a0=L` and real

\[
g=4L^2h-2Lh'-h''.
\]

The two inherited sources are `h=B` and `h=(eta+4)B`, where
`B=exp(-(eta+4)^2/(1-(eta+4)^2))` inside its support. Their definitions,
analytic bounds and degree-24 Taylor-tail proof are inherited pinned inputs,
not reevaluated in this work.

Define normalized amplitudes `U=u/epsilon`, `W=w/epsilon`; hence exact flow
is `U'=W`, `W'=2ikW-g`. Write `X=Re U`, `Z=Re W`, `T=Im W`. At `k>0` the
independently action-derived mode stresses, in the inherited `a0^4/epsilon`
normalization, are

\[
R_m=\frac{(2k^2+3L^2)X-kT-LZ}{2k},\qquad
P_m=\frac{(2k^2/3-L^2)X-kT-LZ}{2k}.
\]

Use one fixed positive measure constant `Pi` throughout:

\[
d\mu=\frac{k^2\,dk}{2\mathrm{Pi}^2},\quad
M_K=\int_0^K\frac{d\mu}{2k}=\frac{K^2}{8\mathrm{Pi}^2}.
\]

`Pi` may be the exact represented producer constant
`14488038916154245685/4611686018427387904`; it is not silently replaced by a
different rounded value or by a newly rounded mathematical pi. The identities
hold for either declared convention if the same constant is used everywhere.
The rational upper bounds below need only `Pi>=3`.

Full stresses are `R=R_m+C_R`, `P=P_m+C_P`, with the full action-derived
metric, mass, volume and fixed-reference W0/W2/W4 contact inventory. Let
`B0=R0+P0` be the matching per-mode baseline. The inherited arbitrary-jet
inventory proves

\[
C_R'=L(C_R-3C_P)-3h'B_0-\frac{Lg}{2k}.                 \tag{1}
\]

Pressure and contacts must continue to be computed directly from that
inventory. Equation (1) is a verification premise, not a definition of them.
This study does not independently rederive every WKB coefficient.

## 2. Full off-shell identity, including the invariant omitted by G alone

Let `f` be a real reconstructed forcing and define normalized residuals

\[
r_u=\widehat U'-\widehat W,\qquad
r_w=\widehat W'-2ik\widehat W+f.
\]

For arbitrary contact approximations define their independent defect

\[
\zeta_C=\widehat C_R'-L(\widehat C_R-3\widehat C_P)
       +3h'\widehat B_0+\frac{Lg}{2k}.
\]

Then, pointwise at `k>0`, the complete reconstructed direct stress obeys

\[
\widehat R'-\widehat F
=\frac{(2k^2+3L^2)\operatorname{Re}r_u
           -k\operatorname{Im}r_w-L\operatorname{Re}r_w}{2k}
 +\frac{L(f-g)}{2k}+\zeta_C,                          \tag{2}
\]
\[
\widehat F=L(\widehat R-3\widehat P)-3h'\widehat B_0.
\]

For raw amplitudes instead, use `e_u=u'-w`,
`e_w=w'-2ik w+epsilon f`; the first numerator in (2) is divided by
`epsilon` and the forcing/contact terms remain as displayed. Mixing these
normalizations would introduce a wrong factor of epsilon.

Indeed, with

\[
c=\operatorname{Re}\widehat U-\frac{\operatorname{Im}\widehat W}{2k},
\quad G=\frac{3L^2\operatorname{Re}\widehat U-L\operatorname{Re}\widehat W}{2k},
\]

`R_m=kc+G`, while

\[
c'=\operatorname{Re}r_u-\frac{\operatorname{Im}r_w}{2k},
\]
\[
G'=L(R_m-3P_m)+\frac{Lf}{2k}
       +\frac{3L^2\operatorname{Re}r_u-L\operatorname{Re}r_w}{2k}.
\]

Thus differencing `G+C_R` alone misses `Delta(kc)` when reconstructed
flow has residuals. The complete equation (2) preserves that contribution.

## 3. Rigorous finite-band conservation bound

Suppose the stated functions are absolutely continuous in time, the
nonnegative residual/source/contact bounds are measurable and their displayed
majorants are integrable in `(k,eta)`. Assume
`|r_u|<=q_u`, `|r_w|<=q_w`, `|f-g|<=q_g` and `|zeta_C|<=q_C`.
Real forcing gives the exact bound

\[
\mathcal B_{\rm Ward}=\int_a^b\!d\eta\int_0^K\!d\mu
\left[\left(k+\frac{3L^2}{2k}\right)q_u
 +\frac{\sqrt{k^2+L^2}}{2k}q_w
 +\frac{L}{2k}q_g+q_C\right].                       \tag{3}
\]

The square root uses Cauchy--Schwarz on the jointly bounded real and imaginary
parts of `r_w`. If only separately bounded components are supplied, use
`q_w,I/2+L q_w,R/(2k)` instead; a component bound is not a complex-modulus bound.

For an independently computed ledger `Qhat`, with a proved quadrature bound
`|Qhat-integral Fhat_K|<=e_Q`, and serialized endpoint errors `e_a,e_b`,

\[
|\widetilde R_K(b)-\widetilde R_K(a)-\widehat Q|
\le \mathcal B_{\rm Ward}+e_Q+e_a+e_b.               \tag{4}
\]

This is a composable certificate once all its inputs have genuine enclosures.
GL-order agreement, accumulator precision and signed residual cancellation
are not substitutes for those inputs.

If the envelopes are independent of momentum, (3) becomes the explicit
finite expression

\[
\int_a^b\!d\eta\left[
 \frac{K^4+3L^2K^2}{8\mathrm{Pi}^2}q_u
 +\frac{(K^2+L^2)^{3/2}-L^3}{12\mathrm{Pi}^2}q_w
 +M_K L q_g+\int_0^Kq_C\,d\mu\right].             \tag{5}
\]

No infrared cutoff is needed for bounded normalized residuals: the weighted
`1/k` terms are `O(k)` near zero. More generally the integrability premise
must be checked, not inferred from having finite momentum samples. The exact
`k=0` mode normalization is not evaluated; the stress integrals are improper
integrals with the displayed finite limit behavior.

For discrete mode measure `nu` and analytic continuum contacts, let
`M_nu=integral dnu/(2k)` and use the direct mixed-target stresses. Their
full defect is

\[
\int d\nu\,\mathcal P(r_u,r_w)
 +L M_\nu(f-g)+L(M_\nu-M_K)g+\zeta_C^K.             \tag{6}
\]

Here `Pcal` is the first fraction in (2). The extra signed momentum moment
term is required. A pressure shift `(M_K-M_nu)g/3` removes it by explicitly
changing the target; it cannot be silently applied. The single moment error
does not bound all other continuum momentum errors.

## 4. Initial-state and trajectory errors: a separate theorem

Compare reconstructed amplitudes with exact flow driven by `g`. Let
`delta U_a,delta W_a` be the actual initial-state errors. Variation of
constants gives, with `E_k(d)=exp(2ikd)` and
`Phi_k(d)=(E_k(d)-1)/(2ik)`, `Phi_0(d)=d`,

\[
\delta W(t)=E_k(t-a)\delta W_a
 +\int_a^tE_k(t-s)[r_w(s)-(f-g)(s)]\,ds,
\]
\[
\delta U(t)=\delta U_a+\Phi_k(t-a)\delta W_a
 +\int_a^t r_u(s)\,ds
 +\int_a^t\Phi_k(t-s)[r_w(s)-(f-g)(s)]\,ds.         \tag{7}
\]

Since `|E|=1` and `|Phi_k(d)|<=min(d,1/k)` for `d>=0`, rigorous envelopes are

\[
b_W=|\delta W_a|+\int_a^t(q_w+q_g),
\]
\[
b_U=|\delta U_a|+(t-a)|\delta W_a|+\int_a^t q_u
       +\int_a^t(t-s)(q_w+q_g).                     \tag{8}
\]

The pointwise direct stress error is consequently bounded by

\[
e_R(k,t)\le\left(k+\frac{3L^2}{2k}\right)b_U
        +\frac{\sqrt{k^2+L^2}}{2k}b_W+e_{C_R},
\]
\[
e_P(k,t)\le\left(\frac{k}{3}+\frac{L^2}{2k}\right)b_U
        +\frac{\sqrt{k^2+L^2}}{2k}b_W+e_{C_P}.       \tag{9}
\]

For the pressure coefficient this uses
`|2k^2/3-L^2|<=2k^2/3+L^2`. Integrating (9) with the declared fixed
measure supplies actual finite-band stress error bounds if all inputs are
certified. It does not erase their initial-state terms.
The corresponding direct work error is bounded separately by

\[
|I_K-\widehat Q|
\le e_Q+\int_a^b[L(e_{R,K}+3e_{P,K})+3|h'|e_{B_0,K}]\,d\eta. \tag{10}
\]

These formulas assume exact fixed geometry and exact `h'`. Geometry and
source-jet uncertainty require additional product-error terms; they are not
included in an amplitude-only error bar.

## 5. Entire retarded kernels: remove momentum quadrature for the forcing error

For the important subcase with identical initial data, no trajectory
residuals and a real force error `delta=f-g`,

\[
X_\delta(t)=-\int_a^t\frac{\sin(2k(t-s))}{2k}\delta(s)\,ds,
\quad Z_\delta(t)=-\int_a^t\cos(2k(t-s))\delta(s)\,ds,
\quad T_\delta=2kX_\delta.
\]

Hence the forcing error in the separately defined density and pressure can
be integrated over continuous momentum exactly. For lag `d>=0`, put

\[
H_0(d)=\int_0^K\sin(2kd)\,dk
=\frac{1-\cos(2Kd)}{2d},
\]
\[
H_1(d)=\int_0^K k\cos(2kd)\,dk
=\frac{K\sin(2Kd)}{2d}+\frac{\cos(2Kd)-1}{4d^2},
\]
\[
H_2(d)=\int_0^K k^2\sin(2kd)\,dk
=-\frac{K^2\cos(2Kd)}{2d}
  +\frac{K\sin(2Kd)}{2d^2}+\frac{\cos(2Kd)-1}{4d^3}.
\]

These are entire functions; use their series or integral definitions near
zero. Their zero-lag values are `H0(0)=0`, `H1(0)=K^2/2`, `H2(0)=0`.
The removable apparent singularities are not an infrared cutoff.

\[
\Psi_R(L,d)=\frac{L H_1(d)}{4\mathrm{Pi}^2}
             -\frac{3L^2H_0(d)}{8\mathrm{Pi}^2},
\]
\[
\Psi_P(L,d)=\frac{H_2(d)}{6\mathrm{Pi}^2}
             +\frac{L H_1(d)}{4\mathrm{Pi}^2}
             +\frac{L^2H_0(d)}{8\mathrm{Pi}^2},
\]
\[
\delta R_K(t)=\int_a^t\Psi_R(L(t),t-s)\delta(s)\,ds,
\quad
\delta P_K(t)=\int_a^t\Psi_P(L(t),t-s)\delta(s)\,ds. \tag{11}
\]

This pressure kernel follows from the direct pressure operator, not from
solving a Ward equation. The kernel check

\[
(L^2\partial_L+\partial_d)\Psi_R=L(\Psi_R-3\Psi_P),
\quad \Psi_R(L,0)=L M_K
\]

is an independent algebra control on the force-mismatch Ward term.
These kernels cover the forcing difference with **unchanged contacts**.
They are not the complete physical change between two different metric
profiles, whose direct contacts would also change.

For nonnegative force-error envelopes let
`alpha(t)=integral_a^t q_g(s) ds`,
`beta(t)=integral_a^t (t-s)q_g(s) ds`. Elementary integral bounds give

\[
|\delta R_K(t)|\le L M_K\alpha
 +\frac{3L^2}{8\mathrm{Pi}^2}\min(K\alpha,K^2\beta)
\le M_K[L\alpha+3L^2\beta].                         \tag{12}
\]
\[
|\delta P_K(t)|\le L M_K\alpha
 +\min\left(\frac{K^3\alpha}{18\mathrm{Pi}^2},
             \frac{K^4\beta}{12\mathrm{Pi}^2}\right)
 +\frac{L^2}{8\mathrm{Pi}^2}\min(K\alpha,K^2\beta). \tag{13}
\]

For example `|H1|<=K^2/2`,
`|H0|<=min(K,d K^2)` and
`|H2|<=min(K^3/3,d K^4/2)`. The lag-integrated inequality is applied before
each minimum in (13). It is valid because both alternatives independently
bound the same integral. These bounds hold over the entire continuous
momentum band; nine narrow probe enclosures are not being reinterpreted as
a completed continuum quadrature.

## 6. Explicit analytic Taylor-tail consequence on the inherited interval

Let `P24` consist of the **exact** degree-24 Taylor polynomial of `g` on each
of the 64 inherited panels. Cauchy radius `R=1/8`, half width `H=1/128`
and uniform complex bound `|g|<=64` give the inherited analytic tail envelope

\[
E(s)=64\sum_{n=25}^{\infty}|(s-c_j)/R|^n
\quad(s\hbox{ in panel }j).
\]

The inherited integrated proof gives

\[
\int_a^b|g-P_{24}|\,ds\le\int_a^b E(s)\,ds
\le\alpha_*:=\frac1{390\,2^{90}}.                   \tag{14}
\]

The envelope is the same even function about each panel center. Its first
moment on that panel is therefore `c_j integral E`; the average of all 64
centers is `(a+b)/2`. Consequently

\[
\int_a^b(b-s)E(s)\,ds
=\frac12\int_a^b E(s)\,ds\le\beta_*:=\alpha_*/2.
\]

For every `t in [a,b]`, positivity and `t-s<=b-s` imply that its prefix
`alpha(t)<=alpha_*` and `beta(t)<=beta_*`. This controls intermediate
times as well as the final endpoint without evaluating any source value.

Using `L<=2/7` and `Pi>=3`, the all-time source-truncation contribution
has the following exact rational upper bounds for `K=64,128,256`:

\[
B_R(K)=\min\left(\frac{K^2}{72}\frac{20}{49},
                 \frac{K^2}{252}+\frac{K}{294}\right)\alpha_*,
\]
\[
B_P(K)=\left[\frac{K^3}{162}+\frac{K^2}{252}
                              +\frac{K}{882}\right]\alpha_*,
\qquad
B_{\rm source\,Ward}(K)=\frac{K^2}{252}\alpha_* .    \tag{15}
\]

At `K=256`, exact rational arithmetic proves

\[
B_R<6\times10^{-28},\qquad
B_P<2.2\times10^{-25},\qquad
B_{\rm source\,Ward}<6\times10^{-28}.
\]

The last bound controls `|M_K integral L(P24-g)|`, the integrated forcing
mismatch with unchanged exact metric contacts. No observed response is an
input. These are analytic **truncation-only** bounds in the inherited scaled
stress/ledger model units. They exclude coefficient enclosure widths,
rounded source polynomial coefficients, propagation/phase arithmetic,
incoming-state error, full contact evaluation, momentum quadrature and time
quadrature. They are not new acceptance tolerances or the achieved errors of
any run. The pressure bound exceeds `1e-26`; the source-moment gate must not
be applied to a distinct pressure quantity without a new registration.

## 7. What conservation cannot establish

An arbitrary incorrect incoming state propagated on shell still satisfies
the same Ward equation. For example `U=A` real constant and `W=0` is a
homogeneous perturbation with
`R_m=A(2k^2+3L^2)/(2k)` and
`P_m=A(2k^2/3-L^2)/(2k)`; it changes direct stress but has zero Ward defect.
Its magnitude is unbounded when `A` is unspecified. Thus (4) never replaces
the state/error theorem (7)--(10).

Similarly, values at finitely many momentum probes do not determine a
momentum integral. A smooth incoming perturbation supported between the
probes vanishes at every probe, has finite nonzero physical integrated stress,
and continues to satisfy Ward. State and spectral regularity bounds are
indispensable for a genuine continuum certificate.

The source kernels and Cauchy corollary apply only to this plane-wave
linearization and the interior source interval. They do not certify a shifted
shell root, full support endpoints, nonlinear geometry, thermalization, a
hot radiation era, a coupled Einstein solution or a higher-dimensional cause
of the Big Bang. All historical metric failures remain unchanged.

## 8. Finite verification and a publishable next scope

The local symbolic verifier checks the full residual identity, canonical
drift, independently defined density and pressure kernel reductions,
finite-K primitive derivatives, zero-lag limits and kernel Ward control.
It rejects symbolic omissions of force work, baseline work, canonical drift,
contact work and pressure terms. It separately checks the rational bounds
in (15). It loads no project numerical module or scientific array.

An independent standard-library verifier can represent exact sparse
polynomials for (2), use formal derivatives of sine/cosine for the finite-K
primitive checks, and prove all remaining comparisons using `Fraction`.
Source-envelope symmetry is a calculus premise documented in (14), not a
claim that a finite algebra checker formalizes Cauchy's theorem.

The meaningful short checkpoint is: **an exact finite-band forcing-error
kernel, a complete off-shell Ward budget, and explicit all-time analytic
source-truncation stress/work bounds**. Numerical full-ledger certification
should only follow a separately frozen implementation that supplies every
missing state, coefficient, residual, contact and quadrature enclosure.
