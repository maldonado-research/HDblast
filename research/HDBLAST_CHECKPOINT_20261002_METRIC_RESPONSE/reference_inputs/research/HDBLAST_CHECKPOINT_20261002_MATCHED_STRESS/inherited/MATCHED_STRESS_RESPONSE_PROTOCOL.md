# Prospective matched linear stress and current response

Prepared 2 October 2026. This is analytic follow-up work, developed after the
variance calibration. It registers or executes no new numerical source,
response, mode, stress, or spectral experiment. Future numerical settings and
gates must be frozen separately. The formulas below apply to the inherited
minimal scalar and fixed positive-reference prescription, per species.

## 1. Geometry, action, and state

Keep de Sitter geometry fixed, with signature (-+++), constant H>0,
`a=-1/(H eta)`, `L=a'/a=-1/eta`, `L'=L^2`, and `x0=r=2H^2`.
The primary calibration uses H=1 as a unit choice. Let `d=delta x`,
`s=a^2 d`, with a smooth compact source and a lower integration boundary in
its source-free past. Keep the incoming BD density matrix fixed. The exact
reference modes are `v_k=exp(-ik eta)/sqrt(2k)`. No finite-time vacuum is reset.
Dots denote cosmic-time derivatives; primes denote conformal-time derivatives.

The physical stress is the minimally coupled stress. Its canonical frequency
being k at this reference does not replace it with an improved conformal
stress. The inherited continuum reference values are

    Q0=H^2/(12 pi^2), rho0=11 H^4/(960 pi^2), p0=-rho0.

The finite prescription is more than the action restricted to S4. Explicitly
inherit `COMMON_ACTION_SMOOTH_FRW.md`, sections 1–5: fourth-order stress and
second-order field-square subtraction at fixed r, and the common local action

    Gamma_loc = integral sqrt(-g)[-V(x)+F(x)R+alpha R^2+beta C^2+gamma E4].

Its reference matching is `V(r)=V_x(r)=V_xx(r)=F(r)=F_x(r)=0`, with the
declared zero reference curvature-squared coefficients and the stated
no-boundary-variation convention for total derivatives. Section 2 of
`DESITTER_COMMON_ACTION_DERIVATION.md` identifies the same finite local
convention. No new finite stress contact or derivative operator is selected.
Sphere-restricted W alone would not supply this off-stationary information.

The inherited PV coefficient bridge is algebraic. Its separate heavy-field
decoupling argument has explicit static-prepared regulator-state hypotheses;
it is not a new proof about arbitrary BD regulator preparations. The local
adiabatic subtraction and trace algebra used here apply directly on the
present smooth compact de Sitter interval.

## 2. Direct physical-mode response, including explicit contacts

Define `D v=v'-L v` and use momentum measure `dmu=k^2 dk/(2 pi^2)`.
The physical per-mode observables are

    e_k = [|D v|^2+(k^2+a^2 x)|v|^2]/(2a^4),
    p_k = [|D v|^2-(k^2/3+a^2 x)|v|^2]/(2a^4),
    Q_k = |v|^2/a^2.

The forced variation satisfies

    delta v''+k^2 delta v=-s v,
    delta v(eta_i)=delta v'(eta_i)=0.

Let

    I_k=integral_eta_i^eta s(t) sin[2k(eta-t)] dt,
    J_k=integral_eta_i^eta s(t) cos[2k(eta-t)] dt,
    A_k=delta|v_k|^2=-I_k/(2k^2).

Then `A_k'=-J_k/k`, `A_k''+4k^2 A_k=-s/k`, and
`delta|v_k'|^2=-k^2 A_k`. Consequently the exact bare variations are

    delta e_k = [3L^2 A_k-L A_k'+s/(2k)]/(2a^4),
    delta p_k = [-(4k^2/3+L^2)A_k-L A_k'-s/(2k)]/(2a^4).

The `+s/(2k)` and `-s/(2k)` pieces are explicit variations of the mass in
the stress operator. They are respectively `+Q0,bare,k d/2` and
`-Q0,bare,k d/2`, before combining them with mode response. Omitting them
changes the observable even if the forced modes are correct.

Equivalently the brackets above are

    delta e: -3L^2 I_k/(2k^2)+L J_k/k+s/(2k),
    delta p: (2/3+L^2/(2k^2))I_k+L J_k/k-s/(2k).

All expressions use the same fixed comoving cutoff before subtraction and
cutoff removal. The apparent singular low-k mode factors have integrable
limits after the momentum measure is included.

## 3. Differentiate the complete inherited subtraction

Write `w=sqrt(k^2+a^2 r)`, and let U be the reference second-order WKB term:

    U=-a''/(2aw)-w''/(4w^2)+3w'^2/(8w^3).

The inherited general fourth-order term is

    W4=-W2^2/(2w)-W2''/(4w^2)+w'' W2/(4w^3)
       +3w' W2'/(4w^3)-3w'^2 W2/(4w^4).

At fixed geometry and fixed r define

    u=delta W2=s/(2w),
    v=delta W4=-Uu/w-u''/(4w^2)+w''u/(4w^3)
               +3w'u'/(4w^3)-3w'^2u/(4w^4),
    b=-k^2/3-a^2 r, c=1-b/w^2,
    j4=L u'/w^2-2Lw'u/w^3+w'u'/(2w^3)-3w'^2u/(4w^4).

Here b is only a subtraction abbreviation, unrelated to the mass-law
coupling used below. The complete variations at `x0=r` reduce to

    delta e_sub = [s/w-L^2 u/w^2+j4]/(4a^4),
    delta p_sub = [c u-s/w+c v+2bUu/w^3+sU/w^2-L^2u/w^2+j4]/(4a^4),
    delta Q_sub = -s/(4a^2 w^3).

The cancellation in delta e_sub follows from `2Uu/w=sU/w^2`; it is not
permission to drop fourth-order subtraction before variation. Pressure
retains source second derivatives through v. The positive reference remains
in w; derivatives of the perturbation belong in u and v, never in w.

The primary direct definitions for a future computation are

    delta rho_K=integral_0^K dmu (delta e_k-delta e_sub),
    delta p_K=integral_0^K dmu (delta p_k-delta p_sub),

together with the already derived combined delta Q_K. Combine the terms
before removing K. Separate divergent integrals do not define contacts.

## 4. Exact closed stress response from the matched variance

Direct reduction of both stress subtractions gives a useful analytic result:

    delta rho = (L^2 delta Q-L delta Q')/(2a^2)+Q0 d/2
                +(L^2 d-L d')/(96 pi^2 a^2),

    delta p = (delta Q''+L delta Q'-3L^2 delta Q)/(6a^2)-Q0 d/6
              +(d''+L d'-11L^2 d)/(288 pi^2 a^2).

Here delta Q is the fully matched retarded response, including its inherited
finite local term. These formulas hold for `x0=r=2H^2` on the fixed de Sitter
geometry and stated fixed state. They are not a general root Hessian or a
metric-response formula. All terms have physical mass dimension four.

For density, the bare bilinear identities already give the first two terms
with bare Q and Q0. The subtraction mismatch is exactly

    [(L^2 delta Q_sub-L delta Q_sub')/(2a^2)
      +Q0,sub,k d/2-delta e_sub]
      = M^2 L(3Ls-s')/(16a^4 w^5),  M^2=a^2 r.

Using `integral_0^infinity k^2 dk/(k^2+M^2)^(5/2)=1/(3M^2)` gives the
displayed density contact without choosing a new finite coefficient.

Pressure can likewise be reduced directly before using a Ward or trace
identity. At a unit observation scale factor, its subtraction mismatch
after the bare Q reduction is

    H^2[70H^6 d-50H^4 d w^2-10H^3 dot(d) w^2
        -5H^2 d w^4+4H dot(d) w^4+ddot(d) w^4]/(24w^9).

Its exact `k^2 dk` integral is
`[ddot(d)+2H dot(d)-11H^2 d]/144`. Including `1/(2pi^2)` and restoring
conformal derivatives gives the pressure contact above. The trace identity
below is an additional derivation and check, not the definition of either
direct stress integral. A pure symbolic development attempt omitted the
specialization `r=2H^2` in this pressure comparison and failed; its actual
output and explicitly reconstructed source are retained in `development/`.
No physical or numerical source evaluation occurred in either attempt.

In cosmic time the density expression is
`H^2 delta Q/2-H dot(delta Q)/2+Q0 d/2+[H^2 d-H dot(d)]/(96pi^2)`.
For the separate stationary past-infinite mass variation, it yields

    rho_x=H^2[5/3-2 gamma_E-ln 2]/(32pi^2),  p_x=-rho_x,

agreeing with independent fixed-r common-action derivatives. That state and
switching limit is an analytic control, not a compact-pulse numerical run.

For a time-domain implementation define
`F[f]=integral_eta_i^eta f(t)[ln(M(eta)(eta-t))+gamma_E+1]dt` and
`q=a^2 delta Q=-F[s']/(8pi^2)`. Smooth zero initial data give the regular
logarithmic-integral derivative formulas

    q'=-(F[s'']+L s)/(8pi^2),
    q''=-(F[s''']+2L s'+L' s)/(8pi^2),
    delta Q'=a^-2(q'-2Lq),
    delta Q''=a^-2[q''-4Lq'+(4L^2-2L')q].

These retain the derivatives of the time-dependent subtraction scale M and
avoid an invalid raw endpoint differentiation. Here `L'=L^2`.

## 5. Trace, Ward identity, and finite-cutoff terms

The inherited local trace remainder is fixed by the common prescription:

    T=-rho+3p=-x Q-(1/2) box Q+A_r,
    A_r=(16pi^2)^-1{[(x-r)-R/6]^2/2
          +(Riem^2-Ric^2)/180+box R/30-box x/6}.

On fixed de Sitter, `box f=-(f''+2L f')/a^2`. Therefore

    delta A_r=-H^2 d/(8pi^2)+(d''+2L d')/(96pi^2 a^2),
    delta T=-2H^2 delta Q-Q0 d
             +(delta Q''+2L delta Q')/(2a^2)+delta A_r.

The differentiated local trace term, including its source derivatives, is
not an adjustable anomaly convention in this inherited model. The two
direct stress responses obey

    delta rho'+3L(delta rho+delta p)=Q0 d'/2.

Both Q-derivative terms and source derivative contacts cancel correctly;
the remaining cancellation uses `Q0=H^2/(12pi^2)`. Ward consistency alone
could not have fixed the stress. With an independently specified delta T it
also permits the reconstruction

    delta rho=a^-4 integral_eta_i^eta a(t)^4[Q0 d'(t)/2-L(t)delta T(t)]dt,
    delta p=(delta T+delta rho)/3.

The homogeneous addition `C/a^4` is fixed to zero by the unchanged initial
state and source-free initial interval. Defining stress by this integral
would make a Ward test tautological; a future test must compare the direct
stress integrals or closed formulas with a separately evaluated ledger.

At finite comoving K, use the matched baseline
`Q0,K=integral_0^K dmu [1/(2ka^2)-1/(2a^2w)+U/(2a^2w^2)]`, which is
generally time dependent. The exact finite-cutoff Ward source is
`Q0,K d'/2`, not the continuum Q0. The density formula remains exact with
Q and Q0 replaced by Q_K and Q0,K and its last local term multiplied by
`vK^3`, where `vK=K/sqrt(K^2+M^2)`.

The finite-cutoff trace uses Q_K, Q0,K and the following explicit remainder:

    delta A_K=(2pi^2)^-1{
      r[-6H^2 d+5H dot(d)+ddot(d)] J5(K/a)/16
      -r^2[50H^2 d+10H dot(d)] J7(K/a)/32
      +70r^3 H^2 d J9(K/a)/64},

    Jn(P)=integral_0^P p^2 dp/(p^2+r)^(n/2)
      =r^((3-n)/2) sum_(ell=0)^((n-5)/2)
        (-1)^ell binomial((n-5)/2,ell) v^(2ell+3)/(2ell+3),
    v=P/sqrt(P^2+r).

It tends to delta A_r above. The exact finite-K trace companion is

    delta T_K=-2H^2 delta Q_K-Q0,K d
              +(delta Q_K''+2L delta Q_K')/(2a^2)+delta A_K.

This supplies an independent finite-band comparator for direct pressure.
Replacing these finite terms by continuum contacts silently adds cutoff
error to a supposed exact identity.

## 6. Current and interpretation

For the inherited quadratic mass law, now denoting its coupling by b,
`delta x=br delta phi` and

    delta j=(br/2)delta Q+(b^2 r Q0/4)delta phi.

The second term is the explicit mass-law contact, distinct from all
renormalization contacts. The common shell loop amplitude multiplies the
entire translated current and stress once. For a scalar perturbation delta phi,
b=0 makes d=0 and the translated response vanishes.

The free canonical energy variation satisfies
`delta(|v'|^2+k^2|v|^2)=0` at first order, but physical minimal stress
contains linear coherent terms from L and the mass. After a compact pulse,
all local source contacts vanish while these retarded stress terms can
remain. Canonical excitation energy starts at quadratic order. The linear
Ward source is `Q0 delta x'/2`; `delta Q delta x'/2` enters only at second
order. No first-order stress result establishes net heating, particle yield,
positive radiation energy, dissipation, or thermalization. The inherited
spectral work identity also contains its fixed-reference drift; it does not
turn canonical oscillator energy into the full physical stress.

## 7. Requirements for a separately registered bounded test

Freeze a compact source, observation domain, precision, cutoff and quadrature
ladders, derivative method, absolute and relative allowances, and a finite
stopping rule before any new stress/source evaluations. The existing pulse
may be reused only by explicit registration; this document selects no new
numerical input or acceptance tolerance.

- Compare direct forced-mode stress with the combined sine/cosine integrals
  and independently evaluated closed formulas. Freeze stable analytic
  derivatives of the logarithmic memory response; do not differentiate its
  singular endpoint by a naive Leibniz boundary term.
- Enforce zero response before the source, unchanged initial data and
  Wronskian controls, and invariance under moving the lower limit within the
  same source-free past.
- Check independently evaluated finite-K Ward and trace identities with
  Q0,K and delta A_K. Separate exact algebra controls from numerical ledger
  and derivative accuracy; no diagnostic may be made true by definition.
- Derive stress tail bounds and derivative-response error bounds before
  claiming continuum accuracy. The variance tail bound alone does not bound
  stress, which contains extra k weights and time derivatives. Separate
  cutoff, finite-k quadrature, time quadrature, mode, and cancellation error.
- Detect wrong minimal/improved stress, omitted explicit mass contacts,
  missing fourth-order pressure subtraction, changed trace-derivative sign,
  moving r, missing mass-law contact, and a varied/reset initial state.
  Preserve failed runs and keep all checks active under Python optimization.
- Check dimensions and the distinct analytic stationary derivative above.
  Do not replace a compact causal history by an equilibrium susceptibility.

This derives the homogeneous scalar-to-stress and scalar-to-current response
around the exact reference on a prescribed geometry. Metric-metric and
metric-scalar kernels, bulk and shell perturbation equations and boundary
conditions, response about an actual shifted stationary root, physical EFT
matching, and any coupled evolution remain outside this bounded task.

## Provenance

`INPUT_SHA256.json` pins the exact inherited files read for this derivation.
They include the prospective variance protocol, full smooth-FRW common
action, its two general symbolic trace/action derivations, the de Sitter
matching derivation and direct exact-mode stress bridge, and the stationary
common-action model. All were read without edits. The closed stress reduction
and its pure symbolic audit are new post-calibration analytic work.
