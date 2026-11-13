# From a linear scalar zero to a conditional nonlinear moment zero

4 November 2026, America/Los_Angeles. A source audit and a conditional
selection/continuity lemma. This note is not a new nonlinear Einstein–scalar
existence theorem, a numerical parameter-radius certificate, or evidence for
a cosmological origin. The methods used below are ordinary symmetry,
implicit-function uniqueness, Taylor's formula and the intermediate value
theorem.

## 1. What the retained project theorems already say

The [full finite-angular contact construction][formal], §§1 and 8, supplies a
conditional complete second Taylor coefficient on a future interaction slab.
Its convention is `g_e=g0+e*q+e^2*k/2` and similarly for the scalar. This
polynomial solves the equations through order two; that result alone is not
an exact family at nonzero parameter.

There is also a later and stronger retained result. The
[finite-regularity quasilinear shell theorem][exact], §§1–5, asserts a
conditional **exact** family on a prescribed future slab `[0,T]`, for all
parameters in one sufficiently small interval. Its solution is in `X_14`,
and its parameter curve is `C^2` into `X_10` and hence a lower piecewise
spacetime `C^4` topology. It assumes the exact localized constraint curve,
the state and hyperbolic margins, compatible boundary construction and a
strict causal support margin. The parameter interval is not numerically
evaluated. The [variation completion][variation] records the derivative
loss and actual-limit argument; formal differentiated equations alone do
not establish this regularity.

The [registered exact-family binding][binding], §§2–4, identifies the
`ell=10^-6` first seed and one selected radial weighted initial curve. It
uses correction support `ell/2<z<5ell/2`, `c_*=2`, `Z=6ell`, but its stated
future slab ends at `T=3ell/2`. That precise slab does not contain the left
endpoint of the newer central past-zero bracket. It cannot silently be
enlarged just because an older linear certificate has a wider domain.

Earlier [intrinsic shell calculations][intrinsic], §§1 and 4–6, and the
[exact-family TT binding][binding], §7, already use `C^2` dependence to infer
family-wise qualitative finite-parameter signs. Thus the missing step is not
an absence of every nonlinear existence or qualitative prediction theorem.
The important remaining tasks are an appropriate registered past-domain
application, a specified nonlinear completion and observer, and effective
constants/remainders for a quantitative amplitude regime.

## 2. A sufficient past slab, with its assumptions kept visible

Retain the exact localized constraint curve and all analytic premises above.
At the time-symmetric background section, background canonical momenta are
zero. Time reversal of exact constraint data sends both canonical momentum
densities to their negatives and leaves the spatial metric and scalar
unchanged. The Hamiltonian constraint is unchanged; the momentum constraint
changes its overall sign. Hence exact constraint satisfaction, radial
support and the initial open background shell collar are preserved.

The quasilinear argument must be applied to this reflected exact data curve,
with its complete corresponding gauge data. Equivalently one applies its
time-oriented proof to negative coordinate time. The scalar potential and
physical spatial-normal shell laws are time-reversal invariant, and the
background `cosh(t)` coefficients are even. The initial-data theorem permits
arbitrary admissible radial displacement and velocity profiles; the
quasilinear proof uses their norms rather than a favorable velocity sign.
These observations justify the application only under the same retained
state-flow, inverse, regularity, support and subsidiary-system hypotheses.
No new independent proof of those analytic theorems is supplied here.

Time reversal is a tensor pullback. In a reflected coordinate `s=-t`,
`g_ss(s)=g_tt(-s)`, `g_sI(s)=-g_tI(-s)`, and `g_IJ(s)=g_IJ(-s)`.
Consequently the initial mixed-component time velocity is even under this
pullback. Negating every coordinate metric time derivative would be wrong.
The natural harmonic completion or the fully pulled-back gauge jets must
be used.

For the central zero, a convenient sufficient choice is

```
L=5ell/2, c_*=2, Z=6ell, T=13ell/8.
L+c_*T=23ell/4 < 6ell,          causal margin = ell/4,
T+Z=61ell/8 < 8ell,             first-field triangle margin = 3ell/8,
T-6145ell/4096 = 511ell/4096 > 0.
```

The first identity uses the support of the **whole** constraint correction,
not merely the first packet. The last identity provides room for the
proper-clock construction below. These exact arithmetic margins do not
themselves certify the state-flow interval, nonlinear characteristic-speed
bound, weighted inverse, or permitted parameter radius. Those remain
retained analytic hypotheses, applied on this slightly longer slab. This
slab concerns the central zero only; it does not cover the later tail or
second late zero near `|t|=2ell`.

## 3. Specify the nonlinear observable before using Taylor's formula

Use the intrinsic shell geodesics normal to the common initial physical
shell section, as in [the intrinsic shell budget][intrinsic], §1. Label
their initial points by the round sphere and transport
`Y=n1^2-n2^2` with those observers. Before their first crossing, their
proper-time cuts have a positive induced three-metric and area measure
`dA_e(tau)`. Retain a common noncrossing tube after family-dependent
parameter shrinkage. A different physical section or observer choice can
give a different observable; joint diffeomorphism transport of fields and
anchors preserves this construction.

Let `Phi_e(tau,Omega)` be the full shell scalar pulled back to these cuts,
including the physical shell embedding, and let `phi_b` be the constant
background shell scalar. Define the actual-area quadrupole moment

```
M_e(tau) = integral Y*(Phi_e-phi_b) dA_e / integral Y^2 dA_e.
```

The denominator stays positive. On the background, `tau=rho_b*t`, and the
specified first physical field is `Phi_1=S(t)Y`. Because the background
shell scalar is constant in time and angle,

```
M_0=0,                    d_e M_e|0 = S(tau/rho_b).
```

First clock, angular-label and area variations multiply either a zero
background derivative or `Phi_0-phi_b=0`; they add no first-order amplitude
term. This is not permission to discard clock or measure terms at arbitrary
parameter or in other observables.

The retained `C^2` parameter curve in a lower `C^4` spacetime norm supplies
two parameter derivatives of the intrinsic geodesic flow, cuts, scalar
pullback and normalized moment on a smaller common compact tube. The
geodesic coefficients use one metric derivative; the time inverse and
positive induced determinants remain smooth there. One must shrink the
parameter interval so the cuts in the following fixed proper-time interval
remain inside the constructed coordinate slab:

```
I=[-rho_b*6145ell/4096, -rho_b*6143ell/4096].
```

Here `rho_b` is the exact background shell radius. At nonzero parameter
these proper-time cuts need not be constant-coordinate-time cuts. The
nonlinear conclusion below is therefore stated in `tau`, not as an
unchanged nonlinear coordinate-`t` bracket.

## 4. A sufficient nonlinear zero-persistence inequality

Assume the preceding application and observable construction have supplied
`M:(-a0,a0)->C(I)` of class `C^2`. Fix `0<a<a0` and set

```
K(a)=sup_{|e|<=a, tau in I} |d_e^2 M_e(tau)| < infinity.
```

For either parameter sign, integral Taylor gives

```
M_e/e - S = e*integral_0^1 (1-s)*d_e^2 M_(s*e) ds,
sup_I |M_e/e-S| <= |e| K(a)/2.
```

The retained [linear cell refinement][linear] gives, at the two background
clock endpoints, `S(t_L)>63` and `S(t_R)<-53`. Hence the sufficient guard

```
0<|e|<=a,                    |e| K(a)<106
```

preserves opposite signs of `M_e/e` at the endpoints. For positive
parameter the signs of `M_e` are positive then negative; for negative
parameter they reverse. Continuity gives at least one zero of `M_e` in
the **open** interval `I`. If `K(a)>0`, it is sufficient to choose
`0<|e|<min(a,106/K(a))`; if `K(a)=0`, every admitted nonzero
`|e|<=a` passes this guard.

This is a conditional nonlinear zero of a specified background-subtracted
quadrupole moment. It is not a zero of the full scalar field at every
angle, a zero of every harmonic, a uniqueness theorem, or a measured event.
No numerical value of `K(a)`, `a0`, the clock-tube radius or an effective
nonlinear amplitude window has been established.

For comparison, an anchored **initial-reference-measure** moment can use
the transported initial round measure instead of `dA_e`. These are distinct
observables. The fixed-round quadrupole projector has `L^infinity` operator
norm `6/pi`, because `integral Y^2=pi^2/3` and `integral |Y|=2pi`.
If the full pulled-back scalar's second parameter derivative is bounded
by `M_raw`, that reference-measure moment has normalized error at most
`3*M_raw*|e|/pi < M_raw*|e|`; thus `|e|M_raw<53` is sufficient. This
raw bound cannot be used for `K(a)` of the actual-area moment without
controlling its area and denominator derivatives as well.

At the single expansion point zero, there is nevertheless an exact useful
cancellation. Write `D_e=integral Y^2 dA_e` and
`N_e=integral Y(Phi_e-phi_b)dA_e`. Then
`N'_0=S D_0` and `N''_0=integral Y Phi_2 dA_0+2S D'_0`.
The quotient rule gives `M''_0=integral Y Phi_2 dA_0/D_0`.
That cancellation does not give a uniform bound on `M''_e` away from zero.

## 5. A selection-dependent oddness lemma

Let `R` be the spatial quarter-turn
`(X1,X2,X3,X4)->(-X2,X1,X3,X4)`. It is in `SO(4)`, normalizes the
coordinate-reflection group, and sends `Y` and the whole registered first
canonical seed `v` to their negatives.

Consider the radial weighted initial construction of [the exact initial
theorem][initial]. Require each actual choice to commute with this `R`:
the annulus, weights and defining function, background density,
constraint and adjoint operators, domain/kernel complement, and the local
inverse. These properties follow for the stated radial geometric choices
and their unique invariant inverse; an unspecified unrelated inverse must
not be substituted. In the canonical affine-seed case

```
x_e=x0+e*v+A u(e),              F(e,u(e))=0,
R*F(e,u)=F(-e,R*u).
```

The quarter-turn maps the invariant weighted spaces to themselves. Local
IFT uniqueness on a common smaller ball therefore gives
`u(-e)=R*u(e)` and `R*x_e=x_(-e)`. More generally the same argument works
for a trial curve satisfying `R*x_trial(e)=x_trial(-e)`. In particular an
`e^2 B/2` addition with `R*B=B` preserves this symmetry: the retained
radial `h`-only mean preparation is such an even direction. General
same-first-tangent second data need not satisfy this condition.

Require the gauge completion, coefficient extensions and exact evolution
to be chosen equivariantly as well. The relative-harmonic reference is
radial, the chosen initial lapse `rho(1+e*uY)` transforms appropriately,
and natural unique evolution/anchored geodesics preserve the relation.
Changing variables by `R` in the moment integral then gives

```
M_(-e)(tau)=-M_e(tau),             M''_0(tau)=0.
```

Only this selected quadrupole moment is asserted to have zero quadratic
parameter coefficient. Other second-order scalar/tensor sectors can be
nonzero: in particular the scalar mean, degree four, and the quadrupole
along the trace-free part of `A^2` need not vanish. The latter is
orthogonal to `A=diag(1,-1,0,0)` because `Tr(A^3)=0`.

Parameter normalization is part of this statement. Replacing the same
curve by `x_(f(e))`, with `f(e)=e+c e^2/2`, preserves its first tangent
but changes the second moment derivative to `c S` when the original
quadratic coefficient is zero. Such a reparametrization does not retain
the selected affine scalar initial sectors. Oddness is not a universal
property of arbitrary nonlinear completions or arbitrary parameters.

With only `C^2` regularity, this lemma gives
`M_e=e*S+o(e^2)`, uniformly on `I`, not a cubic remainder. The odd function
`e*|e|^(3/2)` is `C^2` but is not `O(e^3)`. An explicit third-order
claim requires extra parameter regularity and bounds. The modulus
`omega(a)=sup_{|e|<=a,tau in I}|M''_e|` tends to zero in this selected
case; the zero-persistence guard becomes `|e|omega(a)<106`, still without
a computed numerical radius.

## 6. The concrete missing estimate

For this one chosen observable and completion, validate a bound for
`K(a)` (or the odd-case modulus `omega(a)`) on a **certified** common
parameter/state/clock domain. This requires effective initial inverse and
nonlinear-map constants, evolution/trace and observer-flow bounds, and
the actual completion, not just a norm of a homogeneous difference.
The retained [preparation robustness theorem][class], especially its last
scope paragraph before the numerical certificate, explicitly supplies no
completion-uniform nonlinear radius. Its `h`-class coefficient bounds do
not bound the particular past observable or all higher remainders.

The separate scalar **mean** needed for volume acceleration remains
important even if this selected `Y`-moment's second coefficient vanishes.
[The intrinsic budget][intrinsic], §5, displays its nonzero coupling to
that mean. Neither parity nor a moment zero closes the nonlinear volume
equation, supplies radiation or a physical incoming preparation, calibrates
SI units, or establishes a higher-dimensional origin of the Big Bang.

[formal]: https://github.com/maldonado-research/HDblast/blob/f4001e5d94e6388bfff4f82dee63c174f2dfd5eb/research/HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014/FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md
[exact]: https://github.com/maldonado-research/HDblast/blob/f4001e5d94e6388bfff4f82dee63c174f2dfd5eb/research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md
[variation]: https://github.com/maldonado-research/HDblast/blob/f4001e5d94e6388bfff4f82dee63c174f2dfd5eb/research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/gauge-parameter/FINITE_REGULARITY_VARIATION_COMPLETION.md
[initial]: https://github.com/maldonado-research/HDblast/blob/f4001e5d94e6388bfff4f82dee63c174f2dfd5eb/research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/inputs/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md
[binding]: https://github.com/maldonado-research/HDblast/blob/f4001e5d94e6388bfff4f82dee63c174f2dfd5eb/research/HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/DERIVATION.md
[intrinsic]: https://github.com/maldonado-research/HDblast/blob/f4001e5d94e6388bfff4f82dee63c174f2dfd5eb/research/HDBLAST_INTRINSIC_SHELL_BUDGET_20261101/DERIVATION.md
[class]: https://github.com/maldonado-research/HDblast/blob/f4001e5d94e6388bfff4f82dee63c174f2dfd5eb/research/HDBLAST_INTRINSIC_PREPARATION_ROBUSTNESS_20261104/DERIVATION.md
[linear]: https://github.com/maldonado-research/HDblast/blob/4bfc15d0dee151aab288f33f34cf1ffd32e4c870/research/HDBLAST_PAST_FUTURE_SCALAR_TRANSFER_20261104/refinement/PROOF.md
