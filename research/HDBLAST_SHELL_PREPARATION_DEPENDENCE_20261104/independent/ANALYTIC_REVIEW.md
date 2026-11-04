# Independent preparation-dependence review

4 November 2026. Completed private independent review of the frozen P1
calculation and its elementary sufficient certificate. This is a second AI
derivation and arithmetic check, not external specialist review.

## Inputs and conventions

The references are the 1 November preparation-dependence design and the
published prepared-shell proof, particularly `MEAN_INITIAL_GAUGE.md`,
`CONDITIONED_GAUGE.md`, and the response derivation. All use outward proper
radial coordinate y, `dy=-rho dz`, `r=z/ell`, `x=t/ell`, and
`partial_r=-ell rho partial_y`. The expansions use second parameter
derivatives, `X_e=X_0+e X_1+e^2 X_2/2+o(e^2)`.

I independently differentiated the local flux, the constrained initial map,
the conformal-coordinate expressions and the volume Green kernel. The numerical
review consumes frozen whole-cell boxes as data, without importing either
producer or evaluating background, bump or first-field oracles.

## Exact compact directions and constraints

Let chi be smooth and flat at both ends of the packet and put
`h=rho^-2 chi_y`. Then `integral rho^2 h dy=[chi]=0` exactly. Conversely, for
any smooth flat h with zero weighted moment, its primitive
`chi(y)=integral_(yL)^y rho^2 h ds` is flat and zero on both exterior sides.
This describes the entire compact zero-moment vector space; no numerical
subtraction of two weighted integrals is needed.

Expand the derivative in the flux identity:

```
[rho^4(h_y-a h)]_y
 =rho^4[h_yy+3a h_y-(a_y+4a^2)h].
```

Using the retained background equation
`a_y+4a^2=3/rho^2-(2/3)V`, this gives exactly

```
rho^4 C[h]=-(3/2)[rho^4(h_y-a h)]_y-(9/2)rho^2 h,
C[h]=-(3/2)(h_yy+3a h_y)-V h.
```

The primitive from the lower outward-y endpoint therefore equals
`-(3/2)rho^4(h_y-a h)-(9/2)chi`. Thus

```
zeta=3(h_y-a h)/(2p^2)+9chi/(2rho^4p^2),
s=-p zeta=-3(h_y-a h)/(2p)-9chi/(2rho^4p).
```

Its endpoint constants are zero. Both h and zeta are smooth and compact,
using the nonzero p bound on the support. Their support need only fit in the
whole packet interval, even when h has separated lobes. The homogeneous
canonical mean `(h gamma0+Lie_(zeta d_y)gamma0,0,0,0)` satisfies the full
linearized Hamiltonian and momentum constraints: its Hamiltonian is
`C[h]+rho^-4(rho^4 p^2 zeta)_y=0`; the momentum variation vanishes at the
time-symmetric background. The fixed first constraint Hessian is unchanged.
The absolute new family still satisfies `P B_lambda=-P0 Q`, with the same
entire actual first Hessian Q.

For each fixed finite lambda the weighted IFT has its same derivative in
the correction variable. Its quadratic mean equation is unchanged because
the added datum is in the kernel of P. The retained averaging/injectivity
argument therefore selects precisely the changed canonical mean. This is
conditional on the original weighted-IFT and development premises and is
not a uniform existence-radius statement in lambda.

## Conditioning, all five entries, and r derivatives

The pure second generator removing the spatial Lie derivative has initial
`Z^y=-xi_lambda`, `Z^t=0`, `Z^y_t=0`, `Z^t_t=a xi_lambda`. Its scalar pullback
must be retained. Subtracting the accepted conditioned adapters gives, in
order `(K,L,N,J,w)`,

```
W(0,r)=(h,0,0,h,s),
W_x(0,r)=(0,ell rho(-h_y+a h),0,0,0),
W_r(0,r)=(h_r,0,0,h_r,s_r).
```

The lapse-square value and every u, ut, uy and QM product cancel because the
first field and canonical second momentum are unchanged. They are present
in each absolute adapter. The mixed velocity does not cancel. There is no
rho_b division in `W_x`, and no additional parameter factorial.

For `chi=ell rho_b^3 F(r)`, independently differentiate using
`rho_r=-ell a rho^2`, `a_r=-ell rho a_y`, `p_r=-ell rho p_y`. Write

```
R=(rho_b/rho)^3, b=ell rho a, d=ell rho p_y/p,
c=-3rho_b^3/(2ell rho^4 p), G=F''+4b F'+3ell^2 F.
```

Then the directly usable profile formulas are

```
h=-R F',
h_r=-R(F''+3b F'),
V_L=-R(F''+4b F')=h_r+b h,
s=c G,
s_r=c{F'''+4b F''+3ell^2 F'
        -4ell^2 rho^2(a_y+a^2)F'+(4b+d)G}.
```

The two forms of V_L are an independent algebraic check of both radial
orientation and the harmonic mixed entry. The terms containing a_y and p_y
in s_r must both remain. The Gram matrix is diag(1,2,1,3,1), so cell bounds
can use exactly

```
|W0|_H=sqrt(4h^2+s^2),
|W0_r|_H=sqrt(4h_r^2+s_r^2),
|W0_x|_H=sqrt(2)|V_L|.
```

A preparation variation is not the conditioning diffeomorphism. The raw
scalar variation is zero, so an invariant radial spatial gauge generator
with that scalar variation is zero since p is nonzero. A normal generator
has zero induced spatial/scalar change on this time-symmetric section.
Moreover the equations for a hypothetical vanishing metric variation give
`h+2zeta_y=h+2a zeta=0`; hence `zeta_y=a zeta`, whose compact solution is
zero. Thus a nonzero h changes physical initial data modulo the anchored
gauge. It need not change this particular target observable.

## Difference equation and exact physical observable

Subtract the complete five-component equations and all shell rows for the
two conditioned families. Identical actual first-field bulk and boundary
forcing cancels. The difference retains every internal metric/scalar,
connection, and Robin coupling:

```
W_xx-W_rr=A0 W_x+B0 W_r+D0 W,
L|shell=0, P_N W_r=M0 P_N W.
```

The inherited constraint and harmonic subsidiary conditions are preserved
by the same retained uniqueness theorem. This argument does not replace
the equations with five uncoupled free waves. The free wave is only the
comparison field in the subsequent contraction error estimate.

The background shell scalar is constant in shell time. The proper mean's
clock and volume-measure corrections depend only on the common first jet,
and therefore cancel upon subtraction. Hence `m_h(x)=W_w(x,0)` is exactly
the preparation slope of the anchored proper scalar mean. The absolute
proper-mean correction and every term of the absolute first-only Q1 remain
in each absolute response. Subtraction alone yields

```
nu_tautau+4H nu_tau-2nu/A^2=(F'_b/2)m_h,
nu(0)=nu_tau(0)=0, theta_h=nu_tau.
```

The new coefficients are affine because this entire difference problem is
linear with fixed coefficients, initial data and shell rows. No assertion of
affine dependence of the complete finite-epsilon spacetime is required.

## Positive Green weights and proper-time factors

Write `v=ell T`, `b=ell x`, `u=v-b`. The exact scaled Green function is

```
k(T,x)=sinh(2u)/(2ell)*cosh(b)^2/cosh(v)^2.
```

For `0<=x<=T`, its ratio to T-x lies in `[cosh(v)^-2,1]`: the lower bound
uses `sinh(2u)/(2u)>=1`; the upper bound follows from
`cosh(u+b)>=cosh(u)cosh(b)` and `tanh(u)<=u`.

Direct differentiation gives the useful exact simplification

```
k_T=cosh(b)^2*cosh(v-2b)/cosh(v)^3.
```

Since `|v-2b|<=v` and `b<=v`, `0<k_T<=1`. Also
`k_T>=cosh(v)^-3>=(1+v^2)^-3>=1-3v^2`. With ell=10^-6 and T=5/4,
the relative defects of both k/(T-x) and k_T are less than eta=10^-11.
The bound at x=T is interpreted continuously. The cell weights
`int_l^u (T-x)dx=(u-l)(T-(l+u)/2)` and `u-l` are nonnegative exact
rationals on the frozen history partition.

Background proper time is `tau=rho_b ell x`. Therefore

```
L_M=m_h(T),
L_v=(rho_b ell)^2(F'_b/2) int_1^T k(T,x)m_h(x)dx,
L_Theta=rho_b ell(F'_b/2) int_1^T k_T(T,x)m_h(x)dx.
```

The expansion loses one proper-time factor by differentiation. F'_b is
strictly negative, so it reverses each history's sign. The factor 1/2 is
already in the physical Gauss equation; the slopes vary second parameter
derivatives and acquire no additional 2 or 1/2.

## Complete homogeneous error and scope of an eventual sign

Use T=5/4, Delta=1/4, `q=1201/3200000`. Let D0, Dr, Dv bound the above H
norms on the whole required initial section, including the flat tail. Then

```
N=max((D0+T Dv)/3,Dr+Dv), E=qN/(1-q).
```

The actual initial support and retained finite propagation give the exact
zero wedge `x+r<=1`. The complete coupled resolvent implies
`|m_h(x)-s(x)|<=(x-1)_+ E`. Positive Green weights yield errors
`Delta E`, `Delta^3 E/6`, and `Delta^2 E/2` for the scalar value, weighted
history and unweighted history. Replacing the exact kernels by elementary
weights adds `eta A_v` and `eta A_Theta`, where each A integrates |s| with
the corresponding positive weight. No actual-first-field forcing cap is
added to the difference; no actual forcing is dropped from either absolute
preparation calculation.

Strict slope intervals, if obtained, make the corresponding affine
coefficient nonconstant. A certified finite rational lambda must be checked
by interval addition to the retained base coefficient, not by a nominal
division. Same-sign volume/expansion slopes permit simultaneous reversal of
the two negative base coefficients; nonzero slopes alone do not prove rank
two. An h supported outside the target's initial cone can be nonzero and
physically distinct but invisible at this target.

The review shares the retained analytic constraint/development/gauge/PDE
premises and authenticated background interval boxes. Independent arithmetic
over those boxes tests their aggregation and bookkeeping, not their original
enclosure construction or an independent PDE computation. Any resulting
claim concerns formal second parameter coefficients with a family-dependent
small epsilon interval; no effective interval or nonlinear observed sign is
supplied by this argument.

## Completed arithmetic and independent sufficient bounds

The separately registered whole-cell P1 campaign uses exactly
`chi=ell rho_b^3 Q(r)`, T=5/4, the original N1024 partition restricted to the
causal initial segment, and its registered flat-tail treatment. The independent
`reaggregate.py` reconstructs all 242 cell/target profiles from raw background
and bump jets, verifies all five initial vector entries, checks the background
endpoints against the retained cells, and independently accumulates the
positive history weights and exact dyadic upper Gram square roots. All 41
scientific aggregate fields match the producer's rational endpoints exactly.
Fourteen in-process negative controls reject omitted derivative or adapter
terms, malformed coverage, incorrect Gram/error/history/clock/factorial terms,
and the wrong finite-lambda sign.

The elementary certificate is separately checked by `check_coarse.py`. With
`D=(r-1)(2-r)` its exact second and third derivative numerators are
`P2=1-6D+6D^2` and `P3=1-10D+24D^2-12D^3`. On `0<=D<=3/16`, both decrease,
P2 is at least 11/128, and P3 is between -113/1024 and 1. Continuous extrema
of `w^j exp(-w)`, w>=16/3, plus independent rational Taylor estimates give
the whole-segment caps `(1/3,8,216,6400)` for `(Q,Q',Q'',Q''')`.

Since Q,Q',Q'' and c are positive on the open causal interval, integration
by parts gives the exact positive moment lower bounds `c_min/4` and
`32c_min/9`; the target scalar lower bound is `1408c_min/81`. The checker
reconstructs the broad background inclusions, all seven scalar-derivative
allowances, the complete Gram triangle bounds, the same coupled resolvent,
both proper-time factors and the finite-lambda computation. All 23 scientific
fields match the separate coarse report exactly. Eight additional in-process
mutants reject. This route does not use the fresh bump-jet oracle or quadrature.

Both methods establish positive L_M and strictly negative L_v and L_Theta on
the same direction, and independently select the same finite lambda=-10^21
under the frozen power-of-ten witness rule. With the retained negative base
intervals, both new second coefficients have strictly positive lower bounds.
This is one preparation counterexample to independence of the coefficient
signs from second preparation, at the fixed anchored target. It is not a
rank-two result, a uniform-in-lambda nonlinear theorem, an effective amplitude
bound, or external peer review. The paired normal/optimized reference reports
are byte-identical for each independent checker. In-process mutants are not
counted as separate rejected CLI executions.
