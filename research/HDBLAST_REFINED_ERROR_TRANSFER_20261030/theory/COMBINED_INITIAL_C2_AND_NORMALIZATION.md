# Combining tighter initial estimates and a coordinate seminorm scale

The actual target is the same conditional formal second coefficient `U`.
The reconstruction remains `W1_quadratic` (`W1_dagger`). Its initial
functions, two convolution widths, original bulk source `F/F_x`, physical
shell data `G/M`, first-order fields and finite bulk iteration are fixed.
Only upper estimates are refined. The candidate and its actual error do
not change when an interval majorant becomes smaller.

## Complete initial curvature-error bound

The retained paired comparison expresses the actual-minus-nominal initial
displacement as three telescoping channels. For an exact finite cover by
radial leaves, let `P_ch` be the maximum of that channel's positive full23
`H` norm cap over all leaves. Then

```
sup_r ||actual_u_rr(r)-nominal_u_rr(r)||_H <= sum_ch P_ch.
```

The new paired leaf calculation changes this estimate while preserving
the expression-defined functions. The independent audit recomputes every
new leaf and retains all untouched cells. A parent upper cap and a child
upper cap both bound the same derivative on the child interval; their
minimum is also a valid upper cap. This does not choose a midpoint or
replace the derivative with an interval endpoint.

Separately, the explicit nominal displacement is evaluated through third
radial order. Its construction uses a polynomial nominal background,
rational packet and time profiles, and a polynomial cutoff. The nominal
API supplies order-six jets; the required scaled ordinary third derivative
of the displacement is `6 ell^3 c[3]`. This step requires no new sixth-order
*actual* ODE or packet bound. No unknown Taylor remainder is differentiated.

Let `M3_i` be the new full-cover nominal third-derivative component caps,
`j_i` the retained maximum isolated inner curvature-jump allowance, and
`delta=2^-30`. The positive even kernel has first absolute moment
`35delta/128` and half mass `1/2`. The same nominal displacement is C1, so
its second distributional derivative has no delta masses. A convolution
footprint meets at most one join. The first inner nominal curvature jump
is exactly zero; the outer quadratic continuation matches curvature.
Consequently the same componentwise convolution estimate is

```
s_i=(35delta/128) M3_i+j_i/2.
```

The new complete continuous initial error cap is therefore

```
Eurr=sum_ch P_ch + sqrt(sum_i Gram_i s_i^2).
```

Every Gram weight is positive and the complete 23 slots are retained. The
complete rounded-node/P1 velocity-gradient error remains fixed; a smaller
paired velocity comparison cannot replace its interpolation and smoothing
charges. Recomputing these initial estimates changes no source, candidate
or initial target data.

## Transfer to the same evolved error

Let `e=U-W1_quadratic`, `T=3/2` and `q=3409/8000000`. Write `a,b,d,m`
for the retained actual `H` operator caps of `A,B,D,M`, and let

```
c=a+b+3d/4,
comm=cap(A_x)+cap(B_x)+3cap(D_x)/4.
```

The full absolute candidate derivative/forcing ledger, actual original
`F_x`, C1 error cap `E`, reduced defect cap `F0` and differentiated defect
cap `F1` remain those of the original-source proof. A bound on an error
cannot replace an absolute candidate curvature bound in that ledger.
The original candidate still has `r1(0,r)=0`; its intermediate corrections'
generally nonzero initial accelerations are retained in `F1`.

Only the initial curvature-error cap is replaced. With the other initial
errors `eu,eur,ev,evr` fixed, the differentiated error initial data give

```
w0=Eurr+a ev+b eur+d eu,
J=max((ev+T w0)/3, evr+w0).
```

The same actual shell source cancels from the error law: `g=g_x=M_x=0`.
The physical-side differentiated comparison yields

```
cap(e_xx)=cap(e_xr)=(J+T(F1+comm E))/(1-q),
cap(e_rr)=cap(e_xx)+c E+F0.
```

All coefficients of the replaced upper cap are nonnegative. The checker
reassembles the prior full ledger and reproduces the exact original and
PR103 bounds before substituting the new initial estimate. An independent
audit uses the exact majorant difference `(J_new-J_old)/(1-q)` to check
the transfer. No rounded published bound is subtracted.

The evolved domain is

```
D+={0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x}, x=t/ell, ell=1e-6.
```

The continuous estimates include physical-side one-sided shell and
domain-edge traces. They require the same target regularity, source,
harmonic-neighbourhood and compatible initial/corner premises as the
original proof. No C2 extension across the artificial reflected shell is
claimed.

## A positive actual-target scale

The authenticated actual initial endpoint calculation, independent of any
candidate's affine extension, supplies

```
u0_rr,w0(7/4) in [-3.249080947030276e33,-3.2490809468614435e33].
```

The independently replayed autonomous-ODE route gives a wider overlapping
enclosure whose upper endpoint is still strictly below `-3.249e33`.
The `w:0` Gram weight equals one. Thus the actual target initial curvature
has the positive lower scale

```
sup_[0,7/4] ||U0_rr||_H > L, L=3.249e33.
```

This uses the actual derivative, not the historical affine-candidate error
lower bound `1.6245e33` printed in the same old checkpoint. The input roles
must not be interchanged.

The retained formal target has `U(0,r)=U0(r)`. Under its retained classical
physical-side C2 trace, differentiating this identity tangentially to the
initial surface gives `U_rr(0,r)=U0_rr(r)`, including the one-sided endpoint
trace. The point `(0,7/4)` belongs to the declared closed `D+`. Continuity
gives the same supremum lower bound if one instead takes a supremum over
the interior. Hence

```
A_rr=sup_D+ ||U_rr||_H > L.
```

The target scale is finite under the same continuous-on-compact-domain
regularity premise. Its strictly positive lower bound makes the following
ratios well defined.

## Precisely defined global coordinate ratios

If `Brr` bounds `sup_D+ ||e_rr||_H`, then

```
sup_D+ ||e_rr||_H / sup_D+ ||U_rr||_H < Brr/L.
```

Define the coordinate second-derivative seminorm

```
|V|_(2,H)=max_ab_in_{xx,xr,rr} sup_D+ ||V_ab||_H.
```

It is a seminorm: constant and affine fields can have zero second
derivatives while their values or first derivatives are nonzero. For this
specific target the denominator is positive because `|U|_(2,H)>=A_rr>L`.
Thus the three complete error bounds yield

```
|U-W1_quadratic|_(2,H)/|U|_(2,H) < max(Bxx,Bxr,Brr)/L.
```

The analogous initial ratio is below `Eurr/L`. For a transparent baseline,
the already published PR103 bounds imply initial ratio `<3.648e-6` and
evolved rr/max-second-seminorm ratios `<3.654e-6`. The new calculations
evaluate the improved exact rational quotients and outward readable
numbers, retaining the same conservative scale `L`.

These are ratios of global suprema in the registered coordinates and
positive component weights. They do not bound error divided by the target
at every point. The denominator can be dominated by the large initial
outer-endpoint derivative and can mask large relative errors in smaller
components or local regions. In particular, the ratios are compatible
with the nine nonzero inner-join errors from PR103. Nor does the rr lower
scale separately normalize unknown `sup ||U_xx||_H` or `sup ||U_xr||_H`.

No full C2 norm, invariant curvature normalization, covariant or physical
relative accuracy, observed signal, finite-amplitude nonlinear family or
Big Bang-origin evidence follows. The displayed target remains a formal
second perturbative coefficient under the retained mathematical premises.

## Numerical consequences and a concrete local calibration

The two new initial majorants give a complete initial radial second-error
cap strictly below `6.495e26`. The unchanged error transfer gives each
evolved coordinate second-error cap strictly below `6.643e26`. These upper
certificates are respectively more than 18.23 and 17.86 times tighter than
the PR103 certificates. They refine the proof of the error, not the
underlying reconstruction or its true error.

Using the actual-target scale above gives initial global ratio `<2e-7`
and evolved rr/max-second-seminorm ratios `<2.045e-7`. These are conditional
global coordinate ratios in the stated H norm. Calling this a pointwise
accuracy statement would be false for this very reconstruction.

Indeed, the authenticated PR103 join calculation at initial
`(x,r)=(0,65/64)` supplies complete actual and error intervals for all23
components. Let `A_i` be the maximum absolute endpoint of the actual
component interval. Let `E_lower` be the certified strict lower bound on
the actual pointwise error H norm. Then

```
||U0_rr(65/64)||_H <= sqrt(sum_i Gram_i A_i^2),
||e0_rr(65/64)||_H > E_lower.
```

The exact rational calculation proves

```
E_lower^2 > (999/1000)^2 sum_i Gram_i A_i^2.
```

The target norm is nonzero, since its scalar mean is positive. Therefore
the pointwise initial **same-H-norm** relative error is greater than
`0.999`, or 99.9 percent, at this point. A rational upper square root of
the actual norm yields a sharper verified ratio lower bound of about
`0.9999994986`. The scalar-mean component ratio alone independently exceeds
`0.999` too. The report distinguishes these two quantities.

This is a concrete calibration using the same H norm as the global ratios,
not merely an abstract counterexample. The small ratio of global suprema
and the large local relative error coexist because the former can be
normalized by the much larger outer-endpoint derivative. The absolute
initial-point error here is small relative to that global scale; the
calculation is not a physical model obstruction or a nonlinear result.
It shows why local observables and physically justified norms still need
their own error analysis.

## Reproducibility and controls

All numerical decisions use exact rational arithmetic and explicit runtime
checks that remain active under optimized Python. Source and report roles
are authenticated with SHA256; producer source closures are checked before
scientific imports and on readback. Independent implementations audit the
new derivative leaves, cover, convolution arithmetic, exact transfer and
normalization.

Corrupted input-file controls reject lower-bound or point-value roles,
unsuccessful statuses, changed candidates, incomplete component lists,
altered Gram weights and negative error upper caps. Executed normalization
faults reject an affine-error denominator, a wrong scalar weight, an
omitted classical target trace, and unjustified pointwise, separate-xx,
full-C2 or physical relative claims. These are internal mathematical and
reproducibility checks, not external specialist endorsement or a
proof-assistant certificate.
