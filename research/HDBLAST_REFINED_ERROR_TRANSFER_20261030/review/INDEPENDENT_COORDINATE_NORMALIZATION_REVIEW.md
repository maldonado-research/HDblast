# Independent review of combined C2 transfer and coordinate normalization

30 October 2026. This review concerns the same original-source
`W1_quadratic` reconstruction and formal target as PR102/PR103. The new
initial-data intervals change upper majorants, rather than the functions
or their true error. New derivative enclosures are checked by the separate
initial-refinement audit; this review independently reconstructs their
combined ledger, the C2 transfer and the normalization arithmetic.

## Combined initial ledger and exact C2 delta

All three paired actual-minus-nominal channels remain necessary. If their
uniform displacement-curvature H caps are `B_ch`, the combined bound is

```
Bpaired=sum_ch B_ch,
smoothing_i=(35delta/128) M3_i+inner_jump_i/2,
Inew=Bpaired+outward_sqrt(sum_i Gram_i smoothing_i^2).
```

The refined nominal `M3_i` are bounds on an explicit analytic nominal
function. They do not assume a new sixth derivative of an actual ODE
solution. The unchanged inner jumps, smoothing width and quadratic outer
continuation are retained. The independent checker rebuilds this vector
and positive Gram square directly from authenticated paired/nominal
reports, with no producer code imports.

Let `Iold` denote the PR103 initial curvature-error cap. All other initial
errors, the C1 error cap, actual source derivative, absolute candidate
derivative bounds and reduced-defect bounds are retained. Therefore the
new initial acceleration majorant equals

```
wnew=wold-Iold+Inew,
Jnew=max((v0+(3/2)wnew)/3, v0_r+wnew).
```

The exact difference from the preceding authenticated C2 theorem is

```
new_cap_ab=old_cap_ab+(Jnew-Jold)/(1-q),
q=3409/8000000, ab in {xx,xr,rr}.
```

The independent audit checks these identities with exact rational
arithmetic. It does not subtract rounded advertised bounds. The same
formula holds even if the maximizing branch in `J` changes, provided the
maximum itself is recomputed. Retaining the feedback inverse and all
other initial/source terms is essential.

## A denominator supplied by the actual target

The authenticated actual initial endpoint interval for scalar mean `w:0`
in the PR100 higher-initial report is approximately

```
u0_rr,w0(7/4) in [-3.249080947030276e33,-3.2490809468614435e33].
```

Its full23 positive Gram weight is one. The independent autonomous-ODE
endpoint checker was replayed in isolated normal and optimized Python.
Both reproduce the retained `ENDPOINT.json` byte for byte, SHA256
`e21d2cb172893495b557bbcdcee198d33d21a7ba287c5215d7319eff3d9dc9cc`.
That separate actual-jet route gives a wider overlapping interval,
approximately `[-3.2490825735884044e33,-3.249079962707768e33]`. Both
routes establish the conservative strict lower scale

```
sup_[0,7/4] ||actual u0_rr||_H > L=3.249e33.
```

The actual derivative is the denominator datum. The old affine-candidate
curvature error, the nominal curvature and a numerical interval midpoint
would have different roles and must not be substituted for it.

The original target fixed-point construction uses the actual initial data
in its free wave. Its bulk and boundary lifts have zero Cauchy data, so
`U(0,r)=U0(r)`. Under the retained physical-side C2 regularity, tangential
differentiation gives `U_rr(0,r)=U0_rr(r)`, including the one-sided trace
at `R=7/4`. The declared closed physical domain

```
D+={0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x}
```

contains `(0,7/4)`. Hence `sup_D+ ||U_rr||_H>L`. If only interior
suprema are used, continuity to the same corner gives the same lower
bound. C2 continuity on the compact physical domain makes these suprema
finite; a merely weak initial trace would be insufficient.

## What the resulting ratios mean

Define exactly

```
A0=sup_[0,7/4] ||U0_rr||_H,
Arr=sup_D+ ||U_rr||_H,
|V|_(2,H)=max_ab_in_{xx,xr,rr} sup_D+ ||V_ab||_H.
```

Because each denominator is strictly larger than `L`, positive numerator
caps give the strict bounds

```
sup ||U0_rr-W1q_rr(0,.)||_H / A0 < Inew/L,
sup ||U_rr-W1q_rr||_H / Arr < cap_rr/L,
|U-W1q|_(2,H) / |U|_(2,H) < max_ab cap_ab/L.
```

These are ratios of global coordinate seminorms. A lower bound for `Arr`
does not establish a denominator for `sup||U_xx||_H` or `sup||U_xr||_H`
separately. The maximum-second-derivative expression is a seminorm, not
the full C2 norm: a nonzero affine field has zero second-derivative
seminorm. The particular target has a nonzero denominator, so its
relative seminorm is well-defined.

## Quantitative local calibration using the same H norm

The distinction from pointwise relative accuracy is substantial for this
actual unchanged reconstruction. The pinned PR103 join report at
`(x,r)=(0,65/64)` supplies positive scalar-mean actual and error intervals
both near `0.611423`. Their exact interval quotient lower bound exceeds
`999/1000`.

The same statement also holds for the local complete H norm, rather than
only one component. Let

```
A2=sum_i Gram_i mag(actual_u0_rr_interval_i)^2,
Lerr=the authenticated strict local H-error lower bound.
```

Every actual component lies in its interval, so the actual local target
H norm squared is at most `A2`. The scalar mean's positive lower interval
also ensures that local denominator is nonzero. The independent exact
rational comparison verifies

```
Lerr^2 > (999/1000)^2 A2.
```

Consequently

```
||U0_rr(65/64)-W1q_rr(0,65/64)||_H / ||U0_rr(65/64)||_H > 0.999.
```

The earlier strict error lower bound is approximately `0.6942779566`
and the corresponding complete actual H upper bound is approximately
`0.6942783047`. This is compatible with a tiny ratio of global suprema:
the global denominator is at least `3.249e33`, whereas this local target
is tiny on that coordinate scale. Neither the local mismatch nor the
global ratio is a physical force, model obstruction or observed signal.

The executable checks the denominator role and exact interval arithmetic,
rebuilds all producer ratios and rejects missing trace premises, pointwise
relabeling, separate xx/xr denominator claims, full-C2 relabeling, physical
accuracy claims and the old affine-error denominator role. It authenticates
all supplied reports and reads them back after evaluation. Deterministic
normal/optimized audit reports record final numerical caps and check counts.

These internal mathematical and implementation reviews retain the actual
background, master/gauge, source-map, compatible initial/corner and physical
boundary premises. They are not external specialist endorsement, a proof
of novelty, a physical/covariant relative-accuracy result, a finite-amplitude
nonlinear spacetime or evidence for a higher-dimensional Big Bang origin.
