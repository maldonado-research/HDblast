# Independent internal review of the filtered candidate derivative bounds

30 October 2026. This review concerns the distinct auxiliary candidate
`Z_epsilon`, with `epsilon=1/1024`. It is not an external specialist review,
a novelty assessment, or a statement of second-derivative accuracy to the
unsmoothed target.

## Analytic estimate checked

For a physical-side principal-wave construction

```
W_xx-W_rr=f,
PD W=0,        PN W_r=M PN W+G,
W(0,r)=u_s(r), W_x(0,r)=v_s(r),
```

differentiate in time. With the same constant actual `M`, the differentiated
field `V=W_x` has source `f_x`, shell source `G_x`, and initial data

```
V(0,r)=v_s(r),
V_x(0,r)=u_s,rr(r)+f(0,r).
```

The initial forcing term has coefficient one; it cannot be discarded or
halved. A direct polynomial wave fixture checks this convention. If
`U2,V0,V1` bound respectively `u_s,rr,v_s,v_s,r`, and `f00` bounds `f(0)`,
the free-wave bound is

```
J=max((V0+(3/2)(U2+f00))/3, V1+U2+f00).
```

The principal propagator bounds `X(Kf_x)<=3 sup||f_x||_H/2` and
`X(LG_x)<=sup||G_x||_H` apply. The field and source vanish in the inherited
zero wedge, so the same active shell length `1/4` gives the actual Robin
inverse `1/(1-||M||_H/4)`. This proves the displayed bound for `X(W_x)`.
It controls both `W_xx` and `W_xr`; the principal equation then bounds
`W_rr` by the sum of the `W_xx` and forcing bounds.

For `W0_quadratic`, the frozen forcing is differentiated before estimating:

```
f0_x=A0 Bs_xx+B0 Bs_xr+D0 Bs_x.
```

For `Z_epsilon`, the complete coupled forcing is differentiated before its
output reflection:

```
f1_x=A W0_xx+B W0_xr+D W0_x
     +A_x W0_x+B_x W0_r+D_x W0+(F_epsilon)_x.
```

The checked coefficient expression retains the `D W0_x` term and the time
commutator. It uses `||W0||<=3 X(W0)/4`, rather than an unjustified equality.
The analytic source derivative comes from the smoothing kernel, not from
assuming a bound on the actual `F_x`.

These estimates are **one-sided on the physical domain**. An odd reflected
Dirichlet source can jump at the shell. No assertion of classical C2
regularity across the doubled shell is justified by this argument. Physical
side C2 regularity also retains the stated qualitative source continuity,
initial smoothing and flat-corner compatibility assumptions.

## Initial data and arithmetic checks

The quadratic exterior continuation preserves the authenticated global
componentwise second-derivative bounds. The displacement remains C1,1;
convolution with the even positive kernel therefore has the same uniform
second-derivative bound. Its value and first-derivative smoothing charges
use the exact kernel moments `delta^2/9` and `35 delta/128`.

The smoothed P1 velocity uses the authenticated rounded nodal values
themselves. Convexity bounds its nodal norm and all piecewise slopes in the
positive full23 H norm. The exterior final-slope continuation has the same
Lipschitz cap. The velocity value charge is therefore valid at the outer
endpoint as well as in the interior.

`check_candidate_bound_audit.py` independently reconstructs these initial
caps from the raw profiles and recomputes the two differentiated wave
ledgers from authenticated prior inputs. It imports no scientific producer.
All six candidate coordinate-jet caps match exactly. In particular,

```
sup ||Z_xx||_H, sup ||Z_xr||_H, sup ||Z_rr||_H <1.315e35.
```

The checker also independently charges the source approximation in the C1
inverse and verifies `X(U-Z_epsilon)<3.541e25`. Its 38 named checks include
six deliberately wrong-formula controls. Normal and optimized Python reports
are byte-identical.

## Scope of the result

The first-derivative change from the previous quadratic candidate is bounded
explicitly. The actual target source remains fixed. The new initial target
error forcing is generally `F(0,r)`, even though the filtered source and its
first time derivative vanish initially. This discrepancy is not a shell
corner failure because the actual source vanishes in an open corner collar.

The derivative caps are finite absolute bounds for a numerical candidate.
Their large size provides no claim of a small equation residual, physical
accuracy, natural constraint propagation, exact nonlinear realization,
cosmological prediction or evidence for a higher-dimensional Big Bang.
The actual gauge C4/full `F_x` and quantitative C2 target-error obligations
remain separate. The prior point tables cannot be inherited as exactly
unchanged after this source filtering.

**Assessment:** the candidate absolute second-jet and explicitly charged C1
ledger pass this independent internal review, with the physical-side and
conditional-regularity qualifications stated above. A separate full natural
residual claim requires its own formula and arithmetic review.
