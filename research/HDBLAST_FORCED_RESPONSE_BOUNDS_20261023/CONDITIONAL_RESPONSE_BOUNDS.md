# Quantitative bounds for the auxiliary forced response

This is a standard conditional consequence of the retained causal inverse and
the published continuous source certificate for one five-dimensional
Einstein–scalar shell model. It supplies a lower bound for a full comparison
norm, together with the retained upper bound and prospective Neumann tails.
It does not compute a response waveform or establish a physical force.

## Retained construction

On `D={0<=x<=3, r>=0, x+r<=8}`, retain the complete 23-amplitude system

```text
S=KF+LG,     W=(I-T)^(-1)S,
TV=K(A V_x+B V_r+D0 V)+L(M P_N V|r=0),
||T||_(X->X)<=q=2109/1000000<1,
||V||_X=max(sup_D |V|_H/3, sup_D |V_x|_H, sup_D |V_r|_H).
```

Here H integrates the positive metric Frobenius square plus the scalar square
against probability-normalized sphere measure. The fixed norm is a comparison
norm, not a gravitational energy. The full coupled Robin term, the physical
shell sources and the normal harmonic shell source are retained in T and G.
The sources and compatible C1 class have the registered flat initial wedge;
W has zero initial data. The decomposition is auxiliary: W and the homogeneous
comparison are not separately asserted to solve the physical quadratic
constraint system.

Write `P=partial_x^2-partial_r^2`. The registered operators are

```text
Kf(x,r)=(1/2) integral_(s=0)^x integral_(v=r-x+s)^(r+x-s) f_ext(s,v) dv ds,
Lg(x,r)=-integral_0^max(x-r,0) g(s) ds,
```

where the extension uses the registered odd/even reflection. On the physical
interior, `P(KF)=F` and `P(LG)=0` in distributions. One normalization check is
the retarded fundamental solution `E=(1/2)H(x+r)H(x-r)`: with null coordinates
u=x+r,v=x-r, `P=4 partial_u partial_v` and
`PE=2 delta(u)delta(v)=delta(x)delta(r)`. The shell lift depends on x-r alone,
so its principal wave operator vanishes distributionally. These statements
need no classical second derivative of S or derivative estimate for F.
They do not claim homogeneity of the lift for the full coupled operator.

## Continuous source witness and mean functional

Let lambda take the yy probability sphere mean of the metric block. Its H
dual norm is exactly one: Cauchy–Schwarz gives the upper bound, attained by
the constant pure yy tensor of value one in the selected space. The fixed
functional commutes with x,r derivatives. Consequently

```text
|lambda(S_x)|<=||S||_X,     |lambda(S_r)|<=||S||_X.
```

There is no value-weight factor three or angular point-evaluation factor in
these derivative inequalities. In the retained scalar basis `(1,Z,Fang)`,
the nonconstant entries have probability mean zero. In detail,
`Z=n0^2+n1^2-1/2` and
`Fang=(n0^2-n1^2)^2-(n0^2+n1^2)/2+1/12`, using sphere moments
`mean(ni^2)=1/4`, `mean(ni^4)=1/8`, and
`mean(ni^2 nj^2)=1/24` for i!=j. Thus the registered yy:0 coefficient is
lambda(F), in slot K:0, without a frame rescaling or raised time index.

The published actual-source report encloses the already normalized metric
forcing `F_metric=-2 ell^2 rho^2 QEred` continuously on the square

```text
R=[1/4,1025/4096] x [5/4,5121/4096],     h=1/4096,
lambda(F)>=c_minus
  =44085568717512092120254222950792586054530845
   /2596148429267413814265248164610048
  >16000000000=:c.
```

Its exact report entry is `bulk_witness.sources.metric_forcing.yy[0][0]`.
This is a bound for the actual first-pair source, not a comparison solution
or a finite sample. R is strictly inside the physical interior. Its minimum
r-x is4095/4096, so LG vanishes on a neighborhood of R, and the backward
cones from that neighborhood stay at positive depth. The principal identity
there can alternatively be checked with ordinary unreflected Duhamel
integration, approximating the continuous source uniformly by smooth sources
and passing to distributions. No shell cancellation or new source-derivative
assumption enters the witness.

## Weak test-function lower bound

For every nonnegative chi in `C_c^infinity(R interior)`, the weak equation
PS=F and one integration by parts give

```text
c_minus integral chi <= integral lambda(F) chi
  =-integral lambda(S_x) chi_x+integral lambda(S_r) chi_r
  <=||S||_X (||chi_x||_L1+||chi_r||_L1).
```

Take translated copies of a smooth one-dimensional plateau cutoff: it equals
one on the concentric interval of length h-4epsilon, vanishes outside the
concentric interval of length h-2epsilon, and has monotone increasing and
decreasing transition ramps, where0<epsilon<h/4. Its total variation is two;
its integral m_epsilon lies between h-4epsilon and h-2epsilon. The product
of the two cutoffs has mass m_epsilon^2 and summed first-derivative L1 norm
4m_epsilon. Therefore

```text
||S||_X>=c_minus*m_epsilon/4,
||S||_X>=c_minus*h/4>c*h/4=1953125/2.
```

The first limiting inequality is nonstrict. The strict rounded inequality
uses the separate exact slack c_minus>c after taking the limit.
Since S=(I-T)W, `||S||_X<=(1+q)||W||_X`, and hence

```text
||W||_X >= c_minus*h/[4(1+q)]
        > 976562500000/1002109
        > 974000.
```

The sharper nonstrict endpoint is approximately1034265.5; the rounded strict
endpoint is approximately974507.3. The transfer uses1+q. The factor1-q belongs
to the different resolvent upper estimate below. This proves that the full
auxiliary W is not identically zero, but does not identify a nonzero response
channel, sign, value-only amplitude, point or waveform.

## Retained upper bound and prospective finite-sum errors

The authenticated reports give `F_H<=17387513156573`, `G_H<=3851400` and
`||Kf||_X<=3 sup|f|_H`, `||Lg||_X<=sup|g|_H`. Consequently

```text
S_cap=3 F_H+G_H=52162543321119,
||W||_X<=S_cap/(1-q)
       =52162543321119000000/997891<5.228e13.
```

For the exact true-operator partial sum `V_N=sum_(j=0)^(N-1) T^j S`,

```text
(I-T)V_N-S=-T^N S,
||W-V_N||_X<=q^N S_cap/(1-q).
```

In particular, the four-term tail is

```text
1031965259713929363466110159/997891000000000000000000 < 1035.
```

These are bounds for unevaluated exact partial sums. No numerical candidate,
quadrature, continuous response residual or complete shell-error certificate
is supplied by this arithmetic. Finite partial sums need not satisfy the
exact coupled Robin law. Source approximation, operator approximation,
integration, reconstruction, initial and shell defects must all be included
in any later complete numerical certificate. Tighter source enclosures may
be needed to attain a chosen finite error target.

## Bridge from the separately certified local mean

The separate actual-source cone computation in `compute/EXTERIOR_RESPONSE_PROOF.md`
certifies a positive point interval `[L,U]` for the same global auxiliary W at
`(x,r)=(1/4,5/4)`. Its exact endpoints are recorded in
`compute/ACTUAL_CONE_RESPONSE_NORMAL.json`. Since the X norm weights values by
one third and the mean functional has H dual norm one,

```text
|lambda(W)(1/4,5/4)|<=3||W||_X,
||W||_X>=L/3>8800000.
```

The factor three is essential here, in contrast to the derivative-only weak
test argument above. The local numerical certificate strengthens the global
norm lower bound while the weak-test estimate remains an independent check.
It certifies a specified coordinate mean; it supplies no physical-force or
invariant-observable interpretation. `check_bounds.py` authenticates the exact
accepted numerical report and checks this arithmetic bridge. The separate
producer and reviewer must verify the source-cell and coupling certificate;
the bridge itself does not recalculate it.

## Evidence and limits

`SOURCE_PROVENANCE.json` links every copied dependency to its published path
at commit d110d99f054ad297fdcc0a7cec07499341bae4b8. `INPUT_PINS.json`
authenticates the bytes. `check_bounds.py` verifies those hashes, extracts the
registered witness and q/source caps, checks the scalar basis means and slot
convention, and evaluates the rational bounds. It does not reprove the
upstream geometric/source estimates or the analytic arguments above.

These are standard conditional weak-operator consequences in a fixed gauge,
selected angular space and initial prescription. They give neither a physical
force nor an invariant observable, nonlinear continuation, observational
support for a higher-dimensional origin of the Big Bang, or a novelty claim.
