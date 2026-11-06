# A distinguished weight member and its limits

The pinned historical sources specify a radial defining function and
parameters only conditionally: qweight=d²,
w=d^(2a_w-4) exp(-beta/d), a_w an admissible natural integer and beta>0.
The original theorem and later binding give no numerical shape for d or
beta/d. Neither the annulus nor k=12 resolves that choice. The numerical
member below is newly instantiated for this calculation.

Fix ell=10^-6, A={ell/2<z<5ell/2}xS3 and seed support K=[ell,2ell]. Set

```
r=(z-ell/2)(5ell/2-z)/ell²,
kappa=ell/100, d=kappa r, beta=kappa, a_w=2,
qweight=d², w=exp(-1/r),
mu_H=w²qweight⁴=kappa^8 r^8 exp(-2/r),
mu_M=w²qweight²=kappa^4 r^4 exp(-2/r).
```

Then 0<r<=1 in A and r>=3/4 on K. The square-root weight contrasts are
(4/3)^4 exp(1/3)<4.42 and (4/3)^2 exp(1/3)<2.49. The scalar proof gives
scaled-cover inclusion and comparability estimates. Smooth multiplier,
adjoint API, dotted closure and existence hypotheses of the retained
local theorem continue to apply. This specification does not evaluate a
regularity or evolution constant.

Replacing d,beta by t d,t beta with fixed t>0 multiplies mu_H and mu_M
by positive block constants. Those constants cancel in
mu P*(P mu P*)^-1. After the corresponding constant change of potential
variables, the physical correction range and the local implicit curve
are the same on their common admitted interval, by local uniqueness.
For each fixed t the finite-order weighted norms are equivalent. Thus
kappa can be decreased for a cover normalization without changing this
member's physical correction. Changing r, a_w or beta/d generally selects
a different member; no equivalence with the historical unspecified member
is asserted.

For any fixed admissible radial mu_H, the comparison bound remains

```
||h2||_prob < 5.58e39 C_H,
C_H=sqrt(sup_A mu_H/inf_K mu_H).
```

The analogous momentum statement is
||tau2||_prob <=36 C_M ||s_M||_prob with
C_M=sqrt(sup_A mu_M/inf_K mu_M). Compact K makes these contrasts finite
for each admissible choice, but there is no unique historical numerical
value to report. The scalar/density second data remain zero in the
canonical initial prescription; normal scalar momentum and later scalar
fields need their actual conversions and evolution.

All statements concern amplitude derivatives at the background and the
fixed initial slice. They do not establish nonlinear control on any
specified nonzero amplitude interval.
