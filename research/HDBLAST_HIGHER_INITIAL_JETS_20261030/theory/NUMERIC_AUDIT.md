# Independent audit of the actual time-coefficient ledger

Reviewed `check_time_coefficients.py` and its report with SHA256
`6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d`.
The analytic D_x primitive coefficients agree with the independently
differentiated original full23 assembly. There is no derivative of a frozen
midpoint or of an interval endpoint.

The review uses separate Fraction interval arithmetic, with no import of the
new producer. It authenticates and imports only the original primitive
operator. For the global rectangle and selected shell point, it independently
recomputes all eight coefficient intervals from the reported actual rho, a,
p, tanh and sech ranges. Each exact rational interval lies inside its reported
outward enclosure. It then reconstructs all 529 matrix entries from the eight
primitive coefficients and verifies their containment in the reported D_x
matrix. All ten degree-block weighted Frobenius squares are independently
recomputed exactly and agree with the report. Their square-root caps and the
maximum-over-orthogonal-blocks operator bound are valid.

The reported A_x is a single scalar times the identity. Its exact Fraction
diagonal intervals avoid introducing dyadic rounding outside the inclusive
sharp magnitude cap `3e-12`; no dimension-dependent Frobenius factor is
needed. B_x and M_x are identically zero matrices. The complete global
D_x operator cap is strictly below `2.354e-16`, and

```
||A_x|| + (3/4)||D_x|| < 3.000177e-12.
```

This last expression uses the separately inherited zero-wedge value bound
`||e||<=3X(e)/4`. It is the homogeneous coefficient contribution in the
differentiated error equation, not the complete differentiated forcing.

## Hyperbolic range justification

For nonnegative t, `tanh(t)<=t` because its derivative is at most one. Also

```
d/dt [tanh(t)-t+t^3/3] = t^2-tanh(t)^2 >= 0,
```

so `tanh(t)>=t-t^3/3`. On the declared tiny interval, t is below one and
the lower polynomial is increasing, justifying use of its left endpoint.
Since `sech(t)<=1`, `sech(t)*tanh(t)<=t`, and therefore

```
d/dt [sech(t)-1+t^2/2] = t-sech(t)*tanh(t) >= 0.
```

This proves `sech(t)>=1-t^2/2`. The report's tanh and sech boxes contain
these exact analytic brackets for both domains. Squaring and multiplying
their interval ranges may lose correlation, but only widens the enclosure.

The actual background intervals are inherited from the authenticated
background context; this audit does not independently recertify that
background. Both audit reports are byte-identical in normal and optimized
execution. Neither the ledger nor its audit provides F_x, higher initial
error bounds, a quantitative C2 theorem, a full natural Einstein residual or
physical evidence for a Big Bang mechanism.
