# Independent review of the inner-join curvature diagnostic

30 October 2026. This is an internal mathematical and executable review of a
conditional coordinate initial-data calculation. It is not external specialist
peer review or validation of a physical Big Bang mechanism.

The fixed actual initial target, the PR100/PR102 `W1_quadratic` reconstruction,
its kernel width `delta_u=2^-30`, and the original bulk and shell sources are
preserved. The quadratic repair at the outer endpoint does not enter either
inner kernel footprint.

The diagnostic is mathematically sound. For a C1 nominal displacement with
piecewise C3 restrictions, its distributional second derivative has no delta
term. An even, positive normalized kernel assigns half mass to each side of an
isolated join. Its positive-half first moment is `35 delta_u/256`. Thus

```
u_s,rr(rc)=(kappa_minus+kappa_plus)/2+rem,
abs(rem_i)<=35 delta_u/256*(M3_minus_i+M3_plus_i).
```

The nominal packet cutoff matches through derivative order four. The registered
initial displacement uses packet derivatives through three, giving the C1
matching needed here. The existing cutoff does not match through order five;
its fifth derivative at both ends is 15120 in its normalized cutoff variable.
Exact forced branch evaluations, rather than offset point stencils, supply both
one-sided traces. The local nominal third-derivative bounds cover each entire
half of the kernel footprint. The two inner joins and the repaired outer
endpoint are separated by more than two kernel widths.

The independent checker does not import the new diagnostic producer or its
order-five paired helper. It authenticates the retained original initial
producer and all of that producer's source roots, then reconstructs actual
background jets by an autonomous ODE Taylor recurrence about the authenticated
actual point state. It uses the actual bump derivative polynomials and actual
exponential jets through order five. It separately constructs the forced nominal
branches and their local third-derivative bounds. Reuse of the authenticated
registered full23 reconstruction and interval arithmetic is explicit; this is
an independent derivative route, not wholly independent model implementation.
The actual packet derivative polynomial recurrence, kernel moments, cutoff
matching and ordinary scaled derivative factorials are also checked exactly.

The independent and producer calculations have overlapping actual-curvature and
actual-error intervals in all 46 component/join rows. They match the nominal
one-sided traces, local third-derivative caps and convolution remainders exactly.
They find the same nine sign-excluding actual-error rows:

- At `65/64`, `J:0,J:2,J:4,w:0,w:2,w:4`. All nominal curvature traces at this join
  are exactly zero. The actual packet tail is nonzero, so this result is not
  attributed to a nominal curvature jump. Independent direct-ODE scalar-mean
  error lies between `0.6114227665` and `0.6114230732`.
- At `66/64=33/32`, `w:0,w:2,w:4`. Independent direct-ODE scalar-mean error lies
  between `-5.194525e20` and `-5.192668e20`. The producer's paired route gives a
  tighter compatible interval. These actual-error intervals, rather than the
  nonzero nominal jump allowance alone, establish the mismatch.

For every producer component the checker independently computes the minimum
absolute value of its signed error interval, sums all 23 weighted squares with
the complete positive Gram matrix, and verifies its reported strict H-norm lower
bound. Each of the nine nonzero producer errors also has a zero-width limit
with the same sign. No nonzero error is inferred from an upper bound. Intervals containing
zero remain inconclusive; they do not prove the corresponding error is zero.

The final independent executable has 323 named check groups and four explicit
wrong-formula controls in each execution mode. Normal and optimized reports are
byte identical. It authenticates the producer diagnostic report by SHA256 and
rechecks scientific source closure after evaluation.

These are fixed-width initial radial-curvature error bounds for one auxiliary
coordinate reconstruction. They do not by themselves propagate an error through
time, obstruct the physical model, establish relative physical accuracy, prove
a finite-amplitude nonlinear family or support a cosmological-origin claim.
