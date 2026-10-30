# Refining the initial-curvature error certificate without changing the candidate

This calculation keeps the PR100 `W1_quadratic` / `W1_dagger` initial data,
quadratic outer continuation, smoothing widths, fixed P1 velocity, underlying
actual target and original sources unchanged. It refines an interval proof,
not the initial approximation itself. A smaller upper bound alone does not
prove that the actual error has decreased.

## Dominant cell and exact subdivision

The retained continuous paired comparison has 48 cells on `[1,7/4]`. Its
second telescoping channel on `[65/64,66/64]` gives a displacement-curvature
H upper bound approximately `6.2189008408e33`. Every other cell's same-channel
H cap is below `2.360e26`. The single coarse cutoff cell therefore dominates
the old `6.2189049731e33` complete initial-curvature certificate.

The checker uses the original authenticated order-five `HigherInitial`
algebra. In each tested cell `[a,b]`, it computes the 23 component caps

```
b[ch,i] = mag(2 ell^2 u_i.d[ch].c[2]),
h[ch] = outward_sqrt(sum_i Gram[i] b[ch,i]^2).
```

The channels are the same three telescoping actual-minus-nominal channels
in the retained initial producer. They preserve shared factors before
interval evaluation. The radial derivative scaling follows the authenticated
jet convention; no unknown Taylor remainder is differentiated.

A cell is bisected at its exact rational midpoint whenever any displacement
channel H cap exceeds `1e28`. All accepted leaves and rejected parent cells
are recorded. The resulting ten leaves, with depths from one through seven,
cover exactly `[65/64,66/64]`. All 47 other old cell H caps are retained from
the pinned old report. Thus the new cover has 57 cells, is ordered, has no
gaps or overlaps, and covers exactly the same active footprint. The inactive
inner collar has zero actual and nominal data. No point samples replace a
continuous interval.

For each channel, take the maximum of its cell H caps over the whole cover.
This bounds the uniform H norm of that channel. Sum the three channel bounds
by the triangle inequality. This is valid even if different cells attain the
largest bounds in different channels. It also avoids the extra overestimate
of separately maximizing each component over the entire domain before
packing H.

## Same convolution and outer continuation

Let `delta=2^-30`. The positive even smoothing kernel has
`E|s|=35delta/128` and half mass `1/2`. The nominal displacement is C1, so its
distributional second derivative has no delta term at either inner join.
The first join, `65/64`, has exactly zero nominal curvature jump in all23
slots. Use the retained signed-trace jump *upper allowances* at `66/64`;
they do not themselves establish nonzero actual error. The original global
piecewise nominal third-derivative component caps are retained unchanged.

The same pointwise smoothing allowance is therefore

```
smoothing[i] = (35delta/128) M3[i] + max_inner_jump[i]/2.
```

The footprint contains at most one inner join. At `R=7/4`, the PR100
quadratic continuation has curvature equal to the exact left nominal
curvature, and its third derivative is zero outside. The outer curvature
jump is absent, exactly as in PR100. No midpoint is substituted for the
analytic continuation coefficient. The resulting H smoothing allowance is
approximately `4.1314847171e27`.

Adding this unchanged allowance to the refined paired comparison gives

```
sup_[0,7/4] ||actual u_rr - smoothed nominal u_rr||_H
  < 1.185e28.
```

The exact rational upper cap is approximately `1.1845751478843721e28`, more
than 524990 times smaller than the old upper cap. This is an improvement of
a certificate for the same candidate. It is not an actual-error reduction
factor or a statement of relative physical accuracy.

The complete initial P1-velocity-gradient error certificate is retained
unchanged, below `9.672e24`. In particular, small new paired velocity bounds
are not substituted for its full P1 interpolation and convolution charges.
A separate transfer calculation must insert the new displacement-curvature
cap into the original coordinate C2 proof, retaining all other candidate,
source, initial-data and absolute-correction ledgers.

## Verification and scope

The source core and retained reports are SHA256 authenticated before use and
on readback. The original core authenticates its complete retained source
closure. Normal and Python-optimized executions must agree byte for byte.
Controls reject a missing leaf, duplicated/overlapping leaf, reversed leaf,
omitted smoothing, maximum in place of the sum of telescoping channels, and
the old unresolved coarse cutoff cell masquerading as a refined cell.

The result is a conditional mathematical error certificate for scaled
coordinate initial data in the stated model. It gives no small covariant
Einstein residual, exact nonlinear spacetime, finite-amplitude nonlinear
family, cosmological observable, or evidence for a higher-dimensional origin
of the Big Bang.
