# Scope review: actual causal correction quadrature

This review concerns the conditional registered reduced-system correction
`C=W1-W0`, evaluated at `(x,r)=(9/8,0)`. It does not treat C as the target
error `U-W0`, a nonlinear perturbation family, or a physical force.

## Transfer to the same conditional target

`check_point_transfer.py` adds enclosing correction intervals to the previously
enclosed W0 first jets. It retains the prior whole-domain error certificate
`X(U-W1)<3.541e25`; the numerical quadrature supplies point information and
does not improve that theorem. For each positive Gram weight H_i, the zero
wedge and the selected shell time give component errors at most
`(1/8) E/sqrt(H_i)` for values and `E/sqrt(H_i)` for time derivatives.

The same fixed actual M and G are used in U and W1. Consequently the Robin
normal error is `M PN(U-W1)` and its component cap is `(1/8) E` times the
weighted dual norm of the corresponding M row. Intersecting the W1 radial
enclosure with the independently valid affine law `M PN W1+G` is valid.
Dirichlet values/time derivatives remain exactly zero, while their radial
errors retain the generic derivative cap. The final intersection with the
previous best U boxes is valid because both enclosures target the same U.

A sign of C is a sign of this finite-iteration auxiliary correction. The
transfer script correctly refrains from identifying it with a sign of target
error, an invariant observable, or an actual gravitational force.

## Source normalization and fixed choices

The interior evaluator pins the complete actual geometric-source map and its
continuous source report. Its normalization is already the operator forcing:
`metric=-2 ell^2 rho^2 QEred`, `scalar=ell^2 rho^2 Qphi`. All 93 block Gram
entries, including the mixed-component multiplicities, agree with the
registered 23-component comparison norm after packing. No extra coordinate
scale is applied. The genuine pointwise uniform H cap is 17387513156573;
the response norm `X(KF+LG)` is not used as a pointwise F bound.

F is widened by its genuine component cap multiplied by each exact positive
moment mass, and by the full physical edge measure for the derivative. The
boundary lift retains the fixed actual G through its component caps plus
the higher-order G-history tail. These outward independent interval choices
contain the same fixed actual functions and M throughout; they do not change
the underlying candidate or freeze sources to midpoint/zero values.

## Geometry, signs, and tails

On the physical active triangle, `alpha=s+y`, `theta=(y-s)/(s+y)` has
Jacobian alpha/2. Complete source moments have positive exact masses
`A_j=d^(j+1)/j! [a/(j+1)-d/(j+2)]`, with `a=9/8`, `d=1/8`. Source coupling is
formed before output reflection. For the same time-constant actual M,

```
Z=C|shell=sum_j (-M)^j integral (a-s-y)^j/j! r0(s,y) ds dy,
C_x=PN p-M Z,      C_r=PD p+M Z,
p=integral_0^a r0(s,a-s) ds.
```

The code implements this sign convention and preserves the Dirichlet normal
edge contribution. The edge parameter measure is a, not d. The first omitted
moment tail is multiplied by the geometric factor
`1/(1-m d/(N+2))`, using the complete local residual H cap. The interior
reconstruction's homogeneous and affine-source Robin tails also bound every
omitted term. Zero M rows correctly have no homogeneous Robin normal
correction; their affine source normal value remains inherited from G.

## Inherited support and regularity

Removing alpha<=1 and returning a zero boundary history through time 1 use
the inherited actual zero-wedge/source-support premise. Small uniform source
caps alone would not justify those branches. The module documents this
dependency. Numerical first jets and integral identities do not certify C2
norms or the full Einstein/subsidiary constraint residual. The target transfer
remains conditional on the inherited reduced-system solution, coefficient,
source, initial-data and boundary hypotheses.

The current implementation uses exact rational or outward dyadic interval
arithmetic, whole-panel ranges and positive masses. No quadrature truncation
claim is inferred from sampled source values. Review of the current code
found no substantive source-unit, target-transfer, sign or scope issue.

## Fine numerical result audit

The completed fine quadrature report has SHA-256
`5bd5ce2d40fdf0805d06351160fe7d2e323324983b69aeee5732c1e28013c6f8`.
Its normal and optimized bytes agree; both transfer outputs also agree.
Exact-rational readback verifies 16 value, 16 time-derivative and 17 radial
correction intervals excluding zero. The transfer narrows 53 of the 69
previous best conditional target point boxes.

For the scalar mean at the selected point, the following strictly outward
readable inequalities were checked against the exact correction endpoints:

```
1.596e21 < C_w0     < 6.733e21,
8.355e22 < (C_w0)_x < 1.058e23,
4.147e17 < (C_w0)_r < 1.750e18.
```

All three corrections are positive. The inherited complete scalar-mean
residual at this endpoint has a strictly negative upper endpoint. These facts
are consistent: the causal correction integrates the physical residual over
the past source triangle and includes coupled Robin history. A one-point
residual sign does not determine the sign of that accumulated correction.
Neither sign determines a physical force or the sign of U-W1.

The scalar target-box width ratios relative to the prior PR97 boxes are
approximately 1.12524 for value/radial and 1.33031 for the time derivative.
These are local enclosure gains under the same global error theorem. No
scope or numerical-result inconsistency was found in this audit.
