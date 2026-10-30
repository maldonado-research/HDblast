# A complete reduced-bulk forcing defect at the retained shell sample

This calculation concerns the previous integral-defined Robin reconstruction
`W0` at `(x,r)=(9/8,0)`. The conventions are scaled time `x=t/ell`, inward
scaled radial coordinate `r`, and `ell=1e-6`. Its principal bulk forcing is still frozen at

```
f0=L0 B_s=A0 (B_s)_x+B0mat (B_s)_r+D0 B_s.
```

The actual reduced evolution uses

```
La V=A(x,r)V_x+B(x,r)V_r+D(x,r)V,
U_xx-U_rr=La U+F.
```

The forcing defect for an actual-bulk update is therefore

```
r0=La W0+F-L0 B_s.
```

This sign convention is the actual right-hand side minus the old frozen
principal forcing. A residual defined as principal-minus-right-hand-side
would have the opposite sign. The result below is a registered coordinate
modal forcing defect. It is not a complete pointwise Einstein or constraint
residual, a physical force, or the first actual-bulk correction itself.

## Point-assembled coefficients and prior field enclosures

The program pins the existing `compute/full23_response.py` driver. Its
`Response` constructor authenticates its local initial/kernel dependencies;
the initial `Context` authenticates all seven complete scientific roots before
scientific imports. The actual background Taylor enclosure is then evaluated
at `r=0`, and the driver's full coefficient assembly is evaluated at `x=9/8`.
The retained rational inequalities for `tanh(ell*x)` and `sech(ell*x)` enclose
their actual values. These are point-specific coefficient boxes, not the much
wider coefficient boxes over the full retained cone.

`B_s` and its two first derivatives come from the pinned PR94 shell table.
`W0` and its first derivatives come from the pinned PR95 reconstruction table.
The latter enclose the integral-defined reconstruction, including the
unresolved actual shell source convolution. No old finite partial sum is
substituted for `W0`. The frozen matrices are the same pinned matrices as in
the preceding construction.

The two arithmetic forms

```
La W0-L0 B_s,
(La-L0)B_s+La(W0-B_s)
```

are evaluated independently as rational intervals and intersected. The second
form avoids subtracting enormous free-wave source terms. Both enclose the same
quantity for the same actual coefficient functions. No commutation of the
coupled matrices with parity or independent realization of interval endpoints
is assumed.

## Genuine affine bulk-source bounds

The number previously denoted `SOURCE` bounds the response `X(KF+LG)` and
cannot be used as a pointwise bound on `F`. This checker instead authenticates
`HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_SOURCES_NORMAL.json` and
the complete geometric source map to which that report refers. The source
report's domain

```
0<=x<=3, r>=0, x+r<=8
```

contains the sample, including the physical-side shell boundary. The report
states that it encloses actual evolved common-gauge sources. Its normalization
is exactly

```
metric_forcing=-2 ell^2 rho^2 QEred,
scalar_forcing=ell^2 rho^2 Qphi.
```

Those factors are already present in the saved forcing-component caps and
are not applied a second time. The complete `QEred` source projection remains
an inherited premise, identified by
`PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION`; the operator helper's
older warning about a caller supplying matter terms only does not provide
this premise by itself.

For the seven physical blocks `yy,y0,00,ya,0a,ab,scalar`, the checker compares
all 93 exact entries of the source Gram matrices, including block
multiplicities, against the current diagonal full23 Gram weights under
`pack_fields`. They agree exactly. It also reconstructs the recorded total
source norm square from the component caps. This establishes that the
source cap and present registered H norm have the same normalization. The
uniform bound is

```
||F||_H < 17387513156573.
```

Sharper individual source-component caps are packed directly into all 23
current slots. Each actual `F_i` is enclosed by the corresponding symmetric
interval. This is a uniform magnitude enclosure covering the point, not a
new evaluation of `F` at the point or a determination of its sign.

## Result

Seventeen of the 23 complete forcing-defect intervals exclude zero. At the
primary scalar-mean slot, exact rational endpoints imply the strict readable
bound

```
-9.417e22 < (r0)_w:0 < -9.416e22.
```

The exact interval is approximately
`[-9.416751844940e22,-9.416550797141e22]`. Other examples are positive `K:0`,
negative `N:0` and positive `J:0`. The `u:2`, `u:4`, `H:2`, `H:4` and both tensor slots
remain undetermined. In particular, the tensor slots retain their affine
source enclosures; they are not set to zero merely because the old free-wave
candidate had zero tensor coefficients.

The report includes every exact component interval, both lower-order
arithmetic forms, the individual affine source caps, and lower/upper squared
H-norm bounds. Its whole-defect norm is strictly larger than `9e22` under
these enclosures. This supplies a concrete nonzero forcing defect in the
old auxiliary reconstruction. The sign at one point does not determine the
sign of its spacetime-integrated correction or a physical observable.

## Reproduction and limits

Run `check_bulk_point_residual.py --repo-root PATH`; optional
`--robin-capsule` and `--contact-capsule` select separately located pinned
copies. Normal and optimized Python produce identical JSON. All seven
scientific roots, local driver dependencies and the pinned driver/source/map
files are reauthenticated after evaluation.

This calculation uses only first-derivative field intervals and the known
frozen forcing identity of the earlier reconstruction. It does not numerically
differentiate the field twice or establish a quantitative C2 error. Any
classical principal-wave interpretation retains the earlier qualitative
regularity and corner assumptions. The first-field, background, complete
source-map and operator certificates remain inherited conditional premises.
No full coupled second-order solution, nonlinear family, finite-amplitude
remainder, invariant detection or higher-dimensional Big Bang origin follows.

Pinned new scientific driver/source inputs:

* `full23_response.py`: `efa6a6c2a9c7db596fe024246684f760b807e727fd2ddac8e8126182ba941a73`.
* `ACTUAL_SOURCES_NORMAL.json`: `45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f`.
* `GEOMETRIC_COEFFICIENT_MAPS.json`: `34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615`.
