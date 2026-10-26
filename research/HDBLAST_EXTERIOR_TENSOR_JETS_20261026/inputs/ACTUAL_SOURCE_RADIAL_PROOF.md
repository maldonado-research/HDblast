# Actual full TT source radial derivative

`actual_third_jets.py` encloses the complete projected source `F` and its
scaled source-coordinate derivative `Fr = partial_r F` on whole rectangles
in `0 <= x <= 3`, `r >= x`, `x+r <= 8`.  In particular it covers the
rectangles containing the clipped cells of
`0 <= x <= 1/4`, `1+x <= r <= 3/2-x`.

The new input is a derivative certificate for the actual retained smooth
first solution.  Existence of that solution, the predecessor differentiated
integral and reflection identities, its error certificates, and the formal
all-point geometric projection remain authenticated conditional scientific
premises.  This calculation does not replace the actual solution by the
comparison packet.

## Required derivative order

Every one of the 78 monomials in `metric_forcing/ab/5` is evaluated.  Row 6 is
identical.  Scanning all powers shows the only second physical derivatives
are `Myy`, `Myt`, and `Mtt`, appearing in the three terms

```
-(4/3) A Myy ell^2 rho^2
+(8/3) B Myt ell^2 rho
-(4/3) D Mtt ell^2.
```

The registered frame amplitude `M=2d` contains the gauge field without a
derivative.  All other used metric amplitudes contain at most one gauge
derivative and occur with at most one physical derivative.  The master is
undifferentiated in the used metric amplitudes; the scalar `s`, used only
undifferentiated, contains one master derivative.  After one additional
`r` derivative, the complete source therefore requires gauge derivatives
through total order 3 and master derivatives through total order 2.
Existing authenticated C3 actual error estimates suffice.  No C4 actual
error estimate is inferred from a C3 interval.

`c4_need_audit.py` is a separate exact dependency scan.  Its normal and
optimized reports agree.  The oracle itself also rejects a changed source
dependency signature and rejects any raw variable whose reconstruction
order plus physical derivative order plus one exceeds 3.

## Actual first-solution coefficients

For each mixed scaled derivative of total order at most 3, the evaluator
forms the normalized coefficient of the explicit comparison function and
adds the symmetric actual-minus-comparison error interval divided by the
appropriate factorials.  The errors are the authenticated exterior C3
certificate produced by `exterior_error_caps.py`.  That producer replays
all inherited whole-domain master and gauge error formulas exactly before
restricting to the backward-closed exterior domain.

These coefficients enclose derivatives of one actual smooth solution.
They are not produced by differentiating selected interval values.  All
nonzero master and gauge error caps are retained, as are the actual shell
background uncertainty and differentiated background ODE enclosures.

The degree-3 master/gauge Taylor algebra reconstructs the first metric
amplitudes and scalar.  Amplitudes containing a derivative are certified
only through total order 2; `M=2d` remains certified through total order 3.
Only the source-selected physical derivative expressions are evaluated.
No unused third derivative of a differentiated metric amplitude is exposed
as an actual derivative.

## Physical derivatives and background chain rules

Write `L=1/(ell rho)`.  The physical conversions are

```
j_y  = -L j_r
j_t  = j_x/ell
j_yy = L^2 j_rr + a L j_r
j_yt = -L j_xr/ell
j_tt = j_xx/ell^2.
```

Each is formed as a Taylor-algebra expression, then its value and first
scaled radial coefficient are taken.  Thus the derivative of every `rho`
factor and of the physical conversion itself is included.  Equivalently,
for example,

```
(j_yy)_r = L^2 j_rrr + (2L L_r+a L) j_rr
           +(a_r L+a L_r) j_r.
```

Every background variable actually used by the complete TT map is supplied
as a radial Taylor expression: `rho`, `a`, `ay=-rho^-2-p^2/3`, and
`Vp=V'(eta)`.  The temporal factors `c`, `theta` and the constant `ell` have
zero radial derivative.  The background derivative recurrence is the
authenticated actual ODE recurrence through order 4; no point background
is substituted.

## Full algebraic source differentiation

Each source input is a dual number `(value, radial derivative)`.  Products,
integer powers, sums and inverses obey the exact product and chain rules.
The authenticated full collected source map is evaluated on this ring, so
the second component is `partial_r F`.  Both its value and derivative are
outward intervals on the same grid of denominator `2^112`; all decisions
are exact rational decisions.  There are no finite differences or
floating-point acceptance decisions.

Before importing frozen project code, the oracle authenticates the full
frozen TT tree against its pinned manifest and source hash.  The frozen TT
context then authenticates the seven complete predecessor closures before
its first project import.  All eight closures and the new oracle source
are reauthenticated at completion.

## Bounded verification

The declared verification budget is three points and two whole cells,
with no mesh search.  At every probe:

1. Each source input value agrees exactly in both endpoints with the old
   authenticated actual first-jet oracle.
2. The complete source value agrees exactly with the old full TT source.
3. The entire source value and radial derivative have identical endpoints
   when evaluated with Fraction interval arithmetic and the equivalent
   integer dyadic interval implementation.

`actual_source_radial_normal.json` and
`actual_source_radial_optimized.json` agree byte for byte.  The apex radial
source interval is approximately `[-1.29122e-6,-8.46760e-7]`; these decimal
numbers are display values only.  Rational endpoints in the reports are
the certificate.  The oracle is the new scientific input; the bounded
five-probe verification does not replace its whole-rectangle algebraic
proof.
