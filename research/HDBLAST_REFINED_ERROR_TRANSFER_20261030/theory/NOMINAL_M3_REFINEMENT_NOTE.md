# Refining the explicit nominal third-derivative allowance

This calculation keeps the existing `W1_quadratic` reconstruction, initial
displacement and velocity, smoothing widths, outer quadratic continuation,
actual target, first-order fields, original sources and physical shell law
unchanged. It replaces some intervals used to enclose a derivative of an
explicit nominal function. A smaller upper certificate is not a change in
the numerical candidate or its actual error.

## Derivative being bounded

The retained nominal initial displacement is expression-defined from the
registered polynomial background, explicit rational packet and time-profile
functions, and a polynomial inner cutoff. Its formula uses packet jets
through order three. The existing `HigherInitial.nominal` API evaluates
order-six input jets; each displacement output retains three derivatives.
With its normalized physical-radial Taylor convention, the ordinary third
scaled radial derivative is exactly

```
u_nom,rrr = 3! ell^3 c[3],  ell=1e-6.
```

This is a piecewise ordinary derivative, with one-sided values at joins.
The calculation does not differentiate an unknown remainder, import an
actual sixth-order ODE/packet premise or give an actual sixth derivative.
The nominal expressions are explicit, so their sixth input jets are
available algebraically. Existing source closure and nominal-jet identities
remain authenticated premises.

## A finite exact interval refinement

The authenticated original nominal report contains 189 interpolation cells
covering `[0,7/4]`, including zero, cutoff and plain branches. Its stored
third-derivative component vector is exactly the componentwise maximum of
those cells' `second_derivative_component_caps[1]` vectors. The same vector
is recorded in the higher-initial and preceding refined-initial reports.

The declared work list consists of the 35 old cells on which at least one
of the scalar components `w:0,w:2,w:4` exceeds one tenth of its old global
component bound. Their union is `[273/256,77/64]`, entirely within the plain
packet branch. All 154 other old cells and their component caps are retained.

Each selected cell is evaluated on exact rational dyadic intervals using
its forced `plain` branch. A leaf's raw cap for component i is

```
raw_i = mag(6 ell^3 output_i.c[3]).
```

Because the leaf is a subset of its authenticated old parent, the old parent
cap and new raw cap both bound the same derivative everywhere on that leaf.
Their minimum is therefore a valid upper cap. This minimum is not a midpoint
substitution or an assumption that interval extrema occur together.

A leaf is bisected at its exact midpoint if any scalar effective cap exceeds
the declared tenth-old threshold. Maximum depth is four. Every leaf is
retained at that budget even if the threshold remains unmet, and all misses
are explicitly listed in the report. This finite target is a calculation
budget, not a physical accuracy threshold. Coverage and validity do not
require every target to pass.

Every selected parent's leaves are checked for positive widths, matching
endpoints and exact adjacent boundaries. Together with all retained parents
they cover the full nominal interval without gaps or overlaps. New global
component caps are maxima over this complete cover. Parent minima ensure
no component's certificate increases; observed strict improvement is checked
rather than assumed.

## Same smoothing and jump estimate

The positive even kernel has unchanged width `delta=2^-30` and first absolute
moment `35 delta/128`. The nominal displacement is C1, so its distributional
second derivative has no point-mass term. Its piecewise second derivative
can have a curvature jump. The first inner join at `65/64` has exactly zero
curvature jump; the second inner join's component jump caps are retained.
The outer quadratic continuation matches the exact nominal curvature and
therefore has no curvature jump, while its third derivative is zero outside.
The kernel footprint can encounter at most one inner join.

The same componentwise convolution allowance is

```
smooth_i = (35 delta/128) M3_i + max_inner_jump_i/2.
```

Substitution of the new piecewise nominal M3 caps preserves the original
proof. The complete positive 23-amplitude H norm packs this vector after the
componentwise sum. The old allowance is reconstructed exactly from the pinned
preceding report before computing the new allowance. This report alone does
not replace the paired actual-minus-nominal comparison or the complete initial
error ledger; a separate checker must add those bounds and propagate the
result through the retained C2 comparison theorem.

## Verification and limits

Source and report hashes are checked before scientific imports and on final
readback. Scientific acceptance decisions use explicit runtime checks and
exact rational arithmetic. Normal and optimized Python outputs must agree
byte for byte. Corrupted modes reject missing cover cells, a wrong nominal
branch, omitted factorial or radial scale, dropped inner-jump charges and
maximum in place of the valid parent minimum.

The result concerns coordinate numerical approximation estimates under the
retained model, source and regularity premises. It does not establish a
physical force, pointwise relative accuracy, complete constraints, a small
covariant Einstein residual, a finite-amplitude nonlinear family or a
higher-dimensional origin of the Big Bang.

## Result of the declared finite calculation

All 35 selected parents meet their scalar tenth-old targets within the
maximum depth of four. The 256 accepted leaves and 154 retained old cells
form a complete 410-cell cover. The producer evaluated 477 intervals: 221
split parents and 256 leaves. No threshold-miss leaf is discarded or hidden;
there are zero misses in this execution.

The refined full23 piecewise nominal M3 H cap is approximately
`1.6209832836698667e36`, below `1.621e36`. The complete curvature-smoothing
H allowance is approximately `4.127977880317139e26`, below `4.128e26`,
compared with the former `4.131484717143961e27`. The exact upper-certificate
ratio is approximately `10.0084952898`. These values describe upper bounds
on the same nominal derivative and same smoothing discrepancy; they are
not ratios of actual errors. The complete exact rational quantities and
all component, parent and leaf data are retained in the report.

The dyadic proof boundaries do not introduce new joins of the nominal
function. They divide the enclosure calculation only. Its old branch joins
and their curvature traces remain exactly the same; no additional jump
allowance is required merely because adjacent proof cells have different
upper caps.
