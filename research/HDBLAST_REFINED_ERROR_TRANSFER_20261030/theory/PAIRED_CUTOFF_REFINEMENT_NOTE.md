# Further paired initial-curvature refinement for the same reconstruction

The actual target, registered nominal displacement, P1 velocity, both
convolution widths, outer quadratic continuation, original source and
`W1_quadratic` finite reconstruction remain fixed. This calculation changes
only interval enclosures for the paired initial-data comparison. The same
authenticated order-five algebra evaluates the three telescoping channels
of actual minus nominal data. Its actual-background derivative remainders
are separately proved enclosures; an unknown remainder is not differentiated.

## Authenticated parent cover and bounded refinement

The preceding refinement report, SHA256
`915acf8abf3c3bb2c843eac509315ebdbd568d7c8b10ff2682c900d672e7d4d7`,
supplies ten exact cutoff leaves covering `[65/64,66/64]`. Its remaining
47 cells retain the original complete active interval cover on `[1,7/4]`.
The new calculation retains those 47 rows exactly and subdivides only the
ten authenticated cutoff leaves. Every midpoint and endpoint is rational.

For a new leaf `I` contained in seed leaf `S`, let `b_raw[ch,j,i]` be the
fresh interval cap for the same ordinary scaled derivative. Let `b_S` be
the authenticated component cap from the preceding seed. Both bound the
same function throughout `I`; hence

```
b_eff[ch,j,i] = min(b_raw[ch,j,i],b_S[ch,j,i])
```

is a valid component upper bound on `I`. Every leaf records the raw,
inherited and effective tensors, and its inherited tensor must equal the
authenticated seed tensor exactly. No unreported internal-parent minimum
enters this argument. The retained ordinary derivative conventions are
`u_rr=2 ell^2 c[2]` and `v_r=ell c[1]`.

For each effective component vector, use the same positive full23 Gram
matrix to obtain an outward rational H cap. For each of the three
telescoping channels, take the maximum over the entire initial interval
cover and then add channels by the triangle inequality:

```
paired_u_rr_H <= sum_ch max_leaf H_cap(ch,leaf).
```

The channel maxima may occur on different leaves; summing them remains a
valid upper bound. Neither a maximum in place of the channel sum nor an
incomplete cell cover is accepted.

The declared cutoff stopping threshold is `1e26` for every displacement
channel H cap, with original-cutoff dyadic depth at most 16 and at most
128 final cutoff leaves. A leaf that reaches a budget before meeting its
threshold is still retained with its valid bound and explicitly reported
as a threshold miss. The threshold is a computation choice, not a physical
accuracy criterion. The complete forest, evaluated-node metadata, split
ranges, accepted leaves and stop reasons are saved. Source authentication
is repeated after computation.

## Result of the bounded calculation

The calculation produces 63 cutoff leaves after 116 evaluations, including
53 splits of the ten preceding seeds. Every final cutoff leaf satisfies
the declared `1e26` displacement-channel threshold; neither budget is
exhausted. The complete active cover contains 110 cells: those 63 leaves
and the 47 retained cells, including the original flat-tail interval.

The paired displacement-curvature H upper cap is approximately
`2.3668485863212598e26`, strictly below `2.367e26`. The preceding paired
cap is approximately `7.7142667617e27`, so the new upper certificate is
over 32 times tighter. The remaining ordinary cells now set the global
paired maximum. These ratios compare upper certificates, not actual errors.

Adding even the unchanged preceding convolution charge gives a complete
fallback initial curvature-error cap approximately `4.3681695758e27`,
strictly below `4.369e27`. A separate nominal-M3 calculation can tighten
that allowance further. No such external combination is implicit in this
paired report.

## Convolution and transfer are separate

This report bounds the uniform paired actual-minus-nominal displacement
curvature comparison before its nominal convolution charge. A complete
smoothed initial-error certificate must add a valid smoothing allowance
for the same nominal function, kernel width and curvature jumps. The
report retains the previous allowance as a diagnostic fallback, and
identifies a separate nominal-M3 refinement as an external combination
step. It does not silently drop smoothing or import a changed cutoff.

The full initial P1 velocity-gradient error cap remains exactly the
preceding certified cap, including its secant, rounded-node and
convolution contributions. New paired velocity bounds alone cannot
replace that full error ledger. Likewise, a refined *error* bound cannot
replace an absolute candidate derivative cap in the differentiated bulk
forcing calculation. Evolved C2 transfer and coordinate norm normalization
require their own authenticated arguments.

## Verification and scope

The standalone checker uses Python's standard library and exact rational
arithmetic with explicit runtime checks. Its core and original nominal,
initial, repaired-extension and preceding-refinement reports are pinned.
The core authenticates its complete scientific source closure before
imports and on readback. Normal and optimized execution must give
byte-identical JSON. Structural controls reject missing, overlapping and
reversed leaves, replacing the telescoping sum by a maximum, dropping the
convolution charge and treating an unresolved parent as a successful leaf.
Four injectable corruptions reject an omitted curvature factorial, omitted
radial scaling, changed smoothing width and incomplete authenticated parent
cover before a corrupted certificate can be emitted.

This is a conditional coordinate approximation certificate for the same
formal target and auxiliary finite reconstruction. Refining the certificate
does not change the candidate or its actual error. No pointwise relative
accuracy, covariant or physical accuracy, complete constraints, nonlinear
spacetime family or cosmological-origin conclusion is established here.
