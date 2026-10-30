# Refined coordinate error bounds and their local limitation

The same `W1_quadratic` reconstruction now has conditional continuous
second-derivative error bounds below `6.643e26`, improving the preceding
`1.187e28` upper certificates. Its error in a precisely defined global
coordinate second-derivative seminorm is below `2.045e-7` of the target's
seminorm. Nevertheless, at one initial tail point its local radial second-derivative
relative error in the same angular H norm exceeds **99.9%**.

These statements coexist because a large outer-endpoint derivative can
dominate the global denominator. The calculation makes that limitation
quantitative. The candidate, actual error, original target, first-order
fields, quadratic sources and shell law are unchanged. This checkpoint
improves and interprets bounds; it does not establish physical accuracy.

## Complete-domain refinement

Two independent interval refinements reduce the initial-error majorant:

| Quantity | Previous strict upper | New strict upper |
| --- | ---: | ---: |
| Nominal third-derivative smoothing allowance | 4.132e27 | 4.128e26 |
| Paired actual-minus-nominal initial second-derivative comparison | 7.715e27 | 2.367e26 |
| Complete initial radial second-derivative error | 1.185e28 | 6.495e26 |
| Each evolved second-coordinate error | 1.187e28 | 6.643e26 |

The smoothing calculation replaces 35 selected ordinary-packet cells with
256 exact leaves and retains 154 old cells, giving a complete 410-cell
nominal cover. It uses explicit nominal order-six jets to bound the nominal
third displacement derivative. It requires no new actual sixth-derivative
premise and differentiates no unknown Taylor remainder. All declared
reduction thresholds pass.

The paired comparison replaces the preceding ten cutoff leaves with 63
leaves and retains the other 47 cells, giving a complete 110-cell cover.
Each child component cap may use only its freshly computed enclosure and
its authenticated original seed cap. The unchanged ordinary cells now
control the paired global bound. The same positive full23 angular Gram
weights, second-derivative jumps, parity, outer quadratic continuation and smoothing
widths are retained. The velocity-error certificate is unchanged.

An independent implementation recomputes all 319 new derivative leaves,
both complete covers and the combined initial bound. The evolved C2
estimate replaces only that initial second-derivative error cap in the retained
comparison theorem. All absolute candidate, C1, original source-derivative,
defect, boundary and coefficient ledgers remain fixed. The C1 error upper
bound remains `3.541e25`.

See the [nominal refinement](theory/NOMINAL_M3_REFINEMENT_NOTE.md),
[paired refinement](theory/PAIRED_CUTOFF_REFINEMENT_NOTE.md),
[combined report](reference/COMBINED.json), [C2 report](reference/C2.json),
and [independent initial audit](review/INDEPENDENT_COMBINED_INITIAL_REVIEW.md).

## What the relative statement means

Let `U` be the retained formal second-coefficient target, `W=W1_quadratic`,
`x=t/ell`, `ell=1e-6`, and

```
D+={0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x},
|V|_(2,H)=max_{ab in {xx,xr,rr}} sup_D+ ||V_ab||_H.
```

Here H is the registered positive probability-normalized S3 L2 norm of
the 23 metric/scalar coordinate amplitudes. The actual initial scalar-mean
radial second derivative at `r=7/4`, whose Gram weight is one, is strictly below
`-3.249e33`. The independent actual ODE endpoint calculation retains the
same nonzero scale. The target's retained C2 initial trace and inclusion
of `(0,7/4)` in D+ imply `sup_D+ ||U_rr||_H>3.249e33`. Consequently

```
sup_D+ ||U_rr-W_rr||_H / sup_D+ ||U_rr||_H <2.045e-7,
|U-W|_(2,H) / |U|_(2,H)                 <2.045e-7.
```

The initial radial second-derivative global ratio is also below `2e-7`. These are
ratios of global suprema in the specified coordinates and angular norm.
The seminorm excludes values and first derivatives; it is not the full C2
norm. The rr denominator does not establish separate relative xx or xr
bounds, and none of these ratios is a physical comparison scale.

At the initial point `(0,65/64)`, the previously established signed actual
and error intervals instead give

```
||U0_rr-u_W,rr||_H / ||U0_rr||_H >0.999,
u_W(r)=W(0,r).
```

This is a local radial second-derivative ratio with a positive denominator. The
proof uses all 23 actual component intervals, the same
Gram matrix and the certified local error lower bound. Both this local
error and target are small compared with the outer-endpoint scale. The
existing clipped-tail approximation therefore remains locally inaccurate
even though its normalized global error is small.

See the [normalization derivation](theory/COMBINED_INITIAL_C2_AND_NORMALIZATION.md),
[exact normalization report](reference/NORMALIZATION.json), and
[independent review](review/INDEPENDENT_COORDINATE_NORMALIZATION_REVIEW.md).

## Conditional transfer to the natural equations

A separate result additionally assumes both the actual first harmonic
identity on a neighbourhood and the formal target identity
`DC[U]+QC=0` as a C1 covector identity on the domain. These assumptions
are stated explicitly; this checkpoint does not independently prove the
complete physical initial, shell or subsidiary constraints.

For the trace-reversed Einstein–scalar equations

```
E_ab=Ric_ab-dphi_a*dphi_b-(2/3)V(phi)g_ab,
S=Box_g phi-Vprime(phi),
```

the unchanged-source identity cancels the common quadratic gauge term.
The original reduced defect and the new error jets then bound the natural
second-parameter residuals of W in a positive background-orthonormal-frame
S3 probability L2 comparison norm:

```
||E2[W]|| <1.476e35,
||S2[W]|| <2.118e30.
```

The parameter convention is `g=g0+lambda*q+(lambda^2/2)*W`; the displayed
quantities are second parameter derivatives, not finite-amplitude nonlinear
residuals. They are large absolute bounds and show no physical smallness.
The coordinate target denominator used above does not normalize these
equations. The positive frame norm is not a Lorentz-invariant observable.
The older source-averaged candidate and its extra source charges are not
used here.

See the [natural-equation derivation](theory/NATURAL_W1Q_ERROR_TRANSFER.md),
[exact report](reference/NATURAL.json), and
[independent identity review](review/INDEPENDENT_NATURAL_RESIDUAL_IDENTITY_REVIEW.md).

## Reproduction and scope

From a complete checkout, with Python's standard library:

```
python -I -S -B research/HDBLAST_REFINED_ERROR_TRANSFER_20261030/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-refined-error-fresh
```

Use a fresh output directory outside the capsule. The runner authenticates
inputs before scientific imports and on readback, reproduces 13 reports
including the retained independent endpoint replay in normal and optimized
Python, and compares their bytes with the references. It rejects 43
corrupted formula, scope, source-role or premise executions in each mode.
The natural-transfer flags record the theorem's assumptions; accepting
those flags is not a check that the physical constraints hold.

Independent internal implementations check the fresh derivative leaves,
covering and convolution arithmetic, C2 transfer, coordinate normalization
and natural-residual arithmetic. They retain authenticated upstream model
premises. This is neither external specialist certification nor a formal
proof-assistant result. Established interval methods and the limited
literature comparison remain documented in the preceding checkpoint.

The next useful tasks include improving the numerical reconstruction near
the clipped tail and independently completing the physical constraint
analysis. A changed candidate needs new initial, transfer and point-value
ledgers. No invariant observable, finite-amplitude nonlinear family,
cosmological prediction or evidence of a higher-dimensional Big Bang origin
is established. Frozen manuscripts, earlier candidates and checkpoints are
preserved. No Zenodo version or journal submission accompanies this work.
