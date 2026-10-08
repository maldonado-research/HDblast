# Additive acceptance of dyadic rounding for complete source errors

8 October 2026 UTC. This note supplements, and does not replace or alter, `SOURCE_AND_ENDPOINT_ENGINE_PROOF.md` and the original `MATHEMATICAL_REVIEW_RECEIPT.json` (SHA256 `1763f27b8c878ac835c4e7d47226f3778a771d9cd9fe36d934328d4394d37980`). The original proof, 109 exact checks, 11 negative controls, receipts and reviewed source snapshots remain preserved. No retained array, physical source construction or physical endpoint evaluation was used for this update.

## Exact change and its purpose

The updated `implementation/source/later_source.py` keeps the same exact chosen degree24 source coefficients, the same analytic Taylor tail, and the same exact coefficient discrepancy. It now rounds each cell's **complete** uniform source error upward to the fixed dyadic grid of width `delta=2^-512` before supplying that error to `EndpointFlow`. This prevents accumulation of many distinct odd denominators in the serialized global source-error budget.

The source certificate retains each cell's exact analytic-tail term, coefficient-discrepancy term, and the exact nonnegative increment introduced by this rounding. Its exported complete error is their sum. The physical endpoint engine remains byte-identical to the originally reviewed engine, SHA256 `b1e5c38b79fe1cc7782cf39cb49daa58f311305bbea493a3a720b062f8859c79`. Its existing independent upward rounding of the kernel tail is a separate error contribution and remains counted once.

## Ceiling proof and containment

Let the prior complete cell error be the nonnegative exact rational `q=p/r`, where `p>=0` and `r>0` are integers. Set `S=2^512`. The updated helper implements

```
q_plus = floor((p*S+r-1)/r)/S = ceil(p*S/r)/S.
```

The ceiling identity follows from integer Euclidean division, including the case where `p*S` is divisible by `r`. For every such input,

```
q <= q_plus < q+2^-512.
```

If the old source proof establishes `|g(s)-p_cell(s)|<=q` throughout the cell, the same polynomial therefore satisfies `|g(s)-p_cell(s)|<=q_plus`. There is no new analytic premise and no change to the prescribed real source. The rounding increment is `eta=q_plus-q`, with `0<=eta<2^-512`. Exactly aligned errors, including zero, are unchanged.

The implementation returns `q_plus` as the row's sole `uniform_error`. The engine transports that value directly. `analytic_tail`, `coefficient_error` and `source_error_rounding` are a decomposition of `uniform_error`, not additional independent radii to add again. Explicitly,

```
uniform_error = analytic_tail + coefficient_error + source_error_rounding.
```

Because all returned complete errors are dyadic with denominator dividing `2^512`, summing them with the exact dyadic panel lengths and first-moment factors does not form a least common multiple of their old odd denominators. The retained per-cell explanatory rational literals remain exact; this change specifically controls the aggregate denominator problem.

## Additional endpoint uncertainty, with exact duration factors

For endpoint `t` let `D=t-a` be `1/2` or `1`. The source-cell intervals partition `[a,t]` exactly and have center `c_j` and width `h=1/64`. The engine's source error bounds are

```
B_W = sum_j h q_j,
B_U = sum_j h(t-c_j)q_j.
```

They are the exact integrals of the piecewise constant envelope `q(s)` and of `(t-s)q(s)`. The increments caused by dyadic rounding consequently satisfy

```
0 <= B_W_plus-B_W = integral_a^t eta(s) ds < D*2^-512,
0 <= B_U_plus-B_U = integral_a^t (t-s)eta(s) ds < D^2*2^-513.
```

These are source-response complex-modulus uncertainty increments. When the implementation adds a modulus allowance independently to both real and imaginary interval endpoints, the resulting Cartesian L1 radius can increase by at most twice the corresponding modulus allowance before final grid export. The complete exported L1-radius gate still includes all arithmetic, source, kernel-tail and export contributions.

For `t=-4`, the additional `W` allowance is less than `2^-513` and the `U` allowance is less than `2^-515`. For `t=-7/2`, they are less than `2^-512` and `2^-513`, respectively. No incoming-state uncertainty or earlier BD cap is newly introduced. No force, contact, mode or phase value changes.

## Inspection and independent regression

Static inspection of the updated source confirms that:

1. `outward_error` rejects non-`Fraction` or negative inputs and uses the exact ceiling formula above with the fixed `BITS=512`.
2. `source_cell` first computes the original exact `raw_error=TAIL+discrepancy`, then rounds that complete error once.
3. The emitted `source_error_rounding` is exactly `error-raw_error` and the source row supplies only the rounded `error` to the unchanged engine.
4. Manufactured zero-error source rows preserve zero error and explicitly report zero rounding.
5. Source coefficient construction, source coefficient hashes, declared geometry, derivative inventory and analytic tail are unchanged.

The separate endpoint reviewer exercised 64 fabricated odd-denominator errors without calling the physical source constructor. Their raw summed denominators had 53,865 and 53,866 bits, while rounded aggregate denominators had 513 or 514 bits; the serialized budgets were 898 and 899 bytes. The reviewer checked upward containment, the strict one-cell grid increment, both duration-weighted transport bounds, and target-rectangle containment of an attained zero-momentum source-uncertainty example. The additive receipt binds the inspected control script, PASS receipt and preserved original failure receipt. These checks are not rerun here because the exact mathematical change and its distinctive failure case have already been independently covered.

## Acceptance and unchanged scope

This source-error update is mathematically accepted for the precise updated source bytes and unchanged engine bytes identified in `DYADIC_SOURCE_ERROR_UPDATE_RECEIPT.json`. The original acceptance remains applicable to its original snapshot; this additive acceptance supplies the changed-source proof for the new candidate. Any further source or engine change requires a new inspection.

This acceptance is a static mathematical conclusion, not a public byte GO, execution authorization, or retained-state numerical result. Full frozen-closure review and bounded manufactured pipeline validation remain separate prerequisites. The result still concerns the two retained later endpoints relative to exact evolution from the original saved incoming state. It does not reconstruct the historical unsaved path, certify complete stored stress/contact arithmetic, close continuum momentum quadrature or UV errors, change the failed metric calibration, or establish a higher-dimensional origin.
