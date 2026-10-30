# Transferring a refined initial-data certificate to the unchanged response

This calculation changes an error enclosure for the existing initial data. It
changes neither the expression-defined nominal data, smoothed profiles, free
wave, frozen bulk forcing nor the Robin reconstruction `W` published in
`HDBLAST_ROBIN_RECONSTRUCTION_20261030`. The retained formal target `U` is also
unchanged. A smaller upper certificate therefore does not measure a change in
the unknown true approximation error.

## Independent continuous-coverage and norm audit

The checker authenticates the original initial report and the refinement
report. It verifies all 192 original cells covering `[1,7/4]`; the inherited
initial zero collar covers `[0,1]`. An untouched cell retains its original
three-channel, three-kind, 23-component caps. A refined parent is replaced by
exact dyadic leaves. Every binary path must reproduce its rational endpoints;
sorted leaves must cover the parent without gaps or overlaps, and no leaf may
cross either analytic cutoff join.

For each leaf, the raw interval evaluation and the inherited parent cap both
bound the same fixed difference. The effective magnitude bound is their
componentwise minimum. The checker verifies this clipping identity exactly.
It then takes maxima over every retained leaf and untouched cell before
forming the positive Gram norms. For channel `c` and kind `k`,

```
b[c,k,i] = max_all_cells b_cell[c,k,i],
D[c,k]^2 >= sum_i H[i] b[c,k,i]^2,
delta[k] = sum_c D[c,k].
```

This is a continuous bound on the full initial footprint; it is not an
endpoint sample or an average. The separate refinement producer reruns the
same authenticated paired interval evaluator on every new leaf. This
JSON-only checker independently audits coverage, clipping, maxima and Gram
aggregation without importing that scientific evaluator. Hashes authenticate
the recorded bounds; their analytic meaning and earlier field certificates
remain inherited premises.

## Same-candidate direct C1 ledger

Let the refined actual-minus-nominal bounds be `delta_u`, `delta_ur`,
`delta_v`. All smoothing and interpolation contributions retain their prior
values:

```
a_u = M2 (2^-30)^2/18,
b_u = 35 M2 (2^-30)/128,
a_v = (35/(256*2048)) max_knots ||velocity_slope_jump||_H.
Delta_u  = delta_u + a_u,
Delta_ur = delta_ur + b_u,
Delta_v  = delta_v + eps_v + a_v,
e_initial = max((Delta_u+(3/2)Delta_v)/3, Delta_ur+Delta_v).
```

The checker independently reconstructs the matrix Gram bounds, coefficient
uncertainties, smoothing terms and direct bulk residual from the pinned raw
inputs. It verifies twelve whole-cone quantities against the earlier direct
report: the `B0` size, global Robin correction, actual/frozen contractions,
whole-cone first feedback, coefficient and feedback residuals, and affine
bulk response. Every one remains exactly unchanged.

The direct proof from the previous checkpoint therefore applies with only the
new initial term:

```
r_new = e_initial,new + bulk_residual_response_X_cap,
X(U-W) <= r_new/(1-q_actual)
       <= r_new/(1-3409/8000000).
```

The response norm remains
`X(e)=max(sup||e||_H/3,sup||e_x||_H,sup||e_r||_H)` on the full retained cone.
The old bound on `X(KF+LG)` is still treated in response units, with
`X(KF)<=X(KF+LG)+X(LG)`. It is not replaced by an unproved pointwise bulk
source bound. The actual boundary source cancels in the error equation
because `U` and `W` use the same fixed actual `M` and `G`.

The exact report gives both contraction versions, a strictly outward-rounded
readable bound on a `10^22` grid, and ratios against the preceding direct
certificate. Those ratios concern certified upper bounds rather than true
error measurements.

For the finalized ten-leaf refinement, the conservative inherited-contraction
calculation yields

```
X(U-W) <= 2.021348047763e26 < 2.0214e26.
```

The preceding PR95 direct upper certificate was approximately
`1.734335783353e27`; its ratio to the new readable upper bound exceeds `8.57`.
The initial term now contributes approximately 17.5% of the direct numerator,
while the unchanged bulk residual response contributes approximately 82.5%.
This identifies a change in the conservative budget's dominant term, not a
measurement of the true distribution of error.

## Updated shell intervals

The same authenticated 23-mode reconstruction table at `(x,r)=(9/8,0)` is
reused. All candidate intervals are checked for exact equality with the prior
table. With new global bound `E`, the zero-wedge shell value error is `E/8`.
For each component, the generic errors are at most `E/sqrt(H_i)` for first
derivatives and `E/(8 sqrt(H_i))` for values. Homogeneous Dirichlet values and
time derivatives remain exactly zero.

For Robin rows, the shared actual law gives

```
(U-W)_r = M (U-W),
|(U-W)_r,i| <= (E/8) sqrt(sum_j max|M_ij|^2/H_j).
```

The new intervals are verified to lie inside the older conditional intervals
in all 23 slots. For zero Robin rows the radial error remains exactly zero;
this preserves the actual shell-source interval already present in `W_r` and
does not set the actual tensor source to zero. All signs remain coordinate
and gauge-dependent formal coefficients, with scaled time `x=t/ell` and
inward radial coordinate `r`.

For the registered scalar mean, the new exact interval table implies the
deliberately outward-rounded strict bounds

```
8.29035e30 < U_w:0 < 8.29044e30,
1.974619e32 < partial_x U_w:0 < 1.974625e32,
2.154062e27 < partial_r U_w:0 < 2.154085e27.
```

The last row uses the shared Robin law. The full table contains all 23
components; it reuses the earlier reconstruction enclosures rather than
performing a new independent point evolution.

## Reproduction and limits

`check_refined_response.py` uses only Python's standard library and
authenticated JSON. Its optional `--refined-report`, `--robin-capsule`,
`--transfer-report` and `--source-report` arguments select separately located
pinned copies. Run normally and with `-O`; exact JSON output must agree byte
for byte.

This supplies tighter conditional initial and C1 response certificates for
the same data and candidate. The auxiliary bulk forcing remains frozen. It
does not complete coupled bulk evolution, provide a numerical C2 error,
pointwise Einstein/constraint residual, exact nonlinear spacetime or physical
evidence for a higher-dimensional Big Bang origin.
