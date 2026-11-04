# Frozen admissible-preparation response certificate

4 November 2026. Private numerical result, conditional on the retained
background, operator, gauge, and local-existence premises.

The registered direction is
`chi=ell*rho_b^3*Q(ell*r)`,
`h=-(rho_b/rho)^3*d_r Q(ell*r)`, at `ell=10^-6`, `T=5/4`.
It has exact zero weighted moment because `rho^2 h=chi_y` and `chi` is compact
and flat. The direction, partition, target, tails, coverage, norm allocation,
and sign tests were frozen before the first new numerical evaluation in
`CAMPAIGN_FREEZE.json`. The additive deterministic finite-lambda rule in
`POSTPROCESS_FREEZE.json` was also registered before that evaluation.

## Certified results

The following decimal intervals are deliberately widened outward from the
exact rational endpoint report; their inclusions are checked in
`DISPLAY_AND_CALL_AUDIT.json`.

| Quantity | Outward enclosure |
|---|---:|
| Proper scalar-mean response slope `L_M[h]` | `[1.5393e9,1.7068e9]` |
| Proper volume-second response slope `L_v[h]` | `[-0.00021094,-0.00015795]` |
| Proper expansion-second response slope `L_Theta[h]` | `[-22.9410,-17.4038]` |

All three registered nonzero tests pass. The two physical slopes are strictly
negative. No second direction, rank test, profile tuning, or repeat campaign
with a changed bound was needed.

The registered algebraic rule selects `lambda=-10^21`. For this one fixed,
finite preparation, the exact affine coefficient identities give

| Quantity for `f_lambda=f+lambda*h` | Outward enclosure |
|---|---:|
| Proper volume-second coefficient `v2_lambda` | `[1.0343e17,1.8442e17]` |
| Proper expansion-second coefficient `Theta2_lambda` | `[8.3831e21,1.8496e22]` |

Both coefficients are strictly positive. The rule chose the smallest `10^n`,
`n>=0`, exceeding the exact interval-safe threshold (approximately
`5.18318e20`), then set `lambda=-10^n`. It used the authenticated exact base
intervals rather than displayed decimal approximations. This is postprocessing
of a fixed experiment, not selection of a new direction after seeing a result.

## Computation and evidence

The causal initial interval `[0,5/4]` is exactly zero on `[0,1]`. Its remaining
241 whole cells are the original `N1024` partition restricted to this cone,
with the original monotone flat-tail enclosure on `[1,65/64]`. New `Q` through
`Q_rrr`, `h`, `h_r`, `s=-p*zeta`, `s_r`, and the mixed velocity were enclosed
continuously on every cell. No sampled-value replacement was used.

Only authenticated old **background** fields `rho,p,a`, the global `p_y`, and
shell boxes were reused. Old `f,xi` profiles, preparation norms, and total
response caps were not used as bounds for this new direction. The new profile
and norm boxes therefore constitute new direction data under shared background
uncertainty; they are not an independent PDE observation.

In the conditioned five-mean order `(K,L,N,J,w)` and Gram diagonal
`(1,2,1,3,1)`, the new initial vectors are
`(h,0,0,h,s)`, `(h_r,0,0,h_r,s_r)`, and `(0,V_L,0,0,0)`.
Their componentwise maxima give fresh outward Gram caps and

```
N_h = max((D0+T*Dv)/3,Dr+Dv) < 4.525243e10,
E_h = q*N_h/(1-q) < 1.699018e7,
q   = 1201/3200000.
```

The full coupled homogeneous resolvent remains charged through `E_h`.
Unchanged bulk/shell first forcing, clock and proper-measure corrections, and
`Q1` cancel exactly in the **difference**. They remain part of the absolute
base response. The boundary error is `(x-1)_+ E_h`, so its weighted/unweighted
history charges are `delta^3 E_h/6` and `delta^2 E_h/2`. Whole-cell free moments
are integrated with exact positive weights. The checked kernel comparison
uses `eta=10^-11`; exact `0<=k_T<=1` justifies the expansion error allocation.

The source and local code bytes were pinned in `CODE_AND_INPUT_BINDINGS.json`
before scientific import. The retained source closure was authenticated before
and after each run. The reports retain all exact rational endpoints, raw bump
and background boxes, local-flux overlap checks, fresh profile vectors, caps,
histories, physical factors, sign gates, and finite-lambda calculations.

- `normal/RAW_WITNESS.json`: complete fresh whole-cell/target witnesses.
- `normal/REPORT.json`: exact rational numerical certificate and affine witness.
- `optimized/`: independent `python -O` execution of the same frozen producer.
- `NORMAL_OPTIMIZED_COMPARISON.json`: raw witness, report, and receipt are
  byte-identical between normal and optimized executions.
- `DISPLAY_AND_CALL_AUDIT.json`: outward display inclusion and no-oracle method
  audit.

Each run makes 241 whole-cell bump evaluations and one target bump evaluation.
There are **zero new background and zero evolved actual-first-field oracle
calls**; there are no new actual-oracle point queries. The retained numerical
module is imported only for exact outward bump/interval arithmetic.

## Scope

This proves conditional preparation sensitivity and a simultaneous reversal
of two formal second coefficients in the admissible affine class. It does not
provide a nonlinear epsilon radius uniform in lambda, a quantitative radius
for the displayed finite lambda, total contraction or expansion, a physical
radiation mechanism, or a cosmological prediction. All prior files and public
capsules remain unchanged.
