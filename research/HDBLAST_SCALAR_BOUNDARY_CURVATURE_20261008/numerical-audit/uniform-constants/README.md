# Explicit uniform ranges and curvature suppression

8 October 2026. Exact arithmetic evaluation of the sufficient conditions in
`PROOF_CANDIDATE_USED.md`. The analytic proof and the arithmetic have both passed
separate internal AI-assisted independent reviews. This is not external peer
review, a novelty finding, APS acceptance or observational evidence.

The main output retains its honest status at computation time: its theorem
review was then pending. `REVIEW_ACCEPTANCE.json` binds the subsequently completed
analytic and independently implemented arithmetic reviews to the exact same
proof, producer and report. The recorded pending fields are preserved as
historical provenance.

All models use dimensionless units. Put

    W(η) = 3k + wη²/2 + vη³/6,
    f(η) = f0 + f1η + f2η²/2,
    U = (W′)²/2 − 2W²/3,
    σ = 2W + δf.

The scalar range used to bound polynomial derivatives is |η| ≤ 1/100. The exact
parameter values are:

| Model | k | w | v | f0 | f1 | f2 |
| --- | ---: | ---: | ---: | --- | --- | ---: |
| Registered exact decimal model | 1/9 | 2 | 2 | 1+c | c | 0 |
| Simple polynomial model | 1 | 5 | 0 | 1 | 1 | 0 |

Here **c = 2/(1.0357712571566784) − 4/3**, with the displayed decimal interpreted
as an exact rational number. This specifies the mathematical model and does not
assert that binary64 parameter rounding in the uploaded numerical program is
zero.

For the registered model, the certificate gives b = 1/89600 and

    δ* = 161839258930731 / 1482649651195872051200
       ≈ 1.0915542913337286 × 10⁻⁷.

For the simple model, b = 1/200 and δ* = 1/1000. All twelve sufficient
inequalities for each model hold in exact rational arithmetic. The theorem then
gives existence and uniqueness **inside its specified weighted small-field
tube**, together with uniform remainder bounds. It does not exclude unrelated
global branches.

The following slightly looser constants are easy to quote; their rounding and
range restriction are checked exactly in `READABLE_BOUND_CHECKS.json`:

| Model | Stated range | Kη upper bound | KH upper bound |
| --- | --- | ---: | ---: |
| Registered exact decimal model | 0 < δ ≤ 10⁻⁷ | 4318 | 335 |
| Simple polynomial model | 0 < δ ≤ 10⁻³ | 0.2042 | 0.0236 |

They mean

    |η_b + f1 δ/[4(w−2k)]| ≤ Kη δ²,
    |H² − [k f0 δ/3 + (f0²/36 − k f1²/[24(w−2k)])δ²]| ≤ KH δ³.

These conservative analytic constants were not estimated from sampled solver
residuals. For the registered model, even the smallest uploaded detuning,
0.0005, is more than 4,580 times the exact certified δ*. The simple polynomial
model was not one of the uploaded numerical families. **No uploaded numerical
run gains a certified ODE error bound from this calculation.**

## A strict finite-detuning consequence

Define q = k f1²/[24(w−2k)] and the comparator

    H²_metric = k f0 δ/3 + f0² δ²/36.

The exact checks show q − KH δ* > 0 for both examples, using the tighter exact
KH and δ* in the report. Therefore throughout each certified interval,

    −(q + KH δ*)δ² ≤ H² − H²_metric
                   ≤ −(q − KH δ*)δ² < 0.

With the outward rounded bounds and shorter ranges in the table, the readable
strict inequalities are

    registered: H² − H²_metric < −0.00089 δ²,
    simple:     H² − H²_metric < −0.0138 δ².

The comparator is physically specified: it is the exact constant-scalar branch
of the **separate constant-tension model f(η)=f0**, with the same bulk potential.
In the scalar-dependent tension model with f1 ≠ 0, setting η=0 generally violates
the scalar junction. This comparison concerns the response to changing the
boundary coupling; it does not assert that both branches solve identical
boundary conditions.

## Arithmetic and independent checks

Polynomial derivative bounds use the exact inequality
sup(|P|, |η|≤r) ≤ Σ |p_n|rⁿ. Elementary rational bounds
7 < exp(2) < 8, coth(1) < 4/3 and (1−exp(−2))⁻¹ < 7/6 eliminate
transcendental arithmetic. The positive lower bound on ζ and the lower bound
H0(1/k)² > k²/2 are derived in the report. Every rational output is exact;
floating-point display fields are only approximations.

The producer and its optimized-Python run produce byte-identical results. A
wrong proof hash is rejected before output creation. The independently written
SymPy checker imports no producer source, rederives all five polynomial
expressions and 39 constants per model, verifies the twelve inequalities,
recomputes both remainder bounds and the strict sign, and rejects altered
constant, inequality, sign and parameter inputs. See
`../../independent-math/INDEPENDENT_EXACT_UNIFORM_CONSTANTS.json` and
`../../independent-math/UNIFORM_PROOF_REVIEW.json` in the complete package.

To recompute the producer and readable bounds from this directory:

```sh
python -I -B evaluate_uniform_constants.py --theorem PROOF_CANDIDATE_USED.md --theorem-sha256 544c6dbced51fbc68dccab3bbbe84210545849c6e7447e1b6bdae0fc98e76f51 --out EXACT_UNIFORM_CONSTANTS.json
python -I -O -B evaluate_uniform_constants.py --theorem PROOF_CANDIDATE_USED.md --theorem-sha256 544c6dbced51fbc68dccab3bbbe84210545849c6e7447e1b6bdae0fc98e76f51 --out EXACT_UNIFORM_CONSTANTS_OPTIMIZED.json
python -I -B check_readable_bounds.py --input EXACT_UNIFORM_CONSTANTS.json --input-sha256 95923d9824f3f367385f67392ac10391b489337138851b642aa4553377d11ba1 --out READABLE_BOUND_CHECKS.json
```

These commands use the Python standard library and perform no ODE integrations
or network requests. Existence and remainder estimates do not settle dynamical
stability, cosmological energy transfer, reheating, observational discrimination,
external originality or APS suitability.
