# Results and limits

Algebraic and numerical entries below come from completed JSON receipts in this capsule. Analytic existence and induction statements rest on the displayed proofs and their hypotheses.
Unfinished or failed runs are listed separately at the end.

## Reviewed follow-on results

| Result | Evidence and limits |
|---|---|
| Formal initial-data corrections through every finite order for reflection-invariant commuting quadrupoles | [Reflection theorem](continuation/REFLECTION_FORMAL_INITIAL_DATA.md); local KID characters, covariance and Taylor factors pass 266 exact checks and 4 wrong-claim controls in each normal/optimized mode. The PDE induction is an analytic proof, not a numerical solve. |
| Finite-slab second-order spacetime coefficients before shell contact; arbitrary finite orders in the reflection sector | [Pre-contact proof](continuation/PRECONTACT_SECOND_ORDER_DYNAMICS.md), normally hyperbolic Cauchy theory and gauge/constraint propagation. Uses the union of all correction supports and the actual controlled background interval. No through-contact solution. |
| Exact Gaussian-time initial accelerations and matter-only TT projection `−4𝒯/3` | Independent direct Ricci, ADM conversion, sphere decomposition and causal algebra: 17 checks and 6 wrong-formula controls in each normal/optimized mode. Full quadratic geometric source not evaluated; no claim that the total TT force is nonzero. |
| Validation remains active under Python optimization | Deliberate failures exit nonzero in both modes, with scientific normal/optimized receipts equal after removing the mode marker. See `reviews/reflection_runner_controls.json`, `reviews/PRECONTACT_VERIFICATION_RECEIPT.json` and `reviews/ANGULAR_INJECTED_FAILURE_RECEIPT.json`. |

The two new algebra suites total **283 positive checks and 10 wrong-formula/claim controls per mode**. These counts are separate from the inherited second-order capsule's scripts and do not measure the strength of a PDE existence proof. No convergence or exact nonlinear family follows.

## Claims and their evidence

| Claim | Status | Evidence |
|---|---|---|
| Lemma 1(a): Killing fields of `γ₀` on an annulus are the rotations | Proof (analytic) | Manuscript §2 |
| Lemma 1(b): metric-only adjoint kernel is `span{Ωx_a}` | Proof; `l = 0, 1` reductions symbolic for all `s`; `l = 2, 3, 4` ranks at two sample points | `new/kids_symbolic_l01.json`, `new/kids_tensor_only.json`; referee's own code agrees |
| Lemma 1(c): full KIDs = 4 boosts + 6 rotations | Symbolic | `basis/results/kids.json`, `kid_modes.json` |
| Theorem 2: compactly supported second-order data ⇔ 10 per-copy Taub integrals vanish | Proof, using Delay's Theorems 1.3/5.1; (API) via Lemma 8.1 for the metric operator and §9.3 for divergence; (KRC) via smooth kernel bases (arXiv:1003.0535v4) | Manuscript §3 |
| Explicit commuting correction: no radial derivatives of the coefficients survive | Exact symbolic | `new/fast_verify.json` |
| `P_M[w] = S_M` | Exact symbolic; exact rational (3 random points); independent numerics, max relative residual 5.3×10⁻⁴¹ at 10 points | `new/fast_verify.json`, `new/exact_random_point.json`, `independent/indep_check.json` |
| `P_H[h_T] = S_H − f₀` | Exact symbolic angular reduction for all free source/profile jets; finite exact rational checks (3 seeded points, residual exactly 0); independent numerics: remainder constant on spheres and equal to the engine's `f₀` to 20 digits | `new/angular_exact.json`, `manuscript/ANGULAR_EXACT_VERIFICATION.md`, `new/exact_random_point.json`, `independent/indep_check.json`, `independent/compare_f0.json` |
| `l = 0` identities and solvability with compact support | Symbolic; independent weighted-integral test, relative difference 3.2×10⁻²⁵ | `basis/results/l0_sector.json`, `new/l0_full_sector.json`, `independent/indep_l0.json` |
| Radial gauge annihilated; `ψ` absent from `m` | Symbolic (exact 0); independent: gauge 10⁻³⁶, `ψ`-only weighted integral 7.5×10⁻³⁰ relative | `new/l0_full_sector.json`, `independent/indep_l0.json` |
| Theorem 3: `m[w] = −(2/3)∫s^{3/2}Ω⁵f₀ ds` for `v` in reconstruction gauge | Proof (adjoint form); symbolic reduced densities | Manuscript §5; `new/l0_full_sector.json` |
| Corollary 3.1: velocity density `s^{7/2}[8Φ′²s²K² + 3(2sK′+7K)²]/(8Φ′²Ω) ≥ 0` | Symbolic, hand derivation, independent numerics (20 digits) | `new/wgg_canonical.json`, `independent/indep_check.json` |
| Corollary 3.2: `c₂ = 6s^{13/2}/(Φ′²Ω)`; no `G–K` cross term | Symbolic; referee by hand | `new/wgg_canonical.json` |
| `m` is **not** invariant under first-order gauge changes of `v` | Referee calculation (`m[L_ξ²γ₀] = +6.34×10⁻³` on the test background) | `reviews/REFEREE_REPORT.md` |

Background and first-order sanity checks in the independent code: background
Hamiltonian constraint ≤ 2.6×10⁻³⁹; first-order constraints ≤ 2.9×10⁻³⁹ at
all 10 points. The 5D background relations re-derived from the Ricci tensor
match the original engine's rules exactly (`independent/bg5d.py`).

Wrong-formula controls detected: sign-flipped trace-free metric
coefficient, doubled scalar-momentum coefficient, omitted `l = 0`
subtraction (`new/exact_random_point.json`, `new/fast_verify.json`).
The angular symbolic verifier also detects a sign-flipped `C₁` coefficient
(`new/angular_exact.json`). It uses the regenerated original source and an
independently derived Hamiltonian operator reduction; it does not independently
derive the source Hessian symbolically.

## What these checks do not cover

* Exact checks treat background and profile jets as independent symbols,
  which is stronger than needed. The exact-rational check evaluates three
  seeded sample points without rounding. It is a finite identity test,
  not a symbolic proof. No error probability is claimed: a quantified
  Schwartz–Zippel bound would require a numerator degree bound and an
  explicit analysis of the sampling distribution.
* The independent numerics use one exact test background
  (`1/Ω = 1 + s²/5 + s³/7`), not the registered regular branch. The
  identities are background-independent, so this tests the algebra, not the
  registered numbers.
* Delay's theorem is used as published; its proof was not re-checked here.

## Runs that did not complete

* `basis/code/verify_cc.py`, the companion session's original exact verifier,
  terminated after about 57 CPU minutes with no output (its log is empty).
  It establishes nothing either way. Its job is covered by the checks above.
* A full symbolic expansion of `P_H[h_T] − (S_H − f₀)` in an earlier version of
  `new/fast_verify.py` did not finish in about 45 CPU minutes and was
  stopped; the exact-rational test initially replaced it. The later exact
  angular reduction establishes the same Hamiltonian identity without that
  full coordinate expansion; the historical stalled run is not counted.

## Scope

The recovered localized manuscript concerns initial data for packets vanishing near the shell and the vertex. The separate [reflection proof](continuation/REFLECTION_FORMAL_INITIAL_DATA.md) supplies formal initial coefficients through every finite order in its symmetry sector; the [pre-contact note](continuation/PRECONTACT_SECOND_ORDER_DYNAMICS.md) extends them on a strict pre-contact slab.
Second-order evolution through shell contact, convergence, exact nonlinear
existence, stability, physical energy positivity and any cosmological or
Big Bang inference are not established. No external novelty is certified.
Internal AI review is not external peer review.
