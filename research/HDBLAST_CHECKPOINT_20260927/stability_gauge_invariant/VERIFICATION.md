# Independent audit: linear stability of the static +1 branch (Method A, gauge-invariant)

Audit date: 28 September 2026. Audited folder: `stability_gauge_invariant/`, dated 27 September 2026. All auditor files are in `verify/`. No audited file was modified.

**Overall verdict.** The central result holds up under every check I ran:

- There is no scalar or tensor mode with μ² < 9/4 on the +1 branch, apart from the massless graviton, at all six detunings.
- The calibration tachyon of the original shell is reproduced.

I checked this in two ways:

- **Re-running the package.** The audited scripts reproduce their JSON outputs **bit for bit**; only runtimes differ.
- **Independent work.** I wrote my own derivation (a full nonlinear Einstein tensor, then linearized) and my own solver (a Prüfer angle, own shell polish, Radau/DOP853 variants). These reproduce every key number. The largest difference is 4e-10 in the calibration eigenvalue.

The stability is not marginal. A mode with μ² < 0 would need the shell coefficient B to be ≈ 0 (between −0.66 and +0.008, depending on δ). The actual value is B ≈ 3.48–3.56, so the margin is about 3.5.

The problems I found are minor: one misquoted tolerance, one mislabelled table column, a status label that mixes analytic and numerical evidence, and one post-hoc control criterion. No problem affects the conclusion. The result stays **numerical, not interval-certified**, and it is limited to the scalar and tensor sectors, as the README already says.

## Verdict table

| # | Claim (from the audited summary) | Verdict | Why |
|---|---|---|---|
| 1 | The +1 background is reproduced (φ_b, ρ_b, H², H/H0 match the stated digits), and the independent polish agrees with the package to "≤ 6e-15 relative in ρ_b" | **partially confirmed** | All stated digits match. The re-run is bit-identical. My own integrator and polish give B and ρ_b equal to ≤ 4e-15. However, the audited JSON's own maximum ρ_b difference is **8.5e-14** (δ=0.0003), not 6e-15. This is a harmless misquote. |
| 2 | The sympy derivation gives the gauge-invariant (X,Z) system, the shell condition M = B X + 3λZ/ρ² with B = φ″/φ′+σ″/2, no brane bending in longitudinal gauge, trace junction = momentum constraint, and the tensor equation with Neumann condition | **confirmed** | The re-run gives 55/55 checks and 4/4 controls detected. My independent derivation (`verify/indep_derivation.py`) uses a different method: the full nonlinear G_AB − T_AB, d/dε, a level-set normal on the displaced shell, and background relations solved rather than typed in. It gives 26/26 checks plus 4/4 wrong-formula controls detected. The trace-free junction is exactly −ζ. The scalar junction equals χ′ + σ″χ/2 + 2φ′ψ with ratio exactly 1. The claimed first-order system implies all five field equations. The tensor ratio is −ρ²/2 and the junction is ρ²h′/2. I did not independently re-derive the general-gauge invariance or the ℓ=1 statements; for those I rely on the audited script's re-run. |
| 3 | Calibration: the original shell has one bound state, μ² = −7.717871625260, growth 1.657193631259; δ=0.003 and 0.01 match | **confirmed** | The re-run is identical. My independent solver gives −7.717871625692 (difference 4.3e-10), −7.714023602725 and −7.700561458571, one root each, with root spread across variants ≤ 6e-8 (dominated by the rtol=1e-10 variant). The recorded Chat 9 value −7.717871625176 exists in `D-Blast 3/untitled folder 146/.../independent_GN/RESULT.json`. This is a comparison with the project's own earlier work, not with external literature, as the README states. |
| 4 | +1 branch: no scalar eigenvalue with μ² < 9/4 at any δ; min M̂ ≥ 0.99274; longitudinal cross-check finds no roots; Gauss residual ≤ 3e-15; precision spread ≤ 4.4e-16 | **confirmed** (numerical) | The re-run is identical. With a different formulation (Prüfer angle, 281-point grid, [−400, 2.2499]), my solver finds no root at any δ. Its min M̂ agrees with the audited values to ≤ 3e-15. One caveat: the "precision spread ≤ 4.4e-16" is measured on M̂, which is saturated near 1, so it says little. The independent agreement and the stability margin (below) are the stronger evidence. Also, the README column "M̂(0)" is actually evaluated at the grid point μ² = 0.00243, not at 0. That changes the value by only ~2e-8. |
| 5 | No complex scalar eigenvalues in the two rectangles (winding 0 on +1, 1 on calibration) | **partially confirmed** | I re-ran a subset (calibration, and +1 at δ=0.001 and 0.1, small rectangle). The windings reproduce exactly (1.0000000000000002, 2.7e-17, −1.9e-16), with max phase step 0.29 rad per segment and \|M\| far from 0. I did not re-run the large rectangle or the other δ. Complex modes outside the rectangles are excluded only by the conditional theorem (claim 6). |
| 6 | Analytic bounds: if B>0, the spectrum is real and μ² > −4; B = 3.479–3.555 on +1, and B = −2.68e-4 on the original shell | **confirmed** (conditional on the derived equations) | I re-derived both identities, the boundary terms and the cone vanishing (Re s > −3/2) by hand, and checked them in `verify/identities.py`, including a wrong-sign control. The growth criterion Re μ² < (Im μ²)²/9 ⇔ Re p > 0 has 0 mismatches in 20000 random samples. U″(1) = 28/9, L = 9 and Δ = 18 are exact, so the δ→0 value B → 2 + 14/9 = 32/9 is consistent. The README correctly calls the limit an observation. |
| 7 | ℓ=1 is not physical because Bφ′_b ≠ 0; the ℓ=0 static zero mode is absent (Jacobian condition number 15–17) | **confirmed** (numerical; ℓ=1 logic relies on the audited symbolic checks) | Re-run values: Bφ′_b ranges from −1.39e-4 to −4.46e-2. Jacobian singular values scale ∝ δ with condition number 14.7–17.0 (in the variables log10 e_h and y_b), so the Jacobian is non-singular. I also checked that φ′ has no sign change in the bulk at any δ, so g = φ″/φ′ is regular. |
| 8 | Tensor: only the normalizable massless graviton, no KK bound states in (0, 9/4), no tachyons (Sturm–Liouville) | **partially confirmed** | The Sturm–Liouville bound μ² ≥ 0 and the Neumann condition are correct: I re-derived the tensor junction independently. The absence of bound states in (0, 9/4) is **numerical**, not exact. It is confirmed by the re-run and by my Prüfer tensor solver, which finds roots only at \|μ²\| < 1e-13. The claim's "exact-verified" label overstates the numerical part. |
| 9 | Controls C1–C5 behave as stated | **partially confirmed** | The re-run of `controls.py` is identical in every key. C3, C4 and C5 are sound. My independent intercept for C4 is −7.7197956372 against the closed form −7.7197959184 (difference 2.8e-7). For C1, I confirmed exactly that C₁₄⁽²⁾(cosh y/9) solves the decoupled μ²=0 equation, with C₁₃ as a control. But the C1 docstring expects agreement at O(1e-8), while the observed deviation is O(δ). The pass criterion ("deviation/δ constant") was adapted afterwards. The O(η_b) explanation via U‴(1) = 100/9 is plausible by my order-of-magnitude estimate (≈ −1.7e-4 against −1.2e-4 observed at δ=1e-3). C2 finds its δ=0.001 root at μ² = −28555, far outside the main scan and through a different code path. My B* margin control (below) is a more pertinent test of detection power. |
| 10 | Final verdict: linearly stable in the scalar and tensor sectors at all tested δ; the only marginal mode is the graviton | **confirmed, within the stated scope** | Everything above supports it. The scope limits in the README are adequate: no vector sector, no interval certification, a regular-cone (horizon) normalizability criterion, no nonlinear or tunnelling analysis, no link to a radiation era, and no discovery claim. |
| – | Literature citations (Garriga–Sasaki, Frolov–Kofman, Gen–Sasaki, Garriga–Vilenkin, DeWolfe et al., Hawking–Hertog–Reall, BraneCode, Boucher, Townsend, arXiv:2609.21421) | **unverifiable in this audit** | The auditor's web-search budget was exhausted and paper sites are blocked. The descriptions of the classic references match the auditor's background knowledge, but I did not re-fetch them. The README says arXiv:2609.21421 was "seen only as a search result". Treat it as unverified. |

## Additional auditor result: stability margin (new control)

For real μ² < 0, μ² is an eigenvalue exactly when the shell coefficient equals B*(μ²) = −3λ Z_b/(ρ_b² X_b). I computed B* on the 281-point grid with the independent solver (`verify/indep_spectrum.json`, key `Bstar_control`):

| branch, δ | B* range over μ² ∈ [−400, 0) | actual B | margin B − max B* |
|---|---|---|---|
| original, 0.001 (calibration) | [−0.0283, 2.9e-4] | −2.68e-4 (inside the range, so a root exists, as found) | — |
| +1, 0.0003 | [−0.0021, 2.1e-5] | 3.5553 | 3.5553 |
| +1, 0.001 | [−0.0070, 7.1e-5] | 3.5548 | 3.5547 |
| +1, 0.01 | [−0.0697, 7.2e-4] | 3.5478 | 3.5471 |
| +1, 0.1 | [−0.659, 8.1e-3] | 3.4786 | 3.4705 |

So on the +1 branch, an instability would need the stabilizing shell term to be essentially cancelled. Only a gross error would do that, such as a sign error in σ″/2 (the C2 wrong-sign operator). My first-principles derivation excludes such a sign error, given the registered junctions ρ′/ρ = σ/6 and φ′ = −σ′/2, which it also reproduces. This is the main reason the stability conclusion is robust to numerical detail.

## Checks run (all outputs in `verify/`)

| Check | Script → output | Result |
|---|---|---|
| Re-run backgrounds | `rerun/reproduce_backgrounds.py` → `rerun/BACKGROUNDS.json` | identical to the audited file |
| Re-run symbolic derivation | `rerun/derive_linearized.py` → `rerun/DERIVATION_RESULTS.json` | 55/55, 4/4 controls; identical |
| Re-run full spectrum (2 processes, 4.6 min) | `rerun/spectrum.py` → `rerun/SPECTRUM_RESULTS.json` | identical except `runtime_s` |
| Re-run controls | `rerun/controls.py` → `rerun/CONTROLS_RESULTS.json` | identical |
| Re-run winding subset | `rerun/winding_subset.py` → `rerun/winding_subset.json` | identical windings |
| Independent derivation | `indep_derivation.py` → `indep_derivation.json` | 30/30 booleans true (26 checks + 4 controls) |
| Independent spectrum (Prüfer angle, own polish) | `indep_spectrum.py` → `indep_spectrum.json` | calibration −7.717871625692; +1: no roots; tensor roots only at 0 |
| Analytic identities, growth criterion, Gegenbauer/Δ | `identities.py` → `identities.json` | all true; 0/20000 mismatches |
| Comparison summary | `compare_rerun.py` → `compare_rerun.json` | as quoted above |

One negative note on my own tooling: my first independent background formulation integrated ρ″ as a free variable. It failed, with a relative constraint error of about 6e-7 at the shell. The reason is that the constraint C = H² − 1/ρ² − … is conserved, and round-off near the cone (H² ~ 1/y0²) is carried to the shell. I switched to the first-integral form; this is recorded in the script docstring. This supports the audited code's choice of imposing the constraint.

## Remaining limitations (for the audited result)

- The result is floating-point only, with no interval or validated-numerics certificate.
- μ² < −400 and complex μ² outside the rectangles are covered only by the conditional B>0 theorems.
- The vector sector is not computed. This is expected to be non-dynamical for a single Z2 shell with no bulk vectors, but that is not checked.
- The normalizability-at-the-cone (horizon) criterion is an assumption, standard in the dS-brane literature.
- Linear stability of this configuration says nothing about whether the dynamics reaches it, nor about a radiation era or observational support for the hypothesis.

## Reproduction

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant/verify
export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
(cd rerun && python3 reproduce_backgrounds.py && python3 derive_linearized.py && python3 spectrum.py && python3 controls.py && python3 winding_subset.py)
python3 indep_derivation.py     # ~15 s
python3 identities.py           # ~10 s
python3 indep_spectrum.py       # ~5.5 min, 1 core (reads ../BACKGROUNDS.json only as starting guesses)
python3 compare_rerun.py
```
