# Independent audit: `stability_time_domain` (Method B)

Audit date: 28 September 2026. The auditor did not write the audited code. Audit files are in `verify/` only; nothing else in this folder was changed. All numbers below come from scripts in `verify/` with JSON outputs. This is a floating-point audit; nothing here is an interval certificate.

## Bottom line

The central claim survives the audit: **for homogeneous (FRW-symmetric) perturbations, the static +1 branch shows no growing or slowly decaying physical mode, and the leading physical rate is the continuum edge Re λ = −3/2** (numerical). Three independent lines of evidence support it:

- a bit-for-bit re-run of the audited scripts;
- an independent Chebyshev multi-domain discretisation of the continuum equations, on an independently re-solved background;
- a decoupled-scalar Liouville bound, μ² > 9/4. This bound is conditional on neglecting the metric coupling, which the audit shows to be negligible on this branch.

Several statements in the README and the summary are overstated or wrong in detail. They are listed under "Corrections".

## Verdict table

| # | Claim (audited) | Verdict | Audit evidence |
|---|---|---|---|
| 1 | Calibration: time-domain rate 1.657193368 (linear) and 1.657193408 (nonlinear) vs the Chat 9 reference 1.6571936312 | **confirmed** | Re-runs are bit-identical (`v4…json → rerun_vs_original`). An independent ESPRIT estimator on the saved series gives 1.65719331 at h=5e-5 (−3.2e-7 vs reference) and 1.65719335 (nonlinear). The independent Chebyshev discretisation gives 1.6571938 at N=24 and 1.6571978 at N=32 (`v3_chebyshev.json`). The reference μ² = −7.717871625 maps to λ = 1.6571936312 through λ² + 3λ + μ² = 0 (`v1`). The reference is internal to this project (Chat 9), not external. |
| 2 | Above Re λ = −3/2 there are only gauge modes, far-boundary artifacts and grid-scale oscillations, none with shell support | **partially confirmed** | Reproduced. The Chebyshev check finds only λ ≈ 1 (gauge; f_b ≤ 5e-10) and the far-boundary ladder (0.1035…, f_b ≤ 3e-10) below \|Im\| = 1000. **However**, the semi-discrete operator has **1,166–2,397 eigenvalues with Re λ > 0** per grid. Most are the far-boundary ladder. There are also grid-scale shell-localised modes: Re λ = +0.2477, +0.2459, +0.2440 at h = 4, 3, 2e-4, with f_b ≈ 3.7e-5, labelled constraint-violating. Their real part does **not** shrink under refinement. They are harmless only because RK4 at dt = 0.4h damps them: the fully discrete maximum rate is the gauge value 1.000, confirmed from the saved eigenvalues. The summary headline "the only positive rates in the discrete system are a gauge mode and far-boundary artifacts" is therefore wrong for the semi-discrete system. Item 2 of the README does mention the grid-scale modes. |
| 3 | All shell-supported modes lie on Re λ = −3/2 to ≤1e-8; no bound state with μ² < 9/4 | **confirmed (numerical), with independent support** | Recomputed from the saved eigenvalues: maximum deviation 8.8e-9 to 1.03e-8. The Chebyshev check finds no shell-supported real or low-frequency mode above −3/2 at N = 24 or 32. The decoupled-scalar Liouville potential satisfies V_eff − 9/4 ≥ 0.0055 > 0 on the grid, and the boundary coefficient β = ρ_b(σ″/2 − σ/4) = +238 > 0. So μ² > 9/4 holds for the decoupled operator (conditional). The lowest decoupled μ² is 3.338 (a box mode). |
| 4 | The time-domain rate converges to −3/2 (Richardson −1.50012, second "shell-localised" profile −1.50017) | **partially confirmed; one statement refuted** | The numbers are reproduced exactly. (a) The full h-sequence −1.5370, −1.5393, −1.5062, −1.5011 is **not monotone**. The 3-point Richardson value depends on the assumed order: −1.50012 at the fitted order 2.69, −1.49935 at order 2, −1.50072 at order 4. Only "−1.500 ± 1e-3" is supported. (b) **The "shell-localised" profile is not shell-localised.** Its tuning weight is β = −4.13e9, so after normalisation the near-shell bump has relative weight 2.4e-10. The data is in effect one bump on z ∈ [−0.4, −0.2]. The default profile is likewise dominated by the bump on [−0.16, −0.08]; its near-shell bump has weight 3.8e-5. (c) The evolved signal barely reaches the shell (peak \|f_b\| ≈ 8.5e-9 of the bulk amplitude). Any continuum superposition decays at −3/2 by construction. So the time-domain test is insensitive to a weakly bound, weakly excited mode, and the spectral evidence (item 3) carries the conclusion. |
| 5 | The shell scalar decays promptly, with a floor that shrinks under refinement | **confirmed as a number; weak as evidence** | Max \|f_b\| for t > 5: 1.7e-12, 1.1e-12, 2.9e-13 and 2.5e-14 on the four grids. See 4(c) for why the test is weak. Not reported by the audited work: the near-shell energy (z > −0.1) decays at half-slopes of only −0.85, −0.46 and even +0.20 (bumpB at h=5e-5) over t ∈ [1, 6]. This is energy redistributing toward the shell from initially empty regions. It is not an instability, since it is bounded by the total E, which decays at −3/2. |
| 6 | Wrong-formula control: σ″ → −σ″ gives λ = 167.48925, with partner sum −3 | **confirmed** | Shift-invert re-run bit-identical. ESPRIT on the saved series gives 167.48924533. Chebyshev gives 167.4892458 (N=32), with pair sum −3.0000014. The decoupled scalar alone gives μ² = −28555.12 → λ = 167.4892478, which shows the metric coupling is negligible on this branch. Caveat: the control shows only that a **strong** tachyon is detected. It does not test sensitivity to a weak instability. |
| 7 | Perturbed detuning δ = 0.003 and 0.01: same structure | **confirmed** | The δ = 0.01 shift-invert re-run is bit-identical. The independent background solver confirms that the δ = 0.003 and 0.01 backgrounds are the +1 branch (ρ_b and φ_b agree with the archive to ≤3e-15). |
| 8 | Constraints converge near the shell and fail to converge in the far region (negative result) | **confirmed (by reproduction)** | The values in KEY_RESULTS follow from bit-identical re-runs. The auditor did not re-derive them independently. The linearised constraint expressions are consistent with G_tt − T_tt (factor 2) and G_tz − T_tz (factor 1), derived symbolically (`v1`). |
| 9 | Precision control (η-polynomial vs registered φ-polynomial) | **confirmed** | Saved eigenvalues: gauge 0.99999656 vs 0.99999656; far-boundary mode differs by 6e-9. The η-coefficients of U are verified exactly (`v1`). |
| – | Background: φ_b, ρ_b, H², H/H₀ | **confirmed** | The independent solver (first-integral form, event-located metric junction, 1-D brentq root) gives ρ_b = 129.92476284968 (rel. diff 3.8e-14) and φ_b = 0.9999159473169134. It finds exactly one root in log10(−η_h) ∈ [−40, −3]. H/H₀ = 0.6067217306. |
| – | Equations, junctions, normalisations | **exact-verified (symbolic)** | 32/32 sympy checks, including 3 deliberate wrong-formula controls that are correctly detected. The 1+1 evolution equations for A, B and φ follow exactly from G = T with κ² = 1. The junction normalisation is fixed by the BPS limit (A′ = W/3 = σ/6, φ′ = −W_φ = −σ′/2 at δ = 0). The nonlinear and linearised junction data, σ″ = 4φ, U(φ) and U(1+η), and the regular-cone series all match. |
| – | Literature attributions (Garriga–Sasaki hep-th/9912118; DeWolfe et al. hep-th/9909134) | **unverifiable here** | The auditor's web-search budget was exhausted, so the URLs could not be re-fetched. The same identifiers appear in this checkpoint's `literature/` workstream. The README uses them only for textbook background, not for results. |

## Corrections requested

1. **Summary headline.** Replace "the only positive rates in the discrete system are a gauge mode and far-boundary artifacts" with a statement about the **fully discrete (RK4) system** or the **shell observables**. The semi-discrete operator also has grid-scale, shell-localised, constraint-violating modes with Re λ ≈ +0.245 that do not shrink under refinement.
2. **README item 4.** Change "A second, shell-localised initial profile" to say the second profile is dominated (weight ratio 4e9) by a bulk bump on z ∈ [−0.4, −0.2]. The first profile is dominated by the bump on [−0.16, −0.08]. Say that the two-bump construction effectively yields one bulk bump.
3. **Richardson estimates.** Quote the result as −1.500 ± 1e-3. Note that the h=4e-4 point breaks monotonicity and that the extrapolation is order-sensitive (−1.4994 to −1.5007).
4. **Time-domain caveat.** State explicitly that the time-domain evidence cannot exclude a weakly bound shell mode with small overlap: the peak shell response is ~1e-8. The spectral evidence carries claim 3.
5. **Gauge template match.** Minor: "matches to 1.5e-8" is the best grid (h=1e-4). At h=5e-5 the residual is 3.7e-6; still clearly gauge.

## Checks run (all in `verify/`)

| Script | What it does | Output |
|---|---|---|
| `v1_sympy_equations.py` | Christoffel → Ricci → Einstein for e^{2B}(−dt²+dz²)+e^{2A}dx₃². Evolution equations, constraints, static first integral, BPS/junction normalisation, junction expansions, U coefficients, cone series, dS relation, plus 3 wrong-formula controls. | `v1_sympy_equations.json` (32/32 pass) |
| `v2_background.py` | Independent re-solution of the +1 branch at δ = 1e-3, 3e-3 and 1e-2, with a root scan | `v2_background.json` |
| `v3_chebyshev.py` | Independent 12-subdomain Chebyshev collocation (N = 24, 32) of the linearised continuum system with junction and Dirichlet rows. Original-shell calibration, +1 spectrum, σ″-flip control, and the decoupled-scalar Sturm–Liouville problem with its Liouville bound. | `v3_chebyshev.json`, `v3_chebyshev.log` |
| `rerun_A.sh`, `rerun_B.sh` | Unchanged audited scripts re-run into `verify/rerun/`: 1 dense spectrum, 4 shift-invert runs, 10 evolutions | `verify/rerun/*` |
| `v4_compare_and_extract.py` | Re-run vs original comparison; ESPRIT and log-slope re-extraction; Richardson sensitivity; line deviation and RK4 rates from saved eigenvalues; initial-data composition | `v4_compare_and_extract.json` |

**Tolerances and convergence.**
- Re-runs agree with the originals to 0.0 (bit-identical) in every compared quantity.
- On the original shell, the Chebyshev eigenvalue moves by 4e-6 between N = 24 and N = 32. That sets the accuracy of that cross-check. Its spurious high-frequency modes (\|Im λ\| > 1e4) are discarded by a \|Im\| < 1000 filter.
- A first Chebyshev attempt with N up to 40 was stopped for time; its log is kept in `v3_chebyshev_first_attempt_N24_32_40.log` (N = 40 gave 1.6572153, which shows round-off growth at high N).
- Decoupled-scalar results agree between N = 32 and N = 40 to 1e-7.

**Not checked.** Inhomogeneous (k ≠ 0) and tensor perturbations; nonlinear stability. The audited work states the same limits. The sibling `stability_gauge_invariant` workstream also reports stability; that is a consistency check, not part of this audit.

## Reproduction

```bash
cd research/HDBLAST_CHECKPOINT_20260927/stability_time_domain/verify
python3 v1_sympy_equations.py            # ~1 min
python3 v2_background.py                 # ~1 min
OMP_NUM_THREADS=2 python3 v3_chebyshev.py 24,32   # ~8 min
./rerun_A.sh & ./rerun_B.sh; wait        # ~12 min on 2 cores
python3 v4_compare_and_extract.py
```
