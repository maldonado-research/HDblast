# Independent verification: constraint control in the balanced-disturbance evolutions

This is an adversarial audit of the workstream in this folder (README.md, lab.py, controls.py, calib.py, analyze.py, runs/). The audit changed none of the workstream's files. Every audit script and output is in `verify/`. All checks are floating-point numerical experiments or exact SymPy identities. None of them is an interval certificate.

## Verdict table

| # | Workstream claim | Verdict | Basis (files in `verify/`) |
|---|---|---|---|
| 1 | The baseline and its stall reproduce. H_max(t=0.5) is 0.245, 0.0449, 0.00755, 0.00638, with finest-pair orders 0.24 (t=0.5) and 0.27 (t=1). The lab code matches the wrapper to 7.7e-10. | **confirmed** | `rerun/base_h4e-4` and the short run `rerun/short_base_h5e-5` are bit-identical to the saved runs (`compare_reruns.json`). The table, recomputed from the saved JSON with independent code, matches the README to printed precision (`recompute_saved.json`, no mismatches). The archived checkpoint values 0.0447, 0.00754 and 0.00655 match `src/REFERENCE_MATCHED_TIME_REFINEMENT.json`. |
| 1a | The 0.04 % gap from the archive (0.300596 vs 0.300474) is due to "floating-point/library differences". | **unverifiable** (attribution) | The archived and rerun wrapper JSONs have the same code SHA-256, dt, node count and pinned numpy 2.3.5 / scipy 1.16.3, yet H at t=0 already differs by 2.6 % (2.988e-7 vs 3.067e-7). The difference therefore sits in the sampled initial data and depends on the platform, not on the library version. This fits the noise-floor diagnosis (claim 3), but its precise cause was not found. |
| 2 | Mechanism: the front is seeded near the shell for t≲0.05, and e^{3A}(H+2M) is roughly conserved along z≈−t. The gain from z0=−0.003 is about 2.4e3 at t=0.5 and 5.8e3 at t=1. | **partially-confirmed** | Gains recomputed: 2388 and 5774. At h=5e-5, C₊^w is 2.22e-4 at t=0.05 and 1.16e-4 at t=0.5, so there is no bulk growth. The transport law was re-derived independently (row 5). **Wording error:** the README says the "unweighted, background-normalised H" grows by G. G applies to the *un-normalised* H+2M. The normalised H/S0 also picks up S0(z0)/S0(front) ≈ 9951/27.4 ≈ 364, and the measured growth of max H/S0 is ≈7.4e4. **Caveat:** C₊^w at t=0.05 depends on the time step (CFL 0.2 gives 3.95e-4, CFL 0.4 gives 2.22e-4). It agrees only by t=0.1, so the t=0.05 table is a transient diagnostic. |
| 3 | The stall comes from a non-decreasing initial-data defect floor near the shell, not from xtol and not from Δt. | **confirmed** (numerical; the origin in "sampled profiles" is an inference) | The weighted initial defect peaks at z≈−0.0021. It rises from 2.93e-4 to 4.87e-4 (unmasked: 5.54e-4) between h=1e-4 and h=5e-5. Its odd–even (Nyquist) part grows from 1.6e-4 to 4.1e-4 (`recompute_saved.json: initial_defect_near_shell`), which is the signature of noise amplified by D2 rather than of truncation. The xtol=1e-300 run is identical at t=0.05 and t=0.1, and so is CFL 0.2 at t=0.1. The source of the noise (dense-output interpolation versus cancellation in the profile differences) was not isolated, by the workstream or here. |
| 4 | Remedy A (discrete projection) removes the stall. H_max(t=0.5) is 0.227, 0.0483, 0.00701, 9.35e-4 (orders 2.23, 2.78, 2.91). At t=1 the orders are 2.22, 3.15, 3.13. Causal-window L2 orders are 3.48 and 3.38. | **confirmed** with an extra caveat | The reruns `proj_h4e-4`, `proj_h2e-4` and `short_proj_h5e-5` agree with the saved runs to ≤1e-9 relative at the snapshots. Only the t=0 roundoff values (~1e-14) differ. Newton takes 2–3 iterations. The numbers and L2 orders were recomputed independently. **Missing caveat:** the projection correction does not converge below h=1e-4 (+1.49e-9, then −1.54e-9, with a sign flip). As a result, the self-difference of the warp field *a* converges at only 1.5–2.2 at t=0.5 (magnitude ~1e-8). |
| 5 | Damped transport law: (∂t−∂z)[e^{3A}C₊] = −κe^{3A}C₊ and (∂t+∂z)[e^{3A}C₋] = −κe^{3A}C₊(A_t−A_z)/(A_t+A_z). The wrong sign and the wrong denominator are rejected. | **confirmed** (exact) | `sympy_independent.py` does its own jet-space derivation, with its own total-derivative operators and recursive elimination of second time derivatives. κ may depend on z. It recovers the undamped law as a known limit and rejects both wrong controls. It also derives from the Einstein tensor that H = 2(G−T)_tt and M = (G−T)_tz for the action R/2 − ½(∂φ)² − U, and that EA, EB, EP reproduce (G−T)_zz and (G−T)_xx modulo H. A doubled-kinetic control fails, as it should. The workstream's `controls.py`, rerun into `verify/controls_rerun.json`, gives 17/17 PASS with identical values. The numerical ratios 0.090 vs e^{-2.5} and 3.31 / 11.2 vs e^{1.25} / e^{2.5} match the saved JSON; the agreement is within about 10 %. |
| 5a | Junction and normalisation conventions used by the lab (reused from the frozen solver). | **confirmed** (exact) | The boundary slopes ga and gf equal e^{-B}A_z = σ(φ)/6 and e^{-B}φ_z = −σ_φ/2, expanded exactly in (b, f). gpa and gpf are their exact time derivatives. A σ/3 control fails. The BPS flow φ_y = W_φ, A_y = −W/3 solves the static equations with U = W_φ²/2 − 2W²/3, and σ = 2W reproduces the one-sided Z2 slopes at δ=0. The lab's perturbative H, M, damping denominator and weight equal the continuum expressions. |
| 6 | Remedy B (κ=10, switched off near the shell) is stable to t=1 and lowers H by **13–100×**. B alone stalls at t≤0.5 (order 1.45). A+B gives the smallest residuals. | **partially-confirmed** | The reruns `k10_h4e-4` and `k10p_h4e-4` reproduce, one bit-identically and the other to ~1e-10. It is stable on four grids to t=1. **The factor range "13–100×" does not reproduce from the saved data.** B/baseline is 39–98 at t=0.5, 26–194 at t=1 and only 8–9 at t=0.25; (A+B)/A is 24–42 at t≥0.5. "A+B smallest" holds at t=0.5 (2.28e-5). At t=1, B alone is marginally smaller (1.433e-4 vs 1.444e-4). |
| 7 | Remedy C (order 6 + projection) gives the highest finest-pair orders: 3.85 / 3.67 (max) and 3.97 / 3.91 (L2). | **confirmed** for H, with a caveat | Rerun `o6p_h4e-4` matches to ≤1e-8 relative. The orders were recomputed. **Caveat:** C's field self-differences are *not* the best. Near the shell (z>−0.2) its finest-pair orders are 0.07–0.75 at t=0.5. The magnitudes there, ~1e-10, are at the roundoff floor of the degree-7 closure. Its a and b fields converge at 0.7–1.6. |
| 8 | Smooth, corner-free calibration: H orders 4.00, 3.97 (t=0.25) and 3.99, 3.99 (t=0.5); field orders 3.96–4.01 until a floor of 1e-9 to 1e-10. | **confirmed**, with a scope qualifier | The values in `ANALYSIS.json` match the README. The calibration window excludes the shell, so it tests only the interior scheme and RK4. New supporting evidence: in the *baseline seed* runs, every field converges at order 4.0 at the finest pair within z>−0.2 (`recompute_saved.json`). The shell closure is therefore fourth order too, as far as can be seen. |
| 9 | Asymptotic fourth order for the ε=0.01 seed is not demonstrated; only φ_t self-converges at ~2.1–2.4. | **confirmed as inconclusive**, but the caveat is incomplete | B_t (pb) also converges at only 2.0–2.2 in the baseline and in A, and f at 2.9–3.0. **New localisation:** the largest φ_t self-difference lies on the shell-seeded outgoing front, z = −0.4976 at t=0.5 and −0.9968 at t=1. It does not come from the outer taper (≤1e-11 in the taper's causal future) or from the near-shell region (order ≈4 there). Its L2 size is identical in the baseline, A, B and A+B (3.45e-4 at t=0.5), so no constraint remedy touches it. Pointwise it is 6e-3, against max|φ_t| ≈ 2, which is about 0.3 %. |
| 10 | Negative results: uniform damping becomes unstable (H≈56 at t=1, a grid-scale pb sawtooth growing at ≈2κ); H-only damping hits the guard; KO 0.2, tighter inversion and order 6 without projection give no real gain; the degree-7 closure is roundoff-limited; projection alone at h≥2e-4 gives no gain. | **confirmed** | Saved JSON: H=55.8 at t=1. The Nyquist part dominates pb near the shell, and pb grows by e^{20.6·0.25} between t=0.75 and 1 (≈2κ). H-only stops on the guard at t≈0.050. An independent Fourier analysis gives a high-k growth rate κ(A_z/A_t−1)/2 for H-only damping (80 at r=17, κ=10) and no growth for C₊ damping. The short-run numbers match. The rerun controls reproduce the order-6 closure finding. A vs baseline at h=4e-4 and 2e-4 is 0.227 vs 0.245 and 0.0483 vs 0.0449. |
| 11 | Literature context: Gundlach et al. gr-qc/0504114, Pretorius gr-qc/0407110, GRChombo 1503.03436. No methodological novelty is claimed. | **unverifiable** in this session | The web-search budget was exhausted and arXiv is blocked, so the citations could not be re-checked online. The novelty disclaimer is appropriate. |

## Overall assessment

The workstream's numbers are real and reproducible:

- Every rerun matches the saved output, bit-identically or to ≤1e-8 relative.
- Every table was recomputed with independent code.
- The exact damped-transport identity and the model equations and junctions were re-derived from scratch.

The central conclusions survive:

- The stall is caused by a non-decreasing initial-data defect near the shell.
- Discrete projection (A) restores convergent residuals at about third order.
- Damping (B) lowers residuals strongly.
- Asymptotic fourth order for the seed is not shown.

The corrections needed are modest:

1. The "13–100×" range for Remedy B is wrong. It is about 25–200 at t=0.5–1 and only about 8 at t=0.25.
2. The geometric gain G applies to the un-normalised H+2M, not to the background-normalised H. The latter has an extra factor of about 360.
3. The order caveat should say that B_t, as well as φ_t, converges at about 2. It should also say that the deficit sits on the shell-seeded front z≈−t and is untouched by every remedy.
4. The projected warp and Remedy C's fields do not self-converge well at the finest pair. The magnitudes are ~1e-8 to 1e-10.
5. At t=1, "A+B smallest" is a tie with B.
6. The t=0.05 early-density diagnostic depends on the time step.
7. The archive difference cannot be blamed on library versions, which are identical.
8. The control named "projection change … converges to the continuum data" is true only between h=4e-4 and 2e-4. The README's limitations section already admits this.

None of these changes a label from positive to negative. The README's scope statement is accurate: this is a numerical-method result in a 1+1 model and says nothing about the physical fate of the shell or the Big Bang hypothesis.

## Checks run

| Script (in `verify/`) | What it does | Output |
|---|---|---|
| `rerun_subset.sh` | Reruns h=4e-4 for the baseline, A, B, A+B, C and κ=−5; A at h=2e-4; short (t≤0.1) baseline and A runs at h=5e-5. Single core, ≤20 min each. | `rerun/*.json`, `*.npz`, `*.log`, `rerun/exit_codes.txt` (all exit 0) |
| `compare_reruns.py` | Compares records, snapshots and final states with the saved runs. | `compare_reruns.json` |
| `controls_rerun.py` | Runs the workstream's `controls.py` unmodified, redirected to `verify/`. | `controls_rerun.json` (17/17 PASS, identical values), `controls_rerun.log` |
| `sympy_independent.py` | Einstein tensor from the metric; H, M and evolution equations; a jet-space proof of the damped and undamped transport laws; negative controls; Fourier stability of H-only vs C₊ damping; lab perturbation mapping; junction slopes and their time derivatives; BPS normalisation. | `sympy_independent.json` (28/28 PASS, ≈5 s) |
| `recompute_saved.py` | Independent recomputation of the README tables (H_max, orders, causal L2, early densities); field self-differences in five windows with the argmax location; normalisation (S0) check; odd–even decomposition of the initial defect. | `recompute_saved.json` |

Tolerances:

- A README value "matches" when it equals the recomputed value rounded to three significant figures (±1 in the last digit).
- Orders match to ±0.011.
- ANALYSIS field orders are reproduced to 1e-6.

Precision is double throughout. The environment was Python 3.11, numpy 2.3.5, scipy 1.16.3 and sympy 1.14, with `OMP_NUM_THREADS=1`.

## Reproduction

From `verify/`:

```sh
./rerun_subset.sh                                  # ~15 min on one core
python3 -B compare_reruns.py
OMP_NUM_THREADS=1 python3 -B controls_rerun.py > controls_rerun.log
python3 -B sympy_independent.py
python3 -B recompute_saved.py                      # reads ../runs (saved workstream data)
```
