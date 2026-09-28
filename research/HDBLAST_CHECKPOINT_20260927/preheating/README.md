# Instant-preheating channel on the archived HDBLAST shell roll-off

Research checkpoint workstream, 27–28 September 2026. This folder turns the matter extension's *conditional screening* (checkpoint of 22 September 2026, `matter/MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md` §6–7) into an explicit mode-function particle-production calculation on the actual archived shell trajectory. Nothing outside this folder was modified.

**Bottom line (conditional, negative):** on the fixed archived trajectory, with the χ field's masses kept below the five-dimensional gravity scale M₅=κ₅^(−2/3), this channel does **not** produce a radiation era for any coupling scanned (G=10 to 10⁶, crossings φ*=0.25, 0.5, 0.75, 0.9, decay Yukawa y≤1). The produced energy is at most 1.7×10⁻⁴ of the residual-vacuum radiation threshold for couplings that pass the local screen, and at most 2.8×10⁻³ over the whole scan. Reaching the threshold would require χ masses ≳18 M₅, or ≳6 M₅ if only post-crossing masses are constrained. In any regime where the produced energy *does* reach the threshold, its backreaction on the scalar junction is 10–80%, so the fixed-trajectory calculation stops being self-consistent there. The bulk-coupled feedback, meaning all three modified junctions, is **not** solved here.

## 1. Question

The model has a shell scalar χ with m_χ² = m₀² + ḡ²(φ_b−φ*)². The shell field φ_b(τ) rolls through φ* on the archived registered-detuning trajectory (Chat 14, +1 fate). For which couplings, if any, can the χ quanta created at the crossing:

1. be computed reliably from mode functions, including expansion and the curved trajectory;
2. carry enough energy to exceed the residual-vacuum radiation thresholds of the matter report;
3. leave the φ trajectory unaltered, as the fixed-trajectory calculation assumes;
4. convert into radiation through the optional Yukawa channel χ→ψψ̄?

## 2. Method

**Background.** The background is read (read-only) from the checkpoint's copy of the completed Chat 14 archive, `source_audit/inputs/CHAT14_COMPLETED_REFERENCE.zip` (SHA-256 `c68e1de1…fdf59d`, recorded in PREHEATING_RESULTS.json and CONTROLS.json). The primary run is `v2_t1e3_plus_c4_finecoarse`, which has the finer far grid that the Chat 14 notes recommend. Runs `v2_t1e3_plus_c4` and `v2_t1e3_plus_c4_fine` serve as grid controls, and `v2_t1e3_minus_c4_finecoarse` supplies a secondary crossing at φ*=−0.5.

Data are used only for coordinate t<8.5, before the earliest far-taper signal (s=H₀τ≤6.92). On the minus branch they stop before the reversal (h<0.3, s≤5.83).

φ_b(s) and h(s)=H/H₀ are represented by quintic splines in shell proper time. The scale factor is ln a=∫h ds. It agrees with the independent ∫(1+∂_t a)dt to 3×10⁻⁶.

**Units.** All rates are in H₀=1/ρ_b0, the initial static shell rate, with ρ_b0=78.828 in registered units. The dimensionless coupling is G=ḡ/H₀, identical to the checkpoint screen's G, and q=G|dφ/ds| at the crossing. Gravity enters only through b ≡ κ₅²H₀³ = (H₀/M₅)³, and every density comparison is linear in b.

**Modes.** For X_k=a^{3/2}χ_k:

X″ + Ω²X = 0, Ω² = κ²/a² + μ₀² + G²(φ−φ*)² − (9/4)h² − (3/2)h′

Here κ is the comoving momentum with a=1 at the crossing. The modes start in the second-order adiabatic vacuum at the left edge of an adiabatic window, where the κ=0 adiabaticity max(|Ω′|/Ω², |Ω″|^{1/2}/Ω^{3/2}) falls below 10⁻³. They are integrated with DOP853 at rtol 10⁻¹⁰, and n_k=|β_k|² is read in the same-order adiabatic basis at the right edge. When the edge is the end of the data, the edge adiabaticity is recorded.

Momentum integrals use 48-node Gauss–Legendre quadrature, with the grid extended until the tail is below 10⁻⁹ of the peak. For every case the script also compares against the instant-preheating result n_k=exp(−π(κ²+μ₀²)/q) and its integral N=q^{3/2}/(8π³)e^{−πμ₀²/q} (Felder, Kofman and Linde, as cited in the checkpoint: https://arxiv.org/abs/hep-ph/9812289).

**Densities and subtraction.** After the crossing, densities come from the late-time adiabatic (Parker–Fulling type) particle description with the final n_k: n=∫κ²n_k/(2π²a³), ρ=∫κ²n_kω_k/(2π²a³), p, ⟨χ²⟩=∫κ²n_k/ω_k/(2π²a³), and j=G²(φ−φ*)⟨χ²⟩. These are evaluated from the point where the κ=0 adiabaticity is below 0.05. Interference terms oscillating at 2m are dropped.

Stated renormalization: adiabatic subtraction of the vacuum to the order that removes the UV divergences. The finite local vacuum-polarization (Coleman–Weinberg-type) terms are **assumed absorbed** into the renormalized registered tension σ(φ), which is the matching condition the matter report requires. Their size is reported separately (§4, R8) because it is not small. For the adiabatic method see Parker and Fulling, Phys. Rev. D 9, 341 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.9.341).

**Energy comparisons.** These use the matter report's exact benchmark threshold κ₅²ρ_crit=√(σ_f²+36H_vac²)−σ_f, with σ_f and H_vac from the +1 static branch (`static_branch/PLUS_BRANCH_RESULTS.json`). In H₀ units this gives σ_f/H₀=52.678, H_vac/H₀=0.60672 and κ₅²ρ_crit/H₀=0.12563. The deceleration threshold is 0.12534. The same formula evaluated with σ(φ*) and H(s*) gives a crossing-time comparison.

The ratio R=bρ̂_χ(end)/0.12563 is the largest radiation/vacuum ratio obtainable if *all* χ energy at the data end became radiation instantly.

**Effective-theory condition.** The M₅ cutoff requires every χ mass and √q to lie below λ_c M₅. That means b≤(λ_c/Λ̂)³, with Λ̂ the largest χ mass over the analysed interval ("full"), or over the post-crossing part only ("post", a generous relaxation). The results quote R at λ_c=1 and the λ_c needed for R=1 (N=0) or for one radiation-dominated e-fold (R=e⁴).

**Backreaction.** R_j = κ₅²|j|/|σ′(φ)| is the fractional change of the scalar Neumann datum nφ=−(σ′+κ₅²j)/2. By the budget identity it is also the ratio of the work done on χ to the tension energy released.

**Decay.** Γ=y²m_χ/(8π) for χ→ψψ̄ with massless Dirac ψ, time-dilated per mode. The ledger is ρ̇_r+4Hρ_r=Q. Beyond the data it is extrapolated with φ frozen and H=H_vac, which is a stated extrapolation. Thermalization is not addressed.

## 3. Results

Status labels follow the checkpoint convention.

| # | Result | Status | Evidence |
|---|---|---|---|
| R1 | Flat linear-crossing calibration reproduces n_k=exp(−π(κ²+μ₀²)/q) to ≤8.4×10⁻⁷ relative where n_k>10⁻⁶ (absolute ≤3×10⁻⁹), and N to ≤2.1×10⁻⁹, for q=10, 100, 1160 with μ₀=0 and μ₀=√q/2. The wrong formulas exp(−πκ²/2q) and exp(−2πκ²/q) miss N by −65% and +183% (μ₀=0). | numerical (control) | `CONTROLS.json` → `calibration_flat_linear`, `calibration_pass: true` |
| R2 | On the archived trajectory the produced number obeys N=q^{3/2}/(8π³)(1+D/q+O(q⁻²)) with D=10.648, 7.0034, 3.0032, 1.3416 at φ*=0.25, 0.5, 0.75, 0.9. G·(relative deviation) is constant to 4 digits for G=10³…10⁶ (log–log slope −1.0005, −1.0001). Consequently the instant formula is accurate to 1.2% (φ*=0.5) and 0.5% (φ*=0.75) at G=10³. It is off by 13% and 5% at G=100. At G=10–30 the deviation reaches 5% to a factor of about 10, and the instant approximation is not usable there. | numerical | `CONTROLS.json` → `asymptotic_scaling`; `PREHEATING_RESULTS.json` → `primary.*.analytic` |
| R3 | Closed form of the O(1/q) coefficient for a curved, expanding crossing (see §5): D = (9π/4+15/(2π))h² + (3π/2−15/(8π))ḣ + (45/(32π)−π/8)(φ̈/φ̇)² + (π/8−45/(96π))φ⃛/φ̇ + (45/(8π))hφ̈/φ̇. It matches 11 single- and two-parameter toy backgrounds to ≤3.0×10⁻⁵ absolute, which is the measured numerical floor. It predicts the trajectory D to 9×10⁻⁵ (φ*=0.5) and 2×10⁻⁴ (φ*=0.75) relative, and to ≤3×10⁻⁵ at 0.25 and 0.9. The exponent-only (turning-point) coefficients fail the toys by 0.39. | numerical; closed forms of c_AA, c_B, c_hA and part of c_hh are identified from fits, not derived | `CORRECTION_FIT.json` |
| R4 | Produced energy relative to the residual-vacuum threshold at the M₅ cutoff (λ_c=1): R ≤ 1.7×10⁻⁴ for all screen-passing couplings (maximum at φ*=0.5, G=300) and ≤ 2.8×10⁻³ over the whole scan (G=10, outside validity). The "post" variant gives ≤ 4.0×10⁻³ for screen-passing cases (φ*=0.9, G=3000). Against the crossing-time vacuum the ratios are ≤ 4.4×10⁻⁴ for screen-passing couplings and ≤ 2.1×10⁻³ overall. Relative to the tension, κ₅²ρ/σ ≤ 6.6×10⁻⁶, so the low-energy (ρ≪λ) regime holds. | conditional (fixed trajectory, M₅ cutoff, instantaneous full conversion) | `SUMMARY.json`, table below |
| R5 | Large-G law: R_cut ≃ K/√G with K = \|φ̇*\|^{3/2}\|φ_e−φ*\|/(8π³a_e³Δ_max³ ĉ). K = 1.44×10⁻⁴, 2.79×10⁻³, 1.13×10⁻³, 2.56×10⁻⁴ for φ*=0.25, 0.5, 0.75, 0.9. The numerics agree to ≤3×10⁻⁵ at G=10⁶, with an O(1/q) approach consistent with R2. Stronger coupling makes the channel *less* able to dominate once the cutoff is imposed. | conditional; exact within the instant limit | `SUMMARY.json` → `rows[*].K_closed_form`, `closed_form_rel_dev` |
| R6 | Mass scale needed for R=1: λ_c ≥ 18.1 (full) or ≥ 6.3 (post) over screen-passing cases. One radiation-dominated e-fold needs a further factor e^{4/3}=3.8 (for example λ_c≥68.8 at φ*=0.5, G=300). | conditional | table below |
| R7 | One-shot backreaction at the cutoff: R_j ≤ 2.5×10⁻⁵, so the trajectory is not altered at this order. If instead b is raised until R=1, then R_j = 0.10–0.79, largest for late crossings because σ′→δc near φ=1. Any population large enough to matter back-reacts at O(10–80%) on the scalar junction, and the fixed-trajectory result is not self-consistent there. | conditional | `SUMMARY.json` → `max_Rj_at_cutoff_full_screen_passing`, `range_Rj_at_equality_b_screen_passing` |
| R8 | Vacuum polarization: the unsubtracted one-loop tension shift from χ, \|j_CW\|≈G²\|Δφ\|m²/(16π²) without the log, is 1.1–3.7 times \|σ′\| at the data end at the cutoff for φ*=0.5, G=300–1000. It grows ∝G (3.7×10³ at G=10⁶). Keeping the registered σ(φ) therefore requires tuned φ-dependent counterterms for most of the screened range. | conditional (order-of-magnitude, renormalization scale unspecified) | `PREHEATING_RESULTS.json` → `post.cw_ratio_to_sigma_p_at_cutoff_full_end` |
| R9 | Decay channel: for y=0.01, 0.1, 1 the maximum radiation/vacuum ratio at the cutoff is ≤1.3×10⁻⁴ over screen-passing cases. By the data end y=1 converts 20–54% of the χ energy at φ*=0.5 and 0.75 with G=300–1000, and 9% at φ*=0.25 with G=10³. Decay during the nonadiabatic interval is ~y²/(8π)≤0.04. | conditional | `post.decay` |
| R10 | Reliability: tol 3×10⁻⁴, rtol 10⁻¹², 80 nodes and a wider κ range change N by ≤1.4×10⁻⁷. Order-0 versus order-2 basis and cubic versus quintic spline change it by ≤2.2×10⁻⁴ at G=100 and ≤2.2×10⁻⁶ at G≥10³. Continuing to the data end changes it by ≤6×10⁻⁹. The Wronskian is conserved to ≤1.5×10⁻⁸ for G≥100 (1.3×10⁻⁷ at G=10). The three archived grids agree to ≤2.0×10⁻⁴ in N for G≥10³ (1.2×10⁻³ at φ*=0.75, G=100). φ*=0.49/0.51 gives smooth changes, and μ₀=√q/2 suppresses N by 0.45663 against the instant 0.45594. The Ward identity ρ̇+3h(ρ+p)=jv holds to ≤8.4×10⁻⁵ relative (finite differences). | numerical (convergence evidence, not certified bounds) | `CONTROLS.json`, `PREHEATING_RESULTS.json → grid_spread` |
| R11 | Minus branch: the crossing at φ*=−0.5 (s*=5.274, φ̇=−0.985) produces N matching the instant formula to 14.7%, 1.41% and 0.14% at G=10², 10³, 10⁴. No radiation-era test is possible: that trajectory reverses and collapses at s≈5.89. | numerical; no conclusion on radiation | `minus_branch` |

### Energy, cutoff and backreaction table (primary run, μ₀=0)

The "screen" column marks whether G ≥ the checkpoint's local-eligibility minimum from the finer far grid: 397.6, 136.2, 110.7 and 1488.5 at φ*=0.25, 0.5, 0.75 and 0.9. Column definitions:

- **Deviation from instant N:** relative deviation of the solved N from the instant formula.
- **ρ̂(end):** χ energy density at s=6.92 in H₀⁴ units.
- **R at M₅ cutoff:** R at λ_c=1, "full" variant.
- **λ_c for R=1 / N=1:** mass-to-M₅ factor needed.
- **R_j at cutoff / if R=1:** scalar-junction backreaction.
- **decay y=1:** maximum radiation/vacuum ratio at the cutoff.

| φ* | G | q/H₀² | screen | deviation from instant N | ρ̂(end) | R at M₅ cutoff | λ_c for R=1 | λ_c for N=1 | R_j at cutoff | R_j if R=1 | decay y=1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 100 | 36 | fail | 0.33 | 0.995 | 1.9e-05 | 37.3 | 142 | 1.9e-06 | 0.10 | 5.9e-06 |
| 0.25 | 1e3 | 360 | pass | 0.030 | 244 | 4.7e-06 | 59.7 | 226 | 4.7e-07 | 0.10 | 5.2e-06 |
| 0.25 | 1e6 | 3.6e5 | pass | 3e-05 | 7.48e9 | 1.4e-07 | 191 | 723 | 1.5e-08 | 0.10 | 1.7e-06 |
| 0.5 | 10 | 5.79 | fail | 2.0 | 0.0435 | 2.8e-03 | 7.12 | 27 | 3.9e-04 | 0.14 | 9.0e-05 |
| 0.5 | 100 | 57.9 | fail | 0.13 | 4.94 | 3.1e-04 | 14.7 | 55.8 | 4.8e-05 | 0.15 | 7.4e-05 |
| 0.5 | 300 | 174 | pass | 0.041 | 71.1 | 1.7e-04 | 18.1 | 68.8 | 2.5e-05 | 0.15 | 7.5e-05 |
| 0.5 | 1e3 | 579 | pass | 0.012 | 1.40e3 | 8.9e-05 | 22.4 | 84.9 | 1.4e-05 | 0.15 | 6.2e-05 |
| 0.5 | 1e4 | 5.79e3 | pass | 0.0012 | 4.39e5 | 2.8e-05 | 33.0 | 125 | 4.2e-06 | 0.15 | 1.1e-04 |
| 0.5 | 1e6 | 5.79e5 | pass | 1.2e-05 | 4.38e10 | 2.8e-06 | 71.0 | 269 | 4.2e-07 | 0.15 | 1.3e-05 |
| 0.75 | 100 | 56.4 | fail | 0.054 | 6.36 | 1.2e-04 | 20.3 | 76.9 | 3.6e-05 | 0.30 | 1.6e-05 |
| 0.75 | 300 | 169 | pass | 0.018 | 95.4 | 6.7e-05 | 24.7 | 93.5 | 2.0e-05 | 0.30 | 2.0e-05 |
| 0.75 | 1e3 | 564 | pass | 0.0053 | 1.91e3 | 3.6e-05 | 30.3 | 115 | 1.1e-05 | 0.31 | 2.0e-05 |
| 0.75 | 1e4 | 5.64e3 | pass | 5.3e-04 | 6.01e5 | 1.1e-05 | 44.5 | 169 | 3.5e-06 | 0.31 | 2.2e-05 |
| 0.75 | 1e6 | 5.64e5 | pass | 5.3e-06 | 6.00e10 | 1.1e-06 | 95.9 | 364 | 3.5e-07 | 0.31 | 2.4e-06 |
| 0.9 | 3e3 | 1.05e3 | pass | 0.0013 | 1.16e4 | 4.7e-06 | 59.7 | 227 | 3.7e-06 | 0.79 | 2.8e-06 |
| 0.9 | 1e6 | 3.51e5 | pass | 3.8e-06 | 2.35e10 | 2.6e-07 | 157 | 597 | 2.0e-07 | 0.79 | 3.0e-07 |

The full 36-row table is printed by `build_summary.py` and stored in `SUMMARY.json`. The decay column can exceed "R" because it takes the maximum over time of the radiation density, whereas R uses the χ energy at the data end.

## 4. Conditional conclusion

The following holds within the fixed archived +1 trajectory, the M₅ effective-theory condition on χ masses, the adiabatic vacuum before the crossing, perturbative decay with y≤1, no rescattering or thermalization, and no bulk feedback.

Particle production at the favoured mid-roll crossings (φ*≈0.5, 0.75) is real and is computed reliably for G≳100. It is accurately instant-preheating-like for G≳10³. However, the energy it can deliver is 10⁻⁶–10⁻⁴ of what is needed merely to equal the residual +1-branch vacuum, and a radiation era needs at least e⁴ times that threshold. **In no scanned coupling range does this channel plausibly give a radiation era in this model.**

The only ways around this result lie outside the present assumptions:

- χ masses of tens of M₅ or more;
- a different energy source that depletes the residual vacuum;
- a population so large that its O(10–80%) backreaction requires solving the coupled bulk problem.

The coupled bulk problem is precisely the unsolved step. The result is consistent with, and independent of, the checkpoint's toy 4D budget (N≤0.135 e-folds), which it does not rely on.

**Not solved here:** the modified metric junctions nA=(σ+κ₅²ρ)/6 and nB=(σ−κ₅²(2ρ+3p))/6, the modified scalar junction nφ=−(σ′+κ₅²j)/2, the Weyl term and bulk energy exchange, and hence any self-consistent change of the trajectory. Also not solved: rescattering, nonlinear χ dynamics, fermion back-reaction and Pauli blocking, thermalization, and the actual physical value of M₅ or H₀.

## 5. The O(1/q) correction (R3), status and derivation notes

- The h² and ḣ terms have two parts. The effective mass shift −(9/4)h²−(3/2)ḣ in Ω² gives 9π/4 and 3π/2 exactly, as in the μ₀ term of the instant formula. A complex-turning-point (Dykhne–Davis–Pechukas-type) expansion of the exponent for the redshifting κ²/a² term gives +15/(4π)h² and −15/(8π)ḣ; this project derived it in the course of this work. The remaining +15/(4π)h² was identified numerically. For the method see, for example, https://arxiv.org/html/1402.5669.
- The flat-space φ̈ and φ⃛ terms from the exponent alone are 45/(32π) and −45/(96π). The measured coefficients differ from these by exactly ∓π/8 to fit precision, so the combination (π/8)·d(φ̈/φ̇)/ds is a non-exponent correction. It was identified numerically and is not derived.
- The formula is checked only for the stated five regressors, at q≈2×10³–10⁴ with Richardson extrapolation. It is a result for this project. No literature search was completed to establish whether it is known elsewhere, so no novelty is claimed.

## 6. Literature consulted (search snippets only; paper sites are blocked for download)

- Felder, Kofman and Linde, instant preheating: https://arxiv.org/abs/hep-ph/9812289 (cited through the checkpoint).
- Instant preheating in brane-world and D-brane settings: https://arxiv.org/html/hep-th/0405016v2, https://arxiv.org/html/0905.2284 and https://arxiv.org/html/hep-ph/0205240. These are relevant precedents for brane preheating. None uses this model's junction-coupled scalar, so none fixes G, M₅ or φ* here.
- Parker and Fulling, adiabatic regularization: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.9.341.

## 7. Limitations

- **Trajectory.** The trajectory is the archived pure-tension run, which carries its own known initial junction mismatch and constraint limitations (see the checkpoint audit). Data end at s=6.92, before the trajectory has settled on the +1 branch (h=0.64 there against 0.607). The extrapolation beyond that point freezes φ and H.
- **Cutoff.** The EFT cutoff choice (M₅, factor λ_c) is an assumption. All ratios scale as λ_c³, so the conclusion reverses only for λ_c≳6–20 at the most favourable points.
- **Particle description.** Densities use the late-time particle description. The renormalized stress *during* the nonadiabatic interval is not computed; only the instantaneous zeroth-order-basis occupation is stored for illustration (`primary["0.5|10000"].history`). Coleman–Weinberg terms are assumed absorbed, and R8 shows that this requires tuning.
- **Low couplings.** For G≲30 the particle interpretation is marginal: the order-2 basis sometimes fails at the data end and the code falls back to order 0 (flagged as `order_used`). These rows are shown but marked as failing the screen.
- **Decay.** The rate is the perturbative rest-frame rate with time dilation. Production and decay are treated as sequential.
- **Precision.** All numbers are double-precision numerics with the convergence evidence above. None is a certified error bound.

## 8. Reproduction

Requirements: python3 with numpy ≥2, scipy ≥1.16. Everything is single-threaded and runs in under 15 minutes in total on one core.

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/preheating
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 run_preheating.py      > run_preheating.log      # ~7 min  -> PREHEATING_RESULTS.json
python3 run_controls.py        > run_controls.log        # ~2 min  -> CONTROLS.json
python3 run_correction_fit.py  > run_correction_fit.log  # ~1 min  -> CORRECTION_FIT.json (reads CONTROLS.json)
python3 build_summary.py       > build_summary.log       # seconds -> SUMMARY.json and the tables
```

Files:

- `ph_lib.py`: trajectory loader, mode solver, window finder, densities.
- `run_preheating.py`: production scan, energy, backreaction and decay.
- `run_controls.py`: calibration, wrong-formula, convergence and perturbation controls.
- `run_correction_fit.py`: O(1/q) correction toys and closed forms.
- `build_summary.py`: tables and the K/√G law check.
