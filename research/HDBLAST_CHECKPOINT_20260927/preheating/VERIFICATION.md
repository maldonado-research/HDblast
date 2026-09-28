# Independent audit of the preheating workstream

Audit date: 28 September 2026. The audited files were not modified: their SHA-256 hashes before the audit are in `verify/original_sha256.txt`, and they were re-checked after the audit. Every audit script and output is in `verify/`.

## Overall verdict

The main numerical content holds up.

- **Reproduction.** Re-running all four scripts in a copy (`verify/rerun/`) reproduces `PREHEATING_RESULTS.json`, `CONTROLS.json`, `CORRECTION_FIT.json` and `SUMMARY.json` bit for bit: 18223, 500, 272 and 1076 numbers, with zero differences.
- **Independent solver.** A solver written separately for this audit agrees with the audited N on the archived trajectory to at most 2.2×10⁻⁸ relative. It also confirms the O(1/q) closed form on toy backgrounds the workstream never fitted, including one with all five structures switched on at once.
- **Exact algebra.** Every formula the conclusions rest on checks exactly in sympy: the mode equation, the instant-preheating number, the vacuum and deceleration thresholds, the Ward identity, the K/√G law, the Coleman–Weinberg source and σ′.
- **Headline.** The conditional negative headline is supported within its stated assumptions.

The audit found four problems:

1. **Stale decay numbers in the README.** The decay column of the README table and claim R9 do not match the saved JSON. 15 of the 16 tabulated decay entries differ by more than 10%. The screen-passing maximum is 2.34×10⁻⁴, not "≤1.3×10⁻⁴".
2. **Wrong derivation note (README §5).** A numerical complex-turning-point calculation gives an exponent-only h₀² coefficient of +15/(2π), not +15/(4π). The whole c_hh then follows from the mass shift plus the exponent. The c_hA coefficient 45/(8π) is also exactly the exponent value. Only the ±π/8 pieces of c_AA and c_B remain unexplained. The numbers are unaffected; the attribution is wrong.
3. **Headline depends on when the energy is evaluated.** R uses the χ energy at the end of the data. Taking the maximum over the saved post-window profile at the same cutoff gives up to 5.9×10⁻⁴ (φ*=0.5, G=300), 3.5 times the headline 1.7×10⁻⁴. The conclusion (≪1) is unchanged, but "at most 1.7×10⁻⁴" holds only for the data-end definition.
4. **Minor overstatements.**
   - "G·(relative deviation) constant to 4 digits": the actual spread is 0.39% (φ*=0.5) and 0.10% (φ*=0.75) over G=10³…10⁶, which is 2–3 digits. The drift is the expected O(1/q²) term.
   - "D predicted to ≤3×10⁻⁵ at φ*=0.9": the best agreement is 7.4×10⁻⁵ (G=10⁵) and 4.2×10⁻⁴ at G=10⁶.
   - The two "wrong-formula" calibration controls are fixed by Gaussian algebra (N ratios 2^{±3/2}). They test the comparison metric, not the solver.

## Verdict table

| # | Claim (short) | Verdict | Check performed |
|---|---|---|---|
| R1 | Flat linear-crossing calibration; wrong formulas fail by −65%/+183% | **confirmed** (weak control) | Re-run: identical. The audit's own solver also reproduces n_k=exp(−π(κ²+μ₀²)/q) to 1.1×10⁻⁹ and N to 1.1×10⁻⁹. The wrong-formula deviations equal 2^{−3/2}−1 and 2^{3/2}−1 exactly, so they carry no information about the solver. |
| R2 | N=q^{3/2}/(8π³)(1+D/q); D values; instant formula accuracy vs G | **partially confirmed** | Numbers re-derived from the JSON (`CHECK_CLAIMS.json`). The accuracies (1.2% and 0.5% at G=10³; 13% and 5% at G=100) are correct. D at φ*=0.25 and 0.9 is quoted from the closed form, not measured. The measured values are 10.6476 (G=10⁶) and 1.3415 (G=10⁵), consistent with it. "Constant to 4 digits" is really 2–3 digits over G=10³–10⁶ (0.39% and 0.10% spread). |
| R3 | Closed-form O(1/q) coefficient D | **confirmed numerically; derivation note refuted in part** | Independent solver, Richardson in q: in-sample toys agree to ≤2×10⁻⁵. Out-of-sample toys with all five terms at once and v=0.4 or 0.8 agree to 1.0×10⁻⁵ and 2.2×10⁻⁵. A trajectory-like toy agrees to 1.5×10⁻⁵ (`INDEP_MODES.json`, `INDEP_MIX3.json`). Dimension counting shows the five regressors are the complete set at O(1/q). The mpmath turning-point evaluation (`DDP_EXPONENT.json`) confirms the exponent-only 45/(32π), −45/(96π) and −15/(8π) to 10⁻¹¹. It gives **15/(2π)** h₀² (README: 15/(4π)) and 45/(8π) for the h₀A cross term. So only the π/8·d(φ̈/φ̇)/ds piece is non-exponent. No novelty is claimed, and none is supported. |
| R4 | R ≤1.7×10⁻⁴ (screen-passing) and ≤2.8×10⁻³ (all) at the M₅ cutoff; post ≤4.0×10⁻³; crossing ≤4.4×10⁻⁴ and ≤2.1×10⁻³; κ₅²ρ/σ ≤6.6×10⁻⁶ | **confirmed as defined (data-end energy)**, with caveat | All values recomputed from the JSON. ρ_crit is recomputed from the +1 branch (σ_f/H₀=52.678, h_vac=0.60672, crit=0.125633, exact match). The threshold algebra holds in sympy. Caveat: the maximum over the stored profile is 5.9×10⁻⁴. The result is conditional on the M₅ cutoff, on b being a free parameter, and on the fixed trajectory. |
| R5 | R_cut ≃ K/√G and the K values | **confirmed** | Law exact in sympy. K recomputed: 1.444×10⁻⁴, 2.791×10⁻³, 1.133×10⁻³, 2.562×10⁻⁴. Numerics within ≤3.0×10⁻⁵ at G=10⁶. |
| R6 | λ_c ≥18.1 (full) and ≥6.3 (post); ×e^{4/3} for one e-fold | **confirmed** | Recomputed: 18.130 at (0.5, 300) and 6.275 at (0.9, 3000); 68.8 for one e-fold. |
| R7 | R_j ≤2.5×10⁻⁵ at the cutoff; 0.10–0.79 if R=1 | **confirmed** (conditional) | Recomputed from the JSON: 2.54×10⁻⁵, and 0.1005–0.7917. The junction nφ=−(σ′+κ₅²j)/2 matches the matter report (§2). σ′ is checked in sympy. |
| R8 | CW shift 1.1–3.7×\|σ′\| at φ*=0.5, G=300–1000; 3.7×10³ at G=10⁶ | **confirmed** (order of magnitude) | Values 1.095, 3.651 and 3651. dV_CW/dφ without the log equals G²Δm²/(16π²) exactly (sympy). 15 of 21 screen-passing rows exceed 1, which supports "most of the screened range". |
| R9 | Decay: max radiation/vacuum ≤1.3×10⁻⁴ (screen-passing); y=1 converts 20–54%; ~y²/(8π) | **partially confirmed; the headline number is refuted** | JSON maximum = **2.34×10⁻⁴** at (0.5, 300, y=1). The README table's decay column mismatches the JSON in 15 of 16 rows (for example 0.5\|300: README 7.5×10⁻⁵, JSON 2.34×10⁻⁴). The re-run is bit-identical, so the README text is stale. The conversion fractions 0.20–0.54, 0.087 at φ*=0.25, and Γ=y²m/(8π) are correct. The conclusion (≪1) is unchanged. |
| R10 | Convergence and reliability numbers | **confirmed** | All quoted bounds match `CONTROLS.json` and the grid spread. The one exception is minor: the maximum Wronskian 1.3×10⁻⁷ occurs at (0.75, G=30), not G=10. The audit's independent solver (different basis construction, windows, quadrature and complex state) matches N on the trajectory to 1.4×10⁻⁹–2.2×10⁻⁸. These remain double-precision convergence evidence, not certified bounds. |
| R11 | Minus branch: φ*=−0.5 at s*=5.274; 14.7%, 1.41%, 0.14%; reversal at s≈5.89 | **confirmed** | JSON values match. The archive shows h first negative at s=5.901, and the cut at h<0.3 falls at s=5.839 (`MINUS_BRANCH_CHECK.json`). At G=100 the right-edge adiabaticity is 0.04, so the 14.7% figure is marginal. |
| — | Mode equation Ω²=κ²/a²+m²−(9/4)h²−(3/2)h′ | **confirmed (exact)** | sympy residual 0. Dropping the h′ term leaves a nonzero residual (control). |
| — | Ward identity ρ̇+3h(ρ+p)=jv for the particle description | **confirmed (exact)** | sympy residual 0, per mode. |
| — | Literature: FKL hep-ph/9812289 and Parker–Fulling PRD 9, 341 | **confirmed as standard references** | Their formulas are consistent with the code (sympy check of N). |
| — | Brane-preheating URLs hep-th/0405016, 0905.2284, hep-ph/0205240, and 1402.5669 as a turning-point method reference | **unverifiable** | The web-search budget was exhausted during this audit, so their content could not be confirmed. |

## Checks run (in `verify/`)

| Script | Output | What it does |
|---|---|---|
| `rerun/*.py` (verbatim copies) | `rerun/*.json`, `rerun/*.log` | Full re-run of the four audited scripts (about 13 min on one core) |
| `compare_rerun.py` | `RERUN_COMPARISON.json` | Compares every number with the audited JSON: all equal |
| `check_claims.py` | `CHECK_CLAIMS.json` | Recomputes every README aggregate from the JSON and the checkpoint inputs. Also covers the README decay column, the max over the profile and the D values. |
| `sympy_checks.py` | `SYMPY_CHECKS.json` | Exact algebra (see the table) |
| `indep_modes.py` | `INDEP_MODES.json` | Independent solver: calibration, in-sample and out-of-sample toys, archived trajectory |
| `indep_mix3.py` | `INDEP_MIX3.json` | Repeats the trajectory-like toy. With span 1.4, the audit's toy had a second mass zero at s=−1.46; that was a defect of the audit toy, fixed with span 1.0 and 0.8. |
| `ddp_exponent.py` | `DDP_EXPONENT.json` | mpmath complex-turning-point exponent: coefficients of the exponent-only terms |
| `minus_branch_check.py` | `MINUS_BRANCH_CHECK.json` | Reversal time on the minus-branch archive |

Tolerances: the independent solver uses DOP853 at rtol 10⁻¹¹, a 64-node Gauss–Legendre grid, adiabatic windows at ε=10⁻⁴ (toys) and 10⁻³ (trajectory), and Richardson extrapolation over two or three q values. Its flat-toy floor for D is 1.3×10⁻⁵, and changing ε shifts D by 4×10⁻⁶. The turning-point calculation uses mpmath at 30 digits, and its flat control gives D=−2.5×10⁻¹².

Reproduce:

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/preheating/verify
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
(cd rerun && python3 run_preheating.py > run_preheating.log && python3 run_controls.py > run_controls.log \
  && python3 run_correction_fit.py > run_correction_fit.log && python3 build_summary.py > build_summary.log)
python3 compare_rerun.py; python3 check_claims.py; python3 sympy_checks.py
python3 indep_modes.py; python3 indep_mix3.py; python3 ddp_exponent.py; python3 minus_branch_check.py
```

## Scope remarks

These are not errors, but readers should keep them in mind.

- **Conditional result.** The whole negative result is conditional on three things: the fixed archived trajectory, taking M₅ as the cutoff for χ masses, and b=κ₅²H₀³ being a free parameter bounded only by that cutoff. A lower brane cutoff would strengthen the negative. The workstream states the λ_c³ scaling honestly.
- **Where the calculation stops being self-consistent.** Wherever the channel could matter (R≈1), the workstream's own numbers show 10–80% backreaction and tension shifts larger than σ′. The unsolved coupled bulk problem is therefore the real open question, as the README says.
- **Status labels.** The labels are appropriate. No discovery or novelty claim was found.

## Suggested corrections to the README

1. R9 and the table's decay column: replace them with the JSON values. The screen-passing maximum is 2.34×10⁻⁴ at (φ*=0.5, G=300, y=1).
2. §5 and R3: the exponent-only h₀² term is 15/(2π), as the numerical turning-point calculation shows. c_hh and c_hA then follow from the mass shift plus the exponent; only the ±π/8 terms are fitted.
3. R2: change "constant to 4 digits" to "constant to 0.4% (φ*=0.5) and 0.1% (φ*=0.75)". State that D at φ*=0.25 and 0.9 is the closed-form value, and give the measured values.
4. R3: change the φ*=0.9 agreement to 7×10⁻⁵ (G=10⁵).
5. R4: add that R is evaluated with the data-end energy, and that the maximum over the post-window interval is 5.9×10⁻⁴.
