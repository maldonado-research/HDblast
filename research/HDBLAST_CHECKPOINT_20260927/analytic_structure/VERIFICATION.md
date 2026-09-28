# Independent audit of the `analytic_structure` workstream

Audit date: 28 September 2026. The auditor wrote new code in `verify/` and re-ran the audited scripts on copies in `verify/rerun/`. No audited result file was edited, with one accidental exception: `MANIFEST.sha256.json` was overwritten. See "Process incident" at the end.

**Overall verdict.** The mathematical content is sound and reproducible. Every exact coefficient through O(δ⁸) was re-derived with independent code and matches exactly. A second, independent multiprecision solver agrees with the audited boundary-value solutions to 10⁻³⁹–10⁻⁴¹ relative. The Gegenbauer and leading-order spectrum statements hold under direct ODE tests.

**One claim is refuted.** The workstream says the δ⁹ coefficient of the shell observables "depends on the cone" and "would need global second-order matching". That is wrong. The δ⁹ coefficient is fixed by the local shell expansion. The auditor computed it exactly in c. It does not depend on the free resonance constant ν and has no ln δ. It matches the workstream's own fitted values and the auditor's independent solver.

Minor issues:
- a few quoted ranges are slightly off;
- one control is weak;
- three of the eight literature links could not be traced to a search record.

## Verdict table

| # | Claim (audited) | Verdict | Evidence (this audit) |
|---|---|---|---|
| 1 | The O(δ) and O(δ²) terms of the 22 Sept package are reproduced exactly | **confirmed** | Independent series code (`verify/v1_series_independent.py`) gives −9c/64, (1+c)/27 and (29c²+64c+32)/1152 = (1+c)²/36 − c²/384. Consistent with the package's own float remainders in `static_branch/INDEPENDENT_BRANCH_CHECKS.json`: −0.0158935 and −0.0023243 at δ = 3×10⁻⁴. |
| 2 | New exact coefficients through δ⁸ (φ_b, H², X_b, T_b; ρ_b, y_b, η_h derived from them) | **confirmed** | The auditor's code was written from scratch. It re-derives V and F from W, uses Fractions and a direct division by the linear operators, and does not use `sympy.solve`. All 9×4 coefficients of η_b, H², X_b and T_b match symbolically. Re-running `series_expansion.py` in a copy gives a JSON identical to the original apart from runtime. The first ρ_b and y_b corrections were checked by hand from h₂/h₁ and t₂/t₁. |
| 3 | 30–40-digit BVP solutions confirm the coefficients | **confirmed** | Independent solver `verify/v2_bvp_independent.py`: `mpmath.odefun` with the first integral *enforced*; the audited code uses the second-order R equation and only monitors that integral. Relative agreement with the audited BVP: η_b 2.4×10⁻⁴¹, H² 7.4×10⁻⁴¹, ρ_b 3.0×10⁻⁴¹ at δ = 10⁻³; 1.5×10⁻³⁹ to 8.3×10⁻³⁹ at δ = 10⁻². At δ = 10⁻³: φ_b = 0.99991594731691335387578407751, H² = 5.9240147943288795996237552780×10⁻⁵, the same digits as the audit. Own convergence (dps 40/u₀ = 10⁻³ vs dps 50/u₀ = 5×10⁻⁴): ≤ 1.3×10⁻⁴⁰ for η_b and H². Known limit c = 0, δ = 0.1: H² − (δ/27 + δ²/36) = 4.8×10⁻⁴³. |
| 3a | Quoted ranges ("q₃ rel. err 1e-26…1e-24", "c-perturbation 10⁷–10⁹") | **partially confirmed** (cosmetic) | From `SERIES_VS_BVP.json`: q₃ errors reach 2.2×10⁻³⁰ (c = −0.4, H²), better than quoted. The perturbed-c deviation reaches 6.0×10⁹ (c = −0.4), above the README's "10⁷–10⁹"; the summary's "1e7 to 6e9" is correct. The rest matches the JSON: solver agreement 3.6×10⁻⁴² / 1.8×10⁻³¹, max J2 3.6×10⁻⁵⁰, max first-integral residual 2.8×10⁻³² absolute. |
| 3b | Scaled remainders deviate ∝ δ at every order | **confirmed** | The ratio dev(0.0016)/dev(0.0008) lies between 1.990 and 2.002 for all five quantities, all orders and all three nonzero c (`V4_COMPARE.json`). |
| 4 | Controls fail as they should | **confirmed**, with a note | Values match the JSON: 7.8 / 152 / 5.1 for the wrong h₃, against 1×10⁻⁵ to 3×10⁻⁵ when correct. Note: the "c → c(1+10⁻⁶)" control only shows sensitivity of a δ⁴-scaled remainder to an O(δ) shift, so it is weak as a discriminator. Auditor's own controls: dropping the secular term leaves a nonzero residual of −24576/17 at X²T⁷. A wrong Legendre degree misses the direct ODE by ≥ 4.8×10⁻³. f₁₅ does not solve the m²/k² = 252 equation. |
| 5a | Resonance at (m,j) = (2,7) with S = −24576/17 and κ = 768/17 | **confirmed** (exact) | Auditor's independent recursion gives the same S and κ. With the secular term, all residuals vanish through degree 9. |
| 5b | No ln δ in η_b or H² at δ⁹ | **confirmed**, now exact rather than fitted | Junctions solved with u_b = U and ν kept symbolic: the δ⁹ coefficients of η_b and H² are independent of U, so no ln δ, and independent of ν. X_b at δ⁹ does depend on both, as expected. |
| 5c | "The δ⁹ coefficient depends on cone regularity / is not given by the local recursion; exact value would need global matching" | **refuted** | The free constant ν at the resonance only relabels the growing-mode amplitude, α → α + να², at this order, so it drops out of every shell observable. The δ⁹ coefficient is therefore local, and the auditor obtained it exactly in c (`V1_SERIES_INDEPENDENT.json`, keys `eta_b_delta9_exact`, `H2_delta9_exact`). Registered c: η_b: 3.249147254269×10⁻⁶; H²: −2.759694789×10⁻⁷. The audit's fits gave 3.24914729×10⁻⁶ and −2.7596954×10⁻⁷ (relative 1.2×10⁻⁸ and 2×10⁻⁷). c = 1.3: −2.874978716×10⁻⁴ vs fit −2.874978716×10⁻⁴. Independent BVP: after subtracting the series through δ⁹, the remainder divided by δ¹⁰ is constant, 1.80–1.84×10⁻⁶ for η_b and 1.21–1.23×10⁻⁷ for H² over δ = 10⁻³ to 0.04 (and 1.17×10⁻⁴ / 2.5×10⁻⁵ for c = 1.3). So the exact δ⁹ term is confirmed and the next term is O(δ¹⁰). The README's statement "cone regularity drops out through δ⁸" is true but understated: it drops out at least through δ⁹ for the shell observables. The numbers the workstream reported are correct; only the interpretation is wrong. |
| 5d | η_h acquires a δ¹⁶ ln δ term with B = (81/4)c²[3(1+c)/4]¹⁴ | **confirmed** (algebra) / numerical (fit) | By hand: B = (136/3)(κ/2)X₁²T₁¹⁴ = 1024 X₁²T₁¹⁴ = (81/4)c²[3(1+c)/4]¹⁴. The fit agreement of 2×10⁻⁶ to 5×10⁻⁵ reproduces on re-running `resonance_probe.py` (identical JSON). Unlike (5c), this log is real because η_h measures α itself. The README's wording "δ⁸ · … · δ⁸ ln δ" is awkward but means δ¹⁶ ln δ. |
| 6 | η_h = (136/3)X_bT_b⁷, i.e. −(111537/131072)c(1+c)⁷δ⁸ | **confirmed** | By hand: (136/3)(−9c/64)(3(1+c)/4)⁷ = −111537/131072 c(1+c)⁷. Independent BVP η_h matches the audited 8-term series to a relative 5×10⁻²⁰ at δ = 10⁻³ and 5×10⁻¹² at δ = 10⁻². That is consistent with the relative δ⁸ ln δ correction. |
| 7 | Gegenbauer origin: n(n+4) = 252, Δ = 3W''/W = 18, closed forms, elementary f_ν, no polynomial at φ = −1 | **confirmed** | The audited script re-ran with 28/28 checks and an identical JSON. The auditor's independent sympy checks (`verify/v3_gegenbauer_spectrum.py`) cover: the Σ(j+1)(15−j)e^{(14−2j)u} form; C₁₄⁽²⁾ = ½U′₁₅ with U₁₅(cosh u) = sinh16u/sinh u; the explicit polynomial; C(1) = 680; f_ν and its singular partner (≈ −u⁻³) solving the equation for generic ν; f₁₆ = 2C₁₄⁽²⁾; and b = 18 for the Robin coefficient. Wording nuance: at φ = −1, s = 3W''/W = −18/5 is the *other* root, 4 − Δ, not Δ = 38/5. The README's numbers are right, but the phrase "s … is the conformal dimension Δ" holds only at φ = +1. The claim d ln C/du > 0 holds because C is a positive combination of cosh terms. |
| 8 | Tensor condition (μ − 3/2)P^{−μ}_{1/2}(cosh u_b) = 0; massless graviton plus continuum from m = 3H/2 | **confirmed** (leading order, conditional on the stated setting) | Direct `mpmath.odefun` shooting of the tensor mode equation vs the closed form f′/f = (μ − 3/2)P_{1/2}/(sinh u · P_{3/2}): max error 2.6×10⁻¹⁶ over 9 (u_b, M²) points. The auditor re-checked the Pfaff positivity argument and it is valid: for μ > ½ all terms after the first are negative and the Gauss value is positive; for μ ≤ ½ all terms are ≥ 0. Independent sweep: μ ≤ 12, x ≤ cosh 12, no nonpositive value. Direct-ODE scan: no sign change of f′(u_b) for 0 < M² < 9/4 at 5 radii. The audited spectrum script re-ran with an identical JSON apart from runtime. |
| 9 | Decoupled bulk scalar: no mode below 9H²/4 for H/k < 16.6; bound state M² ≈ 1.004 at H/k = 50 | **confirmed** (conditional: mixing neglected, as the README states) | Independent float shooting gives critical H/k = 16.647 (audit 16.65) and M² = 1.00404 at H/k = 50 (audit 1.00410; the difference comes from the float Frobenius start). F + b is between 31.37 and 31.42 at u_b = 1 and about 31.995 at u_b = 3.364. For M² < 0, f′/f > 0. |
| 10 | Stability of the branch is open (radion mixing omitted) | **confirmed** (honest caveat) | Nothing here tests the coupled scalar–metric sector. This is correctly flagged in the README. |
| 11 | Literature attributions (8 links) | **partially confirmed / unverifiable** | Five links appear with matching titles in the project's literature-workstream search records (`../literature/*_checks/*_sources.json`): Garriga–Sasaki hep-th/9912118, DeWolfe–Freedman–Gubser–Karch hep-th/9909134, Hawking–Hertog–Reall hep-th/0003052, Frolov–Kofman hep-th/0309002, and Garriga–Vilenkin PRD 44, 1007. Karch–Randall hep-th/0011156, Himemoto–Sasaki gr-qc/0010035 and Langlois–Maartens–Sasaki–Wands hep-th/0012044 have no search record in the project. The auditor's web-search budget was exhausted, so they are **unverifiable** here. The identifiers match the auditor's background knowledge, but that is not verification. The attributions are used only as context; no result depends on them. |
| 12 | No novelty beyond the project is claimed; nothing bears on the Big Bang hypothesis | **confirmed** | The README's scope statements are accurate and appropriately modest. |

## Recommended corrections to the audited README and summary

1. Replace "The δ⁹ coefficient itself depends on cone regularity; … would need global second-order matching. That was not attempted." with: the δ⁹ coefficients are locally determined. Their exact values in c are in `verify/V1_SERIES_INDEPENDENT.json`. Registered c: η_b: 3.249147254269×10⁻⁶; H²: −2.759694789×10⁻⁷. The resonance constant ν drops out of the shell observables and affects only η_h.
2. Correct the perturbed-parameter range in the README table: "about 10⁷–10⁹" should read "1.5×10⁷ to 6×10⁹".
3. At φ = −1, say s = 4 − Δ.
4. Label the Karch–Randall, Himemoto–Sasaki and Langlois–Maartens–Sasaki–Wands links as not re-verified, or add their provenance.

## Checks run (all outputs machine-readable)

| Script | Output | What it does | Runtime |
|---|---|---|---|
| `verify/rerun/*` (copies of audited scripts) | `verify/rerun/*.json` vs `*.orig.json` | Re-ran `series_expansion.py 8`, `gegenbauer_structure.py`, `resonance_probe.py`, `spectrum_leading_order.py` and `analyze_series_vs_bvp.py` on copies. All outputs are identical to the originals except the `runtime_s` fields. The audited mpmath BVP scans were *not* re-run; they are replaced by the independent solver. | about 3 min |
| `verify/v1_series_independent.py` | `verify/V1_SERIES_INDEPENDENT.json` | Independent exact series through degree 9 with the secular term and symbolic ν, U and c. Symbolic comparison with the audited coefficients. Exact δ⁹ coefficients. | 28 s |
| `verify/v2_bvp_independent.py` | `verify/runs/V2_BVP_{reg40,reg50,cp13_40,c0}.json`, `verify/logs/v2_*.log` | Independent shooting with `mpmath.odefun` and the first integral enforced. Registered c at δ = 10⁻³ to 0.04; c = 1.3; c = 0. Two precision settings. | 100–370 s per δ, 1 core |
| `verify/v3_gegenbauer_spectrum.py` | `verify/V3_GEGENBAUER_SPECTRUM.json` | Exact superpotential and Gegenbauer identities; direct-ODE tests of the tensor condition; positivity sweep; scalar threshold and bound state; controls. 28/28 checks passed. | 77 s |
| `verify/v4_compare.py` | `verify/V4_COMPARE.json` | Independent BVP vs audited BVP; series through δ⁸ and δ⁹; η_h series; the audited remainder-scaling claim. | a few seconds |

Precision and tolerances:
- **Series:** exact rationals and exact sympy expressions.
- **Solver:** mpmath dps 40 (u₀ = 10⁻³) and dps 50 (u₀ = 5×10⁻⁴). The secant tolerance is 10^(−2·dps+8). J2 residual ≤ 1.1×10⁻⁴² at dps 40.
- **Where the solver's own error dominates:** in η_h. The two settings differ by 3×10⁻²⁴ at δ = 0.01 because of the approximate cone start. The shell observables are insensitive to this (≤ 1.3×10⁻⁴⁰).
- **Spectrum:** scipy DOP853 with rtol 10⁻¹¹; mpmath at dps 25.

## Reproduction

```sh
cd research/HDBLAST_CHECKPOINT_20260927/analytic_structure/verify
python3 v1_series_independent.py
python3 v2_bvp_independent.py 40 0.001 reg 0.001,0.005,0.01,0.02,0.04 reg40
python3 v2_bvp_independent.py 50 0.0005 reg 0.01,0.04 reg50
python3 v2_bvp_independent.py 40 0.001 1.3 0.005,0.01,0.02 cp13_40
python3 v2_bvp_independent.py 40 0.001 0 0.1 c0
python3 v3_gegenbauer_spectrum.py
python3 v4_compare.py
# re-run of the audited scripts on copies (writes only inside verify/rerun):
cd rerun && python3 series_expansion.py 8 && python3 gegenbauer_structure.py && python3 resonance_probe.py \
  && python3 spectrum_leading_order.py && python3 analyze_series_vs_bvp.py
```

## Limitations of this audit

- No interval certification. Agreement between two independent multiprecision solvers is strong evidence, not proof.
- The claim that ν drops out was shown exactly at degree 9 only. Higher resonances such as (3,14) and the decaying-mode keys (1,16) and (2,23) were not examined.
- The spectrum checks cover the same leading-order, decoupled setting as the audit. They cannot address the omitted radion mixing.
- Three literature links could not be verified (row 11).

## Process incident

While inspecting the audited folder, the auditor ran `python3 make_manifest.py --help`. The script ignores arguments, so it regenerated `MANIFEST.sha256.json` in the audited folder at 03:47:48 UTC. The new manifest has 54 entries: the original 29, plus 25 files under `verify/` that existed at that moment.

The audited files themselves were not changed. Their modification times are all at or before 00:17:48 UTC, and the regenerated hashes are of those unchanged files.

An attempt to restore the manifest was blocked by the permission system, so it has **not** been restored. To restore it, delete the 25 `verify/…` keys from `MANIFEST.sha256.json`. Do not re-run `make_manifest.py`: it would again pick up `verify/`. Removing the keys should give the original content, assuming the original was produced by the same script from the same 29 files. Until then, the manifest no longer describes the folder as the workstream published it.
