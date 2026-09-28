# Literature sweep 2: observations (27 September 2026 round)

**Work in progress.** Nothing here is a claim of the programme until it appears in a dated checkpoint or a Zenodo record. This file is separate from `README.md`, which belongs to sweep 1 (theory), so that two concurrent sweeps do not overwrite each other. Merge it into `README.md` when the round is assembled.

## Question

Which 2023–2026 observations could test or constrain the HDBLAST hypothesis? That is, a five-dimensional gravitational event sparking our Big Bang, as represented by the registered 5D Einstein–scalar shell model. The sweep covers:

- pulsar timing arrays (PTAs), including spectral-turnover searches
- DESI dark energy
- CMB limits on ΔN_eff, n_s and tensors
- Big Bang nucleosynthesis (BBN)
- LIGO/Virgo/KAGRA and GW170817 constraints on extra dimensions
- tabletop tests of the inverse-square law
- JWST early galaxies
- the Hubble tension

For each, it states what the observation would constrain in a braneworld or higher-dimensional-origin model. It also asks whether 2025–2026 releases bear on the project's registered PTA "knee" test.

## Method

- **48 distinct web searches**, listed in `observations.md` §7. Three more were refused when the shared search budget ran out. Paper pages could not be opened (arxiv.org, zenodo.org and journal sites are blocked here), so **every source entry is snippet-only** and labelled so. Where a snippet summary merged several pages, the entry has an attribution note.
- Source records live in `observations_checks/sources_observations.py`. The builder `build_observations_md.py` validates them (https URLs, no duplicates, allowed labels, years consistent with arXiv IDs, query indices). It cross-checks **43 numbers quoted in the prose** against the computed JSON, then renders `observations.md` and `observations_sources.json`.
- `observations_checks/knee_and_scale_checks.py` performs every computed number:
  - the registered SBPL knee: slopes, bin placement, amplitude conversions, and ΔN_eff (exact Beta-function closed form plus archived-method reproduction)
  - the braneworld scale-setting translation of the registered static branch
  - dark-radiation bookkeeping
- `observations_checks/project_scan_1yr.py` is a read-only scan of 5,312 project text files for any prior discussion of the 1/yr timing-model degeneracy, with a positive control.
- The project's knee test was read from `new-files/D-Blast 3/GPD-site-deploy/library/answers/checked/what-could-prove-hd-blast-wrong.json`, with its formula from the archived `NEXTLEVEL_DeltaNeff_Audit_PHYSICAL.py`. The later "BK" knee families were read from `Hdblast 5.4 transfers/HDBLAST_ONE_FILE_LATEST_TRANSFER_20260326.md`. All were read only; nothing was modified.

## Results

| Result | Value | Label |
|---|---|---|
| Sources annotated | 51 (constrains 16, challenges 3, method 8, context 24, supports 0) | inventory |
| Do 2025–2026 releases test the registered knee? | No likelihood-level test exists yet. The NANOGrav 20-yr and IPTA DR3 background results were not found. The bearing is indirect: amplitude, index trend, running, piecewise reconstruction | inconclusive |
| Knee position relative to 1/yr | f_k = 0.9972 × (1/yr): within 0.045 bins (NANOGrav 15-yr), 0.029 (EPTA DR2new), 0.013 (MeerKAT). The knee sits on the sky-position/proper-motion fit blind spot. No project file discusses this (0 of 5,312) | numerical; degeneracy from snippet |
| Knee vs NANOGrav 15-yr standard band | Above the 14th bin (27.7 nHz). Inside the band the frozen curve is a γ ≈ 3.02 power law | numerical |
| Frozen height vs published amplitudes at 1/yr (Ω ratio) | NANOGrav 2023: 1.00 (calibration); NANOGrav 2026 chromatic-noise: 1.31; EPTA DR2new: 0.92; PPTA: 1.44; MeerKAT: 0.25 (γ = 13/3) to 0.10 (free) | numerical |
| Trend of 2025–2026 noise-model reanalyses | Toward γ = 13/3 and lower amplitude. Unfavourable to the frozen γ = 3 low-frequency branch; not a rejection | conditional |
| Frozen low-frequency slope f² vs causal f³ tail | A cosmological reading needs a stated mechanism | conditional |
| ΔN_eff of the knee spectrum | 7.09×10⁻⁴: exact closed form Ω_k 2^(1/D) B(s, 1/D−s)/((a1−a2)D) with s = 0.2, matching quadrature to 2e-16 and the archived method to 6e-10 | exact-verified / numerical |
| Ω(25 Hz) of the knee spectrum vs LVK O1–O4a | 2.28×10⁻³⁵ vs ≤ 2.0×10⁻⁹ (26 orders below) | numerical |
| Registered model: ℓ, H·ℓ | ℓ = 9 (sympy exact; σ_c/6 = 1/ℓ). H·ℓ = 0.069271 | exact-verified / numerical |
| Tabletop scale-setting (ℓ ≤ 38.6 µm to 0.1 mm) | M5 ≳ 2.3–3.1×10⁸ GeV; λ^(1/4) ≳ 3.4–5.5 TeV; brane vacuum ρ^(1/4) ≳ 0.76–1.22 TeV; horizon-scale GW today f ≳ 5.2–8.4×10⁻⁵ Hz, 1,600–2,600 × f_knee | conditional |
| Detuning needed to put the static branch's horizon scale at the knee | δ ≈ 2×10⁻¹⁷–1.4×10⁻¹⁶ (or ℓ ≈ 270 m, excluded) | conditional |
| Relaxed RS de Sitter brane as today's dark energy | Needs ℓ ≈ 372 Mpc (3×10²⁹ × the tabletop bound) or δ ≈ 1.1×10⁻⁶² | conditional |
| Dark-radiation (bulk Weyl C/a⁴) allowance | ρ_dr/ρ_γ ≲ 0.016 (ACT DR6), ≲ 0.10 (BBN 2024) (Gaussian approximations) | numerical |

**Controls** (all pass; `knee_and_scale_checks.json → checks`, `builder_selftest.json`):

- Ω(f_k) = Ω_k.
- Slope limits are +2 and −3.
- Analytic vs finite-difference slope agree to 9e-10.
- **Wrong-formula control:** reading a1 = +3 as the low-frequency slope is rejected.
- **Calibration:** A = 2.4×10⁻¹⁵ at γ = 13/3 gives Ω(1/yr) within 0.3% of Ω_k.
- **Perturbed parameters:** f_k ± 10%, and c → 1.01c.
- H·ℓ agrees by two routes to 2e-13.
- The O(δ²) expansion of H² matches the registered value to 3.9e-8.
- **Known limit:** the frequency map gives 1.647×10⁻⁷ Hz at T = 1 GeV, g* = 100, the commonly quoted coefficient.
- The review statement M5 > 10⁸ GeV at ℓ = 0.1 mm is reproduced.
- Builder negative controls: a bad label, a non-https URL, an inconsistent year, a duplicate URL and a corrupted number are all rejected.
- The project scan's positive control finds 4 files.

Precision: IEEE double, mpmath at 30 digits for the Beta function, sympy exact for the model facts. The trapezoid ΔN_eff is identical for n = 5×10⁴ to 4×10⁵. Its 6e-10 offset from the exact value is band truncation at [10⁻¹², 10⁴] Hz, not discretisation error.

## Limitations

- Snippet-only evidence. Numbers transcribed from snippets must be checked against the papers before any external use. Some snippet summaries merged pages, and each such case is flagged.
- The search budget ran out, so absence of IPTA DR3 or NANOGrav 20-yr background results means "not found".
- The knee comparisons are descriptive, based on published power-law summaries rather than free-spectrum posteriors. They are not likelihood tests.
- The scale-setting is order of magnitude. It rests on RS2 low-energy relations with a Z2-doubled bulk (M_Pl² = M5³ℓ, ρ_vac = 3M_Pl²H²), with O((Hℓ)²) corrections ignored, instantaneous thermalisation and f*/H* = 1. The Eöt-Wash 38.6 µm figure is a Yukawa range used as a proxy for the RS power-law correction.
- The registered model has no radiation era and no perturbation spectrum. The CMB, BBN, JWST and H0 entries are therefore requirements, not tests.

## Reproduction

Run from `research/HDBLAST_CHECKPOINT_20260927/literature/`, with python3, numpy 2.x, scipy, sympy and mpmath. Each step uses 1 core and takes under 15 s.

```
python3 observations_checks/knee_and_scale_checks.py        # writes observations_checks/knee_and_scale_checks.json; prints 12 checks
python3 observations_checks/project_scan_1yr.py             # read-only scan of new-files/ text; writes observations_checks/project_scan_1yr.json
python3 observations_checks/build_observations_md.py        # validates records + 43 number checks; writes observations.md, observations_checks/observations_sources.json
python3 observations_checks/build_observations_md.py --selftest   # builder negative controls; writes observations_checks/builder_selftest.json
```

## Files

- `observations.md`: the annotated bibliography and analysis (the main deliverable)
- `observations_checks/sources_observations.py`: the source records and search log
- `observations_checks/knee_and_scale_checks.py` → `knee_and_scale_checks.json`: all computed numbers and controls
- `observations_checks/project_scan_1yr.py` → `project_scan_1yr.json`: the project scan
- `observations_checks/build_observations_md.py` → `observations_sources.json` and `builder_selftest.json`
- `observations_checks/MANIFEST.sha256.json`: hashes of the files above
