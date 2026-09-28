# Literature sweep 2: observations that could test or constrain HDBLAST (2023–2026)

**Work in progress, 28 September 2026.** Nothing here is a claim of the programme until it appears in a dated checkpoint or a Zenodo record.

**How this was done.** There were 48 distinct web searches, listed at the end; three more were refused when the shared search budget ran out. **No paper could be opened.** arxiv.org, zenodo.org and journal sites are blocked for downloads here. Every entry below therefore rests on **search-result snippets only**: titles, URLs and short summaries. Numbers quoted from papers are transcribed from those snippets and must be checked against the papers before reuse. Where a snippet summary did not make clear which page a statement came from, the entry says so. Every number *computed* here comes from `observations_checks/knee_and_scale_checks.py`, which writes `knee_and_scale_checks.json` and passes 12 built-in checks and controls. The builder cross-checks 43 numbers quoted in this text against that JSON.

Labels: **exact-verified**, **numerical**, **conditional**, **negative**, **inconclusive**; for sources, *supports / constrains / challenges / method / context*.

## 1. Bottom line

1. **The registered PTA knee is still untested at likelihood level. No 2025–2026 release performs the test it needs** (**inconclusive**). The new PTA results that bear on it do so only indirectly: amplitudes, index trends, running and piecewise reconstructions. The NANOGrav 20-yr and IPTA DR3 background analyses, which could test it directly, had not appeared in the searches.
2. **The knee sits exactly on a known blind spot** (**numerical** position; the blind spot comes from a snippet). f_k = 3.16×10⁻⁸ Hz is 0.9972 × (1/yr), within 0.045 frequency bins of 1/yr for NANOGrav 15-yr and 0.013 bins for MeerKAT. Fitting each pulsar's sky position and proper motion removes sensitivity around 1/yr. A bend there is therefore maximally degenerate with the timing model. Unless the analysis propagates that absorption, the knee's location is not identifiable. This is the most important methodological finding of this sweep. A scan of 5,312 project text files found **no** discussion of this degeneracy (`observations_checks/project_scan_1yr.json`; the positive control passed). It is standard in the PTA literature but new to this project.
3. **Inside the NANOGrav band the frozen curve is effectively a γ = 3 power law.** Here γ is the timing-residual index, with Ω ∝ f^(5−γ). Its unweighted effective index over the 14 standard bins is γ ≈ 3.02 (**numerical**). Curvature appears only in bins 13–14. The 2025–2026 noise-model reanalyses (EPTA; NANOGrav customized chromatic noise) move the inferred background toward γ = 13/3 and lower amplitude. That trend is **unfavourable** to the frozen low-frequency slope, but it is not a rejection. MeerKAT's free index, γ = 3.60 (+1.31/−0.89), still allows γ = 3.
4. **Amplitude.** The frozen height reproduces the 2023 NANOGrav amplitude at 1/yr (h_c = 2.40×10⁻¹⁵; **numerical**, a calibration). Relative to the 2025–2026 values it has 1.31× the power of NANOGrav's customized-noise amplitude and 1.44× that of PPTA. Against MeerKAT it has only 0.25× (γ = 13/3) to 0.1× (free index). The spread between teams already exceeds a factor of 4 in Ω, so one frozen height cannot match every team.
5. **Scale-setting sharply challenges the idea that the knee is the blast's signature** (**conditional**, order of magnitude). In the registered model H·ℓ = 0.06927 (dimensionless; ℓ = 9 is **exact-verified**). Torsion-balance tests bound ℓ ≲ 38.6–100 µm. That puts the static branch's brane vacuum energy at ρ^(1/4) ≳ 0.76–1.2 TeV. Its horizon-scale GW frequency today would be f ≳ 5.2×10⁻⁵–8.4×10⁻⁵ Hz. That is about 1,600–2,600 times the knee frequency, in the sub-mHz band, not the nHz band. Putting the knee there instead would need δ ≈ 2.0×10⁻¹⁷–1.4×10⁻¹⁶ rather than 10⁻³, or ℓ ≈ 270 m, which is excluded.
6. **Dark energy (DESI).** The relaxing fate is an empty RS de Sitter brane, with w = −1. Identifying it with today's acceleration needs ℓ ≈ 372 Mpc, 3×10²⁹ times the tabletop bound, or δ ≈ 1.1×10⁻⁶² (**conditional**). This is the Randall–Sundrum form of the cosmological-constant problem. HDBLAST currently makes **no** DESI prediction.
7. **Dark radiation is the cleanest direct constraint on a 5D blast remnant.** A bulk black hole or black brane (for example from the collapse fate) appears on the brane as C/a⁴. ACT DR6 (N_eff = 2.86 ± 0.13) implies ρ_dr/ρ_γ ≲ 0.016, and BBN 2024 implies ≲ 0.1 (Gaussian approximations, **numerical**). This applies once a radiation era exists, which the registered model has not yet produced.
8. **Unchanged guardrails** (**numerical**, reproduced). The knee spectrum's ΔN_eff = 7.09×10⁻⁴, which is exact via a Beta-function closed form and matches the archived ≈7×10⁻⁴. Its Ω(25 Hz) = 2.28×10⁻³⁵ is 25.9 orders below the LVK O1–O4a limit. Neither observation tests the knee.
9. **Not yet testable.** The CMB tilt (n_s ≈ 0.968–0.974), tensors (r < 0.034), JWST early galaxies and H0 all require a primordial perturbation spectrum or an expansion history after a radiation era. The registered model has neither yet. These are **open requirements**, not passes or failures.

## 2. What each observation would constrain in a braneworld / higher-dimensional-origin model

| Observation | What it constrains in a braneworld origin model | Current value (snippet) | Status for HDBLAST |
|---|---|---|---|
| PTA background (NANOGrav, EPTA/InPTA, PPTA, CPTA, MPTA) | Spectral shape (knees, breaks, infrared tail), amplitude, isotropy and correlation pattern of any relic GW. Also any radion/KK phase transition at T ~ 0.1–1 GeV | NG15 A = 2.4e-15 (2023), 2.1e-15 with customized noise (2026); EPTA 2.5e-15; PPTA 2.0e-15; MPTA 4.8e-15 | Knee untested at likelihood level; the knee sits on the 1/yr blind spot; index trend unfavourable (conditional) |
| CMB N_eff (Planck, ACT DR6, SPT-3G) | Bulk Weyl dark radiation C/a⁴, KK/radion relics, GW energy density before recombination | N_eff = 2.86 ± 0.13 (ACT DR6) | ρ_dr/ρ_γ ≲ 0.016 (Gaussian approx.); binding once a radiation era exists |
| CMB tensors (BICEP/Keck, SPT-3G) | Energy scale of any inflation-like phase; RS high-energy enhancement of tensors | r < 0.034 (2025 combination) | No HDBLAST tensor prediction yet |
| CMB scalar tilt | Mechanism for adiabatic, nearly scale-invariant perturbations | n_s = 0.968 ± 0.003 (CMB), 0.973–0.974 with BAO | Open requirement |
| BBN | Expansion rate at ~1 MeV: ρ²/2λ term (needs λ^(1/4) ≫ MeV), dark radiation of either sign | ΔN_eff = −0.09 ± 0.28 (one 2024 configuration); braneworld DR −12.1% to +6.2% at 10 MeV (2017) | Automatically satisfied for tabletop-scale λ; DR bound open |
| DESI BAO + SNe | Late-time w(z); a braneworld gives w = −1 on an RS dS brane, phantom-like behaviour only with induced gravity (DGP-type) | w0waCDM preferred at 2.8–4.2σ (DR2); 3.2σ after DES-Dovekie; Lyα DR2 2026 consistent with ΛCDM | Outside the registered model unless a late sector is added |
| LVK stochastic background | High-frequency relic GWs; 25 Hz corresponds to horizon-scale production at T ≈ 1.5×10⁸ GeV (f*/H* = 1) or 1.5×10⁶ GeV (f*/H* = 100) | Ω(25 Hz) ≤ 2.0e-9 (2/3), ≤ 2.8e-9 (flat) | Knee 26 orders below; a TeV-scale blast would land in the sub-mHz band instead |
| GW propagation (GW170817, GWTC-3/4 sirens, graviton mass) | Leakage/damping into non-compact dimensions, bulk shortcuts, KK dispersion | D = 3.95 (+0.09/−0.07) (GWTC-3); D = 4.38 (+1.91/−1.01) (GWTC-4 dark sirens); m_g < 1.92e-23 eV | Bites only if ℓ or a crossover scale is cosmological, which tabletop tests already exclude |
| Tabletop inverse-square law | The AdS radius ℓ, hence M5, λ and the physical unit of the registered model | Yukawa range < 38.6 µm (2020); ℓ < 0.1 mm, M5 > 1e8 GeV (review) | Sets ρ^(1/4) ≳ 0.8–1.2 TeV for the static branch (conditional) |
| JWST early galaxies | Small-scale primordial power, early growth | z ≈ 14.4 galaxies; debated overabundance | Untestable until HDBLAST predicts perturbations |
| Hubble tension | Early expansion (dark radiation), late expansion | Local 73.50 ± 0.81 vs early 67.24 ± 0.35 (7.1σ) | No HDBLAST H0 prediction; ACT DR6 disfavours dark-radiation fixes |

## 3. The project's registered PTA knee test and the 2025–2026 data

**The registered test** (project page `GPD-site-deploy/library/answers/checked/what-could-prove-hd-blast-wrong.json` and Zenodo record 17968738, 17 Dec 2025) is a smooth broken power law. The formula is copied from the archived `NEXTLEVEL_DeltaNeff_Audit_PHYSICAL.py`:

    Omega(f) = Omega_k [((f/f_k)^(a1 D) + (f/f_k)^(a2 D))/2]^(-1/D),  f_k = 3.16e-8 Hz, Omega_k = 8.00e-9, a1 = +3, a2 = -2, D = 2

Failure criteria stated by the project: a straight power law through the knee region; the bend fading under a full noise-aware analysis; teams disagreeing about where it sits; or SMBHBs explaining everything. A later 'BK' (Bessel-K, ν = 2) fingerprint gate, from the 26 Mar 2026 handoff, adds 'LowK' (≈2.96 nHz) and 'HighK' (≈33.1 nHz) knee families. The same handoff records that the IPTA DR2 common-process free spectrum and the archived NANOGrav HD-like export **fail** the frozen BK gate (project-internal; not re-verified here).

### 3.1 Shape facts (exact/numerical, script)

- Local slope d ln Ω/d ln f → +2 far below the knee, -3 far above, and -0.5 at the knee. The equivalent residual indices are γ = 3, 8 and 5.5. Wrong-formula control: reading a1 = +3 as the low-frequency slope, i.e. Ω ∝ f³, is **rejected**. The archived formula gives f² below the knee, matching the book's wording.
- A causal, horizon-limited source in the radiation era has a universal **f³** infrared tail (arXiv 2010.03568). The frozen **f²** therefore needs a stated mechanism if it is to be cosmological.

### 3.2 Where the knee sits relative to each data set (numerical)

| Data set | Span T (yr) | Resolution 1/T (nHz) | Knee in units of 1/T | Knee − 1/yr (bins) |
|---|---|---|---|---|
| NANOGrav_15yr | 16.03 | 1.98 | 15.99 | -0.045 |
| EPTA_DR2new | 10.3 | 3.08 | 10.27 | -0.029 |
| MPTA_4p5yr | 4.5 | 7.04 | 4.487 | -0.013 |
| CPTA_DR1 | 3.0 | 10.6 | 2.992 | -0.0083 |

- NANOGrav 15-yr's 14 standard bins end at 27.7 nHz, so the knee lies **above** the standard band.
- Spans: NANOGrav 16.03 yr is adopted. It is consistent with 'nearly 16 years' in the snippet and with the project page's '14 frequencies from about 2 to 28 nHz'. The EPTA DR2new, MeerKAT and CPTA spans are from snippets.
- Project BK windows: the LowK narrow window is 0.233 nHz wide, 0.12 of NANOGrav's frequency resolution. No current data set can resolve a knee position that finely. Sampling a free spectrum at 0.025 nHz spacing, as the handoff requests, yields strongly correlated values, not independent information. This is a general Fourier-resolution point, not taken from a cited source. The HighK window **contains** 1/yr, so it has the same timing-model degeneracy as the SBPL knee.
- Perturbed-parameter control: moving f_k by −10% / +10% lowers Ω at 1/yr by 12% / 0.4%. The power at 1/yr alone barely pins the knee position. Locating the knee needs the bins on both sides of 1/yr, which is exactly where the timing-model fit absorbs power.

### 3.3 Frozen curve vs a γ = 13/3 power law in the NANOGrav bins (numerical, descriptive)

| Bin | f (nHz) | Ω_SBPL | local γ | Ω_SBPL / Ω_PL(A = 2.4e-15) | Ω_SBPL / Ω_PL(A = 2.1e-15, 2026) |
|---|---|---|---|---|---|
| 1 | 1.98 | 4.43e-11 | 3 | 0.035 | 0.046 |
| 2 | 3.95 | 1.77e-10 | 3 | 0.089 | 0.12 |
| 3 | 5.93 | 3.98e-10 | 3 | 0.15 | 0.2 |
| 5 | 9.88 | 1.11e-09 | 3 | 0.3 | 0.39 |
| 8 | 15.8 | 2.83e-09 | 3 | 0.56 | 0.74 |
| 10 | 19.8 | 4.41e-09 | 3.05 | 0.76 | 0.99 |
| 12 | 23.7 | 6.20e-09 | 3.27 | 0.94 | 1.2 |
| 13 | 25.7 | 7.05e-09 | 3.56 | 1 | 1.3 |
| 14 | 27.7 | 7.71e-09 | 4.05 | 1.1 | 1.4 |

Reading: *if* the true background were the published fixed-index (γ = 13/3) fit, the frozen curve would carry only 3.5–15% of its power in the three lowest bins. This compares a model with a model summary, not with the data. Free-index fits have broad errors: EPTA DR2 γ = 4.19 (+0.73/−0.63), MeerKAT γ = 3.60 (+1.31/−0.89). The curves cross near 25 nHz. PTA sensitivity to red processes is concentrated at the lowest frequencies (general PTA knowledge, not a cited snippet). The frozen curve therefore stands or falls mainly on the low-bin power and the index, not on the knee itself. **This is not a likelihood analysis.** A proper test fits the frozen SBPL, with no free parameters, to the corrected free-spectrum KDEs or the PPL posterior, together with the timing-model transmission.

### 3.4 Do the 2025–2026 releases bear on the knee?

| Result (year) | Bears on | Direction for the frozen knee |
|---|---|---|
| NANOGrav running of the spectral index (2025) | curvature within 2–28 nHz | no curvature required (β consistent with 0, BF 0.69): **inconclusive** |
| NANOGrav piecewise power-law reconstruction (2026) | broken/doubly-broken spectra, BMA | the right comparison target; numbers not visible in snippets: **inconclusive** |
| NANOGrav customized chromatic noise (2026) | amplitude, HD significance | amplitude down to 2.1e-15; frozen height now 1.31× in Ω, within errors: **mildly unfavourable** |
| NANOGrav 15-yr free-spectrum correction / erratum (2026) | the posteriors the project screened | project screens need a rerun on corrected KDEs: **method** |
| EPTA DR2 improved noise models (2025) | index, amplitude | toward γ = 13/3, lower amplitude: **unfavourable to the γ = 3 low branch** |
| MeerKAT PTA 4.5 yr (2024/2025) | amplitude, index at high cadence | amplitude 4–10× the frozen height in Ω; γ = 3 allowed: **cross-team spread** |
| Five-PTA combination (Dec 2025) | amplitude, power-law exponent | power-law only; could host a knee test: **method** |
| NANOGrav 20-yr, IPTA DR3 | direct test | not released in the searches: **pending** |

### 3.5 Frequency–temperature map (numerical, conditional on f*/H* = 1)

A horizon-scale feature at f_k corresponds to production at T ≈ 0.19–0.28 GeV (g* from 106.75 down to 10.75), i.e. the QCD era. The script's map reproduces the commonly quoted coefficient: f ≈ 1.65e-07 Hz at T = 1 GeV, g* = 100. A five-dimensional origin of a feature at that temperature has published precedent: a radion/KK phase transition at MeV–GeV (Phys. Rev. D 108, 095017). That is a different mechanism from the registered HDBLAST shell.

## 4. Scale-setting the registered model with tabletop gravity (conditional)

**Exact-verified** (sympy): W(1) = 1/3, U(1) = -2/27, 1/ℓ = √(−U(1)/6) = 1/9, so ℓ = 9. The critical tension 2W(1) = 2/3 satisfies σ_c/6 = 1/ℓ (verified). This agrees with the theory sweep's cross-check. **Numerical:** H·ℓ = 9/ρ_b = 0.0692709, and 9√(H²) agrees to 2.3e-13. The O(δ²) expansion of H² matches the registered value to 3.9e-08 (relative). ρ/λ ≈ (Hℓ)²/2 = 0.0024, so the static branch is in the low-energy RS regime.

| ℓ bound used | M5 (GeV) | λ^(1/4) (TeV) | ρ_vac^(1/4) (TeV) | T if thermalised (TeV) | f today, f*/H* = 1 (Hz) | f / f_knee |
|---|---|---|---|---|---|---|
| ℓ ≤ 38.6 µm (Eöt-Wash 2020 Yukawa range, proxy) | 3.1×10⁸ | 5.52 | 1.22 | 0.502 | 8.4×10⁻⁵ | 2,640 |
| ℓ ≤ 0.1 mm (review statement) | 2.3×10⁸ | 3.43 | 0.759 | 0.312 | 5.2×10⁻⁵ | 1,640 |

All energies scale as ℓ^(−1/2). A smaller ℓ only raises them, so the table gives **lower bounds** on the energy scales and the frequency. Assumptions: RS2 with a Z2-doubled bulk, M_Pl² = M5³ℓ, ρ_vac = 3M_Pl²H² on the brane, O((Hℓ)²) corrections ignored, full and instantaneous thermalisation, f*/H* = 1, standard radiation era afterwards. The Eöt-Wash 38.6 µm figure is a **Yukawa** range, used only as a proxy. The RS correction is a power law (1 + 2ℓ²/3r²), and the 2026 inverse-square-law review has the proper power-law limits. The review statement 'ℓ < 0.1 mm ⇒ M5 > 10⁸ GeV' is reproduced (control).

**Interpretation (conditional).** If the registered δ = 10⁻³ static branch is the pre-radiation state, and ℓ obeys tabletop bounds, the natural relic-GW frequency lies in the sub-mHz band. A PTA knee would then have to come from later physics, not from the blast's horizon scale. This does not falsify HDBLAST. It shows that the PTA knee and the registered 5D model are currently **disconnected**, and any link needs an explicit mechanism.

## 5. Annotated bibliography

Each entry gives the URL, the year, the finding as reported in the search snippet, its relevance to HDBLAST, and the implication label. **Access: every entry is snippet-only; no full text was read.**

### Pulsar timing arrays (PTA) and the nanohertz background

**The NANOGrav 15-year Data Set: Evidence for a Gravitational-Wave Background** (2023). ApJL 951, L8 (2023).  
URL: https://arxiv.org/abs/2306.16213  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q11

*Finding.* Hellings-Downs-correlated signal in 67 pulsars. For a fiducial f^-2/3 strain spectrum the amplitude is 2.4 (+0.7/-0.6) x 10^-15 (median, 90% interval) at 1/yr. A power-law background is favoured over independent pulsar noise alone with a Bayes factor above 10^14 (snippet wording). Consistent with supermassive black-hole binaries (SMBHBs); exotic cosmological or astrophysical sources not excluded.

*Relevance to HDBLAST.* This amplitude is the calibration anchor of the registered knee: the frozen Omega_k = 8.00e-9 converts to h_c = 2.40e-15 at 1/yr (h = 0.674), reproducing this median to better than 0.1% in strain (0.3% in Omega; script check).

**The NANOGrav 15 yr Data Set: Running of the Spectral Index** (2025). ApJL 978, L29 (2025); arXiv Aug 2024.  
URL: https://arxiv.org/abs/2408.10166  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q1, Q5

*Finding.* Running-power-law fit (logarithmic frequency dependence of the index). 95% credible interval for the running beta in [-0.80, 2.96], consistent with zero. Bayes factor running vs constant power law 0.69 +/- 0.01 (inconclusive). The constant power law still suffices.

*Relevance to HDBLAST.* A knee is an extreme form of spectral curvature. Within the NANOGrav band the frozen curve is almost a pure power law (curvature only in bins 13-14), so this null result neither supports nor excludes it. It does show that current data do not require curvature.

**The NANOGrav 15 yr Data Set: Piecewise Power-Law Reconstruction of the Gravitational-Wave Background** (2026). ApJL (2026), doi 10.3847/2041-8213/ae7086.  
URL: https://arxiv.org/abs/2601.09481  
*Access:* search-snippet-only · *Implication:* method · *Queries:* Q11, Q13, Q16

*Finding.* Piecewise power-law (PPL) spectral reconstruction: constant, broken, doubly broken models combined by Bayesian model averaging. Described as closer to physically realistic (especially cosmological) spectra than the free spectrum. The snippets did not report the fitted break positions or Bayes factors.

*Relevance to HDBLAST.* This is the collaboration's own analogue of the project's knee test. The frozen SBPL should be compared with the PPL posterior (bin-by-bin power and break-frequency posterior), not only with summary points. Result numbers must be read from the paper before any claim.

**The NANOGrav 15 yr Data Set: Impacts of Customized Chromatic Noise Models on Gravitational Wave Analyses** (2026). arXiv June 2026 (companion: arXiv 2606.28571).  
URL: https://arxiv.org/abs/2606.28554  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q18

*Finding.* Customized chromatic (interstellar-medium) noise models for the 15-yr pulsars. Bayes factor for Hellings-Downs correlations over an uncorrelated common red process: 1571 +/- 14 (14 Fourier modes), about 8 times earlier. Power-law amplitude at fixed index reduced to 2.1 (+0.6/-0.5) x 10^-15.

*Relevance to HDBLAST.* Relative to this amplitude, the frozen knee has 1.31 times the power at 1/yr (1.14 times the strain). That is inside the quoted uncertainty. Better noise modelling lowered the amplitude and strengthened the correlation signature. The knee's height was calibrated to the 2023 value, and the direction of change matters for any refreeze.

**KDE Representations of the Gravitational Wave Background Free Spectra Present in the NANOGrav 15-Year Dataset (corrected release; erratum Agazie et al. 2026, ApJL 1006, L67)** (2026). Zenodo data record + ApJL erratum (2026).  
URL: https://zenodo.org/records/21844115  
*Access:* search-snippet-only · *Implication:* method · *Queries:* Q14

*Finding.* A bug in the parallel-tempering routine of PTMCMCSampler biased the free-spectrum posteriors of the 15-yr data set. The record re-releases corrected free-spectrum KDEs. The erratum states SMBHB population inference was affected. The project's own site notes the SMBHB conclusion was unchanged.

*Relevance to HDBLAST.* Any project screen of the knee that used pre-2026 NANOGrav free-spectrum posteriors or 'summary points' derived from them should be rerun on the corrected KDEs. The frozen numbers must stay frozen.

*Attribution note.* Snippet summary combined the Zenodo record text and an erratum citation; the erratum page itself was not seen.

**The NANOGrav 15-year Data Set: Search for Signals from New Physics** (2023). ApJL 951, L11 (2023).  
URL: https://arxiv.org/abs/2306.16219  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q21

*Finding.* Inflation, scalar-induced GWs, first-order phase transitions, cosmic strings and domain walls were tested. All except stable field-theory cosmic strings can reproduce the signal, some with Bayes factors O(10)-O(100) over the SMBHB model. The results are model-sensitive and not conclusive.

*Relevance to HDBLAST.* Sets the bar for any HDBLAST PTA claim. A better fit than a fixed SMBHB template is not evidence for new physics, and the project's BIC screens are weaker than these analyses.

**The NANOGrav 15 yr and 20 yr Datasets: Timing Events and Pulse Shape Changes** (2026). ApJ 1005 (June 2026).  
URL: https://iopscience.iop.org/article/10.3847/1538-4357/ae6db6  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q12

*Finding.* A 2026 paper already describes timing events and pulse-shape changes in the NANOGrav 20-yr dataset. No 20-yr gravitational-wave background result was found in the searches.

*Relevance to HDBLAST.* A 20-yr span gives frequency resolution of about 1.58 nHz, and 1/yr falls near bin 20. The 20-yr background analysis is the next dataset that could test the knee directly. It had not been released as of these searches.

**The second data release from the European Pulsar Timing Array III. Search for gravitational wave signals** (2023). A&A 678, A50 (2023).  
URL: https://arxiv.org/abs/2306.16214  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q7

*Finding.* DR2new (latest 10.3 yr, 25 pulsars, plus about 3.5 yr of InPTA data for 10 of them): Bayes factor 60, false-alarm probability about 0.1% (at least 3 sigma). At fixed index 13/3: A = (2.5 +/- 0.7) x 10^-15 at 1/yr. Full DR2, HD process: log10 A = -14.54 (+0.28/-0.41), gamma = 4.19 (+0.73/-0.63).

*Relevance to HDBLAST.* The frozen knee lies within 0.03 frequency bins of 1/yr for DR2new (resolution 3.08 nHz). The EPTA free-index fit sits near gamma = 13/3, not the gamma = 3 that the frozen curve follows below the knee.

**A clearer view of gravitational-wave signals in pulsar timing arrays (AEI news on improved EPTA DR2 noise models)** (2025). Max Planck Institute for Gravitational Physics news item.  
URL: https://www.aei.mpg.de/1323827/a-a-clearer-view-of-gravitational-wave-signals-in-pulsar-timing-arrays  
*Access:* search-snippet-only · *Implication:* challenges · *Queries:* Q15, Q17

*Finding.* Improved EPTA DR2 noise modelling (pulsar-intrinsic, epoch-correlated and transient noise, with noise-model averaging) makes the background's strain index consistent with -2/3 (gamma = 13/3). It gives a lower amplitude at 1/yr than earlier analyses. Chromatic noise errors have whiter spectra, which push naive fits toward flatter indices.

*Relevance to HDBLAST.* Two independent 2025-2026 reanalyses (EPTA; NANOGrav, entry above) move toward gamma = 13/3 and lower amplitude. The frozen curve's low-frequency branch is Omega ~ f^2 (gamma = 3). This trend is unfavourable to the frozen shape, but it is not a likelihood-level rejection.

*Attribution note.* Year and details from the search summary; the news page date was not visible. The 'whiter spectra' sentence may come from the NANOGrav chromatic-noise work rather than this news item.

**Search for an Isotropic Gravitational-wave Background with the Parkes Pulsar Timing Array** (2023). ApJL (2023), doi 10.3847/2041-8213/acdd02.  
URL: https://arxiv.org/abs/2306.16215  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q9

*Finding.* Common-spectrum process with A = 2.0 +/- 0.2 x 10^-15 (h ~ f^-2/3). Hellings-Downs consistency with false-alarm probability p < about 0.014. The signal strength appeared time-dependent, contrary to expectations for an isotropic background.

*Relevance to HDBLAST.* The frozen knee has 1.44 times PPTA's power at 1/yr. The reported time-dependence is a caution for any fine spectral feature.

**Searching for the nano-Hertz stochastic gravitational wave background with the Chinese Pulsar Timing Array Data Release I** (2023). RAA 23, 075024 (2023).  
URL: https://arxiv.org/abs/2306.16216  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q20

*Finding.* 57 millisecond pulsars, close to 3 yr with FAST. Hellings-Downs evidence at 4.6 sigma around 14 nHz (discrete-frequency method). log10 A_c = -14.4 (+1.0/-2.8) for strain index in [-1.8, 1.5].

*Relevance to HDBLAST.* Resolution about 10.6 nHz, so the knee sits at bin 2.99, essentially at 1/yr. The amplitude is too loose to test the knee.

**The MeerKAT Pulsar Timing Array: The 4.5-year data release and the noise and stochastic signals of the millisecond pulsar population** (2024). arXiv Dec 2024.  
URL: https://arxiv.org/abs/2412.01148  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q3

*Finding.* 83 pulsars, 4.5 yr, high cadence. Common signal log10 A = -14.25 (+0.21/-0.36), gamma = 3.60 (+1.31/-0.89). At gamma = 13/3, log10 A = -14.28 +/- 0.21, ln Bayes factor 4.46. The amplitude is larger than other PTAs report.

*Relevance to HDBLAST.* MeerKAT's free-gamma posterior still allows gamma = 3. Its resolution (7.04 nHz) puts the knee at bin 4.49, so its bins straddle 1/yr more coarsely than NANOGrav's.

**The MeerKAT Pulsar Timing Array: The first search for gravitational waves with the MeerKAT radio telescope** (2024). MNRAS 536, 1489 (2025).  
URL: https://arxiv.org/abs/2412.01153  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q6

*Finding.* Sky-averaged h_c,yr = 7.5 (+0.8/-0.9) x 10^-15 at strain index -0.26, or 4.8 (+0.8/-0.9) x 10^-15 at -2/3. This is inconsistent with other PTAs' common-noise results by at least about 1.4 sigma. Hellings-Downs significance is 3-3.4 sigma depending on noise assumptions.

*Relevance to HDBLAST.* At 1/yr the frozen knee has only 0.25 (gamma = 13/3) to 0.10 (free index) of MeerKAT's power. The cross-PTA amplitude spread already exceeds a factor 4 in Omega. A single frozen height cannot match all teams, which is one of the project's own failure criteria (teams disagreeing).

**The MeerKAT Pulsar Timing Array: Maps of the gravitational-wave sky with the 4.5 year data release** (2024). MNRAS 536, 1501 (2025).  
URL: https://arxiv.org/abs/2412.01214  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q3

*Finding.* Gravitational-wave sky maps from the 4.5-yr data. The snippet reported tentative background evidence but no anisotropy numbers.

*Relevance to HDBLAST.* Anisotropy is an SMBHB discriminator. A cosmological HDBLAST relic would be isotropic to high precision, so detected anisotropy would count against a cosmological reading of the knee.

**Comparing Recent Pulsar Timing Array Results on the Nanohertz Stochastic Gravitational-wave Background** (2024). ApJ 966, 105 (2024); arXiv 2309.00693.  
URL: https://iopscience.iop.org/article/10.3847/1538-4357/ad36be  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q8

*Finding.* EPTA, InPTA, NANOGrav and PPTA results assessed on equal footing: background spectral parameters agree within 1 sigma. A standardized noise model reduces tensions in pulsar noise parameters.

*Relevance to HDBLAST.* Background for the project's cross-team 'overlap' criterion: agreement is at the power-law-parameter level, not at the level of a shared knee.

**Stochastic gravitational-wave background search using data from five pulsar timing arrays** (2025). arXiv Dec 2025 (W.-W. Yu, B. Allen).  
URL: https://arxiv.org/abs/2512.08666  
*Access:* search-snippet-only · *Implication:* method · *Queries:* Q2, Q4

*Finding.* Public pulse arrival times from five PTAs combined into a 121-pulsar data set, about four times larger than any single PTA's, using a 'direct combination' method for shared pulsars. Central result: posterior on amplitude and exponent of a power-law energy-density spectrum.

*Relevance to HDBLAST.* This is the first public multi-PTA combination found (it is not IPTA DR3). It fitted only a power law, so it does not test the knee, but the combined data could host an independent knee test.

**Pulsar timing arrays: the emerging gravitational-wave landscape** (2026). arXiv Mar 2026 (review).  
URL: https://arxiv.org/abs/2603.13643  
*Access:* search-snippet-only · *Implication:* challenges · *Queries:* Q19

*Finding.* Six PTAs (NANOGrav, EPTA, InPTA, PPTA, CPTA, MPTA) report evidence for the background. The review argues that the perceived tension between current amplitudes and standard merger models is largely resolved by new insights into SMBHB populations.

*Relevance to HDBLAST.* If the SMBHB amplitude tension is resolved, a main motivation for an exotic source of the hum weakens.

**International Pulsar Timing Array (home page)** (n/a). IPTA web site.  
URL: https://ipta4gw.org/  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q2, Q8

*Finding.* Appeared in results for IPTA DR3 queries. The snippet summary describes IPTA-DR3 as a reanalysis combining the member PTAs (and CPTA) that is expected to be more sensitive than any single data set. No DR3 background result was found.

*Relevance to HDBLAST.* IPTA DR3 would be the decisive combined test at 1/yr. The project's own site page (checked Sept 2026) also records it as not yet released.

*Attribution note.* The DR3 description could not be tied to a specific page among the results; treat as unverified.

**A sensitivity curve approach to tuning a pulsar timing array in the detection era (with the hasasia package)** (2025). CQG (2025), doi 10.1088/1361-6382/adbbab.  
URL: https://arxiv.org/abs/2409.00336  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q10

*Finding.* Search summary: fitting each pulsar's sky position and proper motion causes an extra loss of sensitivity around f = 1/yr, and fitting parallax causes one around 2/yr. These appear as large spikes in PTA sensitivity curves (computed with hasasia).

*Relevance to HDBLAST.* This is the key methodological fact for the registered knee. The script finds f_k = 0.9972 x (1/yr), within 0.045 bins of 1/yr for NANOGrav 15-yr and 0.013 bins for MeerKAT. A spectral bend exactly where the timing model absorbs power is maximally degenerate with it. A knee test must propagate the timing-model transmission function, or the knee location is not identifiable.

*Attribution note.* The sentence about 1/yr and 2/yr came from the results for this query. The summary did not say whether it came from this paper or from arXiv 2608.00250 (directional anisotropic sensitivity curves).

**Causal gravitational waves as a probe of free streaming particles and the expansion of the Universe** (2020). arXiv Oct 2020.  
URL: https://arxiv.org/abs/2010.03568  
*Access:* search-snippet-only · *Implication:* challenges · *Queries:* Q23

*Finding.* For causal sources with finite correlation length, the infrared tail of a GW spectrum produced in radiation domination scales as f^3. The tail is modified by free-streaming particles and by non-standard expansion.

*Relevance to HDBLAST.* The frozen curve rises as f^2 below the knee, not the universal f^3 causal tail of a horizon-limited radiation-era source. A cosmological HDBLAST reading of the knee therefore needs a named mechanism (non-standard expansion, free-streaming or a specific source) that yields f^2.

**A Universal CMB B-Mode Spectrum from Early Causal Tensor Sources** (2026). arXiv Jan 2026.  
URL: https://arxiv.org/abs/2601.20967  
*Access:* search-snippet-only · *Implication:* method · *Queries:* Q23

*Finding.* Early causal tensor sources share a universal infrared scaling and predict the same B-mode angular distribution.

*Relevance to HDBLAST.* If the blast is a causal (sub-horizon) tensor source, its CMB B-mode imprint is fixed up to amplitude. This is a cross-check between a PTA-band claim and CMB tensor limits.

**Nonlinear Gravitational Wave Memory: Universal Low-Frequency Background** (2025). arXiv Nov 2025.  
URL: https://arxiv.org/abs/2511.08514  
*Access:* search-snippet-only · *Implication:* method · *Queries:* Q23

*Finding.* Nonlinear memory dominates the low-frequency GW spectrum during radiation and kination domination.

*Relevance to HDBLAST.* A second universal infrared contribution that any cosmological spectral template, the knee included, should respect.

**Pulsar timing array stochastic background from light Kaluza-Klein resonances** (2023). Phys. Rev. D 108, 095017 (2023); arXiv 2306.17071.  
URL: https://doi.org/10.1103/PhysRevD.108.095017  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q22

*Finding.* In a warped 5D model (UV brane at the Planck scale, dark brane at the GeV scale), PTA data can be fitted by a first-order confinement phase transition of a radion at the MeV-GeV scale, feebly coupled to the Standard Model. Many existing embeddings are not viable because of radion/graviton phenomenology. A multi-brane set-up is proposed to remain consistent with collider and gravity tests.

*Relevance to HDBLAST.* This is prior art for a five-dimensional origin of the nHz signal, via a phase transition near T ~ 0.1-1 GeV. That matches the script's map of the knee to T ~ 0.19-0.28 GeV for f*/H* = 1. It is not the registered HDBLAST mechanism, and it shows that 5D PTA explanations face strong gravity/collider constraints.

**Periodic Spectral Features in the NANOGrav 15-Year Gravitational Wave Background: A Phenomenological Analysis** (unknown). Zenodo record (not peer reviewed; authorship not seen).  
URL: https://zenodo.org/records/17956111  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q1

*Finding.* Claims a systematic alternation of spectral excess and deficit across frequency bins, with a period of about 4 nHz.

*Relevance to HDBLAST.* A caution: with 14 correlated free-spectrum bins, apparent features are easy to find, and features at the edges of the band (like the knee) are the least constrained.

*Attribution note.* Not a project record (its ID does not appear in the project's files). Not peer reviewed; treat with care.

### CMB: dark radiation, spectral tilt and primordial tensors

**The Atacama Cosmology Telescope: DR6 Constraints on Extended Cosmological Models** (2025). arXiv Mar 2025.  
URL: https://arxiv.org/abs/2503.14454  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q27, Q30, Q32

*Finding.* No evidence for new free-streaming light species: N_eff = 2.86 +/- 0.13. Self-interacting dark radiation N < 0.134. Dark-radiation models are not favoured, because extra radiation increases Silk damping at high multipoles, where ACT DR6 is sensitive.

*Relevance to HDBLAST.* In a braneworld the Weyl ('dark radiation') term C/a^4, set by a bulk black-hole mass, counts as Delta N_eff. A Gaussian approximation from these numbers gives Delta N_eff < 0.071, so rho_dr/rho_gamma < 0.016 at recombination (script). This bounds any bulk remnant of the blast, for example the collapse fate's black brane, and any GW relic, if a radiation era is ever achieved.

**ACT DR6 Insights on alpha-Attractor Inflationary Models and Reheating** (2025). arXiv May 2025.  
URL: https://arxiv.org/abs/2505.01517  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q30

*Finding.* As reported there: Planck 2018 + ACT DR6 gives n_s = 0.9709 +/- 0.0038, rising to 0.9743 +/- 0.0034 with DESI Y1 BAO (the P-ACT-LB combination gives n_s = 0.974 +/- 0.003). Starobinsky-type plateau models become mildly disfavoured at about 2 sigma or more.

*Relevance to HDBLAST.* Any 'blast instead of (or before) inflation' model must produce adiabatic, nearly scale-invariant perturbations with n_s about 0.97. The registered model has no perturbation-spectrum prediction yet, so this is currently an open requirement, not a test it passes or fails.

*Attribution note.* n_s numbers as quoted in the snippet from this and related pages.

**Inflation at the End of 2025: Constraints on r and n_s Using the Latest CMB and BAO Data** (2025). Open Journal of Astrophysics (2026); arXiv Dec 2025.  
URL: https://arxiv.org/abs/2512.10613  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q29, Q31

*Finding.* Planck + SPT + ACT + BICEP/Keck: n_s = 0.9682 +/- 0.0032 and r < 0.034 (95%). Adding DESI BAO leaves r unchanged but shifts n_s to 0.9728 +/- 0.0029, driven by marginal CMB-DESI differences.

*Relevance to HDBLAST.* Primordial tensors: any inflation-like phase on the brane is bounded by r < 0.034. In high-energy RS inflation, tensor and scalar amplitudes are both enhanced (next entries). HDBLAST has no tensor prediction yet.

**SPT-3G D1: Constraints on inflationary gravitational waves with two years of SPT-3G data** (2025). Phys. Rev. D (Dec 2025).  
URL: https://arxiv.org/abs/2505.02827  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q29

*Finding.* SPT-3G alone over the BICEP/Keck field: r < 0.25 (95%), sigma(r) = 0.067, with delensing.

*Relevance to HDBLAST.* An independent tensor channel. It is not competitive yet, but the South Pole Observatory programme targets sigma(r) = 0.001 (snippet).

**Improved limits on the tensor-to-scalar ratio using BICEP and Planck data** (2022). Phys. Rev. D 105, 083524 (2022).  
URL: https://arxiv.org/abs/2112.07961  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q29

*Finding.* Planck + BICEP/Keck 2018 + BAO: r < 0.032 (95%). BK18 alone: r < 0.036.

*Relevance to HDBLAST.* The standard tensor bound, before the 2025 combination above.

**SPT-3G D1: CMB temperature and polarization power spectra and cosmology from 2019 and 2020 observations of the SPT-3G Main field** (2025). Phys. Rev. D (2025).  
URL: https://arxiv.org/abs/2506.20707  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q28

*Finding.* SPT-3G alone: H0 = 66.66 +/- 0.60. SPT + Planck + ACT: H0 = 67.19 +/- 0.38, sigma8 = 0.8137 +/- 0.0037. CMB alone shows no evidence beyond LambdaCDM. There is a 2.8 sigma CMB-vs-DESI DR2 difference within LambdaCDM.

*Relevance to HDBLAST.* The early-universe expansion history HDBLAST must reproduce, if it ever reaches a radiation era.

**Tensor Perturbations from Brane-World Inflation with Curvature Effects** (2014). Phys. Rev. D 89, 063501 (2014).  
URL: https://arxiv.org/abs/1308.5765  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q45

*Finding.* In RS braneworld inflation, tensor (and scalar) amplitudes are enhanced at high energy relative to 4D. Gauss-Bonnet and induced-gravity corrections suppress the RS enhancement; the Gauss-Bonnet term can break the standard consistency relation at high energy.

*Relevance to HDBLAST.* If HDBLAST produced perturbations on the brane at energies comparable to the brane tension, the r and n_t relations would differ from 4D. With the registered H*ell = 0.069 the brane is in the low-energy regime (rho/lambda about 0.0024), so these corrections would be small.

*Attribution note.* The enhancement/consistency statement appeared in the summary for several results; this paper's abstract matches the curvature-effect part.

### Big Bang nucleosynthesis and braneworld dark radiation

**The 2024 BBN baryon abundance update** (2024). arXiv Jan 2024.  
URL: https://arxiv.org/abs/2401.15054  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q32

*Finding.* Delta N_eff from BBN ranges from -0.09 +/- 0.28 (one configuration: BAO+BBN abundances, PArthENoPE v3.0); results depend on the helium data and deuterium-burning rates.

*Relevance to HDBLAST.* Expansion rate at T ~ 1 MeV. For a braneworld this bounds the rho^2/(2 lambda) term, which needs lambda^(1/4) far above MeV. That is automatic if tabletop bounds set lambda^(1/4) of a few TeV. It also bounds dark radiation: a Gaussian approximation gives Delta N_eff < 0.46, so rho_dr/rho_gamma < 0.10 (script).

**New observational limits on dark radiation in braneworld cosmology** (2017). Phys. Rev. D 95, 083516 (2017); arXiv 1706.03630.  
URL: https://doi.org/10.1103/PhysRevD.95.083516  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q33

*Finding.* BBN restricts braneworld and particle dark radiation at 10 MeV to between -12.1% and +6.2% of the total background energy density. An older analysis (arXiv astro-ph/0203272) gave -1.23 <= rho_d/rho_gamma <= 0.11 from BBN, narrowed to -0.41..0.105 with the CMB.

*Relevance to HDBLAST.* A direct bound on the bulk Weyl term. Negative dark radiation, which is allowed in braneworlds, is bounded too. This applies to any HDBLAST scenario in which the blast leaves bulk mass or energy.

*Attribution note.* The snippet summary merged two papers; the -12.1%/+6.2% figure is attributed to the 2017 paper, the -1.23..0.11 figure to astro-ph/0203272.

### DESI and evolving dark energy

**DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints** (2025). arXiv Mar 2025 (+ DESI guide page https://www.desi.lbl.gov/2025/03/19/desi-dr2-results-march-19-guide/).  
URL: https://arxiv.org/abs/2503.14738  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q24

*Finding.* w0waCDM is preferred over LambdaCDM at 3.1 sigma (DESI + CMB) and at 2.8-4.2 sigma when supernovae are added. The best fit has w > -1 today and w < -1 in the past, crossing -1 near z = 0.5.

*Relevance to HDBLAST.* The registered model's relaxing fate is an empty RS de Sitter brane with w = -1 exactly. At the registered detuning, H*ell = 0.069, so identifying that brane with today's dark energy needs ell of about 372 Mpc (excluded by tabletop tests by a factor 3e29) or delta of about 1e-62 (script). HDBLAST therefore makes no DESI prediction. A confirmed evolving w would need an extra late-time sector.

**The Dark Energy Survey Supernova Program: A Reanalysis Of Cosmology Results And Evidence For Evolving Dark Energy With An Updated Type Ia Supernova Calibration** (2025). MNRAS (2026); arXiv Nov 2025.  
URL: https://arxiv.org/abs/2511.07517  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q25

*Finding.* Recalibration (DES-Dovekie) reduces the significance of rejecting LambdaCDM from 4.2 sigma (DES-SN5YR) to 3.2 sigma with CMB + DESI DR2. Only a weak Bayesian preference for w0wa remains.

*Relevance to HDBLAST.* The evolving-dark-energy signal is not settled, so it is not yet a firm constraint on any braneworld late-time sector.

**New DESI DR2 Lyman-alpha Results Shed Light on Dark Energy** (2026). DESI collaboration news, 30 July 2026.  
URL: https://www.desi.lbl.gov/2026/07/30/new-desi-dr2-lyman-alpha-results-shed-light-on-dark-energy/  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q47

*Finding.* The new DR2 Lyman-alpha forest measurements agree with the standard model. DESI notes the hints of evolving dark energy could fade, or a more complex model may be needed.

*Relevance to HDBLAST.* As above: late-time dark energy is currently outside the registered model's reach.

**Braneworld dark energy in light of DESI DR2** (2025). JCAP 11 (2025) 018.  
URL: https://arxiv.org/abs/2507.07193  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q26

*Finding.* Thawing scalar fields on a (4+1)-dimensional ghost-free 'phantom brane' (induced-gravity, DGP-type normal branch) show an effective phantom-divide crossing. A quartic potential explains the DESI DR2 constraints well.

*Relevance to HDBLAST.* A braneworld route to DESI-like behaviour exists, but it uses induced gravity on the brane and a cosmological crossover scale. The registered model is RS2-type (no induced gravity), so this does not carry over without new terms.

### LIGO/Virgo/KAGRA, GW170817 and gravitational-wave propagation

**Upper Limits on the Isotropic Gravitational-Wave Background from the first part of LIGO, Virgo and KAGRA's fourth Observing Run** (2025). Phys. Rev. D, doi 10.1103/wq57-sjt2 (update to April 2025: arXiv 2608.23477).  
URL: https://arxiv.org/abs/2508.20721  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q41

*Finding.* No background detected. Omega_GW(25 Hz) <= 2.0e-9 for a 2/3 power law and <= 2.8e-9 for a flat spectrum (95%). The update through April 2025 quotes the same values.

*Relevance to HDBLAST.* The frozen knee's f^-3 tail gives Omega(25 Hz) = 2.28e-35, 26 orders of magnitude below the limit (script), so the LVK band does not test the knee. If the brane scale is set by tabletop bounds, the static branch's natural horizon-scale frequency is about 5e-5 to 8e-5 Hz (script). That is below the LVK band, in the sub-mHz range targeted by space interferometers.

**GWTC-4.0: Tests of General Relativity. II. Parameterized Tests** (2026). arXiv Mar 2026.  
URL: https://arxiv.org/abs/2603.19020  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q42

*Finding.* 42 confident O4a signals plus 49 earlier events. Graviton-mass bound m_g < 1.92e-23 eV/c^2 (90%), 1.16 times better than GWTC-3.

*Relevance to HDBLAST.* In RS2 the 4D graviton is massless and the massive KK modes are micron-scale, so this bound does not touch the registered model at tabletop-consistent ell. It excludes cosmological-ell variants in which massive-graviton dispersion would appear.

**Limits on the number of spacetime dimensions from GW170817** (2018). JCAP 07 (2018) 048.  
URL: https://arxiv.org/abs/1801.08160  
*Access:* search-snippet-only · *Implication:* method · *Queries:* Q39

*Finding.* Gravitational-wave amplitude damping from leakage into large extra dimensions would bias GW distances; comparing GW and electromagnetic distances for GW170817 limits this and is consistent with D = 4.

*Relevance to HDBLAST.* Foundational method. It is relevant only for non-compact extra dimensions with a crossover scale not far above tens of Mpc.

**Constraining the number of spacetime dimensions from GWTC-3 binary black hole mergers** (2023). Phys. Rev. D 107, 084033 (2023).  
URL: https://arxiv.org/abs/2112.07650  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q40

*Finding.* Using the pair-instability mass gap: screening-independent D = 3.95 (+0.09/-0.07) (68%). With a screening scale: D = 4.23 (+1.50/-0.57) and log10 R_c/Mpc = 4.14 (+0.55/-0.86).

*Relevance to HDBLAST.* Leakage tests constrain HDBLAST only if the brane's effective AdS or crossover length is cosmological. The script's 'dark-energy identification' (ell about 372 Mpc) is in that regime, but tabletop tests already exclude it.

**Searching for Extra Dimensions with Gravitational Waves: Dark-Siren Constraints from GWTC-4** (2026). arXiv June 2026.  
URL: https://arxiv.org/abs/2606.14549  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q40

*Finding.* GWTC-4 dark sirens with a narrow H0 prior give D = 4.38 (+1.91/-1.01), consistent with D = 4. The crossover scale is poorly constrained, and the D limit weakens when the crossover exceeds the distances of the observed population.

*Relevance to HDBLAST.* Same role as the entry above: only large-crossover variants are tested.

**Constraint on the radius of five-dimensional dS spacetime with GW170817 and GRB 170817A** (2020). arXiv Jan 2020.  
URL: https://arxiv.org/abs/2001.06581  
*Access:* search-snippet-only · *Implication:* method · *Queries:* Q39

*Finding.* Title-level only: a constraint on a 5D de Sitter radius from the GW-gamma-ray arrival-time comparison. No number in the snippet.

*Relevance to HDBLAST.* A 'shortcut through the bulk' arrival-time test is the natural GW170817-type test for a de Sitter brane. It should be read before any HDBLAST propagation claim.

### Tabletop tests of Newton's inverse-square law

**New Test of the Gravitational 1/r^2 Law at Separations down to 52 micrometres** (2020). Phys. Rev. Lett. 124, 101101 (2020).  
URL: https://arxiv.org/abs/2002.11761  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q34, Q37

*Finding.* Torsion balance with 18- and 120-fold attractors at 52 um to 3.0 mm separations. Newtonian gravity fits; gravitational-strength Yukawa interactions are limited to ranges below 38.6 um (95%). Motivated by the dark-energy length of about 85 um.

*Relevance to HDBLAST.* Sets the physical unit of the registered model. The RS2 correction is a power law (1 + 2 ell^2/(3 r^2)), not a Yukawa, so 38.6 um is used as a proxy for the upper bound on the AdS radius ell. With ell <= 38.6 um: M5 >= 3.1e8 GeV, lambda^(1/4) >= 5.5 TeV, and the static branch's brane vacuum energy has rho^(1/4) >= 1.2 TeV (script; conditional). The 2026 review's power-law limits (below) are the correct input for a final number.

**Brane-World Gravity (Living Reviews in Relativity)** (unknown). Living Reviews in Relativity (PMC copy); volume/year not shown in the snippet.  
URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC5479361/  
*Access:* search-snippet-only · *Implication:* constrains · *Queries:* Q36

*Finding.* Table-top tests show no deviation from Newton's potential above 0.1 mm, so the AdS5 curvature scale ell < 0.1 mm (ell^-1 > 1e-12 GeV), giving M5 > 1e8 GeV.

*Relevance to HDBLAST.* Control: the script reproduces M5 = 2.3e8 GeV at ell = 0.1 mm with M_Pl^2 = M5^3 ell.

**Short-Range Tests of the Gravitational Inverse-Square Law** (2026). arXiv May 2026 (review).  
URL: https://arxiv.org/abs/2605.18212  
*Access:* search-snippet-only · *Implication:* method · *Queries:* Q35, Q38

*Finding.* A consistent formalism across length scales, updates from the past decade, and a comparison of tabletop and collider results for both Yukawa and power-law potentials, including extra-dimension models.

*Relevance to HDBLAST.* It contains the power-law (RS-type) limits needed to replace the Yukawa proxy for ell.

### The Hubble tension

**The Local Distance Network: a community consensus report on the measurement of the Hubble constant at 1% precision** (2025). A&A (April 2026).  
URL: https://arxiv.org/abs/2510.23823  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q43, Q46

*Finding.* H0 = 73.50 +/- 0.81 km/s/Mpc from a blinded community distance network (Cepheids, TRGB, Miras and more). This is 7.1 sigma from the early-universe LambdaCDM value 67.24 +/- 0.35. A simple local systematic is judged unlikely.

*Relevance to HDBLAST.* Extra dark radiation from a bulk remnant would raise the CMB-inferred H0 only slightly, and ACT DR6 disfavours that route. HDBLAST makes no H0 prediction. One practical effect: converting a strain amplitude into Omega_gw scales as h^-2. With h = 0.735 instead of 0.674, the knee's h_c at 1/yr would be 2.62e-15 instead of 2.40e-15 (script). The frozen Omega_k is tied to h = 0.674.

**JWST Observations Reject Unrecognized Crowding of Cepheid Photometry as an Explanation for the Hubble Tension at 8 sigma Confidence** (2024). ApJL (2024).  
URL: https://arxiv.org/abs/2401.04773  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q43

*Finding.* More than 1000 Cepheids measured with both HST and JWST: mean distance difference -0.01 +/- 0.03 mag. Crowding is rejected as the explanation at 8.2 sigma.

*Relevance to HDBLAST.* Context for the tension's robustness.

**The Chicago Carnegie Hubble Program: Improving the Calibration of SNe Ia with JWST Measurements of the Tip of the Red Giant Branch** (2025). arXiv Mar 2025.  
URL: https://arxiv.org/abs/2503.11769  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q43

*Finding.* The search summary quotes H0 = 69.8 +/- 1.9 from CCHP's JWST Cepheid/TRGB/JAGB work, consistent with the CMB value. It was not verified which CCHP paper this number is from.

*Relevance to HDBLAST.* The local-ladder disagreement is not fully settled, so HDBLAST should not be tuned to either H0.

*Attribution note.* Number attribution uncertain (see finding).

### JWST early galaxies

**Beyond No Tension: JWST z > 10 Galaxies Push Simulations to the Limit** (2025). arXiv Sep 2025.  
URL: https://arxiv.org/abs/2509.07695  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q44

*Finding.* Search summary: spectroscopic galaxies up to z of about 14.4 (MoM-z14) and 14.3 (GS-z14). Hints of an overabundance of bright or massive early galaxies; implied star-formation efficiencies are high. There is an ongoing debate over whether this reflects new physics or systematics (stellar masses, AGN, selection).

*Relevance to HDBLAST.* Tests small-scale primordial power and early structure growth. HDBLAST has no perturbation spectrum, so it is untested here. A future HDBLAST spectrum with enhanced small-scale power would be checked against these counts.

*Attribution note.* Galaxy redshifts came from the summary of several results.

**LambdaCDM not dead yet: massive high-z Balmer break galaxies are less common than previously reported** (2023). arXiv Oct 2023.  
URL: https://arxiv.org/abs/2310.03063  
*Access:* search-snippet-only · *Implication:* context · *Queries:* Q44

*Finding.* Title-level: massive high-z Balmer-break galaxies are less common than earlier claimed.

*Relevance to HDBLAST.* Early 'impossible galaxy' claims have weakened. They are not a motivation for a non-standard origin.

## 6. Limitations

- Snippet-level evidence only. Snippets can garble numbers or merge several papers, and the entries flag the cases noticed. Nothing here substitutes for reading the papers.
- Search coverage is incomplete, and the search budget ran out: three planned queries (IPTA DR3 status, the 'What is the source of the PTA GW signal' paper, EPTA DR2 implications) were refused. Absence of a 2026 IPTA DR3 or NANOGrav 20-yr background result is 'not found', not 'does not exist'.
- The knee comparisons in §3 are descriptive, not likelihood analyses. They use published power-law summaries, not the free-spectrum posteriors.
- The scale-setting in §4 is order of magnitude and conditional on RS2 low-energy relations and on the thermalisation assumptions listed there. The registered model has not produced a radiation era.
- The Gaussian ΔN_eff upper limits are built from central values ± 1σ. They are not the collaborations' own 95% limits.

## 7. Search log

- Q1: `NANOGrav 15-year gravitational wave background spectral turnover low frequency environmental effects`
- Q2: `IPTA third data release DR3 gravitational wave background 2026`
- Q3: `MeerKAT Pulsar Timing Array 4.5 year data gravitational wave background amplitude spectral index`
- Q4: `"Stochastic gravitational-wave background search using data from five pulsar timing arrays"`
- Q5: `NANOGrav 15 yr Data Set Running of the Spectral Index result inconclusive`
- Q6: `MeerKAT pulsar timing array first search gravitational waves common signal amplitude larger than other PTAs frequency`
- Q7: `EPTA DR2new InPTA gravitational wave background amplitude spectral index 10.3 years evidence`
- Q8: `Comparing recent pulsar timing array results on the nanohertz stochastic gravitational-wave background IPTA 2024 consistent`
- Q9: `Parkes Pulsar Timing Array DR3 isotropic gravitational wave background amplitude 13/3 Hellings-Downs significance`
- Q10: `pulsar timing array sensitivity loss at 1/yr frequency pulsar position fit timing model hasasia sensitivity curve`
- Q11: `NANOGrav 15-year evidence gravitational wave background 16.03 yr 14 frequencies amplitude 2.4e-15 gamma 3.2`
- Q12: `NANOGrav 20-year data set release 2026`
- Q13: `"Piecewise Power-Law Reconstruction" NANOGrav 15 yr gravitational-wave background`
- Q14: `NANOGrav 15 yr erratum supermassive black hole binaries 2026 software bug gravitational wave background`
- Q15: `pulsar timing array spectral break knee broken power law gravitational wave background 2025 search`
- Q16: `NANOGrav piecewise power law 2601.09481 broken power law Bayes factor spectral shape result high frequency`
- Q17: `EPTA DR2 improved noise model averaging gravitational wave background lower amplitude strain spectral index -2/3 Max Planck 2025`
- Q18: `"Impacts of Customized Chromatic Noise Models" NANOGrav 15 yr gravitational wave`
- Q19: `"Pulsar timing arrays: the emerging gravitational-wave landscape"`
- Q20: `Chinese Pulsar Timing Array DR1 gravitational wave background 14 nHz 4.6 sigma FAST`
- Q21: `NANOGrav 15 year search for signals from new physics cosmological interpretation Bayes factor 10 to 100 phase transition cosmic strings`
- Q22: `pulsar timing array stochastic background light Kaluza-Klein resonances warped extra dimension phase transition`
- Q23: `causal gravitational wave infrared tail f^3 radiation domination universal low-frequency spectrum cosmological source`
- Q24: `DESI DR2 BAO evolving dark energy w0 wa significance 2025 results`
- Q25: `DES supernova recalibration Dovekie 2025 2026 evolving dark energy significance reduced DESI`
- Q26: `braneworld dark energy DESI DR2 phantom crossing DGP brane fit 2025`
- Q27: `ACT DR6 2025 N_eff constraint extended models dark radiation Hubble tension`
- Q28: `SPT-3G D1 2025 cosmological constraints N_eff Hubble constant CMB-SPA combined`
- Q29: `tensor-to-scalar ratio upper limit 2025 2026 BICEP Keck SPT delensing r < 0.03`
- Q30: `ACT DR6 power spectra cosmological parameters spectral index n_s 0.974 P-ACT-LB inflation implications`
- Q31: `"Inflation at the End of 2025" constraints on r and n_s latest CMB BAO data result`
- Q32: `big bang nucleosynthesis 2024 2025 N_eff constraint deuterium helium primordial abundance ΔN_eff bound`
- Q33: `braneworld dark radiation Weyl term bulk black hole mass constraint BBN CMB Randall-Sundrum brane tension bound`
- Q34: `Eöt-Wash 2020 torsion balance inverse-square law test 52 micrometers extra dimension size bound`
- Q35: `short-range gravity inverse square law test 2024 2025 new limit Yukawa micrometer experiment`
- Q36: `Randall-Sundrum AdS curvature radius bound table-top gravity experiments brane tension TeV fundamental scale 10^8 GeV`
- Q37: `"New Test of the Gravitational 1/r^2 Law at Separations down to 52" Yukawa |alpha|=1 38.6 micrometers dark energy length scale`
- Q38: `"Short-Range Tests of the Gravitational Inverse-Square Law" 2026 review extra dimensions limits`
- Q39: `GW170817 constraint number of spacetime dimensions gravitational wave leakage luminosity distance damping extra dimensions`
- Q40: `gravitational wave leakage extra dimensions constraint dark sirens GWTC-3 screening scale 2023 2024 higher-dimensional`
- Q41: `LIGO Virgo KAGRA O4 upper limits isotropic gravitational-wave background 2025 Omega_GW 25 Hz`
- Q42: `GWTC-4 tests of general relativity 2025 LIGO Virgo KAGRA graviton mass dispersion O4a`
- Q43: `Hubble tension 2026 status SH0ES JWST Cepheid H0 73 CCHP Freedman 70 TRGB JAGB local distance network`
- Q44: `JWST early massive galaxies z>10 overabundance LCDM tension 2025 review bright galaxies star formation efficiency`
- Q45: `Randall-Sundrum high-energy brane inflation tensor amplitude enhancement tensor-to-scalar ratio braneworld consistency relation observational constraint`
- Q46: `Local Distance Network H0DN 2025 1% Hubble constant 73.50 consensus distance ladder paper`
- Q47: `DESI 2026 dark energy update full-shape DR2 evolving dark energy new result`
- Q48: `gravitational wave peak frequency today temperature at production 1.65e-7 Hz T/GeV g* nanohertz QCD epoch`
- refused (budget): `IPTA DR3 combined dataset status 2026 International Pulsar Timing Array combination early results`
- refused (budget): `"What is the source of the PTA GW signal" supermassive black hole binaries environmental effects cosmological models comparison`
- refused (budget): `EPTA second data release implications for massive black holes dark matter and the early Universe`

