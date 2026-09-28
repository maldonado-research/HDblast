# Literature 2022–2026 for HDBLAST: merged, deduplicated and annotated

Checkpoint of 27 September 2026 (searches run 28 September 2026). Built by `synthesis/build_literature.py` from the three sweep records in `literature/`. Every URL below was returned by a web search in this round or is written in the programme's own archive files; none was added by hand.

**Read this first: evidence level.** Paper sites (arxiv.org, zenodo.org, journal sites) were blocked for downloads, so **no paper was opened**. Each entry is based on a search-result snippet (`search-snippet-only`) or on an earlier programme record (`programme-record-only`). Snippets can misattribute authors, venues or numbers. Before any external use, open the paper and check. No independent auditor re-checked the literature sweeps in this round; the stability and mechanism audits could not re-verify citations either, because the shared web-search budget ran out. The methods sweep ran **no** new searches (budget exhausted), so its entries are second-hand.

**Counts.** 188 records in the three sweeps (95 theory, 51 observations, 42 methods) merge to **162 distinct sources** (26 duplicates merged). Of these, **82 are dated 2022–2026**, 3 have no year in the snippet, and 77 are earlier foundational work (listed briefly in the last section).

**Implication labels** (the sweeps' own judgement, relative to HDBLAST): *challenges* = could undercut the hypothesis or the matter extension; *constrains* = sets a bound any HDBLAST version must satisfy; *supports* = consistent with, or provides a mechanism HDBLAST could use (never observational support for HDBLAST itself); *method* = a tool we can use; *context* = prior art or background.

## What the 2022–2026 literature implies for HDBLAST (summary)

1. **Prior art is substantial.** "Our universe as a wall created in a 5D event" is an active published line (dark-bubble cosmology, black-hole-catalysed nucleation, detuned-brane scenarios). HDBLAST must not claim the concept; its own content is the registered model and its computed results.
2. **Every published hot-origin model adds an ingredient** (a bulk black hole, colliding walls, a holographic hot sector, or matter on the wall). A pure-tension wall produces no heat. HDBLAST's own negative results this round (preheating and mechanism screens) agree with that pattern.
3. **Observations constrain, not test, the current model.** The registered 5D model has no radiation era and no perturbation spectrum, so CMB, BBN, JWST, DESI and H0 data are requirements for the future, not tests now. Dark radiation (N_eff) is the sharpest direct bound on any 5D remnant (bulk Weyl term).
4. **The PTA knee is untested at likelihood level** and sits at 0.997×(1/yr), on the pulsar-position/proper-motion fitting blind spot (standard in PTA work, new to this project). 2025–2026 noise reanalyses trend towards γ=13/3, which is unfavourable to the frozen low-frequency slope but is not a rejection.
5. **Most useful new tools:** the Einstein–scalar de Sitter-sliced flow programme (Kiritsis, Nitti and collaborators), dark-bubble junction/perturbation methods, and constraint-damping practice from numerical relativity.

Known inconsistencies between sweep snippets (unresolved, paper not opened): the ACT DR6 N_eff value is quoted as 2.86 ± 0.13 in the observations sweep and 2.89 ± 0.11 in the mechanisms workstream; dark-radiation allowances are quoted in different normalisations (ρ_dr/ρ_γ vs ρ_DR/ρ_SM at production) and are not directly comparable.

## Annotated sources dated 2022–2026 (and undated), grouped by implication

### Challenges (4)

- **Pulsar timing arrays: the emerging gravitational-wave landscape** (2026). arXiv Mar 2026 (review). <https://arxiv.org/abs/2603.13643>
  - *Finding (snippet):* Six PTAs (NANOGrav, EPTA, InPTA, PPTA, CPTA, MPTA) report evidence for the background. The review argues that the perceived tension between current amplitudes and standard merger models is largely resolved by new insights into SMBHB populations.
  - *For HDBLAST:* [observations] If the SMBHB amplitude tension is resolved, a main motivation for an exotic source of the hum weakens.
  - *Labels:* challenges; topic: PTA; access: search-snippet-only; records: observations:PTA_review_2026
- **Dark bubble cosmology and the equivalence principle** — I. Basile, A. Borys, J. Masias (2025). arXiv:2507.03748; Phys. Rev. D (2026). <https://arxiv.org/abs/2507.03748>
  - *Finding (snippet):* Couples the electroweak and strong sectors to the induced braneworld gravity by the dark-bubble mechanism. The electroweak sector is unaffected, but the proton's gravitational and inertial masses differ significantly, severely violating equivalence-principle measurements.
  - *For HDBLAST:* [theory] A sharp negative test that any braneworld with bulk-mediated matter couplings must pass. HDBLAST's proposed matter field chi has a phi-dependent mass, so the bulk scalar couples directly to brane matter; an equivalence-principle and fifth-force check is required before that extension can be called viable.
  - *Labels:* challenges; topic: dark bubbles; access: search-snippet-only; records: theory:A10
- **Dynamical dark energy in 0'B braneworlds** — I. Basile and two co-authors (snippet) (2025). Eur. Phys. J. C (Sept 2025); arXiv:2502.20438. <https://arxiv.org/abs/2502.20438>
  - *Finding (snippet):* Builds a dark-bubble realization in non-supersymmetric type 0'B string theory, the unique scale-separated option among the simplest models. The dark energy varies logarithmically. The blue-shifted spectrum excludes inflation, the late-time predictions conflict with Standard Model couplings, and the model appears not to be phenomenologically viable.
  - *For HDBLAST:* [theory] A concrete example of a 5D-wall cosmology failing when fully specified, published as a negative result. It sets the standard of candour HDBLAST follows.
  - *Labels:* challenges; topic: dark bubbles; access: search-snippet-only; records: theory:A14
- **A clearer view of gravitational-wave signals in pulsar timing arrays (AEI news on improved EPTA DR2 noise models)** (2025). Max Planck Institute for Gravitational Physics news item. <https://www.aei.mpg.de/1323827/a-a-clearer-view-of-gravitational-wave-signals-in-pulsar-timing-arrays>
  - *Finding (snippet):* Improved EPTA DR2 noise modelling (pulsar-intrinsic, epoch-correlated and transient noise, with noise-model averaging) makes the background's strain index consistent with -2/3 (gamma = 13/3). It gives a lower amplitude at 1/yr than earlier analyses. Chromatic noise errors have whiter spectra, which push naive fits toward flatter indices.
  - *For HDBLAST:* [observations] Two independent 2025-2026 reanalyses (EPTA; NANOGrav, entry above) move toward gamma = 13/3 and lower amplitude. The frozen curve's low-frequency branch is Omega ~ f^2 (gamma = 3). This trend is unfavourable to the frozen shape, but it is not a likelihood-level rejection.
  - *Labels:* challenges; topic: PTA; access: search-snippet-only; records: observations:AEI_EPTA_noise

### Constrains (18)

- **The NANOGrav 15 yr Data Set: Impacts of Customized Chromatic Noise Models on Gravitational Wave Analyses** (2026). arXiv June 2026 (companion: arXiv 2606.28571). <https://arxiv.org/abs/2606.28554>
  - *Finding (snippet):* Customized chromatic (interstellar-medium) noise models for the 15-yr pulsars. Bayes factor for Hellings-Downs correlations over an uncorrelated common red process: 1571 +/- 14 (14 Fourier modes), about 8 times earlier. Power-law amplitude at fixed index reduced to 2.1 (+0.6/-0.5) x 10^-15.
  - *For HDBLAST:* [observations] Relative to this amplitude, the frozen knee has 1.31 times the power at 1/yr (1.14 times the strain). That is inside the quoted uncertainty. Better noise modelling lowered the amplitude and strengthened the correlation signature. The knee's height was calibrated to the 2023 value, and the direction of change matters for any refreeze.
  - *Labels:* constrains; topic: PTA; access: search-snippet-only; records: observations:NG15_CNM
- **Searching for Extra Dimensions with Gravitational Waves: Dark-Siren Constraints from GWTC-4** (2026). arXiv June 2026. <https://arxiv.org/abs/2606.14549>
  - *Finding (snippet):* GWTC-4 dark sirens with a narrow H0 prior give D = 4.38 (+1.91/-1.01), consistent with D = 4. The crossover scale is poorly constrained, and the D limit weakens when the crossover exceeds the distances of the observed population.
  - *For HDBLAST:* [observations] Same role as the entry above: only large-crossover variants are tested.
  - *Labels:* constrains; topic: LVK; access: search-snippet-only; records: observations:GWTC4_darksiren_dims
- **Massive Graviton Dark Matter from a Gapped Continuum** — Megias, Prieto and Quiros (2026). arXiv:2607.07295 (v2 3 Aug 2026). <https://arxiv.org/html/2607.07295v2>
  - *Finding (snippet):* Per the programme record: in a linear-dilaton braneworld, brane-to-bulk graviton leakage scales as the eighth power of temperature, and the hidden fluid later scales as matter.
  - *For HDBLAST:* [methods] consistent backreaction must include energy leaving the brane into the bulk. Not every bulk channel acts as dark radiation.
  - *Labels:* constrains; topic: particle production / backreaction; access: programme-record-only; records: methods:D6
- **On the Origin and Fate of Our Universe** — C. Vafa (2025). Gen. Relativ. Gravit. (2025); arXiv:2501.00966. <https://arxiv.org/abs/2501.00966>
  - *Finding (snippet):* A short review of swampland bounds on positive potentials, the de Sitter conjecture, TCC and its relation to the species scale, with implications for inflation and the fate of the universe.
  - *For HDBLAST:* [theory] The current mainstream quantum-gravity view on origins: any HDBLAST claim about an inflating or de Sitter origin should be stated against it.
  - *Labels:* constrains; topic: swampland / strings; access: search-snippet-only; records: theory:G4
- **The NANOGrav 15 yr Data Set: Running of the Spectral Index** (2025). ApJL 978, L29 (2025); arXiv Aug 2024. <https://arxiv.org/abs/2408.10166>
  - *Finding (snippet):* Running-power-law fit (logarithmic frequency dependence of the index). 95% credible interval for the running beta in [-0.80, 2.96], consistent with zero. Bayes factor running vs constant power law 0.69 +/- 0.01 (inconclusive). The constant power law still suffices.
  - *For HDBLAST:* [observations] A knee is an extreme form of spectral curvature. Within the NANOGrav band the frozen curve is almost a pure power law (curvature only in bins 13-14), so this null result neither supports nor excludes it. It does show that current data do not require curvature.
  - *Labels:* constrains; topic: PTA; access: search-snippet-only; records: observations:NG15_running
- **A sensitivity curve approach to tuning a pulsar timing array in the detection era (with the hasasia package)** (2025). CQG (2025), doi 10.1088/1361-6382/adbbab. <https://arxiv.org/abs/2409.00336>
  - *Finding (snippet):* Search summary: fitting each pulsar's sky position and proper motion causes an extra loss of sensitivity around f = 1/yr, and fitting parallax causes one around 2/yr. These appear as large spikes in PTA sensitivity curves (computed with hasasia).
  - *For HDBLAST:* [observations] This is the key methodological fact for the registered knee. The script finds f_k = 0.9972 x (1/yr), within 0.045 bins of 1/yr for NANOGrav 15-yr and 0.013 bins for MeerKAT. A spectral bend exactly where the timing model absorbs power is maximally degenerate with it. A knee test must propagate the timing-model transmission function, or the knee location is not identifiable.
  - *Labels:* constrains; topic: PTA; access: search-snippet-only; records: observations:PTA_sensitivity_1yr
- **The Atacama Cosmology Telescope: DR6 Constraints on Extended Cosmological Models** (2025). arXiv Mar 2025. <https://arxiv.org/abs/2503.14454>
  - *Finding (snippet):* No evidence for new free-streaming light species: N_eff = 2.86 +/- 0.13. Self-interacting dark radiation N < 0.134. Dark-radiation models are not favoured, because extra radiation increases Silk damping at high multipoles, where ACT DR6 is sensitive.
  - *For HDBLAST:* [observations] In a braneworld the Weyl ('dark radiation') term C/a^4, set by a bulk black-hole mass, counts as Delta N_eff. A Gaussian approximation from these numbers gives Delta N_eff < 0.071, so rho_dr/rho_gamma < 0.016 at recombination (script). This bounds any bulk remnant of the blast, for example the collapse fate's black brane, and any GW relic, if a radiation era is ever achieved.
  - *Labels:* constrains; topic: CMB; access: search-snippet-only; records: observations:ACT_DR6_ext
- **ACT DR6 Insights on alpha-Attractor Inflationary Models and Reheating** (2025). arXiv May 2025. <https://arxiv.org/abs/2505.01517>
  - *Finding (snippet):* As reported there: Planck 2018 + ACT DR6 gives n_s = 0.9709 +/- 0.0038, rising to 0.9743 +/- 0.0034 with DESI Y1 BAO (the P-ACT-LB combination gives n_s = 0.974 +/- 0.003). Starobinsky-type plateau models become mildly disfavoured at about 2 sigma or more.
  - *For HDBLAST:* [observations] Any 'blast instead of (or before) inflation' model must produce adiabatic, nearly scale-invariant perturbations with n_s about 0.97. The registered model has no perturbation-spectrum prediction yet, so this is currently an open requirement, not a test it passes or fails.
  - *Labels:* constrains; topic: CMB; access: search-snippet-only; records: observations:ACT_ns
- **Inflation at the End of 2025: Constraints on r and n_s Using the Latest CMB and BAO Data** (2025). Open Journal of Astrophysics (2026); arXiv Dec 2025. <https://arxiv.org/abs/2512.10613>
  - *Finding (snippet):* Planck + SPT + ACT + BICEP/Keck: n_s = 0.9682 +/- 0.0032 and r < 0.034 (95%). Adding DESI BAO leaves r unchanged but shifts n_s to 0.9728 +/- 0.0029, driven by marginal CMB-DESI differences.
  - *For HDBLAST:* [observations] Primordial tensors: any inflation-like phase on the brane is bounded by r < 0.034. In high-energy RS inflation, tensor and scalar amplitudes are both enhanced (next entries). HDBLAST has no tensor prediction yet.
  - *Labels:* constrains; topic: CMB; access: search-snippet-only; records: observations:r_ns_2025
- **Upper Limits on the Isotropic Gravitational-Wave Background from the first part of LIGO, Virgo and KAGRA's fourth Observing Run** (2025). Phys. Rev. D, doi 10.1103/wq57-sjt2 (update to April 2025: arXiv 2608.23477). <https://arxiv.org/abs/2508.20721>
  - *Finding (snippet):* No background detected. Omega_GW(25 Hz) <= 2.0e-9 for a 2/3 power law and <= 2.8e-9 for a flat spectrum (95%). The update through April 2025 quotes the same values.
  - *For HDBLAST:* [observations] The frozen knee's f^-3 tail gives Omega(25 Hz) = 2.28e-35, 26 orders of magnitude below the limit (script), so the LVK band does not test the knee. If the brane scale is set by tabletop bounds, the static branch's natural horizon-scale frequency is about 5e-5 to 8e-5 Hz (script). That is below the LVK band, in the sub-mHz range targeted by space interferometers.
  - *Labels:* constrains; topic: LVK; access: search-snippet-only; records: observations:LVK_O4a_SGWB
- **End-of-the-World Branes and Inflationary Predictions for Rocky and Swampy Landscapes** — B. Hassfeld, A. Hebecker, A. Westphal (2024). JHEP 03 (2025) 196; arXiv:2411.11944. <https://arxiv.org/abs/2411.11944>
  - *Finding (snippet):* A framework for anthropic predictions assuming de Sitter vacua have finite-dimensional Hilbert spaces. Even with eternal inflation, predictions depend on the rates of creating universes from nothing, and these rates are highly sensitive to the existence of end-of-the-world branes. Distinguishes 'swampy' from 'rocky' landscapes.
  - *For HDBLAST:* [theory] Any claim that a 5D event is a typical origin of our universe depends on creation-rate measures. HDBLAST has not addressed measures.
  - *Labels:* constrains; topic: creation / nucleation; access: search-snippet-only; records: theory:B18
- **The MeerKAT Pulsar Timing Array: The 4.5-year data release and the noise and stochastic signals of the millisecond pulsar population** (2024). arXiv Dec 2024. <https://arxiv.org/abs/2412.01148>
  - *Finding (snippet):* 83 pulsars, 4.5 yr, high cadence. Common signal log10 A = -14.25 (+0.21/-0.36), gamma = 3.60 (+1.31/-0.89). At gamma = 13/3, log10 A = -14.28 +/- 0.21, ln Bayes factor 4.46. The amplitude is larger than other PTAs report.
  - *For HDBLAST:* [observations] MeerKAT's free-gamma posterior still allows gamma = 3. Its resolution (7.04 nHz) puts the knee at bin 4.49, so its bins straddle 1/yr more coarsely than NANOGrav's.
  - *Labels:* constrains; topic: PTA; access: search-snippet-only; records: observations:MPTA_data
- **The MeerKAT Pulsar Timing Array: The first search for gravitational waves with the MeerKAT radio telescope** (2024). MNRAS 536, 1489 (2025). <https://arxiv.org/abs/2412.01153>
  - *Finding (snippet):* Sky-averaged h_c,yr = 7.5 (+0.8/-0.9) x 10^-15 at strain index -0.26, or 4.8 (+0.8/-0.9) x 10^-15 at -2/3. This is inconsistent with other PTAs' common-noise results by at least about 1.4 sigma. Hellings-Downs significance is 3-3.4 sigma depending on noise assumptions.
  - *For HDBLAST:* [observations] At 1/yr the frozen knee has only 0.25 (gamma = 13/3) to 0.10 (free index) of MeerKAT's power. The cross-PTA amplitude spread already exceeds a factor 4 in Omega. A single frozen height cannot match all teams, which is one of the project's own failure criteria (teams disagreeing).
  - *Labels:* constrains; topic: PTA; access: search-snippet-only; records: observations:MPTA_search
- **The 2024 BBN baryon abundance update** (2024). arXiv Jan 2024. <https://arxiv.org/abs/2401.15054>
  - *Finding (snippet):* Delta N_eff from BBN ranges from -0.09 +/- 0.28 (one configuration: BAO+BBN abundances, PArthENoPE v3.0); results depend on the helium data and deuterium-burning rates.
  - *For HDBLAST:* [observations] Expansion rate at T ~ 1 MeV. For a braneworld this bounds the rho^2/(2 lambda) term, which needs lambda^(1/4) far above MeV. That is automatic if tabletop bounds set lambda^(1/4) of a few TeV. It also bounds dark radiation: a Gaussian approximation gives Delta N_eff < 0.46, so rho_dr/rho_gamma < 0.10 (script).
  - *Labels:* constrains; topic: BBN; access: search-snippet-only; records: observations:BBN_2024
- **Gravitational reheating at strong coupling** — A. Buchel (2023). JHEP 07 (2023) 159; arXiv:2304.11195. <https://arxiv.org/abs/2304.11195>
  - *Finding (snippet):* Uses gauge/gravity duality to estimate the maximal reheating temperature of strongly coupled theories after a rapid exit from de Sitter; reheating is most efficient when H is much larger than the conformal-breaking scale and the breaking operators are nearly marginal.
  - *For HDBLAST:* [theory] HDBLAST's shell moves between de Sitter rates (H/H0 falling from about 0.64 to 0.607). A bound of this type could say how much gravitational particle production such a mild change can give, probably very little.
  - *For HDBLAST:* [methods] an upper-bound method for gravitational production during a change of de Sitter rate.
  - *Labels:* method, constrains; topic: dS-sliced flows / holographic reheating, particle production / backreaction; access: search-snippet-only; records: theory:F7, methods:D11
- **The second data release from the European Pulsar Timing Array III. Search for gravitational wave signals** (2023). A&A 678, A50 (2023). <https://arxiv.org/abs/2306.16214>
  - *Finding (snippet):* DR2new (latest 10.3 yr, 25 pulsars, plus about 3.5 yr of InPTA data for 10 of them): Bayes factor 60, false-alarm probability about 0.1% (at least 3 sigma). At fixed index 13/3: A = (2.5 +/- 0.7) x 10^-15 at 1/yr. Full DR2, HD process: log10 A = -14.54 (+0.28/-0.41), gamma = 4.19 (+0.73/-0.63).
  - *For HDBLAST:* [observations] The frozen knee lies within 0.03 frequency bins of 1/yr for DR2new (resolution 3.08 nHz). The EPTA free-index fit sits near gamma = 13/3, not the gamma = 3 that the frozen curve follows below the knee.
  - *Labels:* constrains; topic: PTA; access: search-snippet-only; records: observations:EPTA_DR2_III
- **Constraining the number of spacetime dimensions from GWTC-3 binary black hole mergers** (2023). Phys. Rev. D 107, 084033 (2023). <https://arxiv.org/abs/2112.07650>
  - *Finding (snippet):* Using the pair-instability mass gap: screening-independent D = 3.95 (+0.09/-0.07) (68%). With a screening scale: D = 4.23 (+1.50/-0.57) and log10 R_c/Mpc = 4.14 (+0.55/-0.86).
  - *For HDBLAST:* [observations] Leakage tests constrain HDBLAST only if the brane's effective AdS or crossover length is cosmological. The script's 'dark-energy identification' (ell about 372 Mpc) is in that regime, but tabletop tests already exclude it.
  - *Labels:* constrains; topic: LVK; access: search-snippet-only; records: observations:GWTC3_dims
- **Brane-World Gravity (Living Reviews in Relativity)** (year not given). Living Reviews in Relativity (PMC copy); volume/year not shown in the snippet. <https://pmc.ncbi.nlm.nih.gov/articles/PMC5479361/>
  - *Finding (snippet):* Table-top tests show no deviation from Newton's potential above 0.1 mm, so the AdS5 curvature scale ell < 0.1 mm (ell^-1 > 1e-12 GeV), giving M5 > 1e8 GeV.
  - *For HDBLAST:* [observations] Control: the script reproduces M5 = 2.3e8 GeV at ell = 0.1 mm with M_Pl^2 = M5^3 ell.
  - *Labels:* constrains; topic: ISL; access: search-snippet-only; records: observations:LivingReview_RS

### Supports (2)

- **Inflation with a Growing Fifth Dimension** — Harvard group (snippet gives affiliation only) (2025). JHEP 05 (2026) 169; arXiv:2512.04177. <https://arxiv.org/abs/2512.04177>
  - *Finding (snippet):* Inflation with a finite initial time in warped AdS5 with UV and IR branes. The inflaton potential detunes the brane tension, so the fifth dimension grows: a two-field (inflaton plus radion) hyperbolic model with early radion fast-roll and late inflaton slow-roll. The earliest modes are radion-sourced and give a suppressed, blue-tilted scalar spectrum and oscillatory tensors, in principle visible in the CMB. The model links to the dark dimension.
  - *For HDBLAST:* [theory] Structurally the closest 2025-2026 match to HDBLAST's core mechanism: a detuned brane tension drives the dynamics of the extra dimension. It shows how such a mechanism yields observables (large-scale suppression), and is a template for HDBLAST open item 5.
  - *Labels:* supports; topic: swampland / strings; access: search-snippet-only; records: theory:G6
- **A dynamical inflaton coupled to strongly interacting matter** — C. Ecker, E. Kiritsis, W. van der Schee (2023). Phys. Rev. Lett. 130, 251001; arXiv:2302.06618. <https://arxiv.org/abs/2302.06618>
  - *Finding (snippet):* Self-consistently couples the Einstein-inflaton equations to a strongly coupled QFT described holographically, and finds inflation, a reheating phase, and finally a universe dominated by the QFT in thermal equilibrium.
  - *For HDBLAST:* [theory] A controlled demonstration of reheating into a hot (holographic) sector, the step HDBLAST has not achieved. It suggests treating part of the 5D bulk as the hot sector's dual.
  - *For HDBLAST:* [methods] the only item in this session's records that computes reheating with full self-consistent backreaction; the heat bath is a 5D (holographic) sector. Methodologically the closest template for bulk-coupled feedback, which this round's preheating workstream did not solve.
  - *Labels:* supports, method; topic: dS-sliced flows / holographic reheating, particle production / backreaction; access: search-snippet-only; records: theory:F6, methods:D10

### Method (24)

- **Self-gravitating electromagnetic waves in the dark bubble model** — U. Danielsson and one co-author (snippet) (2026). arXiv:2606.16547 (15 June 2026). <https://arxiv.org/abs/2606.16547>
  - *Finding (snippet):* Embeds gravitational and electromagnetic waves using AdS5 pp-wave geometries glued across a three-brane. For localized light beams under mixed AdS5 boundary conditions, the gravitational corrections are consistent with 4D gravity weakening at the 5D AdS scale. On the brane, electromagnetic waves source the Kalb-Ramond B-field in the bulk.
  - *For HDBLAST:* [theory] Exact wall-plus-radiation solutions of this type could serve as analytic benchmarks for any HDBLAST radiation-on-shell calculation.
  - *Labels:* method; topic: dark bubbles; access: search-snippet-only; records: theory:A13
- **Primordial Correlators from a Kaluza-Klein Graviton Continuum** — not captured in snippet (2026). arXiv:2608.01762 (3 Aug 2026). <https://arxiv.org/abs/2608.01762>
  - *Finding (snippet):* A spectral representation for inflationary correlators mediated by a continuum of states, realized in an RS2-like braneworld with the inflaton on a de Sitter brane in AdS5. The tensor sector contains a localized zero mode and a KK graviton continuum with a gap of order H; the query-37 snippet puts the threshold at 3H/2, the principal-series edge.
  - *For HDBLAST:* [theory] HDBLAST's certified 'no scalar bound state below 9H^2/4' (Chat 10) uses exactly this threshold, (3H/2)^2. This paper shows how such a continuum could become an observable (cosmological-collider-like signals), if HDBLAST ever has an inflating stage with perturbations.
  - *Labels:* method; topic: stability of dS walls; access: search-snippet-only; records: theory:E10
- **Complex scalar field thick branes: stability of linear perturbation and evolution of scalar Kaluza-Klein modes coupled with gravity** — not captured in snippet (2026). arXiv:2609.21421 (18 Sept 2026). <https://arxiv.org/abs/2609.21421>
  - *Finding (snippet):* Minkowski, de Sitter and AdS thick branes from a complex scalar show no instability in the scalar and vector sectors. Factorizing the tensor equation excludes tachyonic tensor modes. The scalar zero mode is localized on Minkowski and de Sitter branes, with a gapless continuous scalar KK spectrum.
  - *For HDBLAST:* [theory] A current example of a sector-by-sector stability analysis of de Sitter branes, comparable in scope to HDBLAST Chat 12.
  - *For HDBLAST:* [methods] a current sector-by-sector template, including the vector sector, which this round's stability workstreams did not treat.
  - *Labels:* method; topic: stability of dS walls, gauge-invariant perturbations; access: search-snippet-only; records: theory:E11, methods:C9
- **Effective Dynamics of Inflationary End-of-the-World Branes in AdS3** — K. Fujiki and others (snippet) (2026). arXiv:2609.11643 (Sept 2026). <https://arxiv.org/abs/2609.11643>
  - *Finding (snippet):* Integrates out the bulk to obtain a Liouville-like effective theory for a 2D cosmological EOW brane with a localized scalar; constructs slow-roll inflating trajectories, regular Euclidean geometries continuing to Lorentzian de Sitter, and the semiclassical on-shell action.
  - *For HDBLAST:* [theory] A method for deriving a brane-scalar effective action and a Euclidean continuation in a lower-dimensional toy model, comparable to HDBLAST's closed-form 4D effective theory (Chat 9).
  - *Labels:* method; topic: dS-sliced flows / holographic reheating; access: search-snippet-only; records: theory:F10
- **The NANOGrav 15 yr Data Set: Piecewise Power-Law Reconstruction of the Gravitational-Wave Background** (2026). ApJL (2026), doi 10.3847/2041-8213/ae7086. <https://arxiv.org/abs/2601.09481>
  - *Finding (snippet):* Piecewise power-law (PPL) spectral reconstruction: constant, broken, doubly broken models combined by Bayesian model averaging. Described as closer to physically realistic (especially cosmological) spectra than the free spectrum. The snippets did not report the fitted break positions or Bayes factors.
  - *For HDBLAST:* [observations] This is the collaboration's own analogue of the project's knee test. The frozen SBPL should be compared with the PPL posterior (bin-by-bin power and break-frequency posterior), not only with summary points. Result numbers must be read from the paper before any claim.
  - *Labels:* method; topic: PTA; access: search-snippet-only; records: observations:NG15_PPL
- **KDE Representations of the Gravitational Wave Background Free Spectra Present in the NANOGrav 15-Year Dataset (corrected release; erratum Agazie et al. 2026, ApJL 1006, L67)** (2026). Zenodo data record + ApJL erratum (2026). <https://zenodo.org/records/21844115>
  - *Finding (snippet):* A bug in the parallel-tempering routine of PTMCMCSampler biased the free-spectrum posteriors of the 15-yr data set. The record re-releases corrected free-spectrum KDEs. The erratum states SMBHB population inference was affected. The project's own site notes the SMBHB conclusion was unchanged.
  - *For HDBLAST:* [observations] Any project screen of the knee that used pre-2026 NANOGrav free-spectrum posteriors or 'summary points' derived from them should be rerun on the corrected KDEs. The frozen numbers must stay frozen.
  - *Labels:* method; topic: PTA; access: search-snippet-only; records: observations:NG15_erratum_KDE
- **A Universal CMB B-Mode Spectrum from Early Causal Tensor Sources** (2026). arXiv Jan 2026. <https://arxiv.org/abs/2601.20967>
  - *Finding (snippet):* Early causal tensor sources share a universal infrared scaling and predict the same B-mode angular distribution.
  - *For HDBLAST:* [observations] If the blast is a causal (sub-horizon) tensor source, its CMB B-mode imprint is fixed up to amplitude. This is a cross-check between a PTA-band claim and CMB tensor limits.
  - *Labels:* method; topic: PTA; access: search-snippet-only; records: observations:causal_Bmode
- **Short-Range Tests of the Gravitational Inverse-Square Law** (2026). arXiv May 2026 (review). <https://arxiv.org/abs/2605.18212>
  - *Finding (snippet):* A consistent formalism across length scales, updates from the past decade, and a comparison of tabletop and collider results for both Yukawa and power-law potentials, including extra-dimension models.
  - *For HDBLAST:* [observations] It contains the power-law (RS-type) limits needed to replace the Yukawa proxy for ell.
  - *Labels:* method; topic: ISL; access: search-snippet-only; records: observations:ISL_review_2026
- **Characteristic evolution of conformal scattering: I. Scalar Waves in Minkowski Spacetime** — He, Tian and Zhang (surnames as recorded) (2026). arXiv:2608.15729v1 (16 Aug 2026). <https://arxiv.org/html/2608.15729v1>
  - *Finding (snippet):* Per the programme record: studies compactified double-null scattering with alternative characteristic stencils, and reports loss of convergence near spatial infinity. Its stability, convergence and Richardson-extrapolation evidence are numerical, not interval certificates.
  - *For HDBLAST:* [methods] A compactified or characteristic outer region is one way to stop the far-boundary signals that contaminate the late Chat 14 values (H/H0 near 0.628690 arrives after taper and boundary signals could reach the shell). The reported loss of convergence at the compactified edge is a warning to test any such map with a convergence study at the new boundary.
  - *Labels:* method; topic: constraint control / boundaries; access: programme-record-only; records: methods:A2
- **A posteriori error analysis and adaptivity of a space-time finite element method for the wave equation in second order formulation** — Dong, Georgoulis, Mascotto and Wang (surnames as recorded) (2026). arXiv:2509.08537v2; Numerische Mathematik, DOI 10.1007/s00211-026-01561-3 (per programme record). <https://arxiv.org/html/2509.08537v2>
  - *Finding (snippet):* Per the programme record: explicit a posteriori error bounds from temporal and spatial reconstructions, including mesh-change discontinuities, for a scalar wave equation on a bounded domain.
  - *For HDBLAST:* [methods] a route from 'grid convergence' to computable error estimates for wave evolutions. The HDBLAST system (coupled, nonlinear, with a shell boundary) would need its own constants. At most it can guide an error indicator near the steep wall.
  - *Labels:* method; topic: constraint control / boundaries; access: programme-record-only; records: methods:A3
- **Probing inflationary particle production with the CMB power spectrum** — Jense, Abu El-Haj, Hill and Philcox (2026). arXiv:2606.26823 (25 June 2026). <https://arxiv.org/html/2606.26823v1>
  - *Finding (snippet):* Per the programme record: a particle-production burst from an inflaton-spectator mass crossing is turned into CMB templates normalized by the inverse-covariance norm (their eqs. 26-27, 30); the weak joint-data preference is not an established detection.
  - *For HDBLAST:* [methods] if a mass-crossing burst ever happens during an inflating stage of an HDBLAST variant, this gives the observable-template method. It does not address reheating.
  - *Labels:* method; topic: particle production / backreaction; access: programme-record-only; records: methods:D5
- **Computer-assisted construction of SU(2)-invariant negative Einstein metrics** — Qiu Shi Wang (2026). Ann. Global Anal. Geom. 69, 18 (2026); arXiv:2504.21644v2. <https://arxiv.org/abs/2504.21644v2>
  - *Finding (snippet):* Per the programme record: high-accuracy approximate solutions, rigorous residual bounds and a fixed-point argument yield actual Riemannian Einstein metrics, with explicit constants and reproducible computations.
  - *For HDBLAST:* [methods] the closest structural precedent found. The Euclidean continuation of the +1 branch is a symmetry-reduced (cohomogeneity-one) Riemannian Einstein-scalar metric on a ball with a regular centre and an S^4 shell boundary. The Riemannian approximation -> residual -> fixed-point architecture therefore applies more directly to the static branch than to any Lorentzian evolution. The static Lorentzian problem is the same ODE BVP.
  - *Labels:* method; topic: validated numerics; access: programme-record-only; records: methods:E1
- **Holographic confining theories on space-times with constant positive curvature** — J. Kastikainen, E. Kiritsis, F. Nitti (2025). JHEP 09 (2025) 139; arXiv:2502.04036. <https://arxiv.org/abs/2502.04036>
  - *Finding (snippet):* Two branches of solutions compete, with a phase transition as the curvature varies: the low-curvature phase has the flat-space IR geometry, the high-curvature phase a regular interior; the transition is first- or higher-order depending on the potential's leading exponent.
  - *For HDBLAST:* [theory] HDBLAST also has two static branches (cone near -1 and near +1). The same free-energy (on-shell action) comparison would say which is preferred at the registered H, a cheap, well-defined next calculation.
  - *Labels:* method; topic: dS-sliced flows / holographic reheating; access: search-snippet-only; records: theory:F4
- **Brane Cosmology from AdS/BCFT** — Fujiki, Kanda, Kohara, Takayanagi (surnames only, per snippet) (2025). JHEP 03 (2025) 135; arXiv:2501.05036. <https://arxiv.org/abs/2501.05036>
  - *Finding (snippet):* An end-of-the-world brane in AdS with a brane-localized scalar: the brane equation becomes a Friedmann-like equation, the model can describe creating a universe via a big bang, the near-hyperplane effective action is Liouville gravity with scalar matter, and a timelike g-theorem is proven from the null energy condition (AdS3/BCFT2).
  - *For HDBLAST:* [theory] A recent formal analogue of a shell carrying a scalar in AdS that produces a big-bang-like creation; the null-energy-condition monotonicity is a structural tool HDBLAST could test on its shell trajectories.
  - *Labels:* method; topic: dS-sliced flows / holographic reheating; access: search-snippet-only; records: theory:F9
- **Stochastic gravitational-wave background search using data from five pulsar timing arrays** (2025). arXiv Dec 2025 (W.-W. Yu, B. Allen). <https://arxiv.org/abs/2512.08666>
  - *Finding (snippet):* Public pulse arrival times from five PTAs combined into a 121-pulsar data set, about four times larger than any single PTA's, using a 'direct combination' method for shared pulsars. Central result: posterior on amplitude and exponent of a power-law energy-density spectrum.
  - *For HDBLAST:* [observations] This is the first public multi-PTA combination found (it is not IPTA DR3). It fitted only a power law, so it does not test the knee, but the combined data could host an independent knee test.
  - *Labels:* method; topic: PTA; access: search-snippet-only; records: observations:YuAllen_5PTA
- **Nonlinear Gravitational Wave Memory: Universal Low-Frequency Background** (2025). arXiv Nov 2025. <https://arxiv.org/abs/2511.08514>
  - *Finding (snippet):* Nonlinear memory dominates the low-frequency GW spectrum during radiation and kination domination.
  - *For HDBLAST:* [observations] A second universal infrared contribution that any cosmological spectral template, the knee included, should respect.
  - *Labels:* method; topic: PTA; access: search-snippet-only; records: observations:memory_tail
- **Adaptive FEM with explicit time integration for the wave equation** — M. J. Grote, O. Lakkis, C. S. Santos (2025). arXiv:2507.11193v3; J. Comput. Appl. Math. 481, 117272 (2026) (per programme record). <https://arxiv.org/abs/2507.11193v3>
  - *Finding (snippet):* Per the programme record: combines a posteriori error indicators, evolving spatial meshes and local explicit time steps for the wave equation.
  - *For HDBLAST:* [methods] the constraint workstream of this round found the outgoing constraint front is generated within about 0.02 of the shell in the first ~0.05 time units. Local refinement with local time stepping in exactly that region is the natural use. Not an error certificate.
  - *Labels:* method; topic: constraint control / boundaries; access: programme-record-only; records: methods:A4
- **The dark bubbleography** — S. Banerjee, U. Danielsson, M. Zemsch (2024). JHEP 02 (2024) 102; arXiv:2311.16242. <https://arxiv.org/abs/2311.16242>
  - *Finding (snippet):* Presents the holographic construction of the dark bubble and shows, following holographic renormalization, that non-normalizable modes are essential for a vanishing induced graviton mass in any braneworld model; applies this to the propagator on the wall.
  - *For HDBLAST:* [theory] The claim is stated for 'any braneworld model'. HDBLAST's compact regular-cone bulk has a normalizable zero mode, and its tensor sector was found mode-stable (Chat 12). How these two statements fit together should be checked before HDBLAST claims 4D graviton behaviour.
  - *Labels:* method; topic: dark bubbles; access: search-snippet-only; records: theory:A9
- **Nucleation of de Sitter from the anti de Sitter spacetime in scalar field models** — M. Cadoni, M. Pitzalis, A. P. Sanna (2024). Eur. Phys. J. C (2025); arXiv:2407.10469. <https://arxiv.org/abs/2407.10469>
  - *Finding (snippet):* In Einstein-scalar gravity, gravitational coupling can drive nucleation of de Sitter from AdS through a static, spherically symmetric, metastable scalar lump. Euclidean actions and free energies are compared: the AdS lump is generally less favoured, and the most preferred state is a de Sitter vacuum.
  - *For HDBLAST:* [theory] A worked Euclidean-action comparison in Einstein-scalar gravity. The same comparison between HDBLAST's two static branches (cone near phi=-1 versus near phi=+1) has not been done and is cheap.
  - *Labels:* method; topic: creation / nucleation; access: search-snippet-only; records: theory:B9
- **Crunch from AdS bubble collapse in unbounded potentials** — not captured in snippet (2024). arXiv:2411.07692 (Nov 2024). <https://arxiv.org/abs/2411.07692>
  - *Finding (snippet):* Classical evolution of an AdS bubble nucleated from a Minkowski false vacuum in a potential with an infinitely deep true vacuum along a quartic slope. The interior collapses and forms a spacelike curvature singularity behind an apparent horizon. The scalar's kinetic energy overtakes the negative potential energy, the core density turns positive and trapped surfaces form. This requires no lower bound on the potential.
  - *For HDBLAST:* [theory] HDBLAST's potential is also unbounded below (U ~ -(2/27) phi^6 at large |phi|, checked exactly in theory_checks). Its Chat 13 drafts wrongly blamed the reversal on this; the corrected cause is sub-balanced tension. This paper supplies the right diagnostic for the collapse fate: search for trapped surfaces and an apparent horizon.
  - *For HDBLAST:* [methods] the diagnostic set for the collapse fate: apparent-horizon location, trapped-surface formation, and the sign of the local energy density. The registered U is also unbounded below (sweep 1 checked U ~ -(2/27) phi^6 exactly).
  - *Labels:* method; topic: creation / nucleation, shells and junctions; access: search-snippet-only; records: theory:B17, methods:B5
- **Shedding light on dark bubble cosmology** — I. Basile, U. Danielsson, S. Giri, D. Panizo (2023). JHEP 02 (2024) 112; arXiv:2310.15032. <https://arxiv.org/abs/2310.15032>
  - *Finding (snippet):* Incorporates electromagnetic fields: worldvolume fields backreact on the ambient 5D universe, changing the energy-momentum distribution and the effective gravity induced on the brane, and the resulting 4D cosmology consistently contains electromagnetic waves.
  - *For HDBLAST:* [theory] Shows that radiation on a wall must source bulk fields for consistency. In HDBLAST the only bulk fields are the metric and phi, so a radiation era would need an analogous consistent bulk source (for example a bulk black hole or Weyl term).
  - *Labels:* method; topic: dark bubbles; access: search-snippet-only; records: theory:A7
- **Gravitational waves in dark bubble cosmology** — U. Danielsson and co-authors (snippet) (2022). Phys. Rev. D 106, 024002; arXiv:2202.00545. <https://arxiv.org/abs/2202.00545>
  - *Finding (snippet):* Constructs the 5D uplift of 4D gravitational waves in de Sitter cosmology on a nucleated bubble in AdS5, extends the link between dark bubbles and Vilenkin's quantum cosmology to gravitational perturbations, and explains apparently negative energy contributions in the 4D Einstein equations that distinguish the dark bubble from Randall-Sundrum.
  - *For HDBLAST:* [theory] Method for uplifting tensor perturbations on a nucleated wall. A related snippet (query 61) states that bubble-nucleation boundary conditions select the Vilenkin weighting; that is the kind of weighting any HDBLAST 'creation event' would need.
  - *Labels:* method; topic: dark bubbles; access: search-snippet-only; records: theory:A5
- **Searching for Coleman-de Luccia bubbles in AdS compactifications** — G. Dibitetto, N. Petri (2022). Phys. Rev. D 107, 046020 (2023); arXiv:2207.02172. <https://arxiv.org/abs/2207.02172>
  - *Finding (snippet):* Fully backreacted gravitational instantons are obtained by numerically integrating first-order Hamilton-Jacobi equations: smooth domain walls with de Sitter foliations connecting a supersymmetric and a non-supersymmetric AdS vacuum, showing a nonperturbative instability of the latter.
  - *For HDBLAST:* [theory] The same mathematics as HDBLAST, a superpotential (Hamilton-Jacobi) description of dS-foliated walls, applied in truncations of string theory. It is a numerical method and a possible cross-check for HDBLAST's static-branch solver.
  - *For HDBLAST:* [methods] an independent solver architecture (first-order Hamilton-Jacobi) for the same class of solutions; a cross-check on the +1 branch that shares no code with the checkpoint's shooting solvers.
  - *Labels:* method; topic: creation / nucleation, gauge-invariant perturbations; access: search-snippet-only; records: theory:B11, methods:C11
- **Poschl-Teller non-production points (paper title not recorded)** — Ahmadiniaz et al. (as recorded) (2022 (arXiv number)). arXiv:2205.15946. <https://arxiv.org/pdf/2205.15946>
  - *Finding (snippet):* Per the programme record: section V discusses Poschl-Teller non-production points in an electric-field construction with additional momentum dependence.
  - *For HDBLAST:* [methods] a second exact benchmark family for particle-production codes (reflectionless profiles).
  - *Labels:* method; topic: particle production / backreaction; access: programme-record-only; records: methods:D3

### Context (37)

- **Dark bubbles, dark dimensions and fat gravitons** — U. Danielsson, S. Giri (2026). arXiv:2606.20942 (18 June 2026). <https://arxiv.org/abs/2606.20942>
  - *Finding (snippet):* Uses the instabilities behind the de Sitter swampland conjectures to make accelerated expansion inevitable. The dark bubble is presented as a realization of the dark-dimension proposal and of Sundrum's fat graviton. Predictions: a micron-size dark dimension, gravity fading at micron distances, a string scale of tens of TeV, and a measurable positive spatial curvature.
  - *For HDBLAST:* [theory] Shows the kind of falsifiable output this class of model can produce. HDBLAST has no such prediction yet. Positive spatial curvature is natural for a closed de Sitter slicing, which HDBLAST's static shells also use; this is a possible future observable, but it is not derived for HDBLAST.
  - *Labels:* context; topic: dark bubbles; access: search-snippet-only; records: theory:A12
- **Dark energy from string theory: an introductory review** — D. Andriot (2026). arXiv:2603.25797 (Mar 2026, rev. Apr 2026). <https://arxiv.org/abs/2603.25797>
  - *Finding (snippet):* A review of obtaining dark energy from string theory, either as a cosmological constant (de Sitter solution) or dynamical (quintessence), including historical no-go constraints and attempts to evade them.
  - *For HDBLAST:* [theory] The current reference point for whether any late-time de Sitter phase, such as HDBLAST's +1 plateau, can be UV-complete.
  - *Labels:* context; topic: swampland / strings; access: search-snippet-only; records: theory:G9
- **The NANOGrav 15 yr and 20 yr Datasets: Timing Events and Pulse Shape Changes** (2026). ApJ 1005 (June 2026). <https://iopscience.iop.org/article/10.3847/1538-4357/ae6db6>
  - *Finding (snippet):* A 2026 paper already describes timing events and pulse-shape changes in the NANOGrav 20-yr dataset. No 20-yr gravitational-wave background result was found in the searches.
  - *For HDBLAST:* [observations] A 20-yr span gives frequency resolution of about 1.58 nHz, and 1/yr falls near bin 20. The 20-yr background analysis is the next dataset that could test the knee directly. It had not been released as of these searches.
  - *Labels:* context; topic: PTA; access: search-snippet-only; records: observations:NG20_timing
- **New DESI DR2 Lyman-alpha Results Shed Light on Dark Energy** (2026). DESI collaboration news, 30 July 2026. <https://www.desi.lbl.gov/2026/07/30/new-desi-dr2-lyman-alpha-results-shed-light-on-dark-energy/>
  - *Finding (snippet):* The new DR2 Lyman-alpha forest measurements agree with the standard model. DESI notes the hints of evolving dark energy could fade, or a more complex model may be needed.
  - *For HDBLAST:* [observations] As above: late-time dark energy is currently outside the registered model's reach.
  - *Labels:* context; topic: DESI; access: search-snippet-only; records: observations:DESI_Lya_2026
- **GWTC-4.0: Tests of General Relativity. II. Parameterized Tests** (2026). arXiv Mar 2026. <https://arxiv.org/abs/2603.19020>
  - *Finding (snippet):* 42 confident O4a signals plus 49 earlier events. Graviton-mass bound m_g < 1.92e-23 eV/c^2 (90%), 1.16 times better than GWTC-3.
  - *For HDBLAST:* [observations] In RS2 the 4D graviton is massless and the massive KK modes are micron-scale, so this bound does not touch the registered model at tabletop-consistent ell. It excludes cosmological-ell variants in which massive-graviton dispersion would appear.
  - *Labels:* context; topic: LVK; access: search-snippet-only; records: observations:GWTC4_TGR
- **Particle Creation from Entanglement Entropy** — not stated in the programme record (2026). PTEP 2026, 013A01. <https://academic.oup.com/ptep/article/2026/1/013A01/8400342>
  - *Finding (snippet):* Per the programme record: the main calculation is restricted to low-velocity moving mirrors; the programme treated it as an analogy, not a brane energy source.
  - *For HDBLAST:* [methods] moving-mirror (dynamical Casimir) production is the closest analogue of production by a moving shell. The low-velocity restriction does not cover the Chat 14 roll-off.
  - *Labels:* context; topic: particle production / backreaction; access: programme-record-only; records: methods:D8
- **Weak gravity at micron scales from dark bubble cosmology and its cosmological consequences** — U. Danielsson, S. Giri (2025). arXiv:2511.21362; Phys. Rev. D (2026). <https://arxiv.org/abs/2511.21362>
  - *Finding (snippet):* Predicts that gravity becomes weaker, not stronger, than Newtonian at about 1e-5 m, with explicit table-top predictions. The same effect lowers effective gravity at high energy densities and gives early inflation with nothing beyond radiation. The paper also discusses a quantum origin of the universe in which a 5D black hole of critical size catalyses nucleation of the dark bubble at its horizon, and the black hole's matter becomes the matter on the bubble.
  - *For HDBLAST:* [theory] This is the closest recent published proposal to 'a five-dimensional event created our universe'. It also supplies what HDBLAST lacks: a source of matter, the catalysing black hole. HDBLAST cannot claim the general concept as new and should cite this line of work.
  - *Labels:* context; topic: dark bubbles; access: search-snippet-only; records: theory:A11
- **Self-sustained, out-of-equilibrium inflation** — J. Casalderrey-Solana and colleagues (snippet) (2025). arXiv:2512.18079. <https://arxiv.org/abs/2512.18079>
  - *Finding (snippet):* Holographic de Sitter-invariant states of non-conformal strongly coupled QFTs on dS4: out-of-equilibrium effects can sustain exponential inflation with H far below the QFT scale and the species scale; fine-tuning scales only logarithmically; apparent horizons with growing area signal growing comoving entropy; the regime can be the late-time limit of an FRW start.
  - *For HDBLAST:* [theory] A de Sitter attractor with a coupled hot sector. It contrasts with HDBLAST's attractor, an empty RS de Sitter brane, and suggests one way a coupled matter sector could change the endpoint.
  - *Labels:* context; topic: dS-sliced flows / holographic reheating; access: search-snippet-only; records: theory:F8
- **Stability of non-supersymmetric vacua from calibrations** — not captured in snippet (2025). JHEP 11 (2025) 070; arXiv:2507.02787. <https://arxiv.org/abs/2507.02787>
  - *Finding (snippet):* Calibrations bound D-brane energies and forbid, in the probe approximation, nucleation of D-brane bubbles in several type II AdS4 and AdS5 non-supersymmetric vacua, many of which resisted all decay channels tested. Query 23 snippets summarise the opposite view: the Ooguri-Vafa conjecture that non-SUSY AdS is at best metastable via Brown-Teitelboim brane nucleation.
  - *For HDBLAST:* [theory] The instability of non-supersymmetric AdS, which the dark-bubble origin mechanism needs, is contested. HDBLAST's bulk vacua come from a real superpotential and are perturbatively stable; whether a 5D 'blast' event is available depends on that unresolved UV question.
  - *Labels:* context; topic: swampland / strings; access: search-snippet-only; records: theory:G8
- **SPT-3G D1: Constraints on inflationary gravitational waves with two years of SPT-3G data** (2025). Phys. Rev. D (Dec 2025). <https://arxiv.org/abs/2505.02827>
  - *Finding (snippet):* SPT-3G alone over the BICEP/Keck field: r < 0.25 (95%), sigma(r) = 0.067, with delensing.
  - *For HDBLAST:* [observations] An independent tensor channel. It is not competitive yet, but the South Pole Observatory programme targets sigma(r) = 0.001 (snippet).
  - *Labels:* context; topic: CMB; access: search-snippet-only; records: observations:SPT3G_r
- **SPT-3G D1: CMB temperature and polarization power spectra and cosmology from 2019 and 2020 observations of the SPT-3G Main field** (2025). Phys. Rev. D (2025). <https://arxiv.org/abs/2506.20707>
  - *Finding (snippet):* SPT-3G alone: H0 = 66.66 +/- 0.60. SPT + Planck + ACT: H0 = 67.19 +/- 0.38, sigma8 = 0.8137 +/- 0.0037. CMB alone shows no evidence beyond LambdaCDM. There is a 2.8 sigma CMB-vs-DESI DR2 difference within LambdaCDM.
  - *For HDBLAST:* [observations] The early-universe expansion history HDBLAST must reproduce, if it ever reaches a radiation era.
  - *Labels:* context; topic: CMB; access: search-snippet-only; records: observations:SPT3G_D1
- **DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints** (2025). arXiv Mar 2025 (+ DESI guide page https://www.desi.lbl.gov/2025/03/19/desi-dr2-results-march-19-guide/). <https://arxiv.org/abs/2503.14738>
  - *Finding (snippet):* w0waCDM is preferred over LambdaCDM at 3.1 sigma (DESI + CMB) and at 2.8-4.2 sigma when supernovae are added. The best fit has w > -1 today and w < -1 in the past, crossing -1 near z = 0.5.
  - *For HDBLAST:* [observations] The registered model's relaxing fate is an empty RS de Sitter brane with w = -1 exactly. At the registered detuning, H*ell = 0.069, so identifying that brane with today's dark energy needs ell of about 372 Mpc (excluded by tabletop tests by a factor 3e29) or delta of about 1e-62 (script). HDBLAST therefore makes no DESI prediction. A confirmed evolving w would need an extra late-time sector.
  - *Labels:* context; topic: DESI; access: search-snippet-only; records: observations:DESI_DR2
- **The Dark Energy Survey Supernova Program: A Reanalysis Of Cosmology Results And Evidence For Evolving Dark Energy With An Updated Type Ia Supernova Calibration** (2025). MNRAS (2026); arXiv Nov 2025. <https://arxiv.org/abs/2511.07517>
  - *Finding (snippet):* Recalibration (DES-Dovekie) reduces the significance of rejecting LambdaCDM from 4.2 sigma (DES-SN5YR) to 3.2 sigma with CMB + DESI DR2. Only a weak Bayesian preference for w0wa remains.
  - *For HDBLAST:* [observations] The evolving-dark-energy signal is not settled, so it is not yet a firm constraint on any braneworld late-time sector.
  - *Labels:* context; topic: DESI; access: search-snippet-only; records: observations:DES_Dovekie
- **Braneworld dark energy in light of DESI DR2** (2025). JCAP 11 (2025) 018. <https://arxiv.org/abs/2507.07193>
  - *Finding (snippet):* Thawing scalar fields on a (4+1)-dimensional ghost-free 'phantom brane' (induced-gravity, DGP-type normal branch) show an effective phantom-divide crossing. A quartic potential explains the DESI DR2 constraints well.
  - *For HDBLAST:* [observations] A braneworld route to DESI-like behaviour exists, but it uses induced gravity on the brane and a cosmological crossover scale. The registered model is RS2-type (no induced gravity), so this does not carry over without new terms.
  - *Labels:* context; topic: DESI; access: search-snippet-only; records: observations:brane_DESI
- **The Local Distance Network: a community consensus report on the measurement of the Hubble constant at 1% precision** (2025). A&A (April 2026). <https://arxiv.org/abs/2510.23823>
  - *Finding (snippet):* H0 = 73.50 +/- 0.81 km/s/Mpc from a blinded community distance network (Cepheids, TRGB, Miras and more). This is 7.1 sigma from the early-universe LambdaCDM value 67.24 +/- 0.35. A simple local systematic is judged unlikely.
  - *For HDBLAST:* [observations] Extra dark radiation from a bulk remnant would raise the CMB-inferred H0 only slightly, and ACT DR6 disfavours that route. HDBLAST makes no H0 prediction. One practical effect: converting a strain amplitude into Omega_gw scales as h^-2. With h = 0.735 instead of 0.674, the knee's h_c at 1/yr would be 2.62e-15 instead of 2.40e-15 (script). The frozen Omega_k is tied to h = 0.674.
  - *Labels:* context; topic: H0; access: search-snippet-only; records: observations:H0DN
- **The Chicago Carnegie Hubble Program: Improving the Calibration of SNe Ia with JWST Measurements of the Tip of the Red Giant Branch** (2025). arXiv Mar 2025. <https://arxiv.org/abs/2503.11769>
  - *Finding (snippet):* The search summary quotes H0 = 69.8 +/- 1.9 from CCHP's JWST Cepheid/TRGB/JAGB work, consistent with the CMB value. It was not verified which CCHP paper this number is from.
  - *For HDBLAST:* [observations] The local-ladder disagreement is not fully settled, so HDBLAST should not be tuned to either H0.
  - *Labels:* context; topic: H0; access: search-snippet-only; records: observations:CCHP_2025
- **Beyond No Tension: JWST z > 10 Galaxies Push Simulations to the Limit** (2025). arXiv Sep 2025. <https://arxiv.org/abs/2509.07695>
  - *Finding (snippet):* Search summary: spectroscopic galaxies up to z of about 14.4 (MoM-z14) and 14.3 (GS-z14). Hints of an overabundance of bright or massive early galaxies; implied star-formation efficiencies are high. There is an ongoing debate over whether this reflects new physics or systematics (stellar masses, AGN, selection).
  - *For HDBLAST:* [observations] Tests small-scale primordial power and early structure growth. HDBLAST has no perturbation spectrum, so it is untested here. A future HDBLAST spectrum with enhanced small-scale power would be checked against these counts.
  - *Labels:* context; topic: JWST; access: search-snippet-only; records: observations:JWST_z10
- **Experimental tests of dark bubble cosmology** — U. Danielsson, D. Panizo (snippet) (2024). Phys. Rev. D 109, 026003; arXiv:2311.14589. <https://link.aps.org/doi/10.1103/PhysRevD.109.026003>
  - *Finding (snippet):* The snippet gives only the title, authors and journal; no result is summarized here.
  - *For HDBLAST:* [theory] Listed so that the testability line of the dark bubble programme is not missed; content not assessed.
  - *Labels:* context; topic: dark bubbles; access: search-snippet-only; records: theory:A8
- **De Sitter space constraints on brane tensions and couplings** — S. Hassan, G. Obied, J. March-Russell (2024). arXiv:2411.14529. <https://arxiv.org/abs/2411.14529>
  - *Finding (snippet):* Festina-Lente-type arguments using Nariai de Sitter black holes bound p-brane tensions in de Sitter by the Hubble rate and by Chern-Simons-like worldvolume couplings; D-branes satisfy them; axion domain walls evade them.
  - *For HDBLAST:* [theory] Tangential: HDBLAST's shell has no gauge couplings, so these bounds would matter only if gauge fields were added on the shell.
  - *Labels:* context; topic: swampland / strings; access: search-snippet-only; records: theory:G7
- **The MeerKAT Pulsar Timing Array: Maps of the gravitational-wave sky with the 4.5 year data release** (2024). MNRAS 536, 1501 (2025). <https://arxiv.org/abs/2412.01214>
  - *Finding (snippet):* Gravitational-wave sky maps from the 4.5-yr data. The snippet reported tentative background evidence but no anisotropy numbers.
  - *For HDBLAST:* [observations] Anisotropy is an SMBHB discriminator. A cosmological HDBLAST relic would be isotropic to high precision, so detected anisotropy would count against a cosmological reading of the knee.
  - *Labels:* context; topic: PTA; access: search-snippet-only; records: observations:MPTA_maps
- **Comparing Recent Pulsar Timing Array Results on the Nanohertz Stochastic Gravitational-wave Background** (2024). ApJ 966, 105 (2024); arXiv 2309.00693. <https://iopscience.iop.org/article/10.3847/1538-4357/ad36be>
  - *Finding (snippet):* EPTA, InPTA, NANOGrav and PPTA results assessed on equal footing: background spectral parameters agree within 1 sigma. A standardized noise model reduces tensions in pulsar noise parameters.
  - *For HDBLAST:* [observations] Background for the project's cross-team 'overlap' criterion: agreement is at the power-law-parameter level, not at the level of a shared knee.
  - *Labels:* context; topic: PTA; access: search-snippet-only; records: observations:IPTA_compare
- **JWST Observations Reject Unrecognized Crowding of Cepheid Photometry as an Explanation for the Hubble Tension at 8 sigma Confidence** (2024). ApJL (2024). <https://arxiv.org/abs/2401.04773>
  - *Finding (snippet):* More than 1000 Cepheids measured with both HST and JWST: mean distance difference -0.01 +/- 0.03 mag. Crowding is rejected as the explanation at 8.2 sigma.
  - *For HDBLAST:* [observations] Context for the tension's robustness.
  - *Labels:* context; topic: H0; access: search-snippet-only; records: observations:JWST_crowding
- **Damped energy-norm a posteriori error estimates for fully discrete approximations of the wave equation using C2-reconstructions with the leapfrog scheme** — T. Chaumont-Frelet, A. Ern (2024). arXiv:2403.12954v2 (20 Dec 2024). <https://arxiv.org/abs/2403.12954v2>
  - *Finding (snippet):* Per the programme record: scalar finite-element/leapfrog a posteriori estimates in a damped energy norm; the programme noted that they concern a different formulation and do not supply its constants.
  - *For HDBLAST:* [methods] a second example of energy-norm error control for explicit wave schemes. Relevance is indirect for the method-of-lines finite-difference solvers used in HDBLAST.
  - *Labels:* context; topic: constraint control / boundaries; access: programme-record-only; records: methods:A5
- **Features of a dark energy model in string theory** — S. Banerjee, U. Danielsson, S. Giri (2023). Phys. Rev. D 108, 126009; arXiv:2212.14004. <https://arxiv.org/abs/2212.14004>
  - *Finding (snippet):* Clears up misconceptions about the dark bubble: points out important differences from Randall-Sundrum and explains why gravity neither is, nor needs to be, localized on the dark bubble.
  - *For HDBLAST:* [theory] HDBLAST is Randall-Sundrum-like: a Z2 doubled bulk closed off by a regular cone, so the bulk volume is finite. Its 4D gravity mechanism therefore differs from the dark bubble's, and results cannot be transferred between the two without care.
  - *Labels:* context; topic: dark bubbles; access: search-snippet-only; records: theory:A6
- **Bubbles of cosmology in AdS/CFT** — A. Sahu, P. Simidzija, M. Van Raamsdonk (2023). JHEP (Nov 2023); arXiv:2306.13143. <https://arxiv.org/abs/2306.13143>
  - *Finding (snippet):* The typical big-bang/big-crunch cosmologies of holographic effective theories are not asymptotically AdS, but arbitrarily large spherical bubbles of them can be embedded in asymptotically AdS spacetimes with a Schwarzschild-AdS exterior.
  - *For HDBLAST:* [theory] A holographic framing in which a cosmological region is bounded by a shell with a black-hole exterior. It is a possible interpretation of a future HDBLAST shell with a bulk black-hole mass.
  - *Labels:* context; topic: creation / nucleation; access: search-snippet-only; records: theory:B19
- **The NANOGrav 15-year Data Set: Evidence for a Gravitational-Wave Background** (2023). ApJL 951, L8 (2023). <https://arxiv.org/abs/2306.16213>
  - *Finding (snippet):* Hellings-Downs-correlated signal in 67 pulsars. For a fiducial f^-2/3 strain spectrum the amplitude is 2.4 (+0.7/-0.6) x 10^-15 (median, 90% interval) at 1/yr. A power-law background is favoured over independent pulsar noise alone with a Bayes factor above 10^14 (snippet wording). Consistent with supermassive black-hole binaries (SMBHBs); exotic cosmological or astrophysical sources not excluded.
  - *For HDBLAST:* [observations] This amplitude is the calibration anchor of the registered knee: the frozen Omega_k = 8.00e-9 converts to h_c = 2.40e-15 at 1/yr (h = 0.674), reproducing this median to better than 0.1% in strain (0.3% in Omega; script check).
  - *Labels:* context; topic: PTA; access: search-snippet-only; records: observations:NG15_evidence
- **The NANOGrav 15-year Data Set: Search for Signals from New Physics** (2023). ApJL 951, L11 (2023). <https://arxiv.org/abs/2306.16219>
  - *Finding (snippet):* Inflation, scalar-induced GWs, first-order phase transitions, cosmic strings and domain walls were tested. All except stable field-theory cosmic strings can reproduce the signal, some with Bayes factors O(10)-O(100) over the SMBHB model. The results are model-sensitive and not conclusive.
  - *For HDBLAST:* [observations] Sets the bar for any HDBLAST PTA claim. A better fit than a fixed SMBHB template is not evidence for new physics, and the project's BIC screens are weaker than these analyses.
  - *Labels:* context; topic: PTA; access: search-snippet-only; records: observations:NG15_newphysics
- **Search for an Isotropic Gravitational-wave Background with the Parkes Pulsar Timing Array** (2023). ApJL (2023), doi 10.3847/2041-8213/acdd02. <https://arxiv.org/abs/2306.16215>
  - *Finding (snippet):* Common-spectrum process with A = 2.0 +/- 0.2 x 10^-15 (h ~ f^-2/3). Hellings-Downs consistency with false-alarm probability p < about 0.014. The signal strength appeared time-dependent, contrary to expectations for an isotropic background.
  - *For HDBLAST:* [observations] The frozen knee has 1.44 times PPTA's power at 1/yr. The reported time-dependence is a caution for any fine spectral feature.
  - *Labels:* context; topic: PTA; access: search-snippet-only; records: observations:PPTA_DR3
- **Searching for the nano-Hertz stochastic gravitational wave background with the Chinese Pulsar Timing Array Data Release I** (2023). RAA 23, 075024 (2023). <https://arxiv.org/abs/2306.16216>
  - *Finding (snippet):* 57 millisecond pulsars, close to 3 yr with FAST. Hellings-Downs evidence at 4.6 sigma around 14 nHz (discrete-frequency method). log10 A_c = -14.4 (+1.0/-2.8) for strain index in [-1.8, 1.5].
  - *For HDBLAST:* [observations] Resolution about 10.6 nHz, so the knee sits at bin 2.99, essentially at 1/yr. The amplitude is too loose to test the knee.
  - *Labels:* context; topic: PTA; access: search-snippet-only; records: observations:CPTA_DR1
- **Pulsar timing array stochastic background from light Kaluza-Klein resonances** (2023). Phys. Rev. D 108, 095017 (2023); arXiv 2306.17071. <https://doi.org/10.1103/PhysRevD.108.095017>
  - *Finding (snippet):* In a warped 5D model (UV brane at the Planck scale, dark brane at the GeV scale), PTA data can be fitted by a first-order confinement phase transition of a radion at the MeV-GeV scale, feebly coupled to the Standard Model. Many existing embeddings are not viable because of radion/graviton phenomenology. A multi-brane set-up is proposed to remain consistent with collider and gravity tests.
  - *For HDBLAST:* [observations] This is prior art for a five-dimensional origin of the nHz signal, via a phase transition near T ~ 0.1-1 GeV. That matches the script's map of the knee to T ~ 0.19-0.28 GeV for f*/H* = 1. It is not the registered HDBLAST mechanism, and it shows that 5D PTA explanations face strong gravity/collider constraints.
  - *Labels:* context; topic: PTA; access: search-snippet-only; records: observations:KK_PTA
- **LambdaCDM not dead yet: massive high-z Balmer break galaxies are less common than previously reported** (2023). arXiv Oct 2023. <https://arxiv.org/abs/2310.03063>
  - *Finding (snippet):* Title-level: massive high-z Balmer-break galaxies are less common than earlier claimed.
  - *For HDBLAST:* [observations] Early 'impossible galaxy' claims have weakened. They are not a motivation for a non-standard origin.
  - *Labels:* context; topic: JWST; access: search-snippet-only; records: observations:LCDM_not_dead
- **dS4 universe emergent from Kerr-AdS5 spacetime: bubble nucleation catalyzed by a black hole** — not captured in snippet (2022). JHEP 05 (2023) 107; arXiv:2209.05625. <https://arxiv.org/abs/2209.05625>
  - *Finding (snippet):* Studies nucleation of a vacuum bubble in Kerr-AdS5, sufficient conditions for nucleation with a rotating black hole, and how the black hole changes the transition rate; the bubble carries an emergent dS4 universe.
  - *For HDBLAST:* [theory] Quantitative black-hole-catalysed creation rates for a dS4 wall in AdS5. It is the closest rate calculation to a '5D event creating our universe'.
  - *Labels:* context; topic: creation / nucleation; access: search-snippet-only; records: theory:B5
- **Lectures on the string landscape and the Swampland** — N. B. Agmon, A. Bedroya, M. J. Kang, C. Vafa (2022). arXiv:2212.06187. <https://arxiv.org/pdf/2212.06187>
  - *Finding (snippet):* Lecture notes on the string landscape and on the Swampland programme's constraints for EFTs with a quantum-gravity UV completion.
  - *For HDBLAST:* [theory] Background reference for the swampland criteria cited here.
  - *Labels:* context; topic: swampland / strings; access: search-snippet-only; records: theory:G3
- **The dark dimension scenario (Montero, Vafa, Valenzuela 2022; popular account)** — M. Montero, C. Vafa, I. Valenzuela (per snippet; arXiv:2205.12293 named in snippet, not opened) (2022). Quanta Magazine article (2024) and follow-up EPJC paper arXiv:2309.09330. <https://www.quantamagazine.org/in-a-dark-dimension-physicists-search-for-missing-matter-20240201/>
  - *Finding (snippet):* One micron-size extra dimension with a Kaluza-Klein scale of order meV, motivated by the tiny cosmological constant (about 1e-122 in Planck units) and swampland arguments.
  - *For HDBLAST:* [theory] A concrete single-extra-dimension scenario with laboratory tests; A12 ties the dark bubble to it. HDBLAST's extra-dimension scale (AdS radius ell = 9 in model units at phi=+1) has no physical normalization yet (open item 7).
  - *Labels:* context; topic: swampland / strings; access: search-snippet-only; records: theory:G5
- **Improved limits on the tensor-to-scalar ratio using BICEP and Planck data** (2022). Phys. Rev. D 105, 083524 (2022). <https://arxiv.org/abs/2112.07961>
  - *Finding (snippet):* Planck + BICEP/Keck 2018 + BAO: r < 0.032 (95%). BK18 alone: r < 0.036.
  - *For HDBLAST:* [observations] The standard tensor bound, before the 2025 combination above.
  - *Labels:* context; topic: CMB; access: search-snippet-only; records: observations:BK_Planck_r
- **International Pulsar Timing Array (home page)** (year not given). IPTA web site. <https://ipta4gw.org/>
  - *Finding (snippet):* Appeared in results for IPTA DR3 queries. The snippet summary describes IPTA-DR3 as a reanalysis combining the member PTAs (and CPTA) that is expected to be more sensitive than any single data set. No DR3 background result was found.
  - *For HDBLAST:* [observations] IPTA DR3 would be the decisive combined test at 1/yr. The project's own site page (checked Sept 2026) also records it as not yet released.
  - *Labels:* context; topic: PTA; access: search-snippet-only; records: observations:IPTA_home
- **Periodic Spectral Features in the NANOGrav 15-Year Gravitational Wave Background: A Phenomenological Analysis** (year not given). Zenodo record (not peer reviewed; authorship not seen). <https://zenodo.org/records/17956111>
  - *Finding (snippet):* Claims a systematic alternation of spectral excess and deficit across frequency bins, with a period of about 4 nHz.
  - *For HDBLAST:* [observations] A caution: with 14 correlated free-spectrum bins, apparent features are easy to find, and features at the edges of the band (like the knee) are the least constrained.
  - *Labels:* context; topic: PTA; access: search-snippet-only; records: observations:zenodo_periodic

## Earlier foundational work cited by the sweeps (before 2022)

Listed for completeness; see `literature/theory.md` and `literature/methods.md` for full annotations.

| Year | Source | Label | Why it matters for HDBLAST (first sentence of the sweep note) |
|---|---|---|---|
| 1969 | [Parker (1969), Physical Review 183, 1057 (title not written in the programme record)](https://journals.aps.org/pr/abstract/10.1103/PhysRev.183.1057) | method | a mandatory zero-production control for any gravitational particle-production calculation. |
| 1980 | [Gravitational Effects on and of Vacuum Decay](https://www.semanticscholar.org/paper/Gravitational-Effects-on-and-of-Vacuum-Decay-Coleman-Luccia/a60b863946ca2d85da9d8b48b72257e8ec05e3ed) | context | HDBLAST's balanced tension sigma=2W is the critical, BPS-like tension at which a flat wall exists (checked exactly in theory_checks: 2W = 6/ell at both vacua). |
| 1982 | [Instability of the Kaluza-Klein vacuum](https://ui.adsabs.harvard.edu/abs/1982NuPhB.195..481W/abstract) | context | The foundational example of a higher-dimensional vacuum event with a de Sitter-like expanding wall. |
| 1991 | [Perturbations on domain walls and strings: A covariant theory](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.44.1007) | method | Tachyonic masses on de Sitter walls can partly reflect wall translation, which is not a physical instability. |
| 1998 (arXiv number) | [Instant preheating (title and authors not written in the programme record; see caveat)](https://arxiv.org/abs/hep-ph/9812289v2) | method | the mechanism used by this round's preheating workstream, which found it insufficient on the fixed archived trajectory. |
| 1998 | [Smooth 'creation' of an open universe in five dimensions](https://arxiv.org/abs/hep-th/9804106) | context | Conceptual precedent: a smooth higher-dimensional event can look like a 4D singularity. |
| 1999 | [An Alternative to Compactification](https://ui.adsabs.harvard.edu/abs/1999PhRvL..83.4690R/abstract) | context | The base model that HDBLAST's late-time static branch reduces to. |
| 1999 | [Bent Domain Walls as Braneworlds](https://arxiv.org/pdf/hep-th/9905210) | method | One of the original derivations of the de Sitter brane in AdS5. |
| 1999 | [Dynamics of Anti-de Sitter Domain Walls](https://arxiv.org/abs/hep-th/9910149v1) | method | An exact way to add a bulk black-hole mass behind a shell and obtain an a^-4 ('dark radiation') term. |
| 1999 | [Gravity in the Randall-Sundrum Brane World](https://arxiv.org/abs/hep-th/9911055) | method | Sets the benchmark for gravity on HDBLAST's late-time shell, and indicates what a scalar (radion-like) admixture would do to 4D tests. |
| 1999 | [Modeling the fifth dimension with scalars and gravity](https://arxiv.org/pdf/hep-th/9909134) | method | HDBLAST's registered model uses exactly this structure: U = (1/2)W_phi^2 - (2/3)W^2, phi' = W_phi, A' = -W/3, balanced tension 2W. |
| 1999 | [The Einstein Equations on the 3-Brane World](https://arxiv.org/abs/gr-qc/9910076) | method | The standard way to see the bulk Weyl ('dark radiation') and bulk-scalar contributions to HDBLAST's shell Friedmann equation. |
| 1999 | [Wave function of the radion in a brane world](https://arxiv.org/abs/hep-th/9912160) | method | A reference construction of scalar metric modes; useful when interpreting any HDBLAST shell-scalar coupling to brane matter. |
| 2000 | [Brane New World](https://arxiv.org/abs/hep-th/0003052) | context | An alternative quantum origin of a de Sitter brane in AdS, with a holographic matter sector and a perturbation prediction. |
| 2000 | [Brane-world creation and black holes](https://arxiv.org/abs/hep-th/9912118) | method | HDBLAST's static geometry is this instanton's Lorentzian section with a bulk scalar added: a dS-sliced AdS region closed off by a regular cone, Z2-doubled, with one brane. |
| 2000 | [Radion on the de Sitter brane](https://arxiv.org/pdf/gr-qc/0011078) | supports | Explains why HDBLAST's tachyon cannot be a pure-gravity radion: HDBLAST has one shell, and the regular cone is not a second brane. |
| 2000 | [The big bang as a higher-dimensional shock wave](https://arxiv.org/pdf/gr-qc/0003012) | context | The earlier work whose wording is closest to 'higher-dimensional blast'. |
| 2001 | [A Braneworld Universe From Colliding Bubbles](https://ar5iv.arxiv.org/html/hep-th/0107148) | context | Prior art for a 5D nucleation-and-collision origin of an RS brane universe. |
| 2001 | [Brane Big-Bang Brought by Bulk Bubble](https://arxiv.org/pdf/hep-th/0110286) | supports | Very close to the HDBLAST narrative: a detuned, inflating brane plus a bulk event that ends inflation and reheats it. |
| 2001 | [CFT and Entropy on the Brane](https://www.semanticscholar.org/paper/CFT-and-Entropy-on-the-Brane-Savonije-Verlinde/7e51da58584a6054a76ef7e1d6c126e0c797e242) | supports | Gives the holographic meaning of the a^-4 term from a bulk black hole (item D5). |
| 2001 | [Cosmological Perturbations Generated in the Colliding Bubble Braneworld Universe](https://arxiv.org/abs/hep-th/0111089) | method | A method for primordial perturbations from a 5D origin event. |
| 2001 | [Locally localized gravity](https://iopscience.iop.org/article/10.1088/1126-6708/2001/05/008) | method | The tension classification behind HDBLAST: tension above the critical value gives a de Sitter brane. |
| 2001 | [The Ekpyrotic Universe: Colliding Branes and the Origin of the Hot Big Bang](https://arxiv.org/html/hep-th/0103239v2) | context | Canonical prior art for 'a higher-dimensional event makes the hot Big Bang'. |
| 2001 | [Thick Brane Worlds and Their Stability](https://arxiv.org/abs/hep-th/0107025) | method | A counterpoint: smooth thick de Sitter walls with no thin shell are stable, whereas HDBLAST's thin shell with detuned tension has a tachyon. |
| 2002 | [Born-Again Braneworld](https://arxiv.org/abs/hep-th/0210250) | context | An example where radion (modulus) dynamics of inflating branes, rather than matter, drives the transition, as in HDBLAST's tachyonic roll-off. |
| 2002 | [Cosmological evolution with brane-bulk energy exchange](https://arxiv.org/abs/hep-th/0207060) | method | Classification template for HDBLAST's proposed exchange law rho_dot + 3H(rho+p) = j phi_dot. |
| 2002 | [Observational Constraints on Dark Radiation in Brane Cosmology](https://arxiv.org/abs/astro-ph/0203272) | constrains | Any HDBLAST 'hot' ingredient built from a bulk black hole or Weyl term must satisfy these bounds. |
| 2003 | [Can Inflating Braneworlds be Stabilized?](https://arxiv.org/abs/hep-th/0309002v1) | supports | In the literature, tachyonic moduli of de Sitter branes are generic. |
| 2003 | [Exactly solvable model for cosmological perturbations in dilatonic brane worlds](https://arxiv.org/pdf/hep-th/0307073) | method | An exact benchmark for validating HDBLAST's 5D perturbation and evolution codes on a problem with a known answer, a cheap calibration control. |
| 2003 | [When do colliding bubbles produce an expanding universe?](https://arxiv.org/abs/hep-th/0306151) | method | A published analogue of HDBLAST's two fates, expansion versus reversal and collapse, with an explicit criterion (expansion rate against fifth-dimension momentum transfer) that could be evaluated on the Chat 14 trajectories. |
| 2004 | [BRANECODE: A Program for Simulations of Braneworld Dynamics](https://arxiv.org/abs/hep-ph/0404141v1) | method | Prior numerical art for HDBLAST's rolloff5d evolutions. |
| 2004 | [Collision of Domain Walls and Reheating of the Brane Universe](https://arxiv.org/abs/hep-th/0406235) | method | A calibrated wall-collision reheating estimate against which HDBLAST's formal production proxy (22 Sept matter extension) could be compared, once physical scales are chosen. |
| 2004 | [Cosmological perturbations in a big crunch/big bang space-time](https://arxiv.org/pdf/hep-th/0306109) | method | A matching method for perturbations across a 5D collision event. |
| 2004 | [Fake supergravity and domain wall stability](https://ui.adsabs.harvard.edu/abs/2004PhRvD..69j4027F/abstract) | constrains | Covers HDBLAST's balanced flat wall (delta=0) but, per the snippet, only flat and AdS-sliced walls: no comparable positive-energy theorem is cited for de Sitter-sliced walls. |
| 2004 | [M theory model of a big crunch/big bang transition](https://ui.adsabs.harvard.edu/abs/2004PhRvD..70j6004T/abstract) | context | Shows how a 5D collision can be continued through; relevant if HDBLAST's collapse fate is ever continued past the singularity. |
| 2005 | [Regular collision of dilatonic inflating branes](https://arxiv.org/abs/hep-th/0508145) | supports | Bulk-scalar inflating branes with a radion instability are established; HDBLAST's tachyonic shell mode is of this general type (in a single-shell, regular-cone geometry). |
| 2005 | [Solution of a Braneworld Big Crunch/Big Bang Cosmology](https://arxiv.org/abs/hep-th/0512123) | constrains | A warning for HDBLAST: its closed-form 4D effective theory agrees with 5D to about 1e-6 for the static tachyon, but may fail near fast events (turnaround, collapse). |
| 2006 | [Coupled bulk and brane fields about a de Sitter brane](https://arxiv.org/abs/hep-th/0612202) | challenges | Directly relevant to HDBLAST's matter extension, which couples the brane field chi to the bulk phi. |
| 2006 | [Hidden Supersymmetry of Domain Walls and Cosmologies](https://www.osti.gov/etdeweb/biblio/20777232) | method | A first-order (fake-superpotential) description of curved flows. |
| 2006 | [RS2 gravitational-wave calculation with bulk Kaluza-Klein leakage (title not recorded)](https://arxiv.org/abs/hep-th/0601105v2) | constrains | a numerical precedent for brane-bulk energy bookkeeping in a radiation era; it limits any expansion-only inference. |
| 2007 | [Bouncing Negative-Tension Branes](https://arxiv.org/pdf/0708.0743) | context | Listed for completeness of the brane-collision line. |
| 2007 | [Dynamics of colliding branes and black brane production](https://arxiv.org/abs/gr-qc/0702138) | challenges | Suggests a definite hypothesis for HDBLAST's collapse fate: black-brane formation with the shell trapped. |
| 2007 | [Nonperturbative Instability of AdS5 x S5/Zk](https://arxiv.org/abs/0709.4262) | context | Non-supersymmetric AdS5 backgrounds are generically suspect in string theory. |
| 2008 | [Reheating the Universe in Braneworld Cosmological Models with bulk-brane energy transfer](https://arxiv.org/abs/0805.1792) | method | A direct precedent for the reheating calculation HDBLAST still needs (open item 4). |
| 2010 | [Brane-World Gravity](https://arxiv.org/abs/1004.3962) | context | A standard reference for conventions (junctions, Weyl term, perturbations) when comparing HDBLAST with the literature. |
| 2010 | [Verification Methods: Rigorous Results Using Floating-Point Arithmetic](https://doi.org/10.1017/S096249291000005X) | method | the survey behind interval-Newton/Krawczyk-type existence tests with outward-rounded floating point, the step used by the programme's earlier certified BVP root (M462). |
| 2013 | [Out of the White Hole: A Holographic Origin for the Big Bang](https://arxiv.org/abs/1309.1487) | context | A prominent published model of a 5D gravitational event producing the Big Bang, with its own perturbation mechanism. |
| 2014 | [Simulating the universe(s) II: phenomenology of cosmic bubble collisions in full General Relativity](https://arxiv.org/pdf/1407.2950) | method | Methodology for turning a violent, collision-type origin event into CMB templates, if HDBLAST ever reaches that stage. |
| 2014 | [Tensor Perturbations from Brane-World Inflation with Curvature Effects](https://arxiv.org/abs/1308.5765) | context | If HDBLAST produced perturbations on the brane at energies comparable to the brane tension, the r and n_t relations would differ from 4D. |
| 2017 | [Cosmological Perturbations in the 5D Holographic Big Bang Model](https://arxiv.org/abs/1703.00954) | constrains | Shows how a 5D-origin model is confronted with CMB data, and that doing so can produce tension. |
| 2017 | [Generalized surface tension bounds in vacuum decay](https://arxiv.org/abs/1711.06776) | method | Gives a principled way to define an effective shell tension when the scalar is nonconstant across the wall, as on HDBLAST's +1 branch. |
| 2017 | [Holographic RG flows on curved manifolds and quantum phase transitions](https://arxiv.org/abs/1711.08462) | method | A classification framework for HDBLAST's static de Sitter-sliced solutions. |
| 2017 | [Holographic self-tuning of the cosmological constant; Brane cosmology and the self-tuning of the cosmological constant](https://arxiv.org/abs/1704.05075) | method | The time-dependent PDE-plus-junction problem is the same type as HDBLAST's roll-off; its coordinate choices and solution strategies are directly reusable. |
| 2017 | [New observational limits on dark radiation in braneworld cosmology](https://doi.org/10.1103/PhysRevD.95.083516) | constrains | A direct bound on the bulk Weyl term. |
| 2018 | [De Sitter Space and the Swampland](https://www.semanticscholar.org/paper/De-Sitter-Space-and-the-Swampland-Obied-Ooguri/599c99078a502b7d462d0b8783ff4a2c4436dc19) | constrains | HDBLAST's late-time attractor is an empty de Sitter brane. |
| 2018 | [De Sitter and Anti-de Sitter branes in self-tuning models](https://arxiv.org/abs/1807.09794) | method | The closest published model class: a brane with scalar-dependent tension in an Einstein-dilaton bulk with de Sitter solutions. |
| 2018 | [Emergent de Sitter Cosmology from Decaying Anti-de Sitter Space](https://arxiv.org/abs/1807.01570) | context | The best-developed published precedent for 'our universe is an expanding codimension-one wall in 5D AdS with de Sitter induced on it'. |
| 2018 | [Limits on the number of spacetime dimensions from GW170817](https://arxiv.org/abs/1801.08160) | method | Foundational method. |
| 2018 | [Our universe: An expanding bubble in an extra dimension (press release on A1)](https://www.sciencedaily.com/releases/2018/12/181228164824.htm) | context | Shows that 'our universe is a 5D bubble' has already been widely publicised; HDBLAST communication should cite this line to avoid implying the concept is new. |
| 2019 | [A new kind of cyclic universe](https://arxiv.org/abs/1904.08022) | context | The ekpyrotic programme itself moved from 5D brane collisions to 4D bounces. |
| 2019 | [Exact Bogoliubov coefficients for a sech-profile mass pulse (paper title not recorded)](https://link.springer.com/article/10.1140/epjc/s10052-019-6581-2) | method | an exact calibration for any mode-function particle-production code, alongside the Landau-Zener limit used for instant preheating. |
| 2019 | [Numerical Verification Methods and Computer-Assisted Proofs for Partial Differential Equations](https://doi.org/10.1007/978-981-13-7669-6) | method | the standard reference class for verified solutions of elliptic/ODE BVPs and verified eigenvalue enclosures. |
| 2019 | [de Sitter Cosmology on an expanding bubble](https://arxiv.org/abs/1907.04268) | context | The Unruh-temperature reading applies equally to the HDBLAST static shells (temperature H/2pi of the dS slicing). |
| 2019 | [de Sitter versus Anti de Sitter flows and the (super)gravity landscape (and Part II)](https://arxiv.org/abs/1901.04546) | constrains | Constrains any HDBLAST variant in which a positive-energy 5D region, the 'blast', is joined smoothly to the AdS bulk. |
| 2020 | [Catalytic creation of a bubble universe induced by quintessence in five dimensions](https://arxiv.org/abs/2011.07437) | context | Precedent for catalysed creation of a 4D universe on a 5D bubble. |
| 2020 | [Causal gravitational waves as a probe of free streaming particles and the expansion of the Universe](https://arxiv.org/abs/2010.03568) | challenges | The frozen curve rises as f^2 below the knee, not the universal f^3 causal tail of a horizon-limited radiation-era source. |
| 2020 | [Constraint on the radius of five-dimensional dS spacetime with GW170817 and GRB 170817A](https://arxiv.org/abs/2001.06581) | method | A 'shortcut through the bulk' arrival-time test is the natural GW170817-type test for a de Sitter brane. |
| 2020 | [Cosmology at the end of the world; Cosmology from the vacuum; Accelerating cosmology from a holographic wormhole](https://www.nature.com/articles/s41567-020-0909-6) | context | Microscopic (holographic) control of brane cosmologies inside AdS black holes: a possible long-term UV framing for a shell-plus-bulk-black-hole version of HDBLAST. |
| 2020 | [Dark bubbles: decorating the wall](https://arxiv.org/abs/2001.07433) | method | A worked template for the sector HDBLAST lacks: brane matter that backreacts consistently on a 5D wall geometry. |
| 2020 | [New Test of the Gravitational 1/r^2 Law at Separations down to 52 micrometres](https://arxiv.org/abs/2002.11761) | constrains | Sets the physical unit of the registered model. |
| 2020 | [Nothing really matters](https://ui.adsabs.harvard.edu/abs/2020JHEP...08..040D/abstract) | context | An alternative 'higher-dimensional blast': a bubble of nothing in AdS whose wall is de Sitter. |
| 2020 | [Robustness of slow contraction to cosmic initial conditions](https://arxiv.org/abs/2006.04999) | method | A methodological template for testing whether HDBLAST's two fates are robust to non-perturbative initial disturbances. |
| 2020 | [Verified Bounds for the Determinant of Real or Complex Point or Interval Matrices](https://doi.org/10.1016/j.cam.2019.112610) | method | certifying that a junction Jacobian is nonsingular is the local-uniqueness half of a Krawczyk argument. |
| 2020 | [de Sitter Bubbles and the Swampland](https://arxiv.org/abs/2008.07555) | constrains | TCC-type bounds would limit how long any HDBLAST de Sitter phase (the static shell or the +1 plateau) can last in a UV-complete version. |
| 2021 | [Curing with hemlock: escaping the swampland using instabilities from string theory](https://arxiv.org/abs/2103.17121) | context | States the dark-bubble response to swampland objections. |
| 2021 | [Dark bubbles and black holes](https://link.springer.com/article/10.1007/JHEP09(2021)158) | context | Part of the black-hole-plus-dark-bubble line (see A11, B5). |
| 2021 | [Revisiting Coleman-de Luccia transitions in the AdS regime using holography](https://arxiv.org/abs/2102.11881) | method | HDBLAST's static solutions belong to exactly this class: Einstein-scalar solutions with de Sitter slicing and a regular endpoint where the slice shrinks. |

## Reproduction

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927
python3 synthesis/build_literature.py   # rewrites this file and synthesis/LITERATURE_MERGED.json
```
