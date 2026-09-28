"""Source records for literature sweep 2 (observations), 28 Sept 2026.

Every record is based ONLY on web-search result snippets (titles, URLs, short
summaries).  No paper page could be opened (arxiv.org, zenodo.org and journal
sites are blocked for downloads in this environment).  Numbers quoted in
'finding' are transcribed from those snippets; where the snippet summary did
not make clear which page a statement came from, 'attribution' says so.

Fields: id, title, url, year, venue, access, finding, relevance, implication,
queries (1-based query numbers in QUERIES), attribution (optional).
"""

ACCESS = "search-snippet-only"

QUERIES = [
    "NANOGrav 15-year gravitational wave background spectral turnover low frequency environmental effects",
    "IPTA third data release DR3 gravitational wave background 2026",
    "MeerKAT Pulsar Timing Array 4.5 year data gravitational wave background amplitude spectral index",
    '"Stochastic gravitational-wave background search using data from five pulsar timing arrays"',
    "NANOGrav 15 yr Data Set Running of the Spectral Index result inconclusive",
    "MeerKAT pulsar timing array first search gravitational waves common signal amplitude larger than other PTAs frequency",
    "EPTA DR2new InPTA gravitational wave background amplitude spectral index 10.3 years evidence",
    "Comparing recent pulsar timing array results on the nanohertz stochastic gravitational-wave background IPTA 2024 consistent",
    "Parkes Pulsar Timing Array DR3 isotropic gravitational wave background amplitude 13/3 Hellings-Downs significance",
    "pulsar timing array sensitivity loss at 1/yr frequency pulsar position fit timing model hasasia sensitivity curve",
    "NANOGrav 15-year evidence gravitational wave background 16.03 yr 14 frequencies amplitude 2.4e-15 gamma 3.2",
    "NANOGrav 20-year data set release 2026",
    '"Piecewise Power-Law Reconstruction" NANOGrav 15 yr gravitational-wave background',
    "NANOGrav 15 yr erratum supermassive black hole binaries 2026 software bug gravitational wave background",
    "pulsar timing array spectral break knee broken power law gravitational wave background 2025 search",
    "NANOGrav piecewise power law 2601.09481 broken power law Bayes factor spectral shape result high frequency",
    "EPTA DR2 improved noise model averaging gravitational wave background lower amplitude strain spectral index -2/3 Max Planck 2025",
    '"Impacts of Customized Chromatic Noise Models" NANOGrav 15 yr gravitational wave',
    '"Pulsar timing arrays: the emerging gravitational-wave landscape"',
    "Chinese Pulsar Timing Array DR1 gravitational wave background 14 nHz 4.6 sigma FAST",
    "NANOGrav 15 year search for signals from new physics cosmological interpretation Bayes factor 10 to 100 phase transition cosmic strings",
    "pulsar timing array stochastic background light Kaluza-Klein resonances warped extra dimension phase transition",
    "causal gravitational wave infrared tail f^3 radiation domination universal low-frequency spectrum cosmological source",
    "DESI DR2 BAO evolving dark energy w0 wa significance 2025 results",
    "DES supernova recalibration Dovekie 2025 2026 evolving dark energy significance reduced DESI",
    "braneworld dark energy DESI DR2 phantom crossing DGP brane fit 2025",
    "ACT DR6 2025 N_eff constraint extended models dark radiation Hubble tension",
    "SPT-3G D1 2025 cosmological constraints N_eff Hubble constant CMB-SPA combined",
    "tensor-to-scalar ratio upper limit 2025 2026 BICEP Keck SPT delensing r < 0.03",
    "ACT DR6 power spectra cosmological parameters spectral index n_s 0.974 P-ACT-LB inflation implications",
    '"Inflation at the End of 2025" constraints on r and n_s latest CMB BAO data result',
    "big bang nucleosynthesis 2024 2025 N_eff constraint deuterium helium primordial abundance ΔN_eff bound",
    "braneworld dark radiation Weyl term bulk black hole mass constraint BBN CMB Randall-Sundrum brane tension bound",
    "Eöt-Wash 2020 torsion balance inverse-square law test 52 micrometers extra dimension size bound",
    "short-range gravity inverse square law test 2024 2025 new limit Yukawa micrometer experiment",
    "Randall-Sundrum AdS curvature radius bound table-top gravity experiments brane tension TeV fundamental scale 10^8 GeV",
    '"New Test of the Gravitational 1/r^2 Law at Separations down to 52" Yukawa |alpha|=1 38.6 micrometers dark energy length scale',
    '"Short-Range Tests of the Gravitational Inverse-Square Law" 2026 review extra dimensions limits',
    "GW170817 constraint number of spacetime dimensions gravitational wave leakage luminosity distance damping extra dimensions",
    "gravitational wave leakage extra dimensions constraint dark sirens GWTC-3 screening scale 2023 2024 higher-dimensional",
    "LIGO Virgo KAGRA O4 upper limits isotropic gravitational-wave background 2025 Omega_GW 25 Hz",
    "GWTC-4 tests of general relativity 2025 LIGO Virgo KAGRA graviton mass dispersion O4a",
    "Hubble tension 2026 status SH0ES JWST Cepheid H0 73 CCHP Freedman 70 TRGB JAGB local distance network",
    "JWST early massive galaxies z>10 overabundance LCDM tension 2025 review bright galaxies star formation efficiency",
    "Randall-Sundrum high-energy brane inflation tensor amplitude enhancement tensor-to-scalar ratio braneworld consistency relation observational constraint",
    "Local Distance Network H0DN 2025 1% Hubble constant 73.50 consensus distance ladder paper",
    "DESI 2026 dark energy update full-shape DR2 evolving dark energy new result",
    "gravitational wave peak frequency today temperature at production 1.65e-7 Hz T/GeV g* nanohertz QCD epoch",
]
REFUSED_QUERIES = [
    "IPTA DR3 combined dataset status 2026 International Pulsar Timing Array combination early results",
    '"What is the source of the PTA GW signal" supermassive black hole binaries environmental effects cosmological models comparison',
    "EPTA second data release implications for massive black holes dark matter and the early Universe",
]

S = []


def add(**kw):
    kw.setdefault("access", ACCESS)
    kw.setdefault("attribution", "")
    S.append(kw)


# ---------------------------------------------------------------- PTA -------
add(id="NG15_evidence", section="PTA", year="2023",
    title="The NANOGrav 15-year Data Set: Evidence for a Gravitational-Wave Background",
    url="https://arxiv.org/abs/2306.16213", venue="ApJL 951, L8 (2023)",
    finding=("Hellings-Downs-correlated signal in 67 pulsars. For a fiducial f^-2/3 strain spectrum the amplitude is "
             "2.4 (+0.7/-0.6) x 10^-15 (median, 90% interval) at 1/yr. A power-law background is favoured over independent "
             "pulsar noise alone with a Bayes factor above 10^14 (snippet wording). Consistent with supermassive black-hole "
             "binaries (SMBHBs); exotic cosmological or astrophysical sources not excluded."),
    relevance=("This amplitude is the calibration anchor of the registered knee: the frozen Omega_k = 8.00e-9 converts to "
               "h_c = 2.40e-15 at 1/yr (h = 0.674), reproducing this median to better than 0.1% in strain (0.3% in Omega; script check)."),
    implication="context", queries=[11])
add(id="NG15_running", section="PTA", year="2025",
    title="The NANOGrav 15 yr Data Set: Running of the Spectral Index",
    url="https://arxiv.org/abs/2408.10166", venue="ApJL 978, L29 (2025); arXiv Aug 2024",
    finding=("Running-power-law fit (logarithmic frequency dependence of the index). 95% credible interval for the running "
             "beta in [-0.80, 2.96], consistent with zero. Bayes factor running vs constant power law 0.69 +/- 0.01 "
             "(inconclusive). The constant power law still suffices."),
    relevance=("A knee is an extreme form of spectral curvature. Within the NANOGrav band the frozen curve is almost a pure "
               "power law (curvature only in bins 13-14), so this null result neither supports nor excludes it. It does "
               "show that current data do not require curvature."),
    implication="constrains", queries=[1, 5])
add(id="NG15_PPL", section="PTA", year="2026",
    title="The NANOGrav 15 yr Data Set: Piecewise Power-Law Reconstruction of the Gravitational-Wave Background",
    url="https://arxiv.org/abs/2601.09481", venue="ApJL (2026), doi 10.3847/2041-8213/ae7086",
    finding=("Piecewise power-law (PPL) spectral reconstruction: constant, broken, doubly broken models combined by "
             "Bayesian model averaging. Described as closer to physically realistic (especially cosmological) spectra "
             "than the free spectrum. The snippets did not report the fitted break positions or Bayes factors."),
    relevance=("This is the collaboration's own analogue of the project's knee test. The frozen SBPL should be compared "
               "with the PPL posterior (bin-by-bin power and break-frequency posterior), not only with summary points. "
               "Result numbers must be read from the paper before any claim."),
    implication="method", queries=[11, 13, 16])
add(id="NG15_CNM", section="PTA", year="2026",
    title="The NANOGrav 15 yr Data Set: Impacts of Customized Chromatic Noise Models on Gravitational Wave Analyses",
    url="https://arxiv.org/abs/2606.28554", venue="arXiv June 2026 (companion: arXiv 2606.28571)",
    finding=("Customized chromatic (interstellar-medium) noise models for the 15-yr pulsars. Bayes factor for Hellings-Downs "
             "correlations over an uncorrelated common red process: 1571 +/- 14 (14 Fourier modes), about 8 times "
             "earlier. Power-law amplitude at fixed index reduced to 2.1 (+0.6/-0.5) x 10^-15."),
    relevance=("Relative to this amplitude, the frozen knee has 1.31 times the power at 1/yr (1.14 times the strain). "
               "That is inside the quoted uncertainty. Better noise modelling lowered the amplitude and strengthened the "
               "correlation signature. The knee's height was calibrated to the 2023 value, and the direction of change "
               "matters for any refreeze."),
    implication="constrains", queries=[18])
add(id="NG15_erratum_KDE", section="PTA", year="2026",
    title="KDE Representations of the Gravitational Wave Background Free Spectra Present in the NANOGrav 15-Year Dataset (corrected release; erratum Agazie et al. 2026, ApJL 1006, L67)",
    url="https://zenodo.org/records/21844115", venue="Zenodo data record + ApJL erratum (2026)",
    finding=("A bug in the parallel-tempering routine of PTMCMCSampler biased the free-spectrum posteriors of the 15-yr "
             "data set. The record re-releases corrected free-spectrum KDEs. The erratum states SMBHB population inference "
             "was affected. The project's own site notes the SMBHB conclusion was unchanged."),
    relevance=("Any project screen of the knee that used pre-2026 NANOGrav free-spectrum posteriors or 'summary points' "
               "derived from them should be rerun on the corrected KDEs. The frozen numbers must stay frozen."),
    implication="method", queries=[14],
    attribution="Snippet summary combined the Zenodo record text and an erratum citation; the erratum page itself was not seen.")
add(id="NG15_newphysics", section="PTA", year="2023",
    title="The NANOGrav 15-year Data Set: Search for Signals from New Physics",
    url="https://arxiv.org/abs/2306.16219", venue="ApJL 951, L11 (2023)",
    finding=("Inflation, scalar-induced GWs, first-order phase transitions, cosmic strings and domain walls were tested. "
             "All except stable field-theory cosmic strings can reproduce the signal, some with Bayes factors O(10)-O(100) "
             "over the SMBHB model. The results are model-sensitive and not conclusive."),
    relevance=("Sets the bar for any HDBLAST PTA claim. A better fit than a fixed SMBHB template is not evidence for new "
               "physics, and the project's BIC screens are weaker than these analyses."),
    implication="context", queries=[21])
add(id="NG20_timing", section="PTA", year="2026",
    title="The NANOGrav 15 yr and 20 yr Datasets: Timing Events and Pulse Shape Changes",
    url="https://iopscience.iop.org/article/10.3847/1538-4357/ae6db6", venue="ApJ 1005 (June 2026)",
    finding=("A 2026 paper already describes timing events and pulse-shape changes in the NANOGrav 20-yr dataset. "
             "No 20-yr gravitational-wave background result was found in the searches."),
    relevance=("A 20-yr span gives frequency resolution of about 1.58 nHz, and 1/yr falls near bin 20. The 20-yr "
               "background analysis is the next dataset that could test the knee directly. It had not been released "
               "as of these searches."),
    implication="context", queries=[12])
add(id="EPTA_DR2_III", section="PTA", year="2023",
    title="The second data release from the European Pulsar Timing Array III. Search for gravitational wave signals",
    url="https://arxiv.org/abs/2306.16214", venue="A&A 678, A50 (2023)",
    finding=("DR2new (latest 10.3 yr, 25 pulsars, plus about 3.5 yr of InPTA data for 10 of them): Bayes factor 60, false-alarm "
             "probability about 0.1% (at least 3 sigma). At fixed index 13/3: A = (2.5 +/- 0.7) x 10^-15 at 1/yr. "
             "Full DR2, HD process: log10 A = -14.54 (+0.28/-0.41), gamma = 4.19 (+0.73/-0.63)."),
    relevance=("The frozen knee lies within 0.03 frequency bins of 1/yr for DR2new (resolution 3.08 nHz). The EPTA "
               "free-index fit sits near gamma = 13/3, not the gamma = 3 that the frozen curve follows below the knee."),
    implication="constrains", queries=[7])
add(id="AEI_EPTA_noise", section="PTA", year="2025",
    title="A clearer view of gravitational-wave signals in pulsar timing arrays (AEI news on improved EPTA DR2 noise models)",
    url="https://www.aei.mpg.de/1323827/a-a-clearer-view-of-gravitational-wave-signals-in-pulsar-timing-arrays",
    venue="Max Planck Institute for Gravitational Physics news item",
    finding=("Improved EPTA DR2 noise modelling (pulsar-intrinsic, epoch-correlated and transient noise, with noise-model "
             "averaging) makes the background's strain index consistent with -2/3 (gamma = 13/3). It gives a lower "
             "amplitude at 1/yr than earlier analyses. Chromatic noise errors have whiter spectra, which push naive "
             "fits toward flatter indices."),
    relevance=("Two independent 2025-2026 reanalyses (EPTA; NANOGrav, entry above) move toward gamma = 13/3 and lower "
               "amplitude. The frozen curve's low-frequency branch is Omega ~ f^2 (gamma = 3). This trend is "
               "unfavourable to the frozen shape, but it is not a likelihood-level rejection."),
    implication="challenges", queries=[15, 17],
    attribution=("Year and details from the search summary; the news page date was not visible. The 'whiter spectra' "
                 "sentence may come from the NANOGrav chromatic-noise work rather than this news item."))
add(id="PPTA_DR3", section="PTA", year="2023",
    title="Search for an Isotropic Gravitational-wave Background with the Parkes Pulsar Timing Array",
    url="https://arxiv.org/abs/2306.16215", venue="ApJL (2023), doi 10.3847/2041-8213/acdd02",
    finding=("Common-spectrum process with A = 2.0 +/- 0.2 x 10^-15 (h ~ f^-2/3). Hellings-Downs consistency with "
             "false-alarm probability p < about 0.014. The signal strength appeared time-dependent, contrary to "
             "expectations for an isotropic background."),
    relevance=("The frozen knee has 1.44 times PPTA's power at 1/yr. The reported time-dependence is a caution for any "
               "fine spectral feature."),
    implication="context", queries=[9])
add(id="CPTA_DR1", section="PTA", year="2023",
    title="Searching for the nano-Hertz stochastic gravitational wave background with the Chinese Pulsar Timing Array Data Release I",
    url="https://arxiv.org/abs/2306.16216", venue="RAA 23, 075024 (2023)",
    finding=("57 millisecond pulsars, close to 3 yr with FAST. Hellings-Downs evidence at 4.6 sigma around 14 nHz "
             "(discrete-frequency method). log10 A_c = -14.4 (+1.0/-2.8) for strain index in [-1.8, 1.5]."),
    relevance="Resolution about 10.6 nHz, so the knee sits at bin 2.99, essentially at 1/yr. The amplitude is too loose to test the knee.",
    implication="context", queries=[20])
add(id="MPTA_data", section="PTA", year="2024",
    title="The MeerKAT Pulsar Timing Array: The 4.5-year data release and the noise and stochastic signals of the millisecond pulsar population",
    url="https://arxiv.org/abs/2412.01148", venue="arXiv Dec 2024",
    finding=("83 pulsars, 4.5 yr, high cadence. Common signal log10 A = -14.25 (+0.21/-0.36), gamma = 3.60 (+1.31/-0.89). "
             "At gamma = 13/3, log10 A = -14.28 +/- 0.21, ln Bayes factor 4.46. The amplitude is larger than other PTAs report."),
    relevance=("MeerKAT's free-gamma posterior still allows gamma = 3. Its resolution (7.04 nHz) puts the knee at bin 4.49, "
               "so its bins straddle 1/yr more coarsely than NANOGrav's."),
    implication="constrains", queries=[3])
add(id="MPTA_search", section="PTA", year="2024",
    title="The MeerKAT Pulsar Timing Array: The first search for gravitational waves with the MeerKAT radio telescope",
    url="https://arxiv.org/abs/2412.01153", venue="MNRAS 536, 1489 (2025)",
    finding=("Sky-averaged h_c,yr = 7.5 (+0.8/-0.9) x 10^-15 at strain index -0.26, or 4.8 (+0.8/-0.9) x 10^-15 at -2/3. "
             "This is inconsistent with other PTAs' common-noise results by at least about 1.4 sigma. Hellings-Downs "
             "significance is 3-3.4 sigma depending on noise assumptions."),
    relevance=("At 1/yr the frozen knee has only 0.25 (gamma = 13/3) to 0.10 (free index) of MeerKAT's power. The "
               "cross-PTA amplitude spread already exceeds a factor 4 in Omega. A single frozen height cannot match all "
               "teams, which is one of the project's own failure criteria (teams disagreeing)."),
    implication="constrains", queries=[6])
add(id="MPTA_maps", section="PTA", year="2024",
    title="The MeerKAT Pulsar Timing Array: Maps of the gravitational-wave sky with the 4.5 year data release",
    url="https://arxiv.org/abs/2412.01214", venue="MNRAS 536, 1501 (2025)",
    finding="Gravitational-wave sky maps from the 4.5-yr data. The snippet reported tentative background evidence but no anisotropy numbers.",
    relevance=("Anisotropy is an SMBHB discriminator. A cosmological HDBLAST relic would be isotropic to high precision, "
               "so detected anisotropy would count against a cosmological reading of the knee."),
    implication="context", queries=[3])
add(id="IPTA_compare", section="PTA", year="2024",
    title="Comparing Recent Pulsar Timing Array Results on the Nanohertz Stochastic Gravitational-wave Background",
    url="https://iopscience.iop.org/article/10.3847/1538-4357/ad36be", venue="ApJ 966, 105 (2024); arXiv 2309.00693",
    finding=("EPTA, InPTA, NANOGrav and PPTA results assessed on equal footing: background spectral parameters agree "
             "within 1 sigma. A standardized noise model reduces tensions in pulsar noise parameters."),
    relevance=("Background for the project's cross-team 'overlap' criterion: agreement is at the power-law-parameter "
               "level, not at the level of a shared knee."),
    implication="context", queries=[8])
add(id="YuAllen_5PTA", section="PTA", year="2025",
    title="Stochastic gravitational-wave background search using data from five pulsar timing arrays",
    url="https://arxiv.org/abs/2512.08666", venue="arXiv Dec 2025 (W.-W. Yu, B. Allen)",
    finding=("Public pulse arrival times from five PTAs combined into a 121-pulsar data set, about four times larger "
             "than any single PTA's, using a 'direct combination' method for shared pulsars. Central result: posterior "
             "on amplitude and exponent of a power-law energy-density spectrum."),
    relevance=("This is the first public multi-PTA combination found (it is not IPTA DR3). It fitted only a power law, "
               "so it does not test the knee, but the combined data could host an independent knee test."),
    implication="method", queries=[2, 4])
add(id="PTA_review_2026", section="PTA", year="2026",
    title="Pulsar timing arrays: the emerging gravitational-wave landscape",
    url="https://arxiv.org/abs/2603.13643", venue="arXiv Mar 2026 (review)",
    finding=("Six PTAs (NANOGrav, EPTA, InPTA, PPTA, CPTA, MPTA) report evidence for the background. The review argues "
             "that the perceived tension between current amplitudes and standard merger models is largely resolved by "
             "new insights into SMBHB populations."),
    relevance="If the SMBHB amplitude tension is resolved, a main motivation for an exotic source of the hum weakens.",
    implication="challenges", queries=[19])
add(id="IPTA_home", section="PTA", year="n/a",
    title="International Pulsar Timing Array (home page)",
    url="https://ipta4gw.org/", venue="IPTA web site",
    finding=("Appeared in results for IPTA DR3 queries. The snippet summary describes IPTA-DR3 as a reanalysis combining "
             "the member PTAs (and CPTA) that is expected to be more sensitive than any single data set. No DR3 "
             "background result was found."),
    relevance=("IPTA DR3 would be the decisive combined test at 1/yr. The project's own site page (checked Sept 2026) "
               "also records it as not yet released."),
    implication="context", queries=[2, 8],
    attribution="The DR3 description could not be tied to a specific page among the results; treat as unverified.")
add(id="PTA_sensitivity_1yr", section="PTA", year="2025",
    title="A sensitivity curve approach to tuning a pulsar timing array in the detection era (with the hasasia package)",
    url="https://arxiv.org/abs/2409.00336", venue="CQG (2025), doi 10.1088/1361-6382/adbbab",
    finding=("Search summary: fitting each pulsar's sky position and proper motion causes an extra loss of sensitivity "
             "around f = 1/yr, and fitting parallax causes one around 2/yr. These appear as large spikes in PTA "
             "sensitivity curves (computed with hasasia)."),
    relevance=("This is the key methodological fact for the registered knee. The script finds f_k = 0.9972 x (1/yr), "
               "within 0.045 bins of 1/yr for NANOGrav 15-yr and 0.013 bins for MeerKAT. A spectral bend exactly where "
               "the timing model absorbs power is maximally degenerate with it. A knee test must propagate the "
               "timing-model transmission function, or the knee location is not identifiable."),
    implication="constrains", queries=[10],
    attribution=("The sentence about 1/yr and 2/yr came from the results for this query. The summary did not say whether "
                 "it came from this paper or from arXiv 2608.00250 (directional anisotropic sensitivity curves)."))
add(id="causal_tail", section="PTA", year="2020",
    title="Causal gravitational waves as a probe of free streaming particles and the expansion of the Universe",
    url="https://arxiv.org/abs/2010.03568", venue="arXiv Oct 2020",
    finding=("For causal sources with finite correlation length, the infrared tail of a GW spectrum produced in "
             "radiation domination scales as f^3. The tail is modified by free-streaming particles and by "
             "non-standard expansion."),
    relevance=("The frozen curve rises as f^2 below the knee, not the universal f^3 causal tail of a horizon-limited "
               "radiation-era source. A cosmological HDBLAST reading of the knee therefore needs a named mechanism "
               "(non-standard expansion, free-streaming or a specific source) that yields f^2."),
    implication="challenges", queries=[23])
add(id="causal_Bmode", section="PTA", year="2026",
    title="A Universal CMB B-Mode Spectrum from Early Causal Tensor Sources",
    url="https://arxiv.org/abs/2601.20967", venue="arXiv Jan 2026",
    finding="Early causal tensor sources share a universal infrared scaling and predict the same B-mode angular distribution.",
    relevance=("If the blast is a causal (sub-horizon) tensor source, its CMB B-mode imprint is fixed up to amplitude. "
               "This is a cross-check between a PTA-band claim and CMB tensor limits."),
    implication="method", queries=[23])
add(id="memory_tail", section="PTA", year="2025",
    title="Nonlinear Gravitational Wave Memory: Universal Low-Frequency Background",
    url="https://arxiv.org/abs/2511.08514", venue="arXiv Nov 2025",
    finding="Nonlinear memory dominates the low-frequency GW spectrum during radiation and kination domination.",
    relevance="A second universal infrared contribution that any cosmological spectral template, the knee included, should respect.",
    implication="method", queries=[23])
add(id="KK_PTA", section="PTA", year="2023",
    title="Pulsar timing array stochastic background from light Kaluza-Klein resonances",
    url="https://doi.org/10.1103/PhysRevD.108.095017", venue="Phys. Rev. D 108, 095017 (2023); arXiv 2306.17071",
    finding=("In a warped 5D model (UV brane at the Planck scale, dark brane at the GeV scale), PTA data can be fitted by a "
             "first-order confinement phase transition of a radion at the MeV-GeV scale, feebly coupled to the Standard "
             "Model. Many existing embeddings are not viable because of radion/graviton phenomenology. A multi-brane "
             "set-up is proposed to remain consistent with collider and gravity tests."),
    relevance=("This is prior art for a five-dimensional origin of the nHz signal, via a phase transition near "
               "T ~ 0.1-1 GeV. That matches the script's map of the knee to T ~ 0.19-0.28 GeV for f*/H* = 1. It is not "
               "the registered HDBLAST mechanism, and it shows that 5D PTA explanations face strong gravity/collider "
               "constraints."),
    implication="context", queries=[22])
add(id="zenodo_periodic", section="PTA", year="unknown",
    title="Periodic Spectral Features in the NANOGrav 15-Year Gravitational Wave Background: A Phenomenological Analysis",
    url="https://zenodo.org/records/17956111", venue="Zenodo record (not peer reviewed; authorship not seen)",
    finding="Claims a systematic alternation of spectral excess and deficit across frequency bins, with a period of about 4 nHz.",
    relevance=("A caution: with 14 correlated free-spectrum bins, apparent features are easy to find, and features at "
               "the edges of the band (like the knee) are the least constrained."),
    implication="context", queries=[1],
    attribution="Not a project record (its ID does not appear in the project's files). Not peer reviewed; treat with care.")

# ---------------------------------------------------------------- CMB / BBN -
add(id="ACT_DR6_ext", section="CMB", year="2025",
    title="The Atacama Cosmology Telescope: DR6 Constraints on Extended Cosmological Models",
    url="https://arxiv.org/abs/2503.14454", venue="arXiv Mar 2025",
    finding=("No evidence for new free-streaming light species: N_eff = 2.86 +/- 0.13. Self-interacting dark radiation "
             "N < 0.134. Dark-radiation models are not favoured, because extra radiation increases Silk damping at "
             "high multipoles, where ACT DR6 is sensitive."),
    relevance=("In a braneworld the Weyl ('dark radiation') term C/a^4, set by a bulk black-hole mass, counts as Delta "
               "N_eff. A Gaussian approximation from these numbers gives Delta N_eff < 0.071, so rho_dr/rho_gamma < 0.016 "
               "at recombination (script). This bounds any bulk remnant of the blast, for example the collapse fate's "
               "black brane, and any GW relic, if a radiation era is ever achieved."),
    implication="constrains", queries=[27, 30, 32])
add(id="ACT_ns", section="CMB", year="2025",
    title="ACT DR6 Insights on alpha-Attractor Inflationary Models and Reheating",
    url="https://arxiv.org/abs/2505.01517", venue="arXiv May 2025",
    finding=("As reported there: Planck 2018 + ACT DR6 gives n_s = 0.9709 +/- 0.0038, rising to 0.9743 +/- 0.0034 with "
             "DESI Y1 BAO (the P-ACT-LB combination gives n_s = 0.974 +/- 0.003). Starobinsky-type plateau models become "
             "mildly disfavoured at about 2 sigma or more."),
    relevance=("Any 'blast instead of (or before) inflation' model must produce adiabatic, nearly scale-invariant "
               "perturbations with n_s about 0.97. The registered model has no perturbation-spectrum prediction yet, "
               "so this is currently an open requirement, not a test it passes or fails."),
    implication="constrains", queries=[30],
    attribution="n_s numbers as quoted in the snippet from this and related pages.")
add(id="r_ns_2025", section="CMB", year="2025",
    title="Inflation at the End of 2025: Constraints on r and n_s Using the Latest CMB and BAO Data",
    url="https://arxiv.org/abs/2512.10613", venue="Open Journal of Astrophysics (2026); arXiv Dec 2025",
    finding=("Planck + SPT + ACT + BICEP/Keck: n_s = 0.9682 +/- 0.0032 and r < 0.034 (95%). Adding DESI BAO leaves r "
             "unchanged but shifts n_s to 0.9728 +/- 0.0029, driven by marginal CMB-DESI differences."),
    relevance=("Primordial tensors: any inflation-like phase on the brane is bounded by r < 0.034. In high-energy RS "
               "inflation, tensor and scalar amplitudes are both enhanced (next entries). HDBLAST has no tensor "
               "prediction yet."),
    implication="constrains", queries=[29, 31])
add(id="SPT3G_r", section="CMB", year="2025",
    title="SPT-3G D1: Constraints on inflationary gravitational waves with two years of SPT-3G data",
    url="https://arxiv.org/abs/2505.02827", venue="Phys. Rev. D (Dec 2025)",
    finding="SPT-3G alone over the BICEP/Keck field: r < 0.25 (95%), sigma(r) = 0.067, with delensing.",
    relevance="An independent tensor channel. It is not competitive yet, but the South Pole Observatory programme targets sigma(r) = 0.001 (snippet).",
    implication="context", queries=[29])
add(id="BK_Planck_r", section="CMB", year="2022",
    title="Improved limits on the tensor-to-scalar ratio using BICEP and Planck data",
    url="https://arxiv.org/abs/2112.07961", venue="Phys. Rev. D 105, 083524 (2022)",
    finding="Planck + BICEP/Keck 2018 + BAO: r < 0.032 (95%). BK18 alone: r < 0.036.",
    relevance="The standard tensor bound, before the 2025 combination above.",
    implication="context", queries=[29])
add(id="SPT3G_D1", section="CMB", year="2025",
    title="SPT-3G D1: CMB temperature and polarization power spectra and cosmology from 2019 and 2020 observations of the SPT-3G Main field",
    url="https://arxiv.org/abs/2506.20707", venue="Phys. Rev. D (2025)",
    finding=("SPT-3G alone: H0 = 66.66 +/- 0.60. SPT + Planck + ACT: H0 = 67.19 +/- 0.38, sigma8 = 0.8137 +/- 0.0037. CMB "
             "alone shows no evidence beyond LambdaCDM. There is a 2.8 sigma CMB-vs-DESI DR2 difference within LambdaCDM."),
    relevance="The early-universe expansion history HDBLAST must reproduce, if it ever reaches a radiation era.",
    implication="context", queries=[28])
add(id="RS_tensor", section="CMB", year="2014",
    title="Tensor Perturbations from Brane-World Inflation with Curvature Effects",
    url="https://arxiv.org/abs/1308.5765", venue="Phys. Rev. D 89, 063501 (2014)",
    finding=("In RS braneworld inflation, tensor (and scalar) amplitudes are enhanced at high energy relative to 4D. "
             "Gauss-Bonnet and induced-gravity corrections suppress the RS enhancement; the Gauss-Bonnet term can break "
             "the standard consistency relation at high energy."),
    relevance=("If HDBLAST produced perturbations on the brane at energies comparable to the brane tension, the r and "
               "n_t relations would differ from 4D. With the registered H*ell = 0.069 the brane is in the low-energy "
               "regime (rho/lambda about 0.0024), so these corrections would be small."),
    implication="context", queries=[45],
    attribution="The enhancement/consistency statement appeared in the summary for several results; this paper's abstract matches the curvature-effect part.")
add(id="BBN_2024", section="BBN", year="2024",
    title="The 2024 BBN baryon abundance update",
    url="https://arxiv.org/abs/2401.15054", venue="arXiv Jan 2024",
    finding=("Delta N_eff from BBN ranges from -0.09 +/- 0.28 (one configuration: BAO+BBN abundances, PArthENoPE v3.0); "
             "results depend on the helium data and deuterium-burning rates."),
    relevance=("Expansion rate at T ~ 1 MeV. For a braneworld this bounds the rho^2/(2 lambda) term, which needs "
               "lambda^(1/4) far above MeV. That is automatic if tabletop bounds set lambda^(1/4) of a few TeV. It also "
               "bounds dark radiation: a Gaussian approximation gives Delta N_eff < 0.46, so rho_dr/rho_gamma < 0.10 "
               "(script)."),
    implication="constrains", queries=[32])
add(id="brane_DR_BBN", section="BBN", year="2017",
    title="New observational limits on dark radiation in braneworld cosmology",
    url="https://doi.org/10.1103/PhysRevD.95.083516", venue="Phys. Rev. D 95, 083516 (2017); arXiv 1706.03630",
    finding=("BBN restricts braneworld and particle dark radiation at 10 MeV to between -12.1% and +6.2% of the total "
             "background energy density. An older analysis (arXiv astro-ph/0203272) gave -1.23 <= rho_d/rho_gamma <= "
             "0.11 from BBN, narrowed to -0.41..0.105 with the CMB."),
    relevance=("A direct bound on the bulk Weyl term. Negative dark radiation, which is allowed in braneworlds, is "
               "bounded too. This applies to any HDBLAST scenario in which the blast leaves bulk mass or energy."),
    implication="constrains", queries=[33],
    attribution="The snippet summary merged two papers; the -12.1%/+6.2% figure is attributed to the 2017 paper, the -1.23..0.11 figure to astro-ph/0203272.")

# ---------------------------------------------------------------- DESI ------
add(id="DESI_DR2", section="DESI", year="2025",
    title="DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints",
    url="https://arxiv.org/abs/2503.14738", venue="arXiv Mar 2025 (+ DESI guide page https://www.desi.lbl.gov/2025/03/19/desi-dr2-results-march-19-guide/)",
    finding=("w0waCDM is preferred over LambdaCDM at 3.1 sigma (DESI + CMB) and at 2.8-4.2 sigma when supernovae are added. "
             "The best fit has w > -1 today and w < -1 in the past, crossing -1 near z = 0.5."),
    relevance=("The registered model's relaxing fate is an empty RS de Sitter brane with w = -1 exactly. At the registered "
               "detuning, H*ell = 0.069, so identifying that brane with today's dark energy needs ell of about 372 Mpc "
               "(excluded by tabletop tests by a factor 3e29) or delta of about 1e-62 (script). HDBLAST therefore makes "
               "no DESI prediction. A confirmed evolving w would need an extra late-time sector."),
    implication="context", queries=[24])
add(id="DES_Dovekie", section="DESI", year="2025",
    title="The Dark Energy Survey Supernova Program: A Reanalysis Of Cosmology Results And Evidence For Evolving Dark Energy With An Updated Type Ia Supernova Calibration",
    url="https://arxiv.org/abs/2511.07517", venue="MNRAS (2026); arXiv Nov 2025",
    finding=("Recalibration (DES-Dovekie) reduces the significance of rejecting LambdaCDM from 4.2 sigma (DES-SN5YR) to "
             "3.2 sigma with CMB + DESI DR2. Only a weak Bayesian preference for w0wa remains."),
    relevance="The evolving-dark-energy signal is not settled, so it is not yet a firm constraint on any braneworld late-time sector.",
    implication="context", queries=[25])
add(id="DESI_Lya_2026", section="DESI", year="2026",
    title="New DESI DR2 Lyman-alpha Results Shed Light on Dark Energy",
    url="https://www.desi.lbl.gov/2026/07/30/new-desi-dr2-lyman-alpha-results-shed-light-on-dark-energy/",
    venue="DESI collaboration news, 30 July 2026",
    finding=("The new DR2 Lyman-alpha forest measurements agree with the standard model. DESI notes the hints of "
             "evolving dark energy could fade, or a more complex model may be needed."),
    relevance="As above: late-time dark energy is currently outside the registered model's reach.",
    implication="context", queries=[47])
add(id="brane_DESI", section="DESI", year="2025",
    title="Braneworld dark energy in light of DESI DR2",
    url="https://arxiv.org/abs/2507.07193", venue="JCAP 11 (2025) 018",
    finding=("Thawing scalar fields on a (4+1)-dimensional ghost-free 'phantom brane' (induced-gravity, DGP-type normal "
             "branch) show an effective phantom-divide crossing. A quartic potential explains the DESI DR2 constraints well."),
    relevance=("A braneworld route to DESI-like behaviour exists, but it uses induced gravity on the brane and a "
               "cosmological crossover scale. The registered model is RS2-type (no induced gravity), so this does not "
               "carry over without new terms."),
    implication="context", queries=[26])

# ---------------------------------------------------------------- GW / LVK --
add(id="LVK_O4a_SGWB", section="LVK", year="2025",
    title="Upper Limits on the Isotropic Gravitational-Wave Background from the first part of LIGO, Virgo and KAGRA's fourth Observing Run",
    url="https://arxiv.org/abs/2508.20721", venue="Phys. Rev. D, doi 10.1103/wq57-sjt2 (update to April 2025: arXiv 2608.23477)",
    finding=("No background detected. Omega_GW(25 Hz) <= 2.0e-9 for a 2/3 power law and <= 2.8e-9 for a flat spectrum "
             "(95%). The update through April 2025 quotes the same values."),
    relevance=("The frozen knee's f^-3 tail gives Omega(25 Hz) = 2.28e-35, 26 orders of magnitude below the limit "
               "(script), so the LVK band does not test the knee. If the brane scale is set by tabletop bounds, the "
               "static branch's natural horizon-scale frequency is about 5e-5 to 8e-5 Hz (script). That is below the "
               "LVK band, in the sub-mHz range targeted by space interferometers."),
    implication="constrains", queries=[41])
add(id="GWTC4_TGR", section="LVK", year="2026",
    title="GWTC-4.0: Tests of General Relativity. II. Parameterized Tests",
    url="https://arxiv.org/abs/2603.19020", venue="arXiv Mar 2026",
    finding=("42 confident O4a signals plus 49 earlier events. Graviton-mass bound m_g < 1.92e-23 eV/c^2 (90%), 1.16 times "
             "better than GWTC-3."),
    relevance=("In RS2 the 4D graviton is massless and the massive KK modes are micron-scale, so this bound does not "
               "touch the registered model at tabletop-consistent ell. It excludes cosmological-ell variants in which "
               "massive-graviton dispersion would appear."),
    implication="context", queries=[42])
add(id="Pardo2018", section="LVK", year="2018",
    title="Limits on the number of spacetime dimensions from GW170817",
    url="https://arxiv.org/abs/1801.08160", venue="JCAP 07 (2018) 048",
    finding=("Gravitational-wave amplitude damping from leakage into large extra dimensions would bias GW distances; "
             "comparing GW and electromagnetic distances for GW170817 limits this and is consistent with D = 4."),
    relevance="Foundational method. It is relevant only for non-compact extra dimensions with a crossover scale not far above tens of Mpc.",
    implication="method", queries=[39])
add(id="GWTC3_dims", section="LVK", year="2023",
    title="Constraining the number of spacetime dimensions from GWTC-3 binary black hole mergers",
    url="https://arxiv.org/abs/2112.07650", venue="Phys. Rev. D 107, 084033 (2023)",
    finding=("Using the pair-instability mass gap: screening-independent D = 3.95 (+0.09/-0.07) (68%). With a screening "
             "scale: D = 4.23 (+1.50/-0.57) and log10 R_c/Mpc = 4.14 (+0.55/-0.86)."),
    relevance=("Leakage tests constrain HDBLAST only if the brane's effective AdS or crossover length is cosmological. "
               "The script's 'dark-energy identification' (ell about 372 Mpc) is in that regime, but tabletop tests "
               "already exclude it."),
    implication="constrains", queries=[40])
add(id="GWTC4_darksiren_dims", section="LVK", year="2026",
    title="Searching for Extra Dimensions with Gravitational Waves: Dark-Siren Constraints from GWTC-4",
    url="https://arxiv.org/abs/2606.14549", venue="arXiv June 2026",
    finding=("GWTC-4 dark sirens with a narrow H0 prior give D = 4.38 (+1.91/-1.01), consistent with D = 4. The crossover "
             "scale is poorly constrained, and the D limit weakens when the crossover exceeds the distances of the "
             "observed population."),
    relevance="Same role as the entry above: only large-crossover variants are tested.",
    implication="constrains", queries=[40])
add(id="GW170817_5DdS", section="LVK", year="2020",
    title="Constraint on the radius of five-dimensional dS spacetime with GW170817 and GRB 170817A",
    url="https://arxiv.org/abs/2001.06581", venue="arXiv Jan 2020",
    finding="Title-level only: a constraint on a 5D de Sitter radius from the GW-gamma-ray arrival-time comparison. No number in the snippet.",
    relevance=("A 'shortcut through the bulk' arrival-time test is the natural GW170817-type test for a de Sitter brane. "
               "It should be read before any HDBLAST propagation claim."),
    implication="method", queries=[39])

# ---------------------------------------------------------------- tabletop --
add(id="EotWash2020", section="ISL", year="2020",
    title="New Test of the Gravitational 1/r^2 Law at Separations down to 52 micrometres",
    url="https://arxiv.org/abs/2002.11761", venue="Phys. Rev. Lett. 124, 101101 (2020)",
    finding=("Torsion balance with 18- and 120-fold attractors at 52 um to 3.0 mm separations. Newtonian gravity fits; "
             "gravitational-strength Yukawa interactions are limited to ranges below 38.6 um (95%). Motivated by the "
             "dark-energy length of about 85 um."),
    relevance=("Sets the physical unit of the registered model. The RS2 correction is a power law (1 + 2 ell^2/(3 r^2)), "
               "not a Yukawa, so 38.6 um is used as a proxy for the upper bound on the AdS radius ell. With ell <= 38.6 "
               "um: M5 >= 3.1e8 GeV, lambda^(1/4) >= 5.5 TeV, and the static branch's brane vacuum energy has "
               "rho^(1/4) >= 1.2 TeV (script; conditional). The 2026 review's power-law limits (below) are the correct "
               "input for a final number."),
    implication="constrains", queries=[34, 37])
add(id="LivingReview_RS", section="ISL", year="unknown",
    title="Brane-World Gravity (Living Reviews in Relativity)",
    url="https://pmc.ncbi.nlm.nih.gov/articles/PMC5479361/", venue="Living Reviews in Relativity (PMC copy); volume/year not shown in the snippet",
    finding=("Table-top tests show no deviation from Newton's potential above 0.1 mm, so the AdS5 curvature scale "
             "ell < 0.1 mm (ell^-1 > 1e-12 GeV), giving M5 > 1e8 GeV."),
    relevance="Control: the script reproduces M5 = 2.3e8 GeV at ell = 0.1 mm with M_Pl^2 = M5^3 ell.",
    implication="constrains", queries=[36])
add(id="ISL_review_2026", section="ISL", year="2026",
    title="Short-Range Tests of the Gravitational Inverse-Square Law",
    url="https://arxiv.org/abs/2605.18212", venue="arXiv May 2026 (review)",
    finding=("A consistent formalism across length scales, updates from the past decade, and a comparison of tabletop "
             "and collider results for both Yukawa and power-law potentials, including extra-dimension models."),
    relevance="It contains the power-law (RS-type) limits needed to replace the Yukawa proxy for ell.",
    implication="method", queries=[35, 38])

# ---------------------------------------------------------------- H0 / JWST -
add(id="H0DN", section="H0", year="2025",
    title="The Local Distance Network: a community consensus report on the measurement of the Hubble constant at 1% precision",
    url="https://arxiv.org/abs/2510.23823", venue="A&A (April 2026)",
    finding=("H0 = 73.50 +/- 0.81 km/s/Mpc from a blinded community distance network (Cepheids, TRGB, Miras and more). "
             "This is 7.1 sigma from the early-universe LambdaCDM value 67.24 +/- 0.35. A simple local systematic is "
             "judged unlikely."),
    relevance=("Extra dark radiation from a bulk remnant would raise the CMB-inferred H0 only slightly, and ACT DR6 "
               "disfavours that route. HDBLAST makes no H0 prediction. One practical effect: converting a strain "
               "amplitude into Omega_gw scales as h^-2. With h = 0.735 instead of 0.674, the knee's h_c at 1/yr would be "
               "2.62e-15 instead of 2.40e-15 (script). The frozen Omega_k is tied to h = 0.674."),
    implication="context", queries=[43, 46])
add(id="JWST_crowding", section="H0", year="2024",
    title="JWST Observations Reject Unrecognized Crowding of Cepheid Photometry as an Explanation for the Hubble Tension at 8 sigma Confidence",
    url="https://arxiv.org/abs/2401.04773", venue="ApJL (2024)",
    finding="More than 1000 Cepheids measured with both HST and JWST: mean distance difference -0.01 +/- 0.03 mag. Crowding is rejected as the explanation at 8.2 sigma.",
    relevance="Context for the tension's robustness.",
    implication="context", queries=[43])
add(id="CCHP_2025", section="H0", year="2025",
    title="The Chicago Carnegie Hubble Program: Improving the Calibration of SNe Ia with JWST Measurements of the Tip of the Red Giant Branch",
    url="https://arxiv.org/abs/2503.11769", venue="arXiv Mar 2025",
    finding=("The search summary quotes H0 = 69.8 +/- 1.9 from CCHP's JWST Cepheid/TRGB/JAGB work, consistent with the "
             "CMB value. It was not verified which CCHP paper this number is from."),
    relevance="The local-ladder disagreement is not fully settled, so HDBLAST should not be tuned to either H0.",
    implication="context", queries=[43],
    attribution="Number attribution uncertain (see finding).")
add(id="JWST_z10", section="JWST", year="2025",
    title="Beyond No Tension: JWST z > 10 Galaxies Push Simulations to the Limit",
    url="https://arxiv.org/abs/2509.07695", venue="arXiv Sep 2025",
    finding=("Search summary: spectroscopic galaxies up to z of about 14.4 (MoM-z14) and 14.3 (GS-z14). Hints of an "
             "overabundance of bright or massive early galaxies; implied star-formation efficiencies are high. There is "
             "an ongoing debate over whether this reflects new physics or systematics (stellar masses, AGN, selection)."),
    relevance=("Tests small-scale primordial power and early structure growth. HDBLAST has no perturbation spectrum, so "
               "it is untested here. A future HDBLAST spectrum with enhanced small-scale power would be checked "
               "against these counts."),
    implication="context", queries=[44],
    attribution="Galaxy redshifts came from the summary of several results.")
add(id="LCDM_not_dead", section="JWST", year="2023",
    title="LambdaCDM not dead yet: massive high-z Balmer break galaxies are less common than previously reported",
    url="https://arxiv.org/abs/2310.03063", venue="arXiv Oct 2023",
    finding="Title-level: massive high-z Balmer-break galaxies are less common than earlier claimed.",
    relevance="Early 'impossible galaxy' claims have weakened. They are not a motivation for a non-standard origin.",
    implication="context", queries=[44])
