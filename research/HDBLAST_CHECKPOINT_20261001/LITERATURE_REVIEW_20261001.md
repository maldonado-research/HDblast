# Targeted literature review — 1 October 2026

## Retrieval and scope

Five primary papers were retrieved by the connected Firecrawl research tools: metadata with `research_inspect_paper` and full-text passages with `research_read_paper`. This is a bounded review of particle production, renormalization, backreaction and bulk/KK budgets, not an exhaustive search or external novelty assessment.

Four Consensus search queries for 2024–2026 brane reheating, time-dependent-mass backreaction, stress/current renormalization and brane species cutoffs returned no records before a stalled batch was interrupted after approximately 2,503 seconds. No Consensus paper fetch completed. A later batch of direct arXiv abstract scrapes also returned no result before interruption after approximately 1,114 seconds. Failed calls supply no evidence. The five references below come from completed Firecrawl research retrievals, not those failed batches.

## 1. Higher-dimensional stress bookkeeping

**Anamitra Paul and Sonia Paban (2026), “Stress-Energy Tensor of a Scalar Field on a Product Spacetime with a Time-Dependent Compact Dimension.”**
[arXiv:2603.12444](https://arxiv.org/abs/2603.12444);
[10.1007/JHEP09(2026)188](https://doi.org/10.1007/JHEP09(2026)188).

Fetched text: “The computation we present assumes that the scalars are free and non-minimally coupled to gravity.”

The paper studies FLRW × S¹ with ordinary and compact scale factors a(t), b(t), in d=3,4 spatial dimensions. The mode friction contains (d−1) adot/a + bdot/b. The physical compact spectrum remains discrete while subtraction uses a continuum measure, preserving topology-dependent terms. The authors report a finite stress tensor with vanishing covariant divergence and a massless conformal equal-scale-factor FLRW limit.

**Action for HDBLAST:** check compact-direction pressure, dimensional reduction and the complete stress conservation law when introducing quantum matter. Treat time-dependent compact geometry as part of the mode normalization.

**Restriction:** an unwarped product with a circle is not this warped Z₂ shell model. Their stress approximation substitutes WKB modes where exact modes are unavailable; its special-limit checks do not establish nonadiabatic production accuracy on our trajectory. Its free-field counterterms do not directly solve an interacting phi-dependent mass problem. Cosmological implications are deferred in that paper.

## 2. Derive the force and energy together

**Nathan Herring and Daniel Boyanovsky (2025; retrieved record updated 6 February 2026), “Misalignment dynamics of Scalar Condensates with Yukawa coupling: Particle and Entropy Production.”**
[arXiv:2511.22465](https://arxiv.org/abs/2511.22465).
No DOI was supplied by the retrieved metadata.

Fetched text: “This set of fully renormalized, self-consistent, energy conserving equations are some of the main results of this study.”

Their large-Nf Yukawa treatment jointly renormalizes particle energy and condensate force. Equations III.136–III.150 assemble renormalized energy, mode equations and the condensate equation. The ultraviolet logarithm and field renormalization, together with an initial-time force singularity, are part of the consistent initial-value problem. The reported conserved energy includes condensate kinetic/effective-potential energy and a renormalized particle-energy integral.

**Action for HDBLAST:** derive renormalized rho, pressure and scalar current from the same effective-action/subtraction prescription. Check their work exchange on the numerical history, and specify state preparation rather than treating it as a harmless numerical choice.

**Restriction:** fermionic Yukawa production at large Nf in the displayed nonexpanding condensate system differs from our bosonic spectator on a curved shell. The paper does not present a completed numerical dynamical study. Its counterterms cannot simply be transplanted into m_chi²=m0²+g²(phi−phi_star)².

## 3. Resolve spectra over the actual history

**Ayan Chakraborty, Simon Clery, Md Riajul Haque, Debaprasad Maity and Yann Mambrini (2025), “Generalizing the Bogoliubov vs Boltzmann approaches in gravitational production.”**
[arXiv:2503.21877](https://arxiv.org/abs/2503.21877).
No DOI was supplied by the retrieved metadata.

Fetched abstract: “In the UV regime, oscillations of the inflaton background lead to interference terms that explain the high-frequency oscillations in the spectrum.”

The work computes scalar gravitational production across inflation/reheating with Bogoliubov evolution, analytical infrared/ultraviolet limits and Boltzmann comparisons. The fetched full text discusses ultraviolet spectral agreement and oscillatory interference, including amplitude differences in a numerical comparison.

**Action for HDBLAST:** integrate canonical modes over the actual history with infrared and ultraviolet convergence checks. If multiple crossings occur, retain phases and interference. The existing positive-pulse control tests a solver; it is not a spectrum for the cosmological history.

**Restriction:** their gravitational-production background differs from direct phi-dependent mass production. Their exponential infrared suppression for m_chi/H_e approximately above 3/2, spectral powers and quoted reheating equation-of-state conditions are model-specific; none is a universal HDBLAST threshold.

## 4. Account for the bulk and KK channels

**Luis A. Anchordoqui, Ignatios Antoniadis and Jules Cunat (2026), “Cosmological history after higher dimensional inflation.”**
[arXiv:2606.20486](https://arxiv.org/abs/2606.20486).
Retrieved record dated June 2026; no DOI supplied.

Fetched text: “the enormous multiplicity of accessible KK states can lead to substantial graviton production in the early universe.”

The paper discusses a normalcy temperature for an essentially unpopulated bulk and conventional four-dimensional cosmology, low-temperature reheating and suppressed bulk-graviton channels. For its interval compactification it gives M_star approximately (m_KK/pi)^(d/(2+d)) M_p^(2/(2+d)); its one-extra-dimension thermal estimate gives T_r/GeV approximately no greater than (R_perp/micrometre)^(-1/3).

**Action for HDBLAST:** include a bulk/KK energy budget alongside shell matter. Establish the compactification and scale dictionary before comparing occupied frequencies with an effective gravity cutoff.

**Restriction:** this setup assumes stabilization after higher-dimensional inflation. Its temperature estimate depends on a particular KK spectrum, thermal production and long-lived relic assumptions, and its reheating description is phenomenological. Micron radii or its numerical temperature bound cannot be imposed on the current model without deriving the mapping.

## 5. Separate thermal and nonthermal constraints

**Luis Anchordoqui, Ignatios Antoniadis and Dieter Lüst (2025), “Two Micron-Size Dark Dimensions.”**
[arXiv:2501.11690](https://arxiv.org/abs/2501.11690).
No DOI supplied.

Fetched text: “The non-thermal production of KK modes accompanying the inflaton decay would further constrain” the compactification radius.

The work distinguishes cooling, BBN, relic density, thermal production and nonthermal production. Its two-large-dimension viability depends on KK intratower decays through absent isometries and a constrained radiation-onset temperature.

**Action for HDBLAST:** compute branching into shell/bulk modes, multiplicities and relic evolution. A radiation-like classical component or a particle number does not establish thermalization or cosmological viability.

**Restriction:** this is a particular dark-dimension scenario, with some constraints drawn from earlier literature, not a derived bound on our Einstein–scalar/Z₂-shell construction.

## Resulting priority

The strongest immediate methodological pair is Herring–Boyanovsky for matched source/energy renormalization and Paul–Paban for dimensional stress bookkeeping. Combined with the Bogoliubov study, they support the next test specified in [MATTER_ACTION_AND_CLOCK_CONTRACT.md](MATTER_ACTION_AND_CLOCK_CONTRACT.md): fixed action/state, mapped clocks, modes on the saved shell geometry, matched renormalized stress/current, and a paid-work ledger. Bulk/KK channels follow once the geometry and scale mapping are fixed.

These are useful inputs to a discriminating experiment. None demonstrates HDBLAST reheating, supplies external novelty for the standard canonical transformations, or promotes the sampled 0.31446-e-fold stand-in result to a microscopic hot Big Bang.
