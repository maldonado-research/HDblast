# Einstein Audit — detailed research account and equation register

**Author:** Ricardo Maldonado. **GitHub documentation edition:** 2026-09-30.
**Underlying research snapshot:** 2026-09-08.

This is a public, FTL-only derivative of the author's supplied research handoff. It preserves the detailed account, equation groups, source references, history, claim register, glossary, and outstanding tasks. It is a working-hypothesis/methods publication, with **no empirical FTL detection, beyond-GR result, or external peer review documented**.

The underlying handoff reports earlier inspection of public archives; those inspections were not repeated for this GitHub edition. The private v6.2/v6.3 files remain unavailable. All later test counts, numerical repairs, and named certificates retain their report-only status. “Current,” “latest,” and “current-turn” inside the preserved account refer to its September 8 source snapshot, not a verified September 30 literature or release survey.

Editorial changes: platform placement/payment instructions were omitted; Markdown math delimiters were adapted for GitHub; source equation contents were retained. This derivative does not claim every referenced source or dataset is embedded. See [sources and missing assets](SOURCES_AND_REPRODUCIBILITY.md).

Priority issues: historical **Q** means either occupancy or first severity moment; the SGWB-to-extra-radiation normalization is underspecified; the unit Einstein ball and x-channel coupling are proposed conventions. The private v6.3 mathematics is pending recovery and review.

---


### A.3.0 Scope and evidence tiers

This report isolates **only** the Faster-Than-Light (FTL) / Einstein Audit research program. The author repeatedly instructed that higher-dimensional-blast, Big-Bang-spark, Theory of Everything, dark-matter/dark-energy, black-hole, antimatter, and ε-lattice projects are separate and must not be merged into the FTL program. See `project_sources/09-theory-faster-than-light-instructions.txt`, lines 1–2 and 3510–3517; `project_sources/02-dblast-black-holes-toe-dmde-equations.md`, lines 1281–1285, 1467–1480, and 1966–1979.

This is a research handoff, not a scientific endorsement. The record supports a **conditional mathematical audit framework and a progressively hardened reproducibility protocol**. It does **not** support a detection of FTL propagation, a violation of relativity, an empirical beyond-GR result, or a claim of peer review.

#### Evidence tiers used here

| Tier | Meaning in this report |
|---|---|
| Archived primary record | Public Zenodo metadata and files deposited by the author; establishes what was archived, not scientific validity. |
| Source-text reconstruction | Equations and claims read from the deposited paper, Markdown notes, and code without running the code. |
| Current-turn report | Statements in the visible v6.2/v6.3 continuation; useful for history, but the private v6.2/v6.3 kit itself was not accessible here. |
| Established background | Standard SR, probability, statistics, or GW/cosmology concepts; their particular use and normalization in the Einstein Audit may still be proposed. |
| Proposed construction | Author/conversation-defined quantities such as the Einstein vector, unit Einstein ball, coupling (x=20\Delta N_{\rm eff}J), and phase labels. |
| Internal validation | Tests, synthetic examples, code review, or separate AI-agent review. This is not experimental confirmation or external peer review. |

### A.3.1 Research identity, names, authorship, and status

#### A.3.1.1 Canonical identity

The program's stable short identity is:

> **Einstein Audit for Faster-Than-Light Hidden Sectors (FTL-only)**

Earlier and internal names include:

- **Einstein Audit for Faster-Than-Light Hidden Sectors (MIN)**, where “MIN” means a minimal/analysis-only framework.
- **FTL Einstein Audit**.
- **A1 REALDATA_APPROX**, the approximate-summary-anchor stage.
- **A2**, the intended posterior-level real-data stage.
- **J-first SR envelope**, the v1.2.0 shift from the small-(delta) proxy to the exact severity functional (J[\rho]).
- **Testability Gate** or **bound-existence condition**.
- Internal toolchain labels **v5.8**, **v5.9**, **v6.0**, **v6.1**, **v6.2**, and **v6.3**. These are not the same as the public Zenodo semantic versions.

The v1.2.0 archive contains filenames marked `v1_3` and `v1_4`; the Zenodo description explicitly says these are internal labels and that the archival release version is **v1.2.0**. Do not infer an unpublished public v1.3 or v1.4 from those filenames.

#### A.3.1.2 Author and collaborators

- Author/creator: **Ricardo Maldonado** (also formatted in metadata as `Maldonado, Ricardo` or `Ricardo, Maldonado`).
- ORCID in v1.2.0 and v1.2.1 metadata: **0009-0009-3937-6527**.
- No scientific coauthor or collaborator is explicitly credited in the accessible FTL records.
- NANOGrav, EPTA, IPTA, PPTA, and LVK/GWTC data collaborations are cited or proposed data providers, not collaborators or coauthors of this program.
- The current-turn “independent reviewers” are unnamed. Their review is best described as an internal/separate-agent mathematical and adversarial review, not external independent replication or journal peer review.

#### A.3.1.3 Field, question, scope, and proposed contribution

Fields: special-relativistic kinematics, hypothetical superluminal/hidden-sector phenomenology, early-Universe radiation constraints, stochastic gravitational-wave backgrounds, pulsar-timing-array posterior comparison, black-hole ringdown tests of GR, statistical computing, and information-theory toy models.

Central conditional question:

> If an otherwise hidden sector permits two-way superluminal signaling, how can its special-relativistic causality risk be quantified at population level and jointly constrained by an assumed cosmological radiation budget, PTA posterior consistency, and ringdown GR-deviation summaries?

Proposed contribution:

1. Define an exact antitelephone wedge and a nonnegative mean-square paradox-severity functional (J[\rho]).
2. Compress three heterogeneous channels into an author-defined “Einstein vector” and “Einstein distance.”
3. Derive algebraic envelopes such as (J\le x_{\max}/(20\Delta N_{\rm eff})).
4. Identify a conditional “Testability Gate”: the constructed audit can impose a finite severity bound only if a model's reachable ceiling (J_{\max}) intersects the envelope.
5. Build a fail-closed A2 provenance and inference protocol intended to prevent demo, synthetic, summary-only, or malformed inputs from being promoted as real-data results.

The program is not a dynamical theory of an FTL sector. The accessible record contains no hidden-sector Lagrangian, field content, interaction Hamiltonian, microcausal prescription, preferred-frame dynamics, production mechanism for the proposed SGWB, or demonstrated map from an FTL degree of freedom to a PTA/ringdown observable.

#### A.3.1.4 Public version history

Public metadata was read from the official Zenodo Records API on 2026-09-08. The cited public archive members were also accessed read-only from the official records during the research audit, but they were not among the user-supplied attachments and their downloaded container bytes are not preserved in this task workspace. Any later byte-level audit must re-retrieve and re-hash the official records; the filenames, sizes, MD5 values, page/line references, and extracted conclusions below document this audit pass rather than embedding those remote packages.

| Public version | Date | Record / DOI | Exact archived title | Status |
|---|---:|---|---|---|
| v1.0.0 | 2025-11-26 | [record 17726160](https://zenodo.org/records/17726160), DOI `10.5281/zenodo.17726160` | *Einstein Audit for Faster‑Than‑Light Hidden Sectors (MIN): Causality Indices, Delta N_eff Guardrail, PTA Coherence & Ringdown GR Tests* | Initial analysis-only MIN kit; synthetic/demo inputs; no FTL evidence claim. |
| v1.1.0 | 2025-12-13 | [record 17926499](https://zenodo.org/records/17926499), DOI `10.5281/zenodo.17926499` | *Einstein Audit for Faster‑Than‑Light Hidden Sectors (MIN): Causality, Cosmology, and Multi‑Band Gravitational‑Wave Consistency* | Added activated-cosmology demonstrations, scans, budget fractions, covariance appendix, and toy capacity/flavor extensions. |
| v1.2.0 | 2026-02-05 | [record 18499411](https://zenodo.org/records/18499411), DOI `10.5281/zenodo.18499411` | *Einstein Audit for Faster‑Than‑Light Hidden Sectors (FTL‑only): J‑First SR Envelope, Testability Gate, and A2 Real‑Data Scaffolding* | J-first methods upgrade and A2 scaffolding; A1 values still `REALDATA_APPROX`; no strict real-data result. |
| v1.2.1 | 2026-05-15 | [record 20218030](https://zenodo.org/records/20218030), DOI `10.5281/zenodo.20218030` | *Einstein Audit for Faster-Than-Light Hidden Sectors (FTL-only): A2 Strict Protocol Freeze and Provider-Return Gate* | Latest public release; methods/toolchain freeze only; explicit claim lock; no strict real-data A2 posterior constraints. |

Concept record / all-version DOI: [record 17726159](https://zenodo.org/records/17726159), `10.5281/zenodo.17726159`.

#### A.3.1.5 Post-public private status

The visible continuation reports a private **v6.2 fail-closed integrity layer** followed by a private **v6.3 numerical-and-evidence closure kit**. The final reported v6.3 state was:

- 81/81 tests passed after fresh extraction.
- 534 boolean/null/nonfinite/missing numeric mutations across 14 receipt types produced zero bypasses.
- Two deterministic ZIP builds were byte-identical.
- 85 checksum-covered entries; 86 ZIP members including the checksum file.
- Reported ZIP SHA-256: `f2a14473750a952e18b5c9668888940f60313654d6f2547a512212914a46d9bd`.
- The checked-in real-data candidate failed closed, exit code 2, with 123 unmet machine checks.
- Decision: **GO** for a private methods/preflight release; **NO-GO** for public empirical v1.3.0, an FTL claim, a beyond-GR claim, or a discovery announcement.

These are current-turn reports only. The v6.2/v6.3 ZIP and validation report were not present in the accessible workspace, so their contents and checksum were not independently inspected for this handoff. No new Zenodo version was reported as published. The binding public source therefore remains v1.2.1, while v6.3 supplies a pending correction layer that must be reviewed from its actual files before corpus ingestion.

### A.3.2 Current scientific account

#### A.3.2.1 Established background versus proposed framework

##### Established background

- Lorentz transformations and the colinear velocity-addition rule are standard special relativity.
- In a tachyon-antitelephone setup, reciprocal controllable FTL signaling in relatively moving frames can lead to reception before emission under specified assumptions.
- A stochastic gravitational-wave background is radiation-like, and its integrated energy density can be related to an effective extra-radiation parameter, subject to convention-dependent normalization.
- PTA collaborations infer posteriors for common stochastic-process parameters; ringdown analyses constrain deviations from GR.
- Density overlap, Bayes classification with equal class priors, cross-fitting, Brier scores, covariance quadratic forms, DKW inequalities, Makarov bounds, and Appell (F_1) functions are established mathematical/statistical tools.

##### Proposed by this research

- The particular paradox-strength function (P(\beta,W)) as a severity score.
- The population indices (Q[\rho]), (J[\rho]), and (C[\rho]), including a currently unresolved historical definition of (Q).
- The “causality–cosmology action” (S=\varepsilon_{\rm energy}J) and normalization (S_0=1/6).
- The Einstein vector (X=(x,y,z)), Euclidean “Einstein distance,” and the conventional unit-ball rule (D_E\le1).
- The choice (x=20\Delta N_{\rm eff}J), which multiplies an FTL-severity statistic by a radiation-budget statistic. This is a diagnostic design choice, not a derived physical interaction law.
- The use of a PTA overlap score and a pooled ringdown (z)-score as orthogonal coordinates in the same audit geometry.
- The practical/latent/closed phase terminology.
- The information-channel toy model and its interpretation as an “Einstein capacity bound.”

##### Conditional mathematical consequences

Given the definitions and assumptions above:

- the antitelephone boundary follows algebraically;
- (D_E\le1) yields the no-free-FTL-lunch inequality and the J-first envelope;
- (P\to\beta^2) as (W\to\infty), so (J\to E[\beta^4]) for the stated severity functional and a fixed beta population;
- the testability threshold follows by comparing (J_{\rm bound}) to (J_{\max});
- the classifier-overlap identity follows exactly for equal mixture weights and the true Bayes class probabilities;
- the four phase masses form an exact partition once the draw clouds, denominators, and boundary conventions are fixed.

##### Numerical/computational results

All numerical results recovered here are synthetic, approximate-summary anchors, toy scans, or software tests. None is a measurement of FTL.

##### Speculation and unresolved possibilities

The hidden FTL sector itself, its coupling to gravity/GWs, and any link to actual PTA or ringdown data remain hypothetical. The accessible record does not establish that a common PTA process or a ringdown posterior is generated by an FTL sector.

#### A.3.2.2 Physical setup and mechanism actually specified

The specified physical content is deliberately minimal:

1. Alice and Bob are inertial observers in flat spacetime, with colinear relative speed $v$, $\beta=v/c$, and $|\beta|<1$.
2. Each can send a return signal at the same superluminal factor (W=w/c>1) in their own rest frame.
3. Alice sends at proper time (T>0); Bob returns immediately.
4. The round-trip ratio (t_2/T) defines a paradox-capable region and a proposed severity score.
5. A normalized population (\rho(\beta,W)) represents an ensemble of such configurations.
6. The hidden sector is assumed to have or source an SGWB summarized by (\Omega_{\rm gw}(f)) and (\Delta N_{\rm eff}), but no production mechanism is derived.
7. PTA and ringdown products are treated as consistency diagnostics, not direct FTL detections.

The construction assumes positive colinear $\beta$ in the population integrals. Directional distributions, asymmetric outbound/return speeds, delayed response, preferred frames, curved spacetime, dispersive propagation, stochastic signaling, quantum measurement constraints, and interacting-field dynamics are not derived in the archived core.

### A.3.3 Stable equation and derivation register

All core variables are dimensionless unless units are stated. The EA-* identifiers below are editorial identifiers introduced in the September 8 handoff; they are not original paper equation numbers. The archived paper's equation numbers are included where available.

#### A.3.3.1 Special-relativistic geometry

##### EA-SR-001 — Superluminal factor

$$
W\equiv \frac{w}{c}>1,\qquad u=Wc.
$$

- Symbols: (w) or (u), signal speed; (c), invariant light speed; (W), dimensionless FTL factor.
- Status: definition proposed for this audit; (w>c) is hypothetical.
- Source: `FTL_v1_3_Compiled_Technical_and_Supplements.pdf`, p. 2, Eq. (2); local handoff `02`, lines 1297–1300 and 1994–1998.

##### EA-SR-002 — Lorentz factor and velocity addition

$$
\gamma(\beta)=\frac{1}{\sqrt{1-\beta^2}},
\qquad
u'=\frac{u-v}{1-uv/c^2},
\qquad \beta\equiv \frac vc.
$$

- Status: standard SR for colinear inertial frames.
- Domain: (|\beta|<1); the formal transformation can be applied to a hypothetical spacelike signal, though that does not establish such a signal exists.
- Source: archived PDF, p. 2, Eq. (3).

##### EA-SR-003 — Symmetric antitelephone round-trip time

$$
t_2=T\frac{w^2(1-\beta^2)}{(w-v)^2}
=T\frac{W^2(1-\beta^2)}{(W-\beta)^2}.
$$

- Symbols: (T>0), Alice's send time in the setup; (t_2), return-arrival time in Alice's frame.
- Assumptions: immediate reply, equal FTL speed in each sender's own frame, colinear motion, Minkowski spacetime.
- Status: adapted standard tachyon-antitelephone geometry.
- Source: archived PDF, p. 2, Eq. (4).

##### EA-SR-004 — Exact paradox-wedge boundary

$$
\beta_c(W)=\frac{2W}{1+W^2},
\qquad \text{paradox-capable when }\beta>\beta_c(W).
$$

Derivation from EA-SR-003:

$$
t_2<T
\iff W^2(1-\beta^2)< (W-\beta)^2
\iff \beta\big[(1+W^2)\beta-2W\big]>0.
$$

For the archived domain (\beta>0), this gives the stated boundary.

- Checks: (\beta_c(1)=1), so the wedge vanishes at luminal speed; (\beta_c(W)\to0) as (W\to\infty).
- Source: archived PDF, p. 2, Eq. (5); `02`, lines 2000–2004.

##### EA-SR-005 — Proposed paradox strength

$$
P(\beta,W)=1-\frac{W^2(1-\beta^2)}{(W-\beta)^2}
=1-\frac{t_2}{T}.
$$

For severity integrals the source uses (P_+(\beta,W)=\max(P,0)); (P>0) inside the wedge, (P=0) on the boundary, and (P<0) outside.

- Status: proposed dimensionless severity definition, built from a standard SR time ratio.
- Range: inside the stated domain (0\le P_+\le1); this should be retained as a consistency check.
- Source: archived PDF, p. 2, Eq. (6); `02`, lines 2005–2011.

#### A.3.3.2 Population functionals

##### EA-POP-001 — Population normalization

$$
\int_0^1 d\beta\int_1^\infty dW\,\rho(\beta,W)=1.
$$

- $\rho$ is a dimensionless joint probability density with reciprocal dimensions matching $d\beta\,dW$, both dimensionless.
- The archived support excludes negative relative velocities and angular variables.
- Source: archived PDF, p. 3, Eq. (7).

##### EA-POP-002A — Historical wedge-occupancy definition of (Q)

$$
Q_{\rm occ}[\rho]
=\int_0^1d\beta\int_1^\infty dW\,
\rho(\beta,W)\,\Theta(P(\beta,W)).
$$

- Meaning: probability mass inside the paradox wedge.
- Status: v1.0 metadata and the v1.2.0 compiled paper definition.
- Source: archived PDF, p. 3, Eq. (8); Zenodo v1.0.0 description.

##### EA-POP-002B — Competing first-severity-moment definition of (Q)

The earlier transfer note instead wrote

$$
Q[\rho]=\langle P_+(\beta,W)\rangle.
$$

- Status: conflicting historical definition in `02`, lines 1307–1318. The current-turn v6.3 series for (Q) is consistent with a first positive severity moment, not wedge occupancy.
- Required author decision: preserve distinct symbols, for example (Q_{\rm occ}=E[\mathbf 1_{P>0}]) and (M_1=E[P_+]). Do not merge them silently.

##### EA-POP-003 — Mean-square paradox severity

$$
J[\rho]
=\int_0^1d\beta\int_1^\infty dW\,
\rho(\beta,W)\,[P_+(\beta,W)]^2.
$$

- Meaning: population mean of squared positive severity.
- Range: (0\le J\le1).
- Status: proposed; it is the preferred SR control variable from v1.2.0 onward.
- Source: archived PDF, p. 3, Eq. (9); `02`, lines 2013–2019.

##### EA-POP-004 — Conditional RMS severity

$$
C[\rho]=\sqrt{\frac{J[\rho]}{Q_{\rm occ}[\rho]}}
\qquad(Q_{\rm occ}>0).
$$

- Meaning: RMS of (P_+) conditional on being in the wedge, but only if the denominator is wedge occupancy.
- Revision warning: if (Q) is redefined as (E[P_+]), this conditional-RMS interpretation no longer follows.
- Source: archived PDF, p. 3, Eq. (10).

#### A.3.3.3 Near-luminal limits and the v6.3 correction

##### EA-ASY-001 — Small-FTL variables

$$
W=1+\delta,\qquad 0<\delta\ll1,
\qquad
\varepsilon_{\rm FTL}^2\equiv E[(W-1)^2]=E[\delta^2].
$$

- (\varepsilon_{\rm FTL}^2) is a population second moment, not a speed and not an observed parameter.
- Source: archived PDF, p. 3, Eqs. (11), (14); `02`, lines 2023–2027.

##### EA-ASY-002 — Exact uniform-β wedge fraction

For fixed (W=1+\delta) and uniform (\beta\in[0,1]), wedge occupancy is

$$
f_{\rm wedge}(W)=1-\beta_c(W)
=\frac{(W-1)^2}{1+W^2}
=\frac{\delta^2}{2+2\delta+\delta^2}.
$$

Hence

$$
f_{\rm wedge}
=\frac{\delta^2}{2}-\frac{\delta^3}{2}
+\frac{\delta^4}{4}+O(\delta^5).
$$

- The leading (\delta^2/2) coefficient is correct for **wedge occupancy**.
- Source connection: archived PDF p. 3 Eq. (12) gives the leading term; the exact expression follows directly from archived Eq. (5).

##### EA-ASY-003 — Historical v1.0–v1.2 small-FTL formulas

The archived paper states

$$
Q_{\rm occ}[\rho]\simeq \frac12\varepsilon_{\rm FTL}^2,
\qquad
J[\rho]\simeq \frac16\varepsilon_{\rm FTL}^2,
\qquad
C[\rho]\to\frac1{\sqrt3}.
$$

- Applicability: populations concentrated near (W=1), with the toy uniform-β weighting implicit in the quoted coefficients.
- Source: archived PDF, p. 3, Eqs. (12)–(15), and p. 9 Eq. (47); `02`, lines 2023–2040.
- Revision status: the (J) leading coefficient survives the reported v6.3 correction. The (Q) coefficient depends on which (Q) is meant.

##### EA-ASY-004 — Reported v6.3 corrected first- and second-severity moments

The current-turn v6.3 report gives the uniform-tail onset

$$
Q=\frac{\delta^2}{4}-\frac{\delta^3}{3}
+\frac{9\delta^4}{32}+O(\delta^5),
$$

$$
J=\frac{\delta^2}{6}-\frac{\delta^3}{4}
+\frac{\delta^4}{4}+O(\delta^5).
$$

It explicitly says the earlier $Q\sim\delta^2/2$ coefficient was wrong. This is compatible with $Q=E[P_+]$, but not with wedge occupancy EA-ASY-002. The most specific recoverable interpretation is fixed $W=1+\delta$, uniform $\beta\in[0,1]$, and $Q=E[P_+]$; the unavailable v6.3 source must confirm those conditions. For distributed $\delta$, the displayed powers would require the corresponding moments $E[\delta^2],E[\delta^3],E[\delta^4]$, not substitution of one generic $\delta$. Because the actual v6.3 document was unavailable, the definition and derivation must be verified from that kit before treating EA-ASY-004 as corpus-ready.

##### EA-ASY-005 — v1.2.0 numerical validity bands, now historical

For the toy family (\beta\sim U(0,1)), (\delta\sim U(0,\delta_{\max})), the archived paper reported the largest (\varepsilon_{\rm FTL}^2) retaining a given ratio to the historical leading approximation:

| Accuracy floor | Largest $\varepsilon_{\rm FTL}^2$ meeting the historical $Q_{\rm occ}$ accuracy floor | Largest $\varepsilon_{\rm FTL}^2$ meeting the $J$ accuracy floor |
|---:|---:|---:|
| 99% | (6.0\times10^{-5}) | (2.6\times10^{-5}) |
| 95% | (1.5\times10^{-3}) | (7.0\times10^{-4}) |
| 90% | (6.6\times10^{-3}) | (3.0\times10^{-3}) |
| 80% | (3.0\times10^{-2}) | (1.4\times10^{-2}) |

Source: archived PDF, p. 10, Table 1. The (Q) column is tied to the historical occupancy definition and must not be reused for a first-moment (Q) without recomputation.

#### A.3.3.4 Cosmology module

##### EA-COS-001 — Smooth broken power law

$$
\Omega_{\rm gw}(f)=\Omega_k
\left[
\frac{(f/f_k)^{\alpha_1\Delta}+(f/f_k)^{\alpha_2\Delta}}{2}
\right]^{-1/\Delta}.
$$

- Symbols: (f,f_k) in hertz; (\Omega_{\rm gw}), dimensionless spectral GW energy-density function under an unstated normalization convention; (\Omega_k), value at the knee; (\alpha_1,\alpha_2), slopes; (\Delta), smoothness.
- Archived demo parameters: (\Omega_k=10^{-12}), (f_k=3\times10^{-9}\,\mathrm{Hz}), (\alpha_1=2), (\alpha_2=-1), (\Delta=2), with “conservative” integration limits not specified in the paper text.
- Status: standard-looking SBPL family adapted as a proposed hidden-sector spectrum; no FTL production mechanism derives these parameters.
- Source: archived PDF pp. 3–4, Eq. (16); `02`, lines 2044–2049.

##### EA-COS-002 — Archived SGWB-to-ΔNeff mapping

$$
\Delta N_{\rm eff}
=\mathcal K\int \Omega_{\rm gw}(f)\,d\ln f,
\qquad
\mathcal K=\frac87\left(\frac{11}{4}\right)^{4/3}.
$$

- Dimensionless integral; finite frequency bounds and exterior-tail assumptions are required.
- Status: presented as a cosmology mapping, but the normalization is underspecified. The standard relation applies the prefactor to (\rho_{\rm gw}/\rho_\gamma). If (\Omega_{\rm gw}) is the conventional critical-density fraction, an (\Omega_\gamma)/(h^2) conversion is needed. The archived equation does not say which convention makes the displayed expression complete.
- v6.3 report: a machine PASS was later restricted to an exactly declared piecewise-linear-in-​(\ln f), zero-exterior model with zero uncomputed allowances. Externally asserted cosmology-error bounds remain non-promotable. The exact v6.3 formula was not accessible.
- Source: archived PDF p. 4 Eq. (17); `02`, lines 2050–2057; current-turn v6.3 commentary.

##### EA-COS-003 — Energy-budget fraction

$$
\varepsilon_{\rm energy}
=\frac{\Delta N_{\rm eff}}{\Delta N_{\rm eff}^{\max}},
\qquad \Delta N_{\rm eff}^{\max}=0.3\quad\text{(archived default)}.
$$

- Status: proposed normalization/guardrail, not a timeless universal observational limit.
- Source: archived PDF p. 4 Eq. (18); `02`, lines 2055–2057.

#### A.3.3.5 Observational compression and Einstein geometry

##### EA-OBS-001 — PTA consistency loss

$$
y\equiv1-\mathrm{PCI},\qquad 0\le\mathrm{PCI}\le1.
$$

- PCI is a proposed overlap statistic, not a posterior probability that the FTL hypothesis is true.
- Source: `02`, lines 2061–2069; v1.2.0 patch lines 45–57.

##### EA-OBS-002 — Ringdown standardized deviation

$$
z\equiv\frac{\varepsilon_{\rm pooled}}{\sigma_\varepsilon}.
$$

- (\varepsilon) is intended to be one harmonized, dimensionless, fractional ringdown deviation parameter with GR null value fixed by a manifest.
- (z) is dimensionless and unbounded.
- Source: `02`, lines 2061–2069; v1.2.1 claim-lock and provider contract.

##### EA-EIN-001 — Causality–cosmology action and historical small-δ normalization

$$
S\equiv\varepsilon_{\rm energy}J[\rho],
\qquad S_0\equiv\frac16,
\qquad
\frac{S}{S_0}\simeq\varepsilon_{\rm energy}\varepsilon_{\rm FTL}^2.
$$

- Status: proposed diagnostic construction. Multiplication of (J) and (\Delta N_{\rm eff}), and (S_0=1/6), are definitions chosen to recover the historical small-δ form; they are not derived physical dynamics.
- Source: archived PDF p. 5 Eqs. (20)–(21); `02`, lines 2092–2098.

##### EA-EIN-002 — Einstein vector

Historical form:

$$
X=\left(\frac S{S_0},\;1-\mathrm{PCI},\;
\frac{\varepsilon_{\rm pooled}}{\sigma_\varepsilon}\right)
=(x,y,z).
$$

J-first form:

$$
X=\left(20\Delta N_{\rm eff}J,\;1-\mathrm{PCI},\;
\frac{\varepsilon_{\rm pooled}}{\sigma_\varepsilon}\right).
$$

- Status: proposed compression. “Einstein” is a project label; (X=0) is not a complete mathematical statement of GR.
- Source: archived PDF p. 4 Eq. (19), p. 7 Eq. (29); v1.2.0 patch lines 64–81; `02`, lines 2073–2081 and 2105–2110.

##### EA-EIN-003 — Euclidean Einstein distance and unit ball

$$
D_E=\lVert X\rVert=\sqrt{x^2+y^2+z^2},
\qquad D_E\le1.
$$

- Status: proposed, conventional diagnostic. The threshold 1 has not been calibrated as a frequentist confidence level or Bayesian evidence threshold.
- All coordinates are dimensionless but not naturally commensurate: (y\in[0,1]), (z) is a z-score, and (x) has a chosen normalization.
- Source: archived PDF p. 5 Eqs. (22)–(23); `02`, lines 2073–2087.

##### EA-EIN-004 — Data floor, radicand, and slack

$$
D_{\min}=\sqrt{y^2+z^2},
\qquad
R\equiv1-y^2-z^2
=1-(1-\mathrm{PCI})^2-z^2,
$$

$$
x_{\max}=\sqrt R\quad\text{only for }R\ge0.
$$

- Invalid (R<0) draws must be reported separately and not clipped.
- Source: `02`, lines 2083–2087, 2112–2114; v1.2.1 public description; `FTL_A2_V6_PAPER_INSERT_PROTOCOL.md`, lines 3–17.
- Revision: v6.3 identified the implicit PTA/ringdown draw-coupling assumption in a single empirical (R) posterior and added finite-sample dependence certificates. The exact certificate was unavailable.

##### EA-EIN-005 — Historical small-FTL no-free-lunch inequality

Using (x\simeq\varepsilon_{\rm energy}\varepsilon_{\rm FTL}^2),

$$
(\varepsilon_{\rm energy}\varepsilon_{\rm FTL}^2)^2
\le1-(1-\mathrm{PCI})^2-
\left(\frac{\varepsilon_{\rm pooled}}{\sigma_\varepsilon}\right)^2,
$$

and, when the radicand is nonnegative and (\varepsilon_{\rm energy}>0),

$$
\varepsilon_{\rm FTL}^2\le
\frac{\sqrt{1-(1-\mathrm{PCI})^2-
(\varepsilon_{\rm pooled}/\sigma_\varepsilon)^2}}
{\varepsilon_{\rm energy}}.
$$

- Meaning: algebraic consequence of the proposed Euclidean unit ball and historical small-δ map.
- Limitations: not valid as an empirical FTL bound without real inputs; not valid outside the small-δ or explicit-​(\rho) regime; does not follow from GR alone.
- Source: archived PDF p. 5 Eqs. (24)–(26); `02`, lines 2092–2104.

##### EA-EIN-006 — J-first coupling

$$
x=6\varepsilon_{\rm energy}J
=20\Delta N_{\rm eff}J.
$$

The second equality uses (\Delta N_{\rm eff}^{\max}=0.3).

- Status: preferred v1.2.0 diagnostic map; avoids translating a large-​(W) severity directly into (\varepsilon_{\rm FTL}^2).
- Source: `02`, lines 2105–2110; `FTL_v1_4_NEXT_STEPS.md`, lines 11–31.

##### EA-EIN-007 — J-first SR envelope

$$
J\le J_{\rm bound}(\Delta N_{\rm eff})
=\frac{x_{\max}}{6\varepsilon_{\rm energy}}
=\frac{x_{\max}}{20\Delta N_{\rm eff}}.
$$

Derivation:

$$
D_E\le1\implies x^2\le R\implies
20\Delta N_{\rm eff}J\le\sqrt R=x_{\max},
$$

using (J\ge0), (\Delta N_{\rm eff}>0).

- Status: exact within the proposed audit definitions; it is not a model-independent law of nature.
- Source: `02`, lines 2112–2120; archived patch lines 85–93.

#### A.3.3.6 Testability Gate and beta models

##### EA-GATE-001 — Saturation ceiling

$$
\lim_{W\to\infty}P(\beta,W)=\beta^2,
\qquad
J_{\max}=E[\beta^4].
$$

- Assumptions: fixed beta-population while (W\to\infty); (J) is the mean square of (P_+); enough regularity to exchange the limit and expectation.
- Check: (\beta_c(W)\to0), so all (\beta>0) enter the wedge.
- Source: `02`, lines 2124–2129; archived PDF pp. 17, 21–24.

##### EA-GATE-002 — Bound-existence condition

$$
J_{\rm bound}<J_{\max}
\iff
\Delta N_{\rm eff}>
\Delta N_{\rm eff}^{\rm thr}
\equiv\frac{x_{\max}}{20J_{\max}}.
$$

- Meaning: under the constructed audit, a **finite-$W$** severity crossing can exist only above this strict threshold, assuming $J(W)$ approaches its ceiling without attaining it at finite $W$.
- The non-strict relation $J_{\rm bound}\le J_{\max}$ describes reachability only after admitting the $W\to\infty$ asymptote. Equality does not establish a finite upper bound on $W$.
- Important language: “no finite bound exists” means **this audit channel cannot exclude finite $W$** under the selected beta model. It does not mean FTL is allowed by nature or untestable by all experiments.
- Source: `02`, lines 2131–2138; archived patch lines 95–111.

##### EA-GATE-003 — Beta-family fourth moment

For (\beta\sim\mathrm{Beta}(a,b)), (a,b>0),

$$
E[\beta^4]=
\frac{a(a+1)(a+2)(a+3)}
{(a+b)(a+b+1)(a+b+2)(a+b+3)}.
$$

- Source: `02`, lines 2146–2149; archived PDF p. 24.
- Examples from the v1.2.0 A1 (x_{\max}\simeq0.709) table:

| Beta weighting | (J_{\max}) | (\Delta N_{\rm eff}^{\rm thr}) | Within archived 0.3 guardrail? |
|---|---:|---:|---|
| uniform | 0.200 | 0.177 | yes |
| proportional to (\beta^2) | 0.429 | 0.083 | yes |
| proportional to (\beta^5) | 0.600 | 0.059 | yes |
| proportional to ((1-\beta)^2) | 0.029 | 1.241 | no |
| proportional to (\beta(1-\beta)) | 0.143 | 0.248 | yes |
| proportional to (\beta^2(1-\beta)^2) | 0.119 | 0.298 | yes |

These are model illustrations, not inferred astrophysical beta populations.

##### EA-GATE-004 — Relativistic-tail approximation

$$
J_{\max}\approx f_{\rm rel}\beta_{\rm rel}^4,
\qquad
f_{\rm rel}\beta_{\rm rel}^4
\gtrsim\frac{x_{\max}}{20\Delta N_{\rm eff}},
$$

or

$$
f_{\rm rel}\gtrsim
\frac{x_{\max}}
{20\Delta N_{\rm eff}\beta_{\rm rel}^4}.
$$

- Assumption: a two-component mixture with fraction (f_{\rm rel}) concentrated near (\beta_{\rm rel}), the remainder negligible in fourth moment.
- Status: proposed operational approximation.
- Source: `02`, lines 2153–2159; archived patch lines 127–136.

#### A.3.3.7 Diagnostics, covariance, and information toy model

##### EA-DIAG-001 — Einstein budget fractions

$$
f_{\rm FTL}=\frac{x^2}{D_E^2},\qquad
f_{\rm PTA}=\frac{y^2}{D_E^2},\qquad
f_{\rm RD}=\frac{z^2}{D_E^2},
\qquad
f_{\rm FTL}+f_{\rm PTA}+f_{\rm RD}=1.
$$

- Valid for (D_E>0); descriptive decomposition only.
- Source: `02`, lines 2163–2171; archived PDF p. 8 Eq. (42) gives the standardized diagonal analogue.

##### EA-COV-001 — Covariance matrix

$$
C=
\begin{pmatrix}
\sigma_x^2&\rho_{xy}\sigma_x\sigma_y&\rho_{xz}\sigma_x\sigma_z\\
\rho_{xy}\sigma_x\sigma_y&\sigma_y^2&\rho_{yz}\sigma_y\sigma_z\\
\rho_{xz}\sigma_x\sigma_z&\rho_{yz}\sigma_y\sigma_z&\sigma_z^2
\end{pmatrix}.
$$

- Standard covariance parameterization, proposed for the audit coordinates.
- Source: archived PDF p. 7 Eq. (31).

##### EA-COV-002 — Covariance-weighted quadratic form

$$
\Delta X=X-E_0[X],\qquad E_0[X]=(0,0,0)^T,
$$

$$
\mathcal T_E=\Delta X^TC^{-1}\Delta X,
\qquad
D_E^{({\rm cov})}=\sqrt{\mathcal T_E}.
$$

- Archived code also reports (\mathcal T_\mu=\mu^TC^{-1}\mu) and (\mathcal T_i=X_i^TC^{-1}X_i).
- A numerical ridge was proposed: (C_{\rm reg}=C+10^{-12}\operatorname{tr}(C)I_3).
- Status: standard Mahalanobis structure, proposed audit interpretation.
- Source: archived PDF p. 7 Eqs. (32)–(33); `A2_tools_v1_4/ftl_cov_tension_from_samples.py`, lines 65–84.

##### EA-COV-003 — Approximate null calibration

$$
\mathcal T_E\overset{\rm approx}{\sim}\chi^2_3,
\qquad
p_E=P(\chi^2_3\ge\mathcal T_E),
\qquad
T_E=\Phi^{-1}\left(1-\frac{p_E}{2}\right).
$$

- Assumptions: approximately Gaussian coordinates, correctly specified null mean, known or well-estimated nonsingular covariance, and appropriate independence/effective sample size.
- Caution: (x\ge0), nonlinear transformations, boundary mass at invalid (R), and estimated covariance can violate a simple (\chi^2_3) calibration.
- Source: archived PDF p. 7 Eqs. (34)–(36).

##### EA-COV-004 — Diagonal weighted form and bound

$$
C\simeq\operatorname{diag}(\sigma_x^2,\sigma_y^2,\sigma_z^2),
\quad
\mathcal T_E=\left(\frac{x}{\sigma_x}\right)^2+
\left(\frac{y}{\sigma_y}\right)^2+
\left(\frac{z}{\sigma_z}\right)^2.
$$

For threshold (D_{E,\max}), historical small-δ (x=\varepsilon_{\rm energy}\varepsilon_{\rm FTL}^2) gives

$$
\varepsilon_{\rm FTL}^2\le
\frac{\sigma_x}{\varepsilon_{\rm energy}}
\sqrt{D_{E,\max}^2-
(y/\sigma_y)^2-(z/\sigma_z)^2}.
$$

- Source: archived PDF p. 8 Eqs. (37)–(41).

##### EA-INFO-001 — Toy capacity scaling

With a binary toy channel using (p_0=q), (p_1=2q), clipped below (1/2),

$$
C_{\rm FTL}(q)\sim q\log\frac1q
\qquad(q\to0).
$$

Historical translation:

$$
Q_{\max}\simeq\frac12\varepsilon_{\rm FTL,max}^2,
\qquad
C_{\rm FTL}^{\max}\sim Q_{\max}\log\frac1{Q_{\max}}.
$$

- Status: toy/model-dependent. Logarithm base is not fixed in the displayed derivation even though a demo is quoted in bits per channel use.
- Revision warning: the translation depends on the unresolved (Q) definition and historical coefficient.
- Source: archived PDF p. 6 Eqs. (27)–(28); `02`, lines 1396–1401.

##### EA-FLAVOR-001 — Optional toy flavor–FTL coupling, not core

$$
\varepsilon_\ell=\sqrt{m_\mu/m_\tau},
\qquad
\delta=\kappa\varepsilon_\ell^2,
\qquad
\varepsilon_{\rm FTL}^2\approx\kappa^2\varepsilon_\ell^4,
$$

$$
\kappa\le\frac1{\varepsilon_\ell^2}
\left[
\frac{1-(1-\mathrm{PCI})^2-
(\varepsilon_{\rm pooled}/\sigma_\varepsilon)^2}
{\varepsilon_{\rm energy}^2}
\right]^{1/4}.
$$

- Status: explicitly optional toy extension in the v1.1 handoff, not a derived part of the FTL hidden-sector theory. No mechanism connects lepton masses to superluminality.
- Classification recommendation: exclude from the primary FTL corpus or place in a clearly labeled historical/exploratory appendix.
- Source: `02`, lines 1403–1415.

#### A.3.3.8 A2 statistical estimator equations

##### EA-PCI-001 — K-way density overlap

For normalized PTA posterior densities (p_i(\theta)), (i=1,\dots,K),

$$
m(\theta)=\frac1K\sum_{i=1}^Kp_i(\theta),
\qquad
\mathrm{PCI}_{\rm all}=\int\min_i p_i(\theta)\,d\theta.
$$

- For (K=2), this overlap is (1-\mathrm{TV}(p_1,p_2)).
- The densities must refer to the same parameterization and compatible priors; otherwise overlap confounds prior/model differences with data consistency.
- Source: `FTL_A2_CROSSFIT_PRIMARY_E2E_THEOREM_NOTE_v5_9.md`, lines 3–13.

##### EA-PCI-002 — Equal-prior Bayes identity

$$
\eta_i(\theta)=P(i\mid\theta)
=\frac{p_i(\theta)}{\sum_jp_j(\theta)},
$$

$$
\mathrm{PCI}_{\rm all}
=K\,\mathbb E_{\theta\sim m}\left[\min_i\eta_i(\theta)\right].
$$

Derivation:

$$
K\,m(\theta)\min_i\eta_i(\theta)
=\left(\sum_jp_j\right)
\frac{\min_ip_i}{\sum_jp_j}
=\min_i p_i.
$$

- Exact only for equal class weights and true Bayes probabilities.
- Source: theorem note v5.9, lines 15–25.

##### EA-PCI-003 — Cross-fit classifier estimator

$$
\widehat{\mathrm{PCI}}^{\rm CF}_{\rm all}
=K\frac1n\sum_{r=1}^n
\min_i\widehat\eta_i^{(-s(r))}(\theta_r).
$$

- Fold (s(r)) is withheld from model fitting/calibration for evaluated row (r).
- v6.0 used standardized logistic regression, `CalibratedClassifierCV`, stratified row folds, and balanced bootstrap resampling.
- Later correction: balancing sampled row indices with replacement before row-level folding allowed duplicate source rows across train/test folds. v6.2 added group/provenance protections.
- Crucial limitation: cross-fit removes direct row train/test reuse only when grouping is correct; it does not ensure the classifier family approximates the Bayes posterior.
- Source: theorem note v5.9, lines 27–34; `ftl_a2_crossfit_pci_referee_v5_8.py`, lines 84–98 and 144–176.

##### EA-PCI-004 — Calibration diagnostics used in v6.0

Multiclass Brier score:

$$
B=\frac1n\sum_{r=1}^n\sum_{i=1}^K
(q_{ri}-Y_{ri})^2.
$$

Binned expected calibration error:

$$
\mathrm{ECE}=\sum_b\frac{n_b}{n}
\left|\operatorname{acc}(b)-\operatorname{conf}(b)\right|.
$$

- These diagnostics were recorded but v6.0 did not turn them into an accuracy certificate for PCI.
- Source: `ftl_a2_crossfit_pci_referee_v5_8.py`, lines 113–130.

##### EA-PCI-005 — v6.0 threshold-tax diagnostic

For reference ringdown (z_{\rm ref}),

$$
r(\pi)=2\pi-\pi^2-z_{\rm ref}^2
=1-(1-\pi)^2-z_{\rm ref}^2,
$$

$$
\mathrm{tax}(\pi_{\rm low},\pi_{\rm high})
=\sqrt{\frac{r(\pi_{\rm high})}{r(\pi_{\rm low})}}.
$$

- v6.0 labels: stable at 1% if tax ≤1.01; stable at 5% if ≤1.05; otherwise limited.
- Later result: stability is not accuracy. The v6.2 variance-shift truth case produced stable but badly biased classifier estimates.
- Source: `ftl_a2_crossfit_pci_referee_v5_8.py`, lines 179–207.

##### EA-PCI-006 — Reported v6.3 PCI accuracy inequality

$$
\left|O(q)-\mathrm{PCI}_{\rm all}\right|
\le\sqrt{K(K-1)\rho_B}.
$$

- Current-turn description: (O(q)) is the classifier-derived overlap functional; (K) is the number of PTA classes; (\rho_B) is the raw/excess Brier quantity used by the certificate.
- The exact definitions of (O(q)) and (\rho_B), proof, finite-sample adjustment, and estimator conditions are not present in the accessible sources. v6.3 reportedly permits promotion only when the raw-Brier bound is self-computed; externally asserted bounds are non-promotable.
- Status: reported mathematical advance, awaiting verification from the actual v6.3 note/code.

##### EA-RD-001 — Archived inverse-variance ringdown pooling

For event (i), with posterior-sample mean (\mu_i) and sample variance (s_i^2),

$$
w_i=s_i^{-2},\qquad
\mu_{\rm pool}=\frac{\sum_iw_i\mu_i}{\sum_iw_i},
\qquad
\sigma_{\rm pool}=\left(\sum_iw_i\right)^{-1/2},
\qquad
z=\frac{\mu_{\rm pool}}{\sigma_{\rm pool}}.
$$

- Assumptions: the harmonized event parameters share one null/scale/sign convention; events are sufficiently independent; posterior standard deviations can be used as inverse-variance weights; between-event heterogeneity is negligible.
- v6.0 bootstrap modes: `within_only` primary; `event_and_within` sensitivity.
- Later defect: v6.0 validated manifest fields `scale_to_common`, `shift_to_common`, and `null_value` but the pooling function read raw epsilon arrays and did not apply those transformations. v6.2 corrected this.
- Source: v1.2 A2 runner lines 136–196; v6.0 `ftl_a2_crossfit_primary_e2e_v5_9.py`, lines 124–155; provider referee lines 212–234.

##### EA-PHASE-001 — Common threshold form and phase partition

For any positive SR denominator (K_{\rm SR}),

$$
\Theta_{K_{\rm SR}}=\frac{\sqrt R}{20K_{\rm SR}},\qquad R\ge0.
$$

Let (K_p=J(W_{\rm cap})), (K_e=J_{\max}), and (\Delta=\Delta N_{\rm eff}). Then

$$
P_{\rm invalid}=P(R<0),
$$

$$
P_{\rm practical}(\Delta)
=P\big(0\le R\le(20\Delta K_p)^2\big),
$$

$$
P_{\rm latent}(\Delta)
=P\big((20\Delta K_p)^2<R\le(20\Delta K_e)^2\big),
$$

$$
P_{\rm closed}(\Delta)
=P\big(R>(20\Delta K_e)^2\big).
$$

These sum to one if (0<K_p\le K_e) and the stated inequalities cover all draws.

- Interpretation: “practical” means the boundary can be reached by (W\le W_{\rm cap}); “latent” means only above that cap but below saturation; “closed” means not reachable even at saturation; “invalid” means the observational compression already lies outside the unit ball.
- Status: exact partition of the proposed audit once the joint draw distribution and denominators are defined.
- Source: theorem note v5.9, lines 48–72; v6.0 paper insert lines 13–19.

#### A.3.3.9 Named v6.3 mathematics whose exact equations were inaccessible

The current-turn record reports, but does not expose complete notation/derivations for:

1. A cancellation-resistant fixed-domain/Appell-​(F_1) representation for complete Beta-law severity moments near (W=1).
2. Sharp finite-sample dependence certificates combining Makarov coupling bounds with Dvoretzky–Kiefer–Wolfowitz uncertainty.
3. An exact multiway-overlap/common-component bound.
4. An exact common-covariance Gaussian benchmark.
5. A resolution result that (\alpha=\delta=0.05) requires at least 3,506 genuinely independent units per marginal.

Do not invent these formulas from their names. The v6.3 mathematics note, source code, and tests are required before they can enter the equation corpus.

### A.3.4 Evidence, numerical results, tests, and reproducibility

#### A.3.4.1 What counts as evidence here

The record contains mathematical derivations, synthetic data, approximate anchors, code-path tests, adversarial software tests, and links to real public data products. It contains no completed strict joint analysis of real PTA posterior draws plus harmonized event-level ringdown posteriors. Therefore it contains **no empirical evidence for FTL**.

#### A.3.4.2 Historical demo and A1 values

##### MIN demo in the compiled paper

| Quantity | Archived demo value | Status |
|---|---:|---|
| (\Omega_k) | (10^{-12}) | deliberately tiny synthetic/profile input |
| (f_k) | (3\times10^{-9}\,\mathrm{Hz}) | demo input |
| (\alpha_1,\alpha_2,\Delta) | (2,-1,2) | demo inputs |
| (\Delta N_{\rm eff}) | (\simeq8.7\times10^{-12}) | computed under archived normalization |
| (\varepsilon_{\rm energy}) | (\simeq2.9\times10^{-11}) | derived demo value |
| PCI | (\simeq0.614) | synthetic PTA knee samples |
| ringdown (\varepsilon_{\rm pooled}\pm\sigma_\varepsilon) | (2.6\times10^{-3}\pm2.0\times10^{-2}) | synthetic |
| ringdown (z) | ≈0.13 | synthetic |
| (D_E) | ≈0.41 | demo |
| small-δ bound | (\varepsilon_{\rm FTL}^2\lesssim3\times10^{10}) | non-constraining and outside approximation |
| imposed demo cutoff | (\varepsilon_{\rm FTL}^2\lesssim0.1) | analyst restriction, not observation |
| toy (Q_{\max}), capacity | ≈0.05, ≈(6\times10^{-3}) bits/use | toy only |

Source: archived PDF pp. 4–6.

##### A1 `REALDATA_APPROX` anchors in v1.2.0

| Quantity | Value | Provenance status |
|---|---:|---|
| PCI | ≈0.828 | published-summary-style placeholder, not posterior ingest |
| (y) | ≈0.172 | derived |
| ringdown (\varepsilon_{\rm pooled}) | 0.13 | approximate placeholder |
| (\sigma_\varepsilon) | 0.19 | approximate placeholder |
| (z) | ≈0.684 | derived |
| (D_{\min}) | ≈0.705 | derived |
| (x_{\max}) | ≈0.709 | derived |
| (\Delta N_{\rm eff}\varepsilon_{\rm FTL}^2) | ≤0.213 | historical small-δ product bound |
| (\varepsilon_{\rm FTL,max}^2) at (\Delta N_{\rm eff}=0.30) | 0.709 | outside strict small-δ validity |
| same at 0.10 | 2.126 | invalid as small-δ translation |
| same at 0.05 | 4.252 | invalid as small-δ translation |

The archive repeatedly labels these as approximate anchors, not final evidence. Source: `02`, lines 2185–2193; archived PDF pp. 8–10 and 14–16.

#### A.3.4.3 v6.0/v6.1 public protocol results

The public v6.0 smoke suite executed three synthetic configurations:

- sigmoid calibration + within-event ringdown bootstrap;
- isotonic calibration + within-event bootstrap;
- sigmoid calibration + event-and-within bootstrap.

The v6.0 run index labels all three `CROSSFIT_STABLE_5PCT`, but the overarching promotion state is `SMOKE_TEST_ONLY_SYNTHETIC_ALLOWED`. The primary synthetic report had medians PCI 0.54449, (y=0.45551), (z=0.128561), (R=0.775223), (x_{\max}=0.880468), and zero invalid mass. Its model thresholds were also synthetic. Source: `FTL_A2_V6_PROTOCOL_RUN_INDEX.md`, lines 1–27; `.../primary_sigmoid_within/FTL_A2_CROSSFIT_PRIMARY_E2E_REPORT_v5_9.md`, lines 1–48.

The readiness scan reported `SCREENING_ONLY_UNLOCKED`, with no strict provider-return bundle. The public v1.2.1 claim lock expressly forbids claims of FTL evidence, real-data A2 constraints, strict PTA/ringdown ingestion, or Nobel-level discovery.

#### A.3.4.4 v6.2 corrective audit

The current-turn continuation reports that adversarial review found four material v6.0 defects:

1. A calibrated, cross-fit classifier could return a badly biased PCI for nonlinear posterior differences even when cross-fit looked stable.
2. Balanced sampling with replacement could duplicate source rows across folds, causing group leakage.
3. Explicit synthetic metadata could escape the earlier provenance checks.
4. Ringdown manifest transforms were validated but not applied to the data.

In a known variance-shift case, the exact/known overlap truth was reported as **0.515672**, while three leakage-free classifier-family PCI medians ranged from **0.540 to 0.976**. This demonstrates that “stable cross-fit” is not an accuracy certificate. A classifier-disagreement gate correctly stopped promotion.

The v6.2 mathematical layer was first reported to pass 11 independent unit tests and later to pass 16 baseline tests. These are internal software/mathematical tests, not observations.

#### A.3.4.5 v6.3 corrective audit and final reported state

Additional failures found and reportedly repaired:

- Direct Beta-distribution quadrature could collapse to zero extremely close to (W=1) because the wedge endpoint rounded to one in double precision, even when the true moment was nonzero and potentially enhanced by endpoint-singular Beta tails.
- v6.2's final gate could be satisfied by hand-written `PASS` JSON receipts.
- A hash-and-status graph could still accept a complete set of self-authored minimally shaped receipts.
- JSON `NaN`/`Infinity` values were not uniformly fail-closed.
- An interpolated quantile could claim 0.95 protected empirical mass when only 0.94 was guaranteed for a concrete 100-draw sample.
- Source-lock facts could be trusted from metadata receipts instead of recomputed from staged bytes.
- Replay could compare an arbitrary incomplete subset rather than the exact frozen inventory.
- A deterministic excess-Brier bound and cosmology error allowances were hash-bound but not recomputed.
- Boolean values could enter numeric `DeltaNeff`/`Jmax` fields because of language-level boolean/numeric subtype behavior.
- A legacy v6.2 smoke summary had rewritten itself during replay, invalidating its checksum; the baseline was reportedly restored from an untouched archive.

Reported v6.3 repairs:

- fixed-domain severity-moment evaluation near (W=1);
- DKW–Makarov dependence layer;
- hash-bound evidence DAG tied to exact upstream bytes;
- producer-schema and role-specific semantic validation;
- strict finite numeric parsing;
- exact inventory contract and independent extraction replay;
- recomputation of checksums, TOA overlap, Fourier-grid compatibility, schema extraction, and dual-fetch identity;
- self-computed raw-Brier certificate only;
- exactly specified piecewise-log-frequency cosmology model only;
- source disjointness explicitly provisional until hashed TOA identifiers prove zero overlap;
- separation of byte immutability from redistribution permission.

Final reported test counts rose through 54, 58, 72, 76, 77, and 81 as new adversarial regressions were added. Only the final 81/81 state should be cited, while preserving the earlier counts as development history.

#### A.3.4.6 Falsifiers and contradiction conditions

The framework can be contradicted at multiple levels:

1. **Kinematic derivation failure:** the archived symmetric antitelephone assumptions or algebra do not represent the claimed FTL protocol.
2. **Definition failure:** (P,Q,J), or the (Q) definition cannot be made internally consistent.
3. **Physical-bridge failure:** no defensible model links the hypothetical FTL population to the asserted SGWB/ΔNeff channel.
4. **Normalization failure:** the archived (\Omega_{\rm gw}\to\Delta N_{\rm eff}) convention is incomplete or incompatible with actual spectrum products.
5. **Compression failure:** (P(R<0)) is appreciable; such mass must remain invalid rather than clipped.
6. **Estimator failure:** classifier families disagree beyond the accuracy certificate, known benchmarks fail, or raw Brier bounds are too loose.
7. **Provenance failure:** inputs are synthetic, summary-only, overlapping, transformed incorrectly, unlicensed, mutable, or not hash-bound.
8. **Dependence failure:** conclusions change materially across all joint couplings compatible with the PTA and ringdown marginals.
9. **Testability-gate failure:** (J_{\max}<J_{\rm bound}) within the adopted cosmology domain; this closes this audit channel for that beta model.
10. **Empirical null result:** a properly preregistered, real-data A2 analysis yields only noninformative bounds or no model survives independent reproduction.

The most decisive immediate falsifier is not “observe an FTL signal”; it is to run the fully provenance-closed v6.3 protocol on genuinely independent, compatible posterior products and test whether any claimed bound survives estimator family, dependence, transformation, and source-lock uncertainty.

#### A.3.4.7 Reproducibility record

No research scripts were executed for this export. The following archived methods were inspected only.

Dependencies visible in public code:

- Python 3;
- NumPy;
- pandas;
- Matplotlib;
- SciPy (`quad`) with a NumPy Gauss–Legendre fallback for (J(W));
- scikit-learn: `StandardScaler`, logistic regression, `CalibratedClassifierCV`, `StratifiedKFold`;
- Python standard libraries for JSON, CSV, ZIP, hashing, subprocesses, and paths.

No pinned `requirements.txt`, lockfile, container digest, or exact dependency-version ledger was visible in the public v1.2.0/v1.2.1 top-level file inventories. v6.0 added a compatibility branch for old/new `CalibratedClassifierCV` constructor keywords, but compatibility is not full environment reproducibility.

Archived parameter defaults/examples:

- v1.2 A2 runner: bootstrap 500, histogram bins 160, seed 0, $J_{\max}$ values 0.2, 0.143, 0.119, and a 250-point $\Delta N_{\rm eff}$ grid from $10^{-3}$ to 1.
- v5.9 recommended strict run: observational bootstrap 5000, PCI bootstrap 200, 1500 samples/class, 5 folds, sigmoid calibration, `within_only` ringdown bootstrap, seed as supplied, (W_{\rm cap}=10).
- v6.0 frozen example: observational bootstrap 2000, PCI bootstrap 200, 3000/class, 5 folds, seed 60, (W_{\rm cap}=10), plus isotonic and event-level sensitivity runs.

These commands are historical reproduction instructions, not authorization to execute in this export.

### A.3.5 Sources and provenance

#### A.3.5.1 Local source files

1. `project_sources/02-dblast-black-holes-toe-dmde-equations.md`
   - Lines 1250–1461: v1.1-era FTL new-chat handoff, definitions, equations, diagnostics, numerical-output inventory, and real-data next step.
   - Lines 1960–2246: v1.2.0 canonical handoff with exact wedge, (P,Q,J), J-first envelope, Testability Gate, beta moments, A1 anchors, public files, and A2 mission.
   - Lines 1467–1480 explicitly list separate tracks and warn not to merge them.

2. `project_sources/09-theory-faster-than-light-instructions.txt`
   - Lines 1–2 and 3510–3517: author directive to remain FTL-only and treat other programs as reference-only.
   - Lines 3030–3056: planning notes for formalizing the hidden-sector audit.
   - Lines 3098–3107: explicit synthetic-capacity framing and reported toy limit.
   - Lines 3157–3175: synthetic MIN interpretation and insistence on avoiding physical-discovery claims.
   - Lines 3177–3235: tension-vector, phase, and Monte Carlo planning notes. These are planning/provenance statements, not a final scientific source.

#### A.3.5.2 Public author records and exact top-level files accessed remotely

The following official-record inventories and named archive members were inspected read-only on 2026-09-08. They are external author records, not task attachments. Because their container bytes were not retained with this master, “inspected” here is provenance for this audit pass, not a claim that the research corpus embeds or permanently preserves the packages.

##### v1.0.0 — DOI `10.5281/zenodo.17726160`

- `NEXT_LEVEL_LATEX_MC_20251120_041829.zip`, 159,898 bytes, MD5 `418d9d69b00ee4f98e88bc1e8714385b`.
- `ONE_ZIP_MIN_20251115_013859.zip`, 438,930 bytes, MD5 `2c754fd456800ff74026bda8918d80bb`.

##### v1.1.0 — DOI `10.5281/zenodo.17926499`

- `FTL_Einstein_Audit_SHARE_PACKAGE_v9_ZenodoReady.zip`, MD5 `582f8f65b13e6127c75c045135b70876`.
- `FTL_PUBLIC_EXPLAINER_PACK_v4.zip`, MD5 `93a3211e901b02118d86cef203414e0a`.
- `FTL_Einstein_Audit_PUBLIC_ARTICLE.pdf`, MD5 `b8b42e5e1c3804e40e154b9f899b5272`.
- `FTL_Einstein_Audit_FTLonly_v2.pdf`, MD5 `db2f5418d1b587243b0df7d69956382c`.

##### v1.2.0 — DOI `10.5281/zenodo.18499411`

- `FTL_v1_3_Compiled_Technical_and_Supplements.pdf`, 1,977,645 bytes, MD5 `6730b898fa0b85e52da64369dac8a008`.
- `FTL_v1_3_OnePager_Breakthrough.pdf`, 5,853 bytes, MD5 `dcb1dca54a41b39648a9fde3a098c8fe`.
- `FTL_Zenodo_v1_2_0_release_candidate.zip`, 2,423,492 bytes, MD5 `47b9cc7e1c4718aba80551aef9ed894a`.

Important embedded files inspected:

- `FTL_v1_4_Breakthrough_Summary.md`.
- `FTL_v1_3_PaperPatch_JFirst_TestabilityGate.md`.
- `A2_tools_v1_4/FTL_v1_4_NEXT_STEPS.md`.
- `A2_tools_v1_4/FTL_A2_REALDATA_SOURCES_AND_EXTRACTION.md`.
- `A2_tools_v1_4/ftl_a2_run_end_to_end.py`.
- `A2_tools_v1_4/ftl_cov_tension_from_samples.py`.

The 25-page compiled PDF includes the base technical paper (original equations 1–48), A1 results/figures, the exact-SR model-dependence addendum, the phase diagram, and the beta-fourth-moment gate note.

##### v1.2.1 — DOI `10.5281/zenodo.20218030`

- `FTL_A2_ZENODO_PROTOCOL_RELEASE_DECISION_KIT_v6_1.zip`, 1,463,506 bytes, MD5 `45b332403cfd5b9688e5bccec27f5a24`.
- `FTL_A2_STRICT_PROTOCOL_FREEZE_KIT_v6_0.zip`, 1,502,150 bytes, MD5 `d9fb0abc2ae08c06bf2091171273191f`.

Important embedded files inspected:

- `FTL_A2_STRICT_PROTOCOL_FREEZE_SUMMARY_v6_0.md`.
- `DOCS/FTL_A2_STRICT_PROTOCOL_FREEZE_RUNBOOK_v6_0.md`.
- `DOCS/FTL_A2_V6_PAPER_INSERT_PROTOCOL.md`.
- `CONFIGS/FTL_A2_V6_STRICT_PROVIDER_BUNDLE_CONTRACT.yaml`.
- `VENDOR/.../FTL_A2_CROSSFIT_PRIMARY_E2E_THEOREM_NOTE_v5_9.md`.
- `VENDOR/.../ftl_a2_crossfit_pci_referee_v5_8.py`.
- `VENDOR/.../ftl_a2_crossfit_primary_e2e_v5_9.py`.
- v6.0 readiness and synthetic smoke outputs.
- v6.1 claim-lock, version-strategy, release-decision, and Zenodo-description files.

The precise release relationship is documented in `DOCS/FTL_A2_ZENODO_DECISION_MEMO_v6_1.md`: lines 7–15 recommend publishing the work as the narrow methods/tooling release **v1.2.1** and reserving v1.3 for the first strict posterior-based A2 result; lines 34–53 record `SCREENING_ONLY_UNLOCKED` and enumerate the missing real-data files; lines 55–87 state the allowed/disallowed claims and version reservation. `ZENODO_RELEASE_FILES/README_v1_2_1_PROTOCOL_FREEZE.md`, lines 3–23, likewise states that v1.2.1 is a methods/tooling update without FTL evidence or strict A2 posterior results, and lines 25–40 identify the locked inputs and v6.0 ZIP as the protocol payload. `DOCS/FTL_A2_CLAIM_LOCK_v1_2_1.md`, lines 3–14 and 25–46, separates green-light protocol claims from prohibited empirical claims. Thus record 20218030 is not a rebranding of v6.0 as a physics result: it is a v6.1 release-decision wrapper plus the frozen v6.0 protocol payload, archived under public semantic version v1.2.1.

#### A.3.5.3 External data products named in the record

These are candidate inputs, not incorporated FTL evidence:

| Product | DOI / record | Intended role | Status in FTL analysis |
|---|---|---|---|
| EPTA DR2 III search chains | `10.5281/zenodo.8091568` | PTA posterior (\log_{10}A_{\rm gw},\gamma_{\rm gw}) | Public `chains.zip`; not completed in strict A2. |
| EPTA DR2 II posterior distributions | `10.5281/zenodo.8025019` | lighter derived PTA representation | Candidate only. |
| IPTA DR2 GWB MCMC output | `10.5281/zenodo.5787557` | PTA posterior | Candidate; archived note says discard ~25% burn-in. |
| NANOGrav 15-year free-spectrum KDE | `10.5281/zenodo.10344086` | possible pseudo-draws / sensitivity | Derived KDE is not automatically equivalent to raw joint posterior. |
| NANOGrav 15-Year Data Set | `10.5281/zenodo.16051178` | current v6.3 source-lock target | Timing-data release; full compatible posterior/prior/license/schema remains unresolved. |
| GWTC-3 Tests of GR | `10.5281/zenodo.17461225` | older event-level ringdown products | Candidate only; exact common epsilon not fixed. |
| GWTC-5.0 Tests of GR | `10.5281/zenodo.21454847` | current pSEOB/ringdown target | Public record includes `pSEOB.tar.gz`; not ingested in strict A2. |

The current API lists GWTC-5.0 publication date 2026-07-22, while the v6.3 conversation says a pSEOB source date was corrected to official 2026-07-21 metadata. This discrepancy must be resolved at file/metadata-field level rather than silently choosing one date.

#### A.3.5.4 Missing underlying sources

- The private v6.2 and v6.3 kits, source tree, validation report, exact Appell-​(F_1) formula, dependence-certificate derivation, test fixtures, and checksum ledger were unavailable.
- The public v1.2.1 archive contains the older v6.0 implementation that v6.2/v6.3 reportedly supersede; it must not be presented as the current correct implementation.
- No strict real-data provider-return bundle was available.
- No frozen PTA source bytes, prior files, TOA-identifier disjointness ledger, Fourier-grid compatibility receipt, or final redistribution-license decisions were available.
- No harmonized GWTC-5 pSEOB event mini-CSVs or transformation manifest were available.
- No external referee report, journal decision, independent lab reproduction, or physical experiment is documented.

An audit report, inventory, or checksum is not the underlying science. The actual posterior chains, priors, transformation definitions, source bytes, and implementation are required for scientific review.

### A.3.6 Chronological development history

#### 2025-11-26 — v1.0.0 MIN

- Introduced the antitelephone population, (P,Q,J,C), (\varepsilon_{\rm FTL}^2), SBPL/ΔNeff module, PTA PCI, ringdown (z), Einstein vector/distance, no-free-lunch inequality, and toy capacity.
- Explicitly analysis-only and pre-theory; no real PTA/ringdown data included.
- Historical (Q) was described as wedge occupancy with leading (\tfrac12\varepsilon_{\rm FTL}^2).

#### 2025-12-13 — v1.1.0

- Added non-tiny demo ΔNeff cases 0.1 and 0.3, two-dimensional energy scans, budget fractions, covariance-weighted tension appendix, optional flavor–FTL toy coupling, and public explainer.
- Continued to use the small-FTL (x\simeq\varepsilon_{\rm energy}\varepsilon_{\rm FTL}^2) translation.
- The transfer note at `02`, lines 1311–1318 wrote (Q=\langle P\rangle), conflicting with the archived occupancy definition.

#### 2026-01-29 to 2026-02-05 — J-first upgrade and v1.2.0

- Exact SR calculations showed the small-δ translation could be used far outside its validity when ΔNeff was small.
- Retained (J) as primary and shifted to (x=20\Delta N_{\rm eff}J).
- Derived the saturation ceiling (J_{\max}=E[\beta^4]), the J envelope, and the Testability Gate.
- Added beta-family and relativistic-tail illustrations.
- Added posterior-level A2 scaffolding and covariance tooling.
- Preserved A1 `REALDATA_APPROX` values only as placeholders.

#### February–May 2026 — A2 pipeline versions v5.x

- Histogram-overlap PCI was progressively demoted to a resolution sensitivity.
- Established the K-way overlap identity and cross-fit classifier estimator.
- Added provider-return validation, ringdown harmonization manifest, radicand (R), invalid-mass reporting, finite-​(W) thresholds, and phase partition.

#### 2026-05-06 to 2026-05-15 — v6.0/v6.1 and public v1.2.1

- v6.0 froze the sequence: provider referee → cross-fit calibrated PCI → (R) posterior → (x_{\max}) → SR thresholds/phases.
- Added old/new scikit-learn calibration-constructor compatibility.
- Synthetic smoke runs exercised primary and sensitivity configurations.
- v6.1 made the publication decision: publish only as v1.2.1 methods/tooling; reserve v1.3.0 for genuine posterior-level real data.
- Public claim lock forbade FTL evidence and real-data A2 claims.

#### Private v6.2 continuation — date not stated in accessible record

- Red-team audit invalidated “stable cross-fit” as an accuracy certificate.
- Found duplicate-row fold leakage, synthetic-marker gaps, and unapplied ringdown transforms.
- Added multi-classifier disagreement and fail-closed integrity tests.
- Result: methods correction, no empirical FTL discovery.

#### Private v6.3 continuation — current-turn, session date 2026-09-08

- Corrected the uniform-tail first-moment coefficient and added stable near-​(W=1) moment evaluation.
- Added finite-sample dependence treatment and a PCI accuracy inequality.
- Hardened evidence provenance, receipt semantics, strict numeric parsing, source locks, exact inventory replay, quantile coverage, Brier/cosmology recomputation, and deterministic builds.
- Final private artifact reportedly passed 81 tests and 534 receipt mutations, but the real-data candidate failed 123 checks.
- Explicit final decision remained no-go for empirical v1.3.0.

### A.3.7 Claim register

| Claim ID | Claim | Status | Assumptions / support | Main limitation |
|---|---|---|---|---|
| FTL-CL-001 | The stated symmetric two-leg FTL protocol is paradox-capable for (\beta>2W/(1+W^2)). | Conditional mathematical result | EA-SR-003/004; standard SR algebra | Specific signaling protocol; no physical FTL mechanism. |
| FTL-CL-002 | (P=1-t_2/T) measures paradox severity. | Proposed definition | Dimensionless and sign-matches wedge | Severity choice is conventional, not unique. |
| FTL-CL-003 | (J=E[P_+^2]) is a useful population severity. | Proposed, retained | Bounded, nonnegative, exact in chosen setup | Depends on (\rho(\beta,W)) and weighting. |
| FTL-CL-004 | Historical (Q\simeq\tfrac12\varepsilon_{\rm FTL}^2). | Definition-dependent / partly superseded | Correct for uniform-β wedge occupancy | v6.3 (Q) series uses ¼, implying a different (Q). |
| FTL-CL-005 | (J\simeq\tfrac16\varepsilon_{\rm FTL}^2). | Conditional asymptotic | Uniform-tail/near-​(W=1) assumptions; v6.3 preserves leading term | Higher-order and population dependence matter. |
| FTL-CL-006 | The displayed SBPL integral yields ΔNeff. | Qualified/underspecified | Archived code and formula | Spectrum normalization, frequency bounds, tails, and photon-density conversion unclear. |
| FTL-CL-007 | (D_E\le1) defines “Einstein-allowed.” | Proposed convention | Definition of unit ball | Not a theorem of GR or calibrated hypothesis test. |
| FTL-CL-008 | No-free-lunch inequality bounds (\varepsilon_{\rm FTL}^2). | Conditional algebra | Unit ball plus small-δ map | Most A1 numeric bounds fall outside small-δ validity. |
| FTL-CL-009 | J-first envelope bounds (J). | Conditional algebra, current core | EA-EIN-006/007 | Physical (x\propto\Delta N_{\rm eff}J) bridge remains assumed. |
| FTL-CL-010 | (J\to E[\beta^4]) at (W\to\infty). | Conditional mathematical limit | Exact (P), fixed beta distribution | Does not constrain a physical population absent a beta model. |
| FTL-CL-011 | Testability Gate separates active/inactive audit regions. | Conditional result | Compare reachable (J_{\max}) to (J_{\rm bound}) | Only this constructed audit channel; not universal testability. |
| FTL-CL-012 | (\mathrm{PCI}_{\rm all}=K E_m[\min\eta_i]). | Exact probability identity | Normalized densities, equal mixture, true Bayes probabilities | Estimated classifiers may be inaccurate. |
| FTL-CL-013 | Cross-fitting makes PCI reliable. | Superseded as stated | Removes direct train-on-test leakage if grouped correctly | v6.2 truth case showed severe stable bias. |
| FTL-CL-014 | v6.3 Brier inequality certifies PCI accuracy. | Reported, not file-verified | Current-turn formula and tests | Exact definitions/proof/source unavailable. |
| FTL-CL-015 | Inverse-variance ringdown pooling yields (z). | Conditional method | Common epsilon, transformations, independence, homogeneity | v6.0 failed to apply transforms; physical common parameter unresolved. |
| FTL-CL-016 | Phase masses sum to one. | Exact conditional partition | Valid joint draws and (K_p\le K_e) | Joint PTA/ringdown dependence and finite-sample uncertainty matter. |
| FTL-CL-017 | v6.3 is reproducibly hardened. | Internally reported | 81 tests, deterministic build, mutation sweep | Private artifact unavailable; not external replication. |
| FTL-CL-018 | Observations support FTL. | Unsupported / prohibited | None | No strict real-data A2 analysis and no direct FTL observation. |

#### A.3.8.3 Terminology and symbol glossary

| Term/symbol | Meaning | Status/caution |
|---|---|---|
| FTL | faster than light | Hypothetical in this program. |
| MIN | minimal framework | Historical release label. |
| A1 | approximate-summary-anchor stage | Not posterior-level real data. |
| A2 | strict posterior-level PTA + ringdown stage | Still locked. |
| (W=w/c) | signal-speed factor | (W>1) hypothetical. |
| (\beta=v/c) | observer relative speed | Archived support (0\le\beta<1). |
| (\delta=W-1) | near-luminal excess | Small-δ expansion variable. |
| (P) | proposed paradox severity | (P_+=\max(P,0)) in moments. |
| (Q) | historically either wedge occupancy or first positive moment | Must be disambiguated before ingestion. |
| (J) | (E[P_+^2]), mean-square severity | Preferred SR variable. |
| (C) | conditional RMS severity | Only with occupancy denominator. |
| (\rho(\beta,W)) | population density | Model assumption, not inferred. |
| (\varepsilon_{\rm FTL}^2) | (E[(W-1)^2]) | Population moment, not observation. |
| (\Omega_{\rm gw}(f)) | SGWB spectrum | Normalization convention must be explicit. |
| (\Delta N_{\rm eff}) | extra-radiation parameter | Standard observable; audit bridge proposed. |
| PCI | Probability of Consistency Index | K-way posterior overlap, not FTL probability. |
| (\varepsilon_{\rm pooled}) | common ringdown deviation summary | Exact parameter must be harmonized. |
| (z) | (\varepsilon_{\rm pooled}/\sigma_\varepsilon) | Proposed audit coordinate. |
| (x) | SR–cosmology coordinate | J-first (20\Delta N_{\rm eff}J). |
| (R) | (1-(1-\mathrm{PCI})^2-z^2) | Canonical observational radicand. |
| (D_E) | Euclidean norm of (X) | Author-defined “Einstein distance.” |
| (J_{\max}) | (E[\beta^4]) saturation ceiling | Conditional on exact setup. |
| (W_{\rm cap}) | finite practical speed cap | Archived default 10; analyst choice. |
| DKW–Makarov | finite-sample/marginal dependence certificate | Reported v6.3; exact implementation unavailable. |

### A.3.9 Outstanding decisions and research tasks

#### Mathematical/definition decisions

1. Resolve (Q): wedge occupancy versus first positive severity moment. Preserve both if both are useful.
2. Provide the exact fixed-domain/Appell-​(F_1) moment formula, proof, parameter domain, branch conventions, and limiting tests.
3. Provide the exact PCI Brier-bound definitions and proof.
4. Provide the exact DKW–Makarov certificate, confidence allocation, effective-sample assumptions, and derivation of 3,506 units.
5. Explain whether (D_E\le1) is only a visualization convention or can be calibrated as a statistical decision rule.
6. State whether covariance-weighted (\mathcal T_E) replaces or merely supplements the Euclidean audit, and address boundary/non-Gaussian null calibration.

#### Physical-model tasks

7. Specify a causal/dynamical hidden-sector model: fields, action/Lagrangian, symmetries, preferred frame if any, stability, unitarity, energy positivity, and interactions with the visible sector.
8. Derive, rather than assume, any relation between (J), hidden-sector energy density, ΔNeff, and an SGWB.
9. Specify why PTA and ringdown should jointly constrain one hidden sector and what alternative astrophysical explanations predict.
10. Resolve the ΔNeff normalization convention, integration domain, transfer functions, and exterior spectral tails.
11. Define the beta population physically and show how it can be inferred or bounded.

#### Real-data and reproducibility tasks

12. Freeze exact EPTA/IPTA/NANOGrav posterior products under a common model and compatible priors.
13. Prove independent sampling units and PTA disjointness with hashed TOA identifiers.
14. Obtain the unresolved compatible NANOGrav full posterior, prior, schema, and license rather than substituting timing data or KDE without qualification.
15. Select the GWTC-5 pSEOB parameter and events; apply signed/scale/null transforms; validate the harmonization physically, not only syntactically.
16. Bind source bytes, extraction code, environment, schemas, manifests, decisions, and receipts in a non-circular replay graph.
17. Run the fully frozen protocol on real data; keep v1.3.0 unpublished unless every gate passes.
18. Obtain external human statistical, SR, GW, and cosmology review and independent reproduction.


### A.3.10 Completeness and exclusions

Included in the underlying handoff review (referenced/extracted; underlying artifacts are not embedded here):

- both local FTL-bearing source files with exact line ranges;
- public Zenodo metadata for all four versions;
- the public v1.2.0 compiled paper, notes, A2 source code, and file inventory;
- the public v1.2.1 v6.0/v6.1 claim-lock, protocol, code, and synthetic-output inventory;
- all substantive v6.2/v6.3 corrections and reported final metrics visible in the current turn.

Not included or not recoverable:

- the actual private v6.2/v6.3 kit and validation report;
- exact formulas/derivations for the Appell-​(F_1), DKW–Makarov, multiway overlap, and common-Gaussian results;
- strict real-data posterior inputs and outputs;
- a full bibliography from the compiled paper (the accessible PDF is largely self-contained and does not provide a conventional reference list in the extracted pages);
- essential figures as images. Their captions and mathematical meanings are summarized, but the plots/phase heatmaps should be preserved separately from the public PDFs if visual fidelity is required.

Research assets that Markdown cannot faithfully replace:

- the exact SR-model envelope plots, beta-family phase heatmap, A1 contour plots, and code-generated calibration plots;
- posterior chains/KDE grids and GWTC HDF5/tar products;
- deterministic ZIP manifests, role-specific JSON schemas, adversarial fixtures, and checksum evidence graph from the unavailable v6.3 kit.

