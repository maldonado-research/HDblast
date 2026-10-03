# Einstein Audit for Faster-Than-Light Hidden Sectors (FTL-only)

**Ricardo Maldonado · working hypothesis and methods research**  
ORCID: [0009-0009-3937-6527](https://orcid.org/0009-0009-3937-6527)  
GitHub documentation edition: **2026-09-30** · research snapshot: **2026-09-08**  
First GitHub publication: **2026-10-01**

Could a hypothetical faster-than-light sector be subjected to clear mathematical and observational consistency checks? The Einstein Audit asks that question by examining causal loops, an assumed extra-radiation budget, and selected gravitational-wave consistency measures.

**This research does not demonstrate faster-than-light travel or communication. No completed strict real-data FTL constraint, physical discovery, or external peer review is documented in the available sources.**

This folder is a **separate research program** hosted in the author's existing HDblast repository. HDBLAST's higher-dimensional origin hypothesis is not part of this FTL audit and is not evidence for it. Black-hole ringdown appears here as a proposed observational input; the separate horizon-reflection hypothesis is not incorporated.

## Read the research

| Document | Purpose |
|---|---|
| [Detailed account and equations](docs/RESEARCH_ACCOUNT_AND_EQUATIONS.md) | Definitions, assumptions, derivations, historical results, limitations, sources, claim register, glossary, and research tasks |
| [Core mathematics guide](docs/CORE_MATHEMATICS.md) | A short route through the SR geometry, population severity, J-first envelope, and posterior-overlap identity |
| [Archived equation register](docs/ARCHIVED_V1_2_1_EQUATION_REGISTER.md) | Historical EA-001–EA-097 formulas, preserved separately with correction warnings |
| [Evidence and open problems](docs/EVIDENCE_AND_OPEN_PROBLEMS.md) | What is established, proposed, reported, missing, and potentially falsifiable |
| [Sources and reproducibility](docs/SOURCES_AND_REPRODUCIBILITY.md) | Archival references, source snapshot, missing code/data/proofs, and reproduction requirements |
| [100 FTL Explore questions](docs/100_FTL_EXPLORE_QUESTIONS.md) | Educational questions specifically about this working hypothesis |
| [Publication history](CHANGELOG.md) | Public archival versions versus private toolchain reports |
| [Citation metadata](CITATION.cff) | Attribution for this GitHub documentation edition |

## The central calculation

In a specific two-observer thought experiment, each observer can send a hypothetical signal at speed $u=Wc>c$ in their own frame. The observers move apart with relative speed $v=\beta c<c$. Alice sends at time $T>0$, and Bob replies immediately.

For collinear motion in flat spacetime, the reply time in Alice's frame is

$$
\frac{t_2}{T}=\frac{W^2(1-\beta^2)}{(W-\beta)^2}.
$$

The modeled reply precedes the sending event when $\beta>2W/(1+W^2)$. That is a conditional kinematic result, not evidence that nature permits the assumed communication rule.

The audit defines a positive paradox strength $P_+=\max(1-t_2/T,0)$ and averages its square over a stated population:

$$
J=\mathbb E[P_+^2].
$$

Under a fixed observer-speed population, $J$ approaches $J_{\max}=\mathbb E[\beta^4]$ as $W$ grows. This ceiling determines whether this particular diagnostic can impose a finite speed bound.

## The proposed observational bridge

The author-defined coordinates are

$$
x=20\Delta N_{\rm eff}J,\qquad y=1-\mathrm{PCI},\qquad z=\frac{\varepsilon_{\rm pooled}}{\sigma_\varepsilon}.
$$

The basic audit uses $D_E=\sqrt{x^2+y^2+z^2}\le1$. This unit-ball rule and the coupling to $\Delta N_{\rm eff}$ are **proposed diagnostics**, not consequences of Einstein's field equations or a calibrated discovery threshold. PCI is overlap between compatible posterior distributions; it is not a probability that FTL is true.

For positive $\Delta N_{\rm eff}$ and $R=1-y^2-z^2\ge0$, the construction gives

$$
J\le\frac{\sqrt R}{20\Delta N_{\rm eff}}.
$$

If this allowance is at or above $J_{\max}$, this channel provides no finite speed exclusion. Non-exclusion does not support FTL. Draws with $R<0$ must be reported separately rather than clipped into an admissible result.

## What is still missing?

The record does not supply a physical FTL-sector action, a production mechanism, a derived connection from its degrees of freedom to PTA/ringdown data, or a completed authenticated posterior analysis. It also contains important corrections: overloaded Q notation, an underspecified radiation conversion, classifier accuracy failures, and reported problems with provenance and ringdown transformations.

Public v1.2.1 is the latest documented archival release in the source snapshot. A later private toolchain v6.3 is described in the conversation, but its actual files are unavailable for this publication. Its test counts and new certificate formulas are therefore reported history, not independently verified assets. This documentation edition is **not public empirical v1.3.0**.

## Archival publications

The source record identifies:

- [v1.2.1 — DOI 10.5281/zenodo.20218030](https://doi.org/10.5281/zenodo.20218030): A2 Strict Protocol Freeze and Provider-Return Gate; methods/tooling only.
- [v1.2.0 — DOI 10.5281/zenodo.18499411](https://doi.org/10.5281/zenodo.18499411): J-First SR Envelope, Testability Gate, and A2 Real-Data Scaffolding.
- [All-version concept DOI 10.5281/zenodo.17726159](https://doi.org/10.5281/zenodo.17726159).

These identifiers are preserved from the supplied sources. Live Zenodo retrieval was unavailable during this edition; no claim is made that a complete current release list was verified on September 30. A Zenodo deposit does not constitute peer review.

## Review and attribution

Mathematical, statistical, cosmological, and gravitational-wave critiques are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md). Cite the version actually used and identify whether a result is conditional mathematics, a synthetic example, a reported software test, or a physical observation.

Documentation was prepared with OpenAI ChatGPT/Codex assistance under the author's direction. The review performed for this edition is internal and does not count as independent scientific replication. No research script or new experiment was run for this publication. No raw private transcripts, unrelated hypothesis archives, book manuscript, third-party posterior data, or mixed article ZIP is included.

The repository's existing [licensing policy](../../LICENSE) applies to author-owned documentation. Referenced external datasets and publications retain their own terms; their contents are not redistributed here.
