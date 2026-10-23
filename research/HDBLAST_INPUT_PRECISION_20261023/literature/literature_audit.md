# Bounded methods audit — 23 October 2026

The closest useful recent sources support a familiar proof pattern: keep a structured approximate solution, enclose its defect and propagation rigorously, and then certify an error bound. They do **not** supply an HDBLAST certificate. The present route can use a scalar comparison enclosure for a regular shell root whose existence is retained from the theorem; it need not prove a new invariant-manifold existence theorem.

Scope: three query framings (invariant-manifold residual validation; correlated validated IVP enclosures; Einstein–scalar shell/linearization), including a 2025–23 October 2026 date window, and four retained primary records. Searches were bounded to five hits per query/refinement. This is a methods comparison, not an exhaustive literature or novelty assessment. No empirical result below is evidence for the HDBLAST claims.

## Primary records

| Source | Extent actually checked | Applicability and limit |
| --- | --- | --- |
| **Olivier Hénot and Akitoshi Takayasu (2026), _Computer-Assisted Proofs in Dynamical Systems: A Case Study of a Heteroclinic Orbit in the Shimizu–Morioka System_.** [arXiv:2605.07500v1](https://arxiv.org/abs/2605.07500v1); [full HTML](https://arxiv.org/html/2605.07500v1). | Abstract; full-text §1.1, Theorem 1.1; §3 invariance equation (14), sequence norm (18), residual estimate (35), Proposition 3.3 and error bound (39). Other sections were searched for context, not independently checked line by line. | A current, concrete example of validating a polynomial manifold approximation by a residual **and** operator/tail estimates. Their theorem requires an injective approximate inverse and bounds `Y`, `Z(R)<1`, with `Y/(1−Z)≤R`; a small residual alone is insufficient. The parameterization uses validated equilibrium/eigenpair data and nonresonance assumptions. These hypotheses are not established here for the shell equations. The useful transfer is certificate design, not the paper's existence, uniqueness or transversality conclusion. |
| **Salma Iraky and Andrew Sogokon (2026), _TERA: A Unified Taylor Model Enabled Reachability Analysis Framework_.** [arXiv:2607.01189v1](https://arxiv.org/abs/2607.01189v1); [full HTML](https://arxiv.org/html/2607.01189v1). | Complete short HTML text, with attention to §II and its notes on outward/correct rounding and dependency preservation. The implementation and benchmark certificates were not audited. | Represents trajectories as `P(x0,t)+I`, preserves functional dependence in the polynomial, and uses a rigorous interval remainder; its text explicitly requires correct rounding. This directly motivates carrying one shared uncertainty packet through an actual-minus-nominal calculation and bounding the resulting polynomial difference. It does not prove that the local paired jets, truncation remainder, or background packet used here are valid. Its stochastic probability bounds do not apply to a deterministic shell certificate. |
| **Martin Berz and Kyoko Makino (1998), _Verified Integration of ODEs and Flows Using Differential Algebraic Methods on High-Order Taylor Models_, Reliable Computing 4(4), 361–369.** [DOI](https://doi.org/10.1023/A:1024467732637); [Crossref metadata](https://api.crossref.org/works/10.1023%2FA%3A1024467732637). | Exact-title Crossref record checked: authors, title, journal, date, pages and DOI. No abstract was returned. Full text was **not** read; the listed publisher host was outside the configured destination allowlist and was not requested. TERA §II explicitly identifies this as its foundational validated ODE method. | Established primary provenance for Taylor-model flow enclosure. Retained as a bibliographic pointer, not independent verification of a specific algorithm or theorem. No detailed claim about its hypotheses is imported from inaccessible full text. |
| **Marc Mars (2005), _First and Second Order Perturbations of Hypersurfaces_.** [arXiv:gr-qc/0507005v1](https://arxiv.org/abs/gr-qc/0507005v1); [primary metadata/abstract](https://export.arxiv.org/api/query?id_list=gr-qc/0507005). | arXiv title, author and full abstract only. | Covers induced-metric/extrinsic-curvature perturbations when both the spacetime metric and hypersurface move, with an application to perturbed matching. Relevant if shell matching is perturbed. It is not a validated IVP method, an Einstein–scalar shell existence result, or a precision enclosure for the retained root. No full-text formula from it is used. |

## Transfer to the current precision calculation

Use `z=p/eta`, `H=rho^(-2)` and the exact reduced vector field, including

`z_y = U'(eta)/eta − 4 a z − z^2`.

For a polynomial graph `z=Z(x)`, `x=(eta,H)`, write the reduced system as `x_y=F_x(x,z)`, `z_y=F_z(x,z)` and define the graph defect

`R_Z(x) = DZ(x) F_x(x,Z(x)) − F_z(x,Z(x))`.

Along the exact retained solution, `e=z−Z(x)` obeys the identity

`e_y = F_z(x,Z+e)−F_z(x,Z) − DZ(x)[F_x(x,Z+e)−F_x(x,Z)] − R_Z(x)`.

Thus a scalar majorant must retain the coupling through `eta_y=eta z` and any dependence of `a` on the state. Enclose the defect and the complete error right-hand side with directed rounding, prove the comparison tube is invariant over its stated domain, and apply it to the theorem-retained existent regular root. This is a necessary enclosure for that root; a new existence, uniqueness, stability or shell-matching theorem does not follow.

For a local paired jet, form the actual-minus-nominal background and observable differences in the **same** packet, including their common parameters and all truncation/roundoff remainders. Subtract coefficients/expressions before taking a range. Two independent interval evaluations followed by subtraction discard dependency. Increasing arithmetic precision alone does not restore that dependency when input uncertainty dominates.

## Coverage and retrieval limits

The 2025–2026 primary arXiv queries produced the two directly relevant 2026 records above. No 2025 record is promoted to the core comparison after this bounded screening. That does not establish an absence of relevant 2025 work. Recent Einstein–scalar shell hits concerned additional Maxwell/Gauss–Bonnet/other sectors or unrelated observables; they were screened out rather than treated as support for these equations.

The connected Firecrawl account returned **insufficient credits** for all three initial searches; none ran, and none was retried there. Retrieval continued through allowed arXiv and Crossref APIs over the configured proxy. The current network policy was observed as enforced, with current observations. Exact queries, returned hit identifiers and retrieval extents are in [search_provenance.json](search_provenance.json). No third-party PDF was copied or published. This audit changes no research manuscript or frozen result.
