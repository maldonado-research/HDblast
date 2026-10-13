# Bounded primary-literature comparison for the regular-cone phase-space checkpoint
Public adaptation: executor paths and cache-location details are removed; scientific statements, source attribution, reading scopes and limitations are preserved. Original source bodies are not redistributed.

Reading/search date: 13 October 2026. This is a targeted comparison with selected
original papers, not an exhaustive survey or a priority certificate. Project
inputs were the current `JACOBI_RESPONSE_TO_SPECTRAL_POSITIVITY.md`,
`NEXT_CALCULATIONS.md`, and the prior bounded stability review. This note does
not independently certify the checkpoint's new field calculations.

## Main finding

The literature already supplies the covariant presymplectic current, the
boundary criterion for a diffeomorphism to be a degeneracy, the effect of GHY
and corner terms, and radiative null-boundary charge/flux prescriptions. A
particularly direct de Sitter-brane precedent already separates radial
normalizability from smoothness at the light-cone horizon: Garriga and Sasaki
explicitly say their normalizable graviton zero mode is singular there. A
model-specific reconstruction and flux calculation is therefore needed to
identify which physical completion the HDBLAST radial theorem controls.

The narrow potential contribution is an explicit dictionary for this doubled
Einstein–scalar shell between its existing radial energy domain, reconstructed
covariant perturbations, the chosen cone-limit prescription, and the allowed
zero-charge gauge quotient. Neither the names of the phase-space methods nor
the observation that a normalizable separated graviton can be horizon singular
is new. The present bounded search does not establish whether the exact
HDBLAST dictionary occurs in other work.

## Source register: original passages actually read

Page numbers below are printed manuscript pages, not journal pagination.
"Selected passages" means the sections listed were inspected; it does not
mean a complete-paper reading. Original source bodies are not redistributed in this subset.

| Source and identifier | Reading scope | Exact useful statement | Transfer limit |
|---|---|---|---|
| Vivek Iyer and Robert M. Wald, *Some Properties of Noether Charge and a Proposal for Dynamical Black Hole Entropy*, [gr-qc/9403028v1](https://arxiv.org/abs/gr-qc/9403028v1), DOI [10.1103/PhysRevD.50.846](https://doi.org/10.1103/PhysRevD.50.846) | §3 pp.8–10, §6 pp.14–15; selected passages | Potential ambiguity `Θ→Θ+δμ+dY`, Eq.(43); current and Cauchy integral Eqs.(44)–(45); boundary ambiguity Eq.(47); Noether-current identity Eq.(75), surface variation Eq.(77); Hamiltonian existence requires boundary integrability Eq.(80), and on-shell H is a surface term Eq.(82). | These formulas do not prescribe the HDBLAST field space. Noncompact Cauchy convergence, falloff and surface independence are explicitly extra hypotheses, p.10. Orientations and antisymmetrization must be translated consistently. |
| Robert M. Wald and Andreas Zoupas, *A General Definition of “Conserved Quantities” in General Relativity and Other Theories of Gravity*, [gr-qc/9911095v1](https://arxiv.org/abs/gr-qc/9911095v1), DOI [10.1103/PhysRevD.61.084027](https://doi.org/10.1103/PhysRevD.61.084027) | §§2–4 pp.7–20; selected passages; §5 null-infinity formulas not used as an HDBLAST theorem | On shell `ω(δφ,Lξφ)=d(δQ−ξ·θ)`, Eq.(13). Hamiltonian integrability requires `∫∂Σ ξ·ω=0` for every pair of allowed linearized solutions, Eq.(18), with differentiability/simple-connectivity qualifications. Case I has vanishing boundary pullback of ω; Case II can have radiation and nonconserved charges, pp.14–16. Boundary potential and modified charge are Eqs.(25)–(26); flux is `Fξ=Θ(φ,Lξφ)`, Eq.(32). | §4 assumes the **full** ω extends continuously to B, not only its pullbacks to approaching slices (p.12 footnote6), and smooth compactified cuts (p.13). The boundary potential's existence/uniqueness and stationary reference require proof. An energy completion need not satisfy these smooth-extension assumptions. |
| Daniel Harlow and Jie-qiang Wu, *Covariant phase space with boundaries*, [1906.08616v4](https://arxiv.org/abs/1906.08616v4), DOI [10.1007/JHEP10(2020)146](https://doi.org/10.1007/JHEP10(2020)146) | §§2.1–2.4 pp.7–22; scalar §3.2 pp.24–25; GR §3.5 pp.28–31; §4.3 p.48; selected passages | Boundary variational relation `Θ+δℓ=dC` after restriction to Γ, Eq.(2.20), gives `ω=δ(Θ−dC)`, Eq.(2.23), and `Ω=δ[∫ΣΘ−∫∂ΣC]`, Eq.(2.27). Covariant charge includes C, Eq.(2.51). GR Eqs.(3.43)–(3.45) give the GHY completion and `c^μ=−γ^{μν}n^αδg_{να}/(16πG)`. Physical quotient uses the zero modes of Ω, Eq.(2.8). | p.3 explicitly: “in this paper we do not consider null boundaries.” The GR example fixes the induced boundary metric, Eqs.(3.46)–(3.48). A dynamical common shell with summed junctions and moving embeddings requires a new application; the null limit and its convergence are additional calculations. |
| Venkatesa Chandrasekaran, Éanna E. Flanagan and Kartik Prabhu, *Symmetries and charges of general relativity at null boundaries*, [1807.11499v3](https://arxiv.org/abs/1807.11499v3), DOI [10.1007/JHEP11(2018)125](https://doi.org/10.1007/JHEP11(2018)125) | §2.3 pp.10–11, §5.1 pp.35–38, §6.4–6.5 pp.46–49; selected passages | Null field space fixes an equivalence class of universal structure `p=[ℓ^a,κ,ℓ̂_a]`, Eqs.(5.1)–(5.4). p.36 explicitly distinguishes diffeomorphism freedom from symplectic degeneracy. Their localized potential is `Θ_ijk=(ε_ijk/16π)h^{AB}(σ_AB−q_ABθ/2)`, Eq.(6.26), with charge/flux Eqs.(6.27)–(6.32). | The explicit construction is vacuum GR in **3+1 dimensions** and its chosen field space. Its vanishing-potential criterion is shear **and expansion** free. A 5D light cone has a three-dimensional screen and is expanding; matter and the cone vertex require separate treatment. The source does not establish zero charge for HDBLAST bending. |
| Jaume Garriga and Misao Sasaki, *Brane-world creation and black holes*, [hep-th/9912118v2](https://arxiv.org/abs/hep-th/9912118v2), DOI [10.1103/PhysRevD.62.043523](https://doi.org/10.1103/PhysRevD.62.043523) | §§I–III pp.2–11 and selected later passages; previously retrieved original; detailed notes in `BRANE_PRIMARY_READING.md` | The dS-sliced chart covers only the exterior of the r=0 light cone, pp.3–4. Normalizable tensor zero mode and continuum threshold occur in Eqs.(2.12)–(2.14), pp.6–7. §III p.9 gives `C_(5)^2=C_(γ)^2/a^4(r)` and explicitly calls the zero mode singular on the Rindler horizon. p.10 says small zero-mode fluctuations still have finite action; the small-r GHY term vanishes. | Scalar-free AdS interiors differ from HDBLAST. The authors' Green-function explanation cites regular isolated-source fields in the flat RS setting; it is not a proved arbitrary-packet completeness theorem for the dressed dS cone. Vanishing Weyl contraction alone is not a complete smoothness criterion. |
| Uchida Gen and Misao Sasaki, *Radion on the de Sitter brane*, [gr-qc/0011078v2](https://arxiv.org/abs/gr-qc/0011078v2), DOI [10.1143/PTP.105.591](https://doi.org/10.1143/PTP.105.591) | §§I–IV pp.1–7 and Appendix A pp.7–10; previously retrieved original; detailed notes separate | Projected-Weyl radial orthogonality uses self-adjointness and weight `a^(3−n)dz`, Eqs.(2.19)–(2.21), p.4. Local displacement/relative displacement distinction and one-vacuum-brane radion limit occur pp.6–7. Appendix A retains bulk-domain and surface displacement contributions before Eq.(A18). | Projected-Weyl Green-function norm is not a metric/scalar action norm. Its special constant radial mode differs from the metric graviton zero mode. Local gauge/bending observations do not determine a global cone charge quotient. |
| Lev Kofman, Johannes Martin and Marco Peloso, *Exact identification of the radion and its coupling to the observable sector*, [hep-ph/0401189v1](https://arxiv.org/abs/hep-ph/0401189v1), DOI [10.1103/PhysRevD.70.085015](https://doi.org/10.1103/PhysRevD.70.085015) | pp.2–7 and Appendix A p.16, with delegated PDF verification of Eqs.(26)–(27); previously retrieved original | §III B p.6 explicitly retains y-total derivatives. Reduced action Eq.(19), Wronskian pairing Eq.(22), boundary-subtracted norm Eq.(26), and positive physical coefficient Eq.(27) establish that the linear eigenproblem alone does not normalize physical modes. Eq.(27)'s outer measure is `dy/A`, verified visually in PDF. | Flat slices and a finite two-brane orbifold; it supplies boundary-norm precedent, not a cone symplectic-completion theorem. Its divided variables require separate scalar-free/degenerate analysis. |
| Marco Bertola, Jacques Bros, Vittorio Gorini, Ugo Moschella and Richard Schaeffer, *Decomposing Quantum Fields on Branes*, [hep-th/0003098v1](https://arxiv.org/abs/hep-th/0003098v1), DOI [10.1016/S0550-3213(00)00280-7](https://doi.org/10.1016/S0550-3213(00)00280-7) | §2 pp.5–8, §3 pp.9–10, §4 opening and §4.1 pp.11–14; selected passages | The radial Hilbert weight is `ω^(d−2)dv_X`, Eq.(23). Individual spectral fields are distributions requiring smearing, and inversion is checked case by case, p.6 after Eqs.(25)–(26). Canonical commutator decomposition uses a product Cauchy surface and correctly rescaled time normal, Eqs.(32)–(35). Extra global-invariance, normal-analyticity and convergence hypotheses give diagonal two-point decomposition, Eqs.(46)–(48). §4.1 foliates the spacelike exterior of the Minkowski origin by dS leaves, and Eqs.(49)–(56) give an essentially self-adjoint radial decomposition with continuum threshold `(d−1)^2/4`, hence 9/4 for d=4. | A free positive-mass scalar in the Minkowski vacuum, entire exterior R∈(0,∞), no shell and no coupled gravity. Its threshold does not exclude shell-bound HDBLAST modes. It supplies an explicit ambient spectral/state dictionary under additional hypotheses, not a smooth classical metric/scalar completion or null flux theorem. The PDF's later TeX-generation date is not its arXiv version date. |
| Luca Ciambelli, Laurent Freidel and Robert G. Leigh, *The Phase Space of Gravity on Null Hypersurfaces*, [2608.14449v1](https://arxiv.org/abs/2608.14449v1), submitted 14 August 2026; no DOI supplied in retrieved metadata | §2.1 p.3; §3.1–3.2 pp.10–13; §3.2 endpoint/flux discussion p.15; §4.2 pp.20–21; conclusions pp.29–30; selected passages | With specified corner polarization, canonical null potential Eq.(60) gives diffeomorphism charge/flux Eqs.(62)–(65). Generator rescaling has corner-area charge Eq.(73); Ehresmann-connection shifts satisfy `I_ζΩ=0`, Eq.(76). Endpoint cuts approaching caustics are allowed only when their boundary terms are well defined, p.15. After transport Eqs.(163)–(164), p.20 explicitly requires function-space and boundary choices, fixes/projects possible zero/edge modes, and acknowledges different admissible Green kernels. | Derivation uses 4D vacuum GR on a caustic-free segment and kinematical reduction **before** Raychaudhuri/Damour constraints. p.30 states a higher-dimensional extension but leaves constraint reduction, finite-region corner flux/edge completion, and caustic behavior as future work. It is recent relevant precedent, not a completed 5D dressed-cone+shell phase-space theorem. |

J. Lee and Robert M. Wald, *Local symmetries and constraints*,
*J. Math. Phys.* **31**, 725–743 (1990), DOI
[10.1063/1.528801](https://doi.org/10.1063/1.528801), was located as
INSPIRE record 300075. Its original full text was **not retrieved or read**;
the returned metadata/abstract establish identity only. No exact equation or
theorem is attributed to a new original reading of that paper here. The
directly read Iyer–Wald and Harlow–Wu originals provide the current formalism
comparison.

## What the sources require for the present checkpoint

1. **Declare the physical region and null limit.** A dS radial chart's y→0
   limit with bounded slice coordinates collapses toward the cone vertex;
   reaching a finite nonzero null generator requires a correlated radial/time
   limit in inertial or null coordinates. An integral over shrinking timelike
   cuts at fixed tangential time and a characteristic radiation integral are
   different questions. The geometry observation follows from the HDBLAST
   leading metric `dy²+y²γ`; Garriga–Sasaki supplies the primary precedent for
   the chart/light-cone distinction. No source read permits silently
   exchanging these limits.
2. **State the completion.** Smooth extendible covariant perturbations,
   closures of compactly supported data in the energy norm, and fields with
   externally specified characteristic data need not be the same space.
   Radial self-adjointness controls the Green pairing on the specified radial
   domain. It does not alone prove characteristic flux vanishes or impose
   inertial smoothness on every separated mode. The known singular yet
   normalizable tensor zero mode makes this distinction concrete.
3. **Retain the action and corner information.** For the two interiors, sum
   bulk potentials and GHY boundary completions, then include one common shell
   action and its embedding variation. Shell/junction cancellation must hold
   for the admitted paired variations. A regulator boundary becomes harmless
   only after its flux and corner limits have been evaluated, with declared
   orientation and regularity. Harlow–Wu gives a derivation template, not
   the missing null-limit theorem.
4. **Classify gauge by its pairing.** On shell the contraction with a
   diffeomorphism is a boundary variation. A proposed generator is a pure
   redundancy only if it preserves the specified state space and has zero
   pairing with every admitted variation, including shell/corner terms.
   Merely having `δg=Lξg`, `δη=Lξη` locally does not prove this. A constant
   generator charge may be set by reference normalization; the decisive
   degeneracy statement is that its variation/pairing vanishes. An internal
   artificial cut's contributions can cancel between regions if the full
   extension is retained, while treating that cut as a physical boundary can
   retain charge/flux data.
5. **Separate quadratic energy from charge flux and smoothness.** The
   action-derived positive generic kinetic coefficients already found in
   HDBLAST are relevant to energy and the reduced symplectic bilinear.
   They do not prove the reconstructed tensor/scalar curvature has a smooth
   extension. Conversely, excluding an individually singular separated
   mode from a smooth class does not automatically exclude its profile from
   an energy completion. A characteristic flux obstruction limits the
   physical interpretation without changing an already proved radial
   self-adjoint operator theorem.

These are application requirements and source-scope limits, not a claim that
the HDBLAST completion necessarily fails. The new derivation should make
its admitted data, boundary conditions and result explicit.

For the scalar contribution, Harlow–Wu §3.2 pp.24–25 explicitly gives
`θ_φ^μ=−∇^μφ δφ`, Eq.(3.10), and its fixed-background scalar symplectic
current, Eqs.(3.13)–(3.15). In coupled Einstein–scalar theory, the metric and
volume variation must also be taken when antisymmetrizing the scalar
potential; the fixed-background example is not the entire coupled current.
Wald–Zoupas §6 p.35, additionally inspected here, explicitly says that a
minimally coupled scalar does not add the Einstein–Klein–Gordon Noether
charge Q term, but does add to θ, ω and the flux. Their concrete statement
uses asymptotically flat null infinity and a smooth conformally rescaled
massless scalar; those falloff hypotheses must not be imposed on the HDBLAST
finite cone by analogy.

The brane-source comparison was checked again by a separate primary reader
against its detailed notes, with no material correction. A useful precision
for future wording is that `C²` is quadratic in the linearized Weyl tensor
around conformally flat AdS. The quoted source establishes generic
non-conformally-flat zero-mode horizon singularity; the contraction's
vanishing for a special null field is not a proof of regularity.

## Search evidence and limits

`PHASE_SEARCH_LOG_20261013.json` records five queries with UTC retrieval
times, exact URLs, caps, status codes, response hashes and returned metadata:
three INSPIRE and two arXiv queries, each capped at 20 items. The recent
queries use a 2025–13 October 2026 submission cutoff. Returned totals are
not reading counts. The broad recent INSPIRE query found 98 records, of
which only its first 20 metadata records were inspected; substantial
irrelevant matches demonstrate that it is not an exhaustive topical survey.
The dS-foliation INSPIRE conjunction returned zero but the arXiv query
returned 43, of which the first 20 were inspected. Zero results from one
query do not establish absence of relevant work.

The recent search located *The Phase Space of Gravity on Null Hypersurfaces*,
[2608.14449v1](https://arxiv.org/abs/2608.14449v1), submitted 14 August 2026.
Its original was subsequently fetched and the specified sections read,
including the explicit boundary-kernel qualifications and future-work scope
registered above. It is newly inspected here, not newly released since the
previous checkpoint. The older slicing source was selected by known title
and verified by a further INSPIRE query (`arxiv:hep-th/0003098`, cap5,
one returned record), in addition to the five discovery queries. Recent work
exists; absence of an exact HDBLAST match in this bounded search supplies no
first-ever or global-novelty claim.

Six new original PDFs were retrieved and converted
with `pdftotext -layout`, and the selected passages registered above were
read. `PHASE_PRIMARY_SOURCE_RECEIPTS.json` pins the first four versions,
metadata, URLs, dates, bytes and hashes; `PHASE_FOLLOWUP_SOURCE_RECEIPTS.json`
pins the two direct-PDF readbacks and slicing identity check. Although the
follow-up PDF URLs were unversioned, each PDF has an explicit v1 arXiv header,
and the retained content hash identifies the body actually read. Three earlier original brane sources were
read from previously retrieved originals; selected earlier download receipts
are included in `CACHED_BRANE_SOURCE_RECEIPTS.json`. `BRANE_PRIMARY_READING.md` records the
precise delegated reading scope and the PDF-specific equation checks.

Runtime and networking instructions were read; policy permitted these
destinations and `observations_current` was true, but runtime HTTP policy
state was `unknown`. Specific successful readbacks establish those HTTP
responses, not global policy enforcement or secret readiness. Network
reads retained the inherited managed proxy, configured CA trust and TLS
verification. No credentials were printed. Follow-up arXiv API metadata
retrieval timed out once and returned HTTP429 on a later retry. After that
the API was not retried; selected direct-PDF URLs and the separate INSPIRE
identity query succeeded through the same proxy/TLS path. The API failures
do not establish version metadata; the actual original readbacks and header
versions are recorded separately. No proxy or TLS
bypass, external writing, publication, or repository editing occurred.

The review does not claim full citation-graph coverage, priority,
nonlinear stability, cosmological origin, observational discrimination,
or all-sector physical completeness.
