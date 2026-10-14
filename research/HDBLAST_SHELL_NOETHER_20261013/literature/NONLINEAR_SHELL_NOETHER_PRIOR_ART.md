# Nonlinear shell Noether identities: prior art and transfer limits

Prepared 14 October 2026 UTC; literature cutoff **13 October 2026 inclusive**.
This is a bounded, selective primary-source audit for the HDBLAST continuation.
It is not an exhaustive search, external review, or novelty determination.
Equation and page references below refer to the identified arXiv versions,
except for the separately identified published erratum. Full sources and raw
API responses remain in private caches and are not publication inputs.

## What the sources establish

Covariant phase space with a completed boundary action, moving integration
domains, edge-mode extensions, Hamiltonian thin-shell mechanics, and quadratic
linearization-stability conditions are established subjects. The model-specific
task is to combine them for **two Einstein–scalar regions with exactly matched
induced fields, two outward GHY terms, and one common scalar-dependent tension
action**, with the signs and cut terms fixed by that action.

The papers read here do not, by citation alone, establish the unrestricted
nonlinear degeneracy identity for that coupled 4+1 dimensional interface. In
particular, neither conservation of a restricted linearized current nor a
standard formula for a smooth external boundary supplies the missing internal
shell cancellation. Necessary quadratic conditions also do not establish an
exact nonlinear branch or a completed constraint-range theorem.

## Primary sources actually read

| Source | Exact reading scope and useful formulas | Limit of transfer |
|---|---|---|
| V. Iyer and R. M. Wald, *Some Properties of Noether Charge and a Proposal for Dynamical Black Hole Entropy*, [gr-qc/9403028v1](https://arxiv.org/abs/gr-qc/9403028v1), [DOI](https://doi.org/10.1103/PhysRevD.50.846) | §3 pp.9–10, Eqs.(43)–(45); §6 pp.14–15, Eqs.(73)–(82). Potential ambiguity, symplectic integral, fixed-vector Noether identity and Hamiltonian relation. After Eq.(82), an on-shell Hamiltonian on a compact closed Cauchy slice vanishes. | Smooth bulk formalism, not a derived Israel interface identity. Keep the chosen symplectic sign convention; the printed v1 Eq.(44) orders the variations oppositely to Wald–Zoupas Eq.(4). |
| R. M. Wald and A. Zoupas, *A General Definition of “Conserved Quantities” in General Relativity and Other Theories of Gravity*, [gr-qc/9911095v1](https://arxiv.org/abs/gr-qc/9911095v1), [DOI](https://doi.org/10.1103/PhysRevD.61.084027) | §§2–3 pp.6–11, Eqs.(2)–(17), especially (4), (8), (10), (13)–(16). On a background solution, arbitrary variation gives `ω(δ, Lξ)=ξᵃδCₐ+d(δQξ−iξθ)`. Linearized solutions remove `δC`; boundary integrability remains a separate question. | Does not identify a shell constraint distribution or its adjoint trace domain. The full off-shell identity has additional equation-of-motion terms. |
| D. Harlow and J.-q. Wu, *Covariant phase space with boundaries*, [1906.08616v4](https://arxiv.org/abs/1906.08616v4), [DOI](https://doi.org/10.1007/JHEP10(2020)146) | §2.1 pp.7–11, §2.2 pp.12–14, §2.3 pp.20–22; GR §3.5 pp.28–31. Eqs.(2.20), (2.23), (2.27), (2.47)–(2.51), (3.43)–(3.48), (3.60)–(3.62). | The GR example fixes the induced metric on a fixed timelike external boundary. The paper expressly flags singularities at symmetry fixed points and the heuristic functional-analytic treatment. It is not a smooth-manifold theorem at a symmetric shell background. |
| W. Donnelly and L. Freidel, *Local subsystems in gauge theory and gravity*, [1601.04744v2](https://arxiv.org/abs/1601.04744v2), [DOI](https://doi.org/10.1007/JHEP09(2016)102) | §2.5 pp.15–17, Eqs.(2.44)–(2.53); §3 pp.17–22, Eqs.(3.5), (3.18)–(3.32). Gluing is a constraint reduction by the diagonal boundary group; gravitational embedding transforms inversely to the bulk diffeomorphism. | The displayed detailed fusion example is Yang–Mills. The gravity calculation at (3.8), (3.28) uses vacuum equations. Its artificial subregion edge variables are not automatically physical shell matter. |
| A. J. Speranza, *Local phase space and edge modes for diffeomorphism-invariant theories*, [1706.05061v2](https://arxiv.org/abs/1706.05061v2), [DOI](https://doi.org/10.1007/JHEP02(2018)021) | §§3–4 pp.9–14, Eqs.(3.1)–(3.9), (4.1)–(4.16). Pullback variation, inverse embedding gauge action, extended potential and all moving-cut terms. | Auxiliary edge-mode construction for a subregion. Adding its corner charge to a physical-shell potential needs an action/quotient justification; it cannot be used to assume the desired gluing identity. |
| P. Hájíček and J. Kijowski, *Lagrangian and Hamiltonian Formalism for Discontinuous Fluid and Gravitational Field*, [gr-qc/9707020v1](https://arxiv.org/abs/gr-qc/9707020v1), [DOI](https://doi.org/10.1103/PhysRevD.57.914), **with [published erratum](https://doi.org/10.1103/PhysRevD.61.129901)** | Abstract/introduction pp.1–5; §§3.2–3.4 pp.29–34, Eqs.(92)–(107); §3.5 matching Eqs.(110)–(112); §3.6 pp.41–43; all four pages of the erratum. Two bulk actions, one shell action, induced matching, Israel variation, area/boost-angle corner pair, and differentiability limitations. | 3+1 dimensional fluid shell; adapted-coordinate embeddings are held fixed. This is direct physical-interface precedent but not a 4+1 Einstein–scalar theorem. The erratum must accompany any use of its detailed formulas. |
| E. Altas and B. Tekin, *Second Order Perturbation Theory in General Relativity: Taub Charges as Integral Constraints*, [1903.11982v2](https://arxiv.org/abs/1903.11982v2), [DOI](https://doi.org/10.1103/PhysRevD.99.104078) | §I pp.1–3, Eqs.(1)–(9); beginning §II pp.3–4, Eqs.(10)–(16). Quadratic source plus linear correction, Killing contraction, boundary ADT term, compact closed-slice Taub necessity. | Smooth vacuum/cosmological Einstein backgrounds. No internal Israel jump, scalar matching or shell delta derivative. No shell coefficient or sufficiency statement is imported. |
| E. Altas, E. Kilicarslan and B. Tekin, *Einstein–Yang–Mills Theory: Gauge Invariant Charges and Linearization Instability*, [2105.11744v2](https://arxiv.org/abs/2105.11744v2), [DOI](https://doi.org/10.1140/epjc/s10052-021-09460-7) | §III pp.6–8, Eqs.(34)–(47), conclusion and classical bibliography. A coupled-field example deriving a weaker necessary integral condition from the second-order equation. | Yang–Mills gauge-sector test. The authors explicitly do not analyze the gravitational linearization-instability sector there. It supplies no Einstein–scalar-shell theorem. |
| H. Adami, M. Golshani, M. M. Sheikh-Jabbari, V. Taghiloo and M. H. Vahidinia, *Covariant Phase Space Formalism for Fluctuating Boundaries*, [2407.03259v1](https://arxiv.org/abs/2407.03259v1), [DOI](https://doi.org/10.1007/JHEP09(2024)157) | §§3.1–3.3 pp.5–8, Eqs.(3.1)–(3.17), and opening §3.4. Varying a Heaviside domain produces the displacement term; varying the integrated potential produces a codimension-two corner form. | Fluctuating region boundaries. Its charge construction does not impose the physical shell's Israel/scalar residuals. The abstract's general charge-integrability claim was not audited here. |
| H.-F. Ding, *Conserved charges and the first law of black holes for field dependent symmetry generators in the covariant phase space formalism*, [2609.27454v2](https://arxiv.org/abs/2609.27454v2), 29 September 2026 | Abstract and §2 pp.3–8, especially Eqs.(23)–(30): variable parameters require the subtraction `Qδϵ`; charge density (27), integrability (30). | Recent black-hole/gauge-field application, not a shell or nonlinear existence result. It illustrates why a field-dependent embedding compensator cannot be differentiated as a fixed parameter. |

## The boundary identity and physical gluing

Harlow–Wu's boundary variational relation is

```text
(θ + δℓ)|Γ = dC,
Ω = δ[∫Cauchy θ − ∫cut C].
```

For Einstein–Hilbert plus GHY, their Eq.(3.43) gives, in their conventions,

```text
θ|Γ + δℓ = −(16πG)^−1 εΓ (Kᵃᵇ − Khᵃᵇ) δg_ab + dC,
cᵃ = −(16πG)^−1 hᵃᵇ nᶜ δg_bc,   C = i_c εΓ.
```

For their external Dirichlet problem, the metric term vanishes because the
induced variation is fixed to zero. Their Brown–York Hamiltonian, Eqs.(3.60)–
(3.62), is generally nonzero. For the physical shell, the common induced
metric varies: one instead sums both sides and the **single** shell variation.
Its metric and scalar coefficients must be the coupled junction residuals.
Only after these vanish, and embedding/cut contributions are accounted for,
can the analogous completed relation imply internal cancellation.

Hájíček–Kijowski is a useful independent check. Their Eq.(94) is two smooth EH
bulk integrals minus `(8πG)^−1∫√|γ|[L]` plus one shell action. They use one
directed normal and `[L]=L+−L−`; this converts to the sum with two outward
normals, not equality of the two extrinsic curvatures. Equation (95) makes
`γ=X−*g−=X+*g+` part of the configuration space. Equation (100), **as corrected**,
gives the Israel residual and shell-cut potential
`−(8πG)^−1∫√λ δ[α]` together with the matter cut potential. Equation (107)
retains the corresponding area/angle and matter pairs in the Hamiltonian
variation. A shell-volume term alone would miss this gravitational cut pair.

The erratum is consequential here. Items 13–14 restore square-root density
factors in Eq.(99) and identify the final integral in Eqs.(99)–(100) as an
**external** `Σ+` integral; the uncorrected v1 prints `Σ` there. That term is
not an uncancelled physical-shell residual. Item 16 corrects an angle-label
typo; later items correct signs and densities in component Hamiltonian
formulas. The current audit uses this paper as a structural check, with its
erratum, and does not import its fluid component equations as our model's.

Their pp.42–43 also explain why a singular shell constraint need not be a
differentiable phase-space function when smeared only over the shell. The
total Hamiltonian's differentiability uses the bulk terms and continuity of
bulk/surface smearings. This directly cautions against omitting interface
conditions when identifying the adjoint kernel with Killing initial data.

## Moving embeddings and exact induced matching

Speranza Eqs.(3.2), (3.5) give

```text
δ(X*α) = X*(δα + Lχ α),       Iξ χ = −ξ.
```

Thus a spacetime transformation acts on the embedding by `δξX=−ξ∘X` when
the bulk variation is `Lξ`. Fixing a Killing background including its
parametrized shell requires quotienting the common shell reparametrization
or including the compensating reference-shell vector. A background-tangent
embedding variation is not literally zero before this quotient/compensation.

The same pullback identity gives the following **model application**, rather
than a separately quoted shell theorem. For a differentiable family and a
smooth local extension `χ_ε=∂εX_ε∘X_ε^−1`, second induced variations contain

```text
∂ε²(Xε* gε)|0
 = X0*[g̈ + 2 Lχ ġ + Lχ̇ g + Lχ² g]|0,
```

and the same formula with the scalar. Exact matching requires equality of
these full expressions on the two sides. Equality only of the second-order
bulk fields evaluated at the old graph omits graph/evaluation terms.

Speranza's on-shell extended potential and cut form are, Eqs.(4.9), (4.12),

```text
Θext = ∫Σ(θ + iχL) + ∫∂ΣQχ,
Ωext = ∫Σω + ∫∂Σ[δQχ + LχQχ + iχθ + (1/2)iχiχL].
```

His Eq.(4.7) explains that the potential without `dQχ` need not annihilate
diffeomorphisms; a particular ambiguity choice supplies the extended
degeneracy. Also, `iχL` generally survives in a theory with matter. These
facts make the source useful for auditing moving-domain terms, but they do
not license inserting an auxiliary edge-mode completion into a fixed
physical-shell action and then declaring gluing proved.

Adami et al. provide a complementary distributional derivation: Eqs.(3.3)–
(3.6) retain `L Δ(B)δB`, and (3.8)–(3.11) retain the variation of the domain
and its corner, including a derivative of the delta distribution during the
calculation. Their Eq.(3.16) is a flux balance. Conservation or integrable
subregion charges alone do not imply a common physical-shell gauge generator
has zero charge.

## Quadratic necessity, normalization, and historical scope

With one fixed convention `Ω=δΘ` and `δH_K=Ω(δ,G_K)`, a generator fixing the
whole background gives the Hessian relation
`D²H_K[v,v]=Ω0(v,DG_K[v])`. If the complete constrained Hamiltonian vanishes
on a closed doubled solution and a sufficiently differentiable exact family
exists, differentiating gives the necessary zero. The corresponding Taylor
coefficient is one half of that Hessian when the expansion uses
`Φε=Φ0+εv+(ε²/2)w+…`. Sign and scale belong to the chosen action, orientation,
current convention and generator; they should be derived locally, not copied
from a paper with a different antisymmetrization or perturbation expansion.

Altas–Tekin Eqs.(5)–(9) make the smooth compact-slice obstruction explicit:
the linear correction becomes a boundary term after a Killing contraction,
leaving an integrated quadratic source when the boundary is empty. In a
split shell geometry each bulk slice has a cut. Its cancellation requires
the shell identity; first-order support away from the cut cannot eliminate
an arbitrary second-order correction there.

Historical metadata/reference checks, **not original-paper readings in this
audit**, identify:

- J. Lee and R. M. Wald, *Local symmetries and constraints*, J. Math. Phys.
  31 (1990), 725–743, [DOI](https://doi.org/10.1063/1.528801).
- A. E. Fischer and J. E. Marsden, *Linearization stability of the Einstein
  equations*, Bull. AMS 79 (1973), 997–1003,
  [DOI](https://doi.org/10.1090/S0002-9904-1973-13299-9).
- V. Moncrief, *Spacetime symmetries and linearization stability of the
  Einstein equations. I*, J. Math. Phys. 16 (1975), 493–498,
  [DOI](https://doi.org/10.1063/1.522572).
- J. M. Arms, *Linearization stability of gravitational and gauge fields*,
  J. Math. Phys. 20 (1979), 443–453,
  [DOI](https://doi.org/10.1063/1.524094).
- J. M. Arms and J. E. Marsden, *The absence of Killing fields is necessary
  for linearization stability of Einstein's equations*, Indiana Univ. Math.
  J. 28 (1979), 119–125; and A. E. Fischer, J. E. Marsden and V. Moncrief,
  *The structure of the space of solutions of Einstein's equations. I. One
  Killing field*, Ann. Inst. H. Poincaré A 33 (1980), 147–194. These two were
  identified through the bibliography of the read 2021 original.
- W. Israel, *Singular hypersurfaces and thin shells in general relativity*,
  Nuovo Cimento B 44 (1966), 1–14,
  [DOI](https://doi.org/10.1007/BF02710419), with its separate 1967 correction
  [DOI](https://doi.org/10.1007/BF02712210). The correction is not the original
  1966 article; both appeared in the metadata result.
- A. H. Taub, *Variational principles in general relativity*, the 1970 CIME
  lectures cited in the read 2019/2021 papers. No original text was inspected.

The classical smooth-theory hypotheses and sufficiency/range theorems were
not independently audited here. They must not be restated as a theorem for
piecewise-smooth scalar/tension interfaces. Passing the ten selected
canonical Killing-pairing tests establishes only their necessary vanishing
once the model's full identity is proved; it is not second-order or nonlinear
integrability.

## Search and access limits

Twelve specified searches used INSPIRE (eight) and Crossref (four). Ten
returned results; two Crossref queries returned HTTP 429. Failed queries are
not negative results. One recent INSPIRE query was date filtered through
13 October 2026; the final targeted query covered recent linearization work
and fluctuating-boundary titles. Only the returned first ten INSPIRE or first
five Crossref results were screened, so neither ranking nor empty selections
supports an absence claim. Broad queries produced unrelated uses of “shell,”
“edge modes,” and “linearization stability”; in particular radial thin-shell
stability is not automatically the Taub integrability question.

Seven selected arXiv originals and one APS erratum were freshly retrieved;
three covariant-phase-space originals were reused from a prior private cache.
The 2026 field-dependent-generator paper was selectively read; other recent
hits were title/metadata screened only. All network requests preserved the
configured proxy and TLS verification. Source bodies, hashes, dates, query
text and access failures are retained privately. The public source register
states exactly which sections were read; download and reading counts are
not evidence of exhaustive coverage.
