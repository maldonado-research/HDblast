# Focused high-order method update

Scientific cutoff: **6 November 2026**, continuing the retained review. Retrieval provenance: **12 November 2026 UTC**. The directory name records retrieval, not a later scientific cutoff. No publication after the scientific cutoff is used. This targeted update supplements the eight-paper retained review; it does not repeat that survey or claim exhaustive coverage.

The useful addition is a precise route from rotation symmetry to a first-wave angular estimate, together with a recent primary example of time-first Neumann energy and spatial recovery. These sources support techniques. None supplies the numerical high-order Einstein–scalar shell theorem or a finite nonlinear actual-area observable remainder.

## Primary passages inspected

| Source/version | Relevant passage actually read | Useful content and limits |
| --- | --- | --- |
| Lee Lindblom, Nicholas W. Taylor and Fan Zhang, *Scalar, Vector and Tensor Harmonics on the Three-Sphere*, [arXiv:1709.08020v1](https://arxiv.org/html/1709.08020v1), 23 September 2017; related journal DOI [10.1007/s10714-017-2303-y](https://doi.org/10.1007/s10714-017-2303-y). | §§2–4: scalar eigenvalue (1); scalar-derived vector definitions and eigenvalues (12), (15); tensor constructions (25), (28), and eigenvalues (35)–(40). | Gives the exact rank-dependent rough-Laplacian shifts on the round three-sphere. This is particularly useful for catching a componentwise scalar-eigenvalue substitution. It supplies harmonic identities, not the coupled time/radial evolution or numerical Sobolev conversion. |
| Alan Rios Fukelman, Matías Sempé and Guillermo A. Silva, *Notes on Gauge Fields and Discrete Series representations in de Sitter spacetimes*, [arXiv:2310.14955v3](https://arxiv.org/html/2310.14955v3), revised 27 December 2023; published [3 January 2024, JHEP 2024, 11](https://doi.org/10.1007/JHEP01(2024)011). | Appendix B.2: Casimirs (B.8)–(B.13), Lie-derivative realization (B.27), scalar Laplacian normalization (B.30)–(B.33), and the representation labels following (B.33). | Identifies scalar degree k with the SO(4) representation [k/2,k/2]. Its Casimir is negative on nonconstant harmonics, so its sign must be reversed for a positive angular energy. The main paper concerns free gauge fields; only the geometric representation identities are used here. |
| Soichiro Katayama, Hideo Kubo and Sandra Lucente, *Almost global existence for exterior Neumann problems of semilinear wave equations in 2D*, [arXiv:1208.3531v1](https://arxiv.org/html/1208.3531v1), 17 August 2012. | Lemma 2.2 and its proof; §§4.1–4.2, especially (4.2)–(4.3). | Time differentiation preserves its fixed homogeneous Neumann rows. Spatial derivatives are recovered by elliptic estimates and replacing the Laplacian through the equation, rather than imposing homogeneous Neumann data on arbitrary spatial derivatives. The scalar two-dimensional exterior-domain theorem does not transfer to the four-spatial-dimensional tensor shell or supply its constants. |
| Fei Hou, Huicheng Yin and Meng Yuan, *Long time smooth solutions of 2-D quadratic quasilinear wave equations in exterior domains with Neumann boundary conditions*, [arXiv:2602.14019v1](https://arxiv.org/html/2602.14019v1), 15 February 2026. | Boundary admissibility and compatibility (1.3)–(1.4); elliptic Lemma 2.1; Lemma 4.1; Lemma 4.3 and its recovery recurrence (4.18)–(4.22); selected good-unknown/boundary terms in §5.1. | A new-to-this-review 2026 quasilinear example: first control time energies, then recover spatial derivatives. It explicitly assumes a principal nonlinearity compatible with Neumann flux and treats modified-unknown boundary terms. Its two-dimensional scalar exterior geometry, structural assumptions and decay mechanism differ from HDblast. Its unspecified estimate constants cannot be inserted as model-specific numerical constants. |

The relevant full text was retrieved for every row; the conclusions above do not rest only on abstracts. Exact version metadata was checked at the primary arXiv pages.

## Direct application to the retained first wave

The following is a **model-specific inference/derivation to verify in the new proof**, not a theorem quoted from the papers.

Let X_AB=n_A∂_B−n_B∂_A, 0≤A<B≤3, be the six unit-sphere rotation generators, and let T_AB be their natural Lie actions on every scalar/vector/symmetric-tensor slot. Include their zeroth-order fiber actions in the moving frame. Define the positive group Casimir

    C = −Σ_(A<B) T_AB².

On scalars C=−Δ_κ. The registered initial Y=n_0²−n_1² is the restriction of a homogeneous harmonic polynomial of degree two, so CY=8Y. Equivariant operations such as covariant differentiation, multiplication by κ, trace and contraction preserve its group-representation label. The completed first initial jet must be checked in this sense; individual tensor components are not scalar functions with the same rough-Laplacian equation. The representation and scalar normalizations are supported by [Appendix B.2](https://arxiv.org/html/2310.14955v3#A2.SS2).

For example, on a unit S3 the rough-Laplacian eigenvalues at k=2 are −8 for Y, −6 for its gradient, −8 for Yκ, and −2 for its traceless scalar-derived Hessian. These follow from [equations (1), (15), (35), (38)](https://arxiv.org/html/1709.08020v1). They describe a common scalar-generated group type with different connection-Laplacian shifts.

A rigorous propagation argument is short once its interfaces are checked:

1. Prove the **complete** background linearized operator L_bg commutes with the natural rotations. This includes reference-connection, frame and ordered lower maps, not only the principal scalar symbol.
2. Check the fixed D/N projections, both boundaries, the normal operator and the completed Cauchy data are equivariant. Rotation acts at fixed time and normal coordinate.
3. For the homogeneous first wave, w=(C−8)V1 then solves the same homogeneous linear problem with zero Cauchy and boundary data.
4. Apply uniqueness in the retained class to obtain CV1=8V1 throughout its common existence interval.

This preserves the complete scalar-generated isotypic sector. It does not assert that every component remains a scalar multiple of Y, or that the nonlinear exact family has a finite angular expansion.

There is also a useful exact ordered-word identity. For any section u in this eigenspace, with the invariant positive L2 bundle product, each T_i is skew-adjoint and [C,T_i]=0. Consequently

    Σ_(all ordered words I of length m) ||T_I u||²_(L2)
        = 8^m ||u||²_(L2).

For m=1 this is integration by parts. Induction applies the same identity to each shorter word, which remains in the eigenvalue-eight subspace. The assertion requires the invariant natural bundle action; a bare component derivative does not have this property.

To use this for the retained N14 budget, the proof still needs explicit finite constants converting lifted Lie words to the specified covariant-slot Sobolev norm, including extra derivative slots and all time/normal recovery terms. Rotation symmetry by itself does not compute those constants.

## Boundary-safe route for the high-order estimate

The two Neumann papers above independently exhibit the same useful organization: commute in time while preserving the actual fixed boundary condition, then use the equation and an elliptic/spatial recurrence to recover derivatives that do not preserve it. The 2026 paper's admissibility condition (1.3) also makes clear that quasilinear boundary flux must be checked, rather than inferred from the word “Neumann.” [Hou–Yin–Yuan, §§1.1, 4.1](https://arxiv.org/html/2602.14019v1).

For the retained shell, the appropriate implementation remains its own direct normal-recovery equation: isolate the nonzero c^zz D_z²u term, charge the differentiated coefficient and lower-map products, and induct on normal order. Time and lifted SO(4) words preserve its fixed component D/N split. Normal differentiation generally does not.

There are three separate quantities to keep explicit:

- The boundary-preserving commuted energy and every time coefficient commutator.
- The finite normal-recovery recurrence and its source/initial-time traces.
- The conversion from rotation-word control to the retained covariant Sobolev norm.

This is compatible with an exact angular simplification for V1 and an angular-uniform finite-Sobolev estimate for arbitrary evolved data. The latter is still needed for the actual nonlinear family and must retain variable-coefficient rotation commutators.

## Scope, provenance and remaining limit

Fresh focused searches covered S3 tensor harmonics and SO(4) Casimirs, time-first Neumann/quasilinear energy, spatial recovery, and recent 2024–2026 methods. The four sources above are additions to the retained applicability review; no optional cosmology collection was needed to resolve this mathematical question. An older tensor-Casimir paper was retrieved during discovery but is not relied on, since the selected sources give more direct S3 normalizations.

This update provides method evidence and an explicit conditional Casimir argument. It does not certify that the newly evaluated N14 constants close, an actual-family amplitude interval, a common proper-time observation domain, or the finite nonlinear actual-area quadrupole remainder. See [LITERATURE_SOURCES.json](LITERATURE_SOURCES.json) for exact retrieval metadata and coverage.

