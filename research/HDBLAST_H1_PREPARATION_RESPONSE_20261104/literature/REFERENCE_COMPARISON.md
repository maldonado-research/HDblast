# Preparation selection: focused primary-source comparison

4 November 2026. This is a focused literature audit, not an exhaustive search or a priority claim. Existing research capsules and papers are unchanged. Full PDFs were retrieved into this private research directory; their exact transport hashes are recorded in `PRIMARY_SOURCE_PINS.json`.

## Direct overlap: first-tangent canonical energy

S. Hollands and R. M. Wald, *Stability of Black Holes and Black Branes*, Commun. Math. Phys. **321** (2013) 629–680, DOI [10.1007/s00220-012-1638-1](https://doi.org/10.1007/s00220-012-1638-1), [arXiv:1201.0463v4](https://arxiv.org/abs/1201.0463v4).

Section 2.2, equations (19)–(24), constructs the symplectic current and integrated product from the background and two first perturbations. Section 2.3, page 17, equation (46), defines

`E = W_Sigma(g; gamma, Lie_t gamma)`.

Equation (48) treats this as a bilinear form in first perturbations. The paragraph after equation (47) explicitly distinguishes its first-tangent dependence from second variations of mass, angular momentum and horizon area. It already observes that the difference of two second perturbations is a homogeneous first-order solution and cancels from their first-law combination. Thus shared-first canonical-energy blindness is not a new general principle.

The paper's stability conclusions concern stationary, asymptotically flat vacuum black holes with bifurcate Killing horizons, restricted perturbations, horizon gauge conditions and asymptotic falloff. Compactly supported gauge degeneracy is discussed before Lemma 2 on page 16; allowed noncompact/asymptotic generators require further conditions. Conservation uses the background Killing field and appropriate symplectic boundary conditions. These hypotheses cannot be imported into the finite Einstein–scalar shell by renaming a coordinate or proper clock as `t`.

**Safe consequence for the current model.** A preparation criterion that is genuinely a function only of the fixed first tangent assigns the same value to the base and the admissible opposite-sign preparation. A threshold on that value either admits both or rejects both. This uses the retained HDBLAST witness to diagnose a candidate selection rule, rather than extending the black-hole stability theorem. It neither computes a positive shell canonical energy nor selects a physical initial state.

## General Noether charges are local, but Hamiltonians require conditions

V. Iyer and R. M. Wald, *Some Properties of Noether Charge and a Proposal for Dynamical Black Hole Entropy*, Phys. Rev. D **50** (1994) 846–864, DOI [10.1103/PhysRevD.50.846](https://doi.org/10.1103/PhysRevD.50.846), [arXiv:gr-qc/9403028v1](https://arxiv.org/abs/gr-qc/9403028v1).

The paper covers diffeomorphism-invariant theories with arbitrary tensor matter. Equation (20) defines the symplectic potential via `delta L = E delta Phi + dTheta`; equations (49)–(51) give the Noether current and its on-shell local charge. Proposition 4.1, equation (52), describes the charge as a local expression in fields, generator and finitely many derivatives, with explicitly discussed ambiguities. Equations (78)–(82), pages 14–15, show that a Hamiltonian, **if it exists**, has an on-shell surface representation; equation (80) states the boundary condition needed for the displayed integration. The subsequent asymptotically flat energy definition is not automatically a finite-shell energy.

**Safe consequence for the current model.** If a specified initial charge is a smooth functional of the same finite boundary jets, boundary cut/anchor, generator and reference prescription, agreement of those data through parameter order two gives equal charge coefficients through order two by the chain rule. Compact support of the added initial datum supplies the collar agreement. This is a locality/Taylor statement, not a claim that every physical energy has been constructed or is blind to all orders. Future charges after the perturbation reaches a boundary, field-dependent interior normalizations and unequal higher-order data are outside it.

## Finite-boundary completion is part of the theory

D. Harlow and J.-q. Wu, *Covariant phase space with boundaries*, JHEP **10** (2020) 146, DOI [10.1007/JHEP10(2020)146](https://doi.org/10.1007/JHEP10(2020)146), [arXiv:1906.08616v4](https://arxiv.org/abs/1906.08616v4).

The introduction, pages 3–4, stresses that spatial boundary conditions define the Hamiltonian system; Dirichlet and Neumann scalar cavities are different systems. The sufficient variational condition is `(Theta+delta ell)|Gamma=dC`, equation (2.20). Equations (2.23)–(2.27) include the boundary correction in the symplectic form and use the resulting vanishing boundary flux to prove slice independence. Equations (2.48)–(2.51) give the Hamiltonian and its boundary form for suitable boundary-preserving diffeomorphisms. The correction involving `C` must not be silently discarded.

Section 3.2 gives a two-derivative scalar example; section 3.5 treats general relativity with its boundary action, equations (3.35)–(3.46), and recovers the Brown–York generator, equations (3.61)–(3.62), for the stated boundary-isometry setting. Combining their separate examples is not a derivation of the present dynamical scalar-tension shell's phase space. Its shell action, admissible variations, embedding/clock treatment, symplectic flux and generator still need to be specified and checked.

## Standard wave energy estimates and the proposed weaker-norm route

H.-O. Kreiss, O. Reula, O. Sarbach and J. Winicour, *Well-posed initial-boundary value problem for the harmonic Einstein equations using energy estimates*, Class. Quantum Grav. **24** (2007) 5973–5984, DOI [10.1088/0264-9381/24/23/017](https://doi.org/10.1088/0264-9381/24/23/017), [arXiv:0707.4188v2](https://arxiv.org/abs/0707.4188v2).

Section 2, Lemma 1, equation (4), constructs a positive energy with mixed derivative terms. Theorem 1, equations (5)–(7), derives an estimate with boundary dissipation by integration by parts. Appendix 1, Lemma 2, equations (18)–(19), explains stability under lower-order perturbations and adds a boundary value norm. Appendix 2, pages 14–16, chooses an outward-tilted timelike multiplier `u^a=T^a+delta N^a`, obtains a coercive energy, and absorbs mixed boundary terms by a Young-type inequality, equations (22)–(24). Thus choosing a modified energy to expose usable boundary flux is an established method.

The geometric extension is H.-O. Kreiss, O. Reula, O. Sarbach and J. Winicour, *Boundary conditions for coupled quasilinear wave equations with application to isolated systems*, Commun. Math. Phys. **289** (2009) 1099–1129, DOI [10.1007/s00220-009-0788-2](https://doi.org/10.1007/s00220-009-0788-2), [arXiv:0807.3207v1](https://arxiv.org/abs/0807.3207v1). Its introduction expressly describes a nonstandard energy based on a timelike direction adapted to the boundary conditions. Section II explains that freedom in the evolution/multiplier direction permits boundary as well as bulk control. Its general theorem uses a specified class of boundary derivative couplings, with a smallness/triangular structure condition, and suitable data regularity and compatibility.

**Comparison with the proposed HDBLAST estimate.** The tilted bulk quadratic form, positive boundary value term, flux absorption and continuous Gronwall argument are applications of standard multiplier/energy methods. A contribution here would be the explicit estimate for the retained five-channel operator and its conversion from a weaker initial profile norm to the two registered proper shell observables. The particular positive Gram matrix, nonsymmetric Robin block, Dirichlet row, shrinking characteristic cone, coefficient caps and numerical target bounds must be checked directly; these papers do not supply those constants or that model-specific estimate.

In particular, the 2007 theorem's displayed boundary condition has a strictly positive time-derivative coefficient, and the 2008/2009 general theorem has its own admissible boundary-coupling conditions. Neither theorem should simply be cited as proving this mixed Dirichlet/Robin problem. Estimating a mixed term by `|<V,MU>_H| <= ||M||_H |V|_H |U|_H` does not require `M` to be self-adjoint, but whether it is absorbed with the claimed constants is part of the current direct proof. Replacing that term by the derivative of a Robin quadratic potential would need additional symmetry assumptions and must not be done silently.

This mathematical energy norm is not thereby a positive, conserved physical Noether energy of the Einstein–scalar shell or a physical preparation law. A successful wave estimate can establish a sufficient class-to-observable implication while leaving the physical action/state selecting that class unresolved. The profile norm, its normalization and the retained exact-family and linear-evolution premises must remain explicit.

For bibliographic clarity, [arXiv:gr-qc/0602051v3](https://arxiv.org/abs/gr-qc/0602051v3) is the earlier **two-author** Kreiss–Winicour paper, *Problems which are well-posed in a generalized sense with applications to the Einstein equations*, Class. Quantum Grav. **23** (2006) S405–S420, DOI 10.1088/0264-9381/23/16/S07. It uses a pseudodifferential approach. It is useful background, but is not the four-author integration-by-parts energy paper.

## Recent adjacent work screened

1. A. Bhattacharya and O. Parrikar, *Covariant phase space and the semi-classical Einstein equation*, [arXiv:2510.19939v2](https://arxiv.org/abs/2510.19939v2), revised 14 April 2026. Section 2.5 reviews the classical identity with matter; the abstract and sections 3–4 propose the quantum-state Berry curvature contribution and a subregion construction. No quantum state, renormalized stress tensor or Berry-curvature data are supplied by the current classical HDBLAST preparation. This paper suggests a separate possible state-based research program; it does not select the current second preparation or certify its radius.
2. D. Anninos et al., *The Disordered Universe 2025: Features in a Lambda>0 Universe*, [arXiv:2610.37117v1](https://arxiv.org/abs/2610.37117v1), 27 October 2026. Sections 5.1–5.4 discuss finite timelike boundaries and their boundary-condition dependence. Page 36, Question 7, asks for a portion of pure de Sitter space with a timelike boundary that is thermodynamically and dynamically stable at linear order, or a proof none exists. This is posed as an open problem, not a stability theorem for a scalar shell. The conformal/Dirichlet boundary settings in the discussion differ from the current shell matching problem.

The remaining search hits (a gravitational-second-law paper, parallel-gravity surface charges and canonical quantization of spherical Einstein–scalar spacetimes) were screened at search-excerpt level only and were not used as premises. No comprehensive 2024–2026 literature coverage is claimed.

## Essential boundary of the proposed result

Do not identify a quadratic first-tangent statistic with the unrestricted second derivative of an arbitrary Hamiltonian. In Hollands–Wald equation (44), the background Killing condition is used explicitly. Without such a symmetry, differentiating `delta H=Omega(delta Phi,Lie_xi Phi)` can retain terms involving the second datum and `Lie_xi Phi0`. Therefore the present safe claim concerns statistics known to depend only on the first tangent and, separately, specified local initial-boundary jet functionals. It is not a no-go theorem for every energy, every state condition or every nonlinear preparation mechanism.

The useful next calculation is a complete variational/symplectic specification for one proposed physical preparation rule, including its action and shell boundary terms. A profile Sobolev bound or a chosen numerical threshold alone is not that physical specification.
