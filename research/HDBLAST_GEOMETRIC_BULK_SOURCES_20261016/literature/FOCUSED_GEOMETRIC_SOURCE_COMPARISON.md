# Primary comparison for the geometric quadratic bulk source

Checked 16 October 2026 UTC. The closest general references already supply parameter-derivative formulas for metric curvature and second-order gauge changes. They provide independent algebraic checks for this five-dimensional Einstein–scalar calculation. The remaining model-specific work is to evaluate the registered reduced source in the full angular space and certify its approximation error. None of the papers compared here certifies that calculation or its material-shell response.

## Three primary curvature and gauge references

All three papers are by David Brizuela, José M. Martín-García and Guillermo A. Mena Marugán. The primary arXiv records list only the pinned v1 for each paper; a discovery index's later update date is not an arXiv revision date.

| Primary source | Inspected result and concrete overlap | Boundary of transfer |
| --- | --- | --- |
| [xPert: Computer algebra for metric perturbation theory, arXiv:0807.0824v1](https://arxiv.org/abs/0807.0824v1), 6 July 2008; [Gen. Relativ. Gravit. 41, 2415–2431 (2009)](https://doi.org/10.1007/s10714-009-0773-2) | Section II treats arbitrary dimension and background. Equations (1)–(4) use parameter derivatives with the expansion coefficient `epsilon^n/n!`. Equations (8)–(12) give connection, Riemann and Ricci perturbations. Equations (18)–(19) expand a perturbed covariant derivative over all tensor slots. | Supports the universal Ricci Hessian and an independent symbolic implementation. It supplies no interval enclosure, angular closure, scalar-potential convention, relative-harmonic reduction or shell inverse. The package's reported speed is not an error certificate. |
| [Second and higher-order perturbations of a spherical spacetime, gr-qc/0607025v1](https://arxiv.org/abs/gr-qc/0607025v1), 6 July 2006; [Phys. Rev. D 74, 044039](https://doi.org/10.1103/PhysRevD.74.044039) | Section II, equations (9)–(16), gives the earlier connection/curvature formulas before specializing to a spherical background. Section V.2, equations (66)–(69), retains the second gauge generator, twice the Lie derivative of the first perturbation, and the squared first Lie derivative. | The subsequent evolution equations use a four-dimensional `M2 x S2` background and Regge–Wheeler gauge. Their harmonic coefficients cannot be substituted for a five-dimensional `M2 x S3` relative-harmonic calculation. |
| [High-order gauge-invariant perturbations of a spherical spacetime, gr-qc/0703069v1](https://arxiv.org/abs/gr-qc/0703069v1), 12 March 2007; [Phys. Rev. D 76, 024004](https://doi.org/10.1103/PhysRevD.76.024004) | Section II, equations (8)–(12), gives tensor gauge transformations and their composition. Equation (10) is the complete second-order tensor change. Section V.3 explains how low angular gauge modes obstruct local invariants and enter higher orders through mode coupling. | The `S2` low-mode analysis does not prove the present `S3` closure or invariance of a selected response. It does justify checking the entire generated sector space rather than declaring a coefficient invariant from its coordinate label. |

## Direct checks for the registered source

Let `h` denote the first metric perturbation, with all contractions and derivatives taken in the background. Set

```text
P^C_AB = 1/2 g^CD (nabla_A h_BD + nabla_B h_AD - nabla_D h_AB),
J^C_AB = -2 h^C_D P^D_AB.
```

On the path `g(epsilon)=g+epsilon h`, the second metric parameter derivative vanishes. The cited connection formula then gives `delta Gamma=P` and `delta^2 Gamma=J`. With the registered Ricci convention, its contraction gives

```text
Q_Ric_AB = nabla_C J^C_AB - nabla_B J^C_AC
           + 2 (P^C_CD P^D_AB - P^C_BD P^D_AC).
```

Thus the factor two in the product is a parameter-derivative factor. The coefficient of `epsilon^2` is half this Hessian. A path with a nonzero second metric jet also has its linearized Ricci contribution; that term is separate from this quadratic source.

The covariant derivative of `nabla_C h_AB` acts on its derivative index `C` as well as `A,B`. Likewise `nabla J` acts on every index of the `(1,2)` tensor `J`. Treating frame components or scalar amplitudes as scalars through these operations omits connection terms. The xPert derivative formulas make this an independent check of the derivative-slot implementation, in addition to the explicit curvature formula.

The tensor gauge law is

```text
delta^2 T_new - delta^2 T_old
  = L_xi2 T + 2 L_xi1(delta T_old) + L_xi1^2 T.
```

For a single fixed flow the second independent flow generator is zero, while `L_xi1^2 T` remains. Apply the chosen convention to both metric and scalar. These laws concern geometric tensors. They do not replace the direct calculation of the relative-harmonic source `Q_C` or the registered reduction `Q_Ered = Q_E - nabla_(A Q_C_B)` with one-half symmetrization.

## Shell geometry and recent numerical methods

Marc Mars, [First and Second Order Perturbations of Hypersurfaces, gr-qc/0507005v1](https://arxiv.org/abs/gr-qc/0507005v1), 1 July 2005; [Class. Quantum Grav. 22, 3325–3348](https://doi.org/10.1088/0264-9381/22/16/013), remains the closest geometric shell comparison. Proposition 3, equations (40)–(41), gives second variations of the induced metric and extrinsic curvature. Section 7 separates ambient spacetime gauge from intrinsic hypersurface reparameterization and uses a covariant second displacement convention. Translate that convention before comparing a fixed-flow shell embedding. Its matching theorem equates the first and second fundamental forms; a shell with Israel stress instead has a prescribed extrinsic-curvature jump. The geometric variation formulas can check each side, but the theorem does not supply the stressed matching law or a hyperbolic inverse.

A useful application of the established gauge law is the intrinsic shell scalar `I=gamma^-1(d phi_shell,d phi_shell)`. Here `gamma` is the induced metric and `phi_shell` is the scalar pulled back by the shell embedding. If the background shell scalar is constant, then `I0=I1=0`. Writing `psi=delta phi_shell`, direct differentiation gives `I''=2 gamma0^-1(d psi,d psi)`. The displayed second gauge law makes `I''` invariant under shell reparameterization, because every added term contains `I0` or `I1`; consistent simultaneous transformation of the ambient fields and embedding preserves the pullback. This is the familiar vanishing-lower-order gauge criterion applied to a geometric shell diagnostic. It is a kinematic consequence, not a novelty claim or a theorem for the second response. The first pullback perturbation `psi` must include any embedding contribution.

Two current reconstruction papers remain useful for the numerical stage: Zhaonan Dong, Lorenzo Mascotto and Zuodong Wang, [A priori and a posteriori error estimates of a C0-in-time method for the wave equation in second order formulation, arXiv:2411.03264v3](https://arxiv.org/abs/2411.03264v3), revised 3 May 2026; and Zhaonan Dong, Emmanuil H. Georgoulis, Lorenzo Mascotto and Zuodong Wang, [A posteriori error analysis and adaptivity of a space-time finite element method for the wave equation in second order formulation, arXiv:2509.08537v2](https://arxiv.org/abs/2509.08537v2), revised 23 April 2026. The former's equation (44) and Proposition 3.2 repair temporal derivative jumps with a continuous reconstruction. The latter's sections 4.1–4.4 and Theorem 5.3 keep temporal, spatial and mesh-transfer defects separate. Their scalar Dirichlet-wave displacement estimates in `L-infinity(time;L2(space))` do not certify the required first-field scaled-C3 approximation errors, pointwise shell traces or the registered response norm. Use the reconstruction and accounting ideas alongside the actual differentiated integral equations and continuous residual estimate.

Zhongshan An and Michael T. Anderson, [Well-posed geometric boundary data in General Relativity, II: twisted Dirichlet boundary data, arXiv:2507.15567v3](https://arxiv.org/abs/2507.15567v3), revised 30 May 2026, remains a useful harmonic-gauge corner/subsidiary comparison, with different vacuum boundary data. Their current [Part III, arXiv:2503.12599v3](https://arxiv.org/abs/2503.12599v3), revised 18 June 2026, states a weakened injective/dense-range result after a gap and reports failures of general energy estimates for its conformal-mean-curvature data. Neither paper supplies a numerical inverse for the scalar-dependent material shell.

This is a bounded primary comparison of curvature/gauge, shell-geometry and numerical-method sources. It is not an exhaustive literature search, a complete proof audit, or a novelty or priority assessment. No theorem for the actual coupled material-shell response is imported.
