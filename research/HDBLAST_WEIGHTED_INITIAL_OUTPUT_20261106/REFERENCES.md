# Primary-source comparison and research scope

Sources checked for this 6 November 2026 working checkpoint. This is a
focused comparison; absence of an identical formula in these sources is
not a proof of priority. Internal review is not external peer review.

- **E. Delay, Smooth compactly supported solutions of some underdetermined
  elliptic PDE, with gluing applications**,
  [arXiv:1003.0535v4](https://arxiv.org/html/1003.0535v4).
  Theorem 3.8 is the weighted projected isomorphism used by the retained
  construction, assuming the asymptotic Poincare inequality. Its proof
  is variational; the minimum-energy comparison here is an application
  of that established framework. Proposition 3.6, section 8, sections
  9.3/9.5 and Appendix 10 distinguish coercivity, adjoint hypotheses,
  symmetric divergence, scalar curvature and scaled regularity. The
  matter-dependent zeroth-order term must be retained in our metric
  adjoint. The momentum block requires the Killing/Korn argument.
  These results do not supply our numerical high-order constant.

- **P. T. Chrusciel and E. Delay, On mapping properties of the general
  relativistic constraints operator in weighted function spaces, with
  applications**, [author manuscript](https://homepage.univie.ac.at/piotr.chrusciel/papers/cs/grqc.pdf).
  Section 5 and the appendices develop exponential boundary weights and
  adjoint estimates. Their boundary-collar estimates explain the weight
  machinery, but do not evaluate a coercivity constant for this entire
  four-dimensional Einstein–scalar annulus. Vacuum source normalization
  cannot replace the canonical matter calculation in this checkpoint.

- **Y. Mao, S.-J. Oh and Z. Tao, Initial data gluing in the asymptotically
  flat regime via solution operators with prescribed support properties**,
  [arXiv:2308.13031v1](https://arxiv.org/html/2308.13031v1).
  Explicit support-controlled constraint inverses and quantitative mapping
  estimates have substantial precedent. This work treats flat
  three-dimensional vacuum linearization, which differs from our curved
  scalar-frozen four-dimensional metric operator and selected minimizer.

- **P. T. Chrusciel, A. Cogo and F. Nutzi, A Bogovskii-type operator for
  Corvino–Schoen hyperbolic gluing**,
  [arXiv:2409.07502v1](https://arxiv.org/html/2409.07502v1).
  Proposition 1.1 gives a compact-support scalar-curvature right inverse
  modulo a finite-rank static obstruction, gaining two Sobolev
  derivatives on a bounded hyperbolic domain. It does not identify our
  weighted minimizer or its scalar-gradient correction. Our monopole
  particular uses a source derivative and p^-2; it does not claim that
  optimal regularity property.

- **P. Isett, Y. Mao, S.-J. Oh and Z. Tao, Integral formulas for
  under/overdetermined differential operators via recovery on curves and
  the finite-dimensional cokernel condition I: General theory**,
  [arXiv:2509.04617v1](https://arxiv.org/html/2509.04617v1).
  Recovery identities for adjoints yielding support-controlled right
  inverses by duality are established general methods. Theorems 1.2,
  1.11 and 1.14 impose curve/domain and symbol conditions and provide
  regularity results beyond the specialized particular used here. A
  separate announced companion theorem is not treated as an inspected
  proof or as a numerical constant for our annulus.

- **B. Le Floch and P. G. LeFloch, Optimal localization for the Einstein
  constraints**, [arXiv:2312.17706v2](https://arxiv.org/pdf/2312.17706v2),
  revised in 2026. The variational seed-to-solution projection and
  weighted fourth-order Hamiltonian estimates are relevant comparisons.
  Their asymptotically Euclidean vacuum conical setting, decay and
  stability hypotheses differ from this compact Einstein–scalar annulus.
  No evaluated H14 constant transfers from that theory.

- **F. Nutzi, Solutions of the constraints with controlled decay to Kerr,
  including Schwartz decay**,
  [arXiv:2603.23484v1](https://arxiv.org/html/2603.23484v1).
  This 2026 work provides a current comparison for quadratic nonlinear
  completion and support-sensitive optimal weighted mapping in vacuum
  on R3. It does not evaluate the present scalar-affine weighted output.

- **J. Isenberg and D. Maxwell, A Phase Space Approach to the Conformal
  Construction of Non-Vacuum Initial Data Sets in General Relativity**,
  [arXiv:2106.15027v2](https://arxiv.org/html/2106.15027v2).
  Canonical matter fields and conjugate densities are established
  phase-space data. This motivates preserving the distinction between
  affine scalar density and physical normal scalar momentum. The
  conformal construction is a different method from the localized inverse.

The contribution claimed here is the specialized actual-source
calculation, explicit modal comparison tensors, same-selected-output
energy bounds and evaluated background-dependent L2 ceilings. The
angular inversion also has an internal antecedent in the pinned
[localized constraint manuscript](inputs/CANONICAL_CONSTRAINTS.md).
Neither weighted minimization nor the broad dual-inverse principle is
claimed as new.

For k=12, the potential regularity estimate requires its L2_w norm;
the comparison controls qweight² P_g* N in L2_w. The positive
nonconstant-mode energy term controls N in L2_mu_H, where
mu_H=w²qweight⁴, and does not automatically control N in L2_w. A
quantified coercivity and scaled regularity/extension estimate, or a
direct physical-tensor estimate, is still needed for H14/H13. The
small scalar gradient makes monopole conditioning a particular concern.
Exact kernel absence alone does not quantify it. After that come the
coupled second evolution, physical pullbacks and common-domain bounds.
