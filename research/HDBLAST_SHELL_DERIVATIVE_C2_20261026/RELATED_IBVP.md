# Focused literature: bulk, boundary, regularity, constraints

This is a bounded primary-source comparison, not a complete literature search
or novelty determination. The referenced models are distinct from the HDBLAST
Einstein–scalar problem and cannot by themselves certify its full boundary
system or numerical error.

- C. Deffayet, “A note on the wellposedness of scalar brane world cosmological
  perturbations,” [arXiv:hep-th/0409302](https://arxiv.org/html/hep-th/0409302),
  §§2–3. The 2005 paper studies linear scalar perturbations with a maximally
  symmetric 5D bulk and specified brane matter. It differentiates the master
  PDE to obtain a first-order characteristic system and analyzes initial data
  and compatibility with the boundary. Higher derivatives of the enlarged
  system are constrained by the original PDE, not independently selectable;
  inconsistent corner data can produce derivative discontinuities (§3.2).
  This is a useful comparison for differentiating HDBLAST's coupled modal
  system, not a direct well-posedness theorem for its bulk scalar field,
  second-order source, or mixed shell system.

- H.-O. Kreiss, O. Reula, O. Sarbach and J. Winicour, “Boundary conditions for
  coupled quasilinear wave equations with application to isolated systems,”
  [arXiv:0807.3207](https://arxiv.org/html/0807.3207), §I.1, theorem 1, and
  §III. The authors prove well-posedness for a class of systems using a
  symmetric-hyperbolic reduction and maximally dissipative boundary conditions,
  with smooth compatible initial/boundary data and their stated coefficient
  conditions. Harmonic Einstein equations are an application. Their particular
  boundary hypotheses must be checked and cannot be assumed for HDBLAST's
  Dirichlet/Robin shell; a continuous-dependence theorem is also not itself a
  numerical C² enclosure.

- S. Frittelli and R. Gómez, “Initial-boundary-value problem of the
  self-gravitating scalar field in the Bondi-Sachs gauge,”
  [Physical Review D 75, 044021 (2007)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.75.044021).
  The primary abstract describes a spherically symmetric Einstein–scalar
  symmetric-hyperbolic system and a constraint-preserving boundary condition.
  It concerns 4D spherical symmetry and does not prove constraint propagation
  for the present 5D quadrupole shell; the paper's abstract was inspected.

- J. Martin, G. Felder, A. Frolov, M. Peloso and L. Kofman, “Braneworld
  dynamics with the BraneCode,” [arXiv:hep-th/0309001](https://arxiv.org/html/hep-th/0309001),
  §§2–3 and Appendix A. This earlier work performs numerical nonlinear 5D
  evolution with a bulk scalar and two orbifold branes and explicitly enforces
  initial constraints and boundary conditions. Its model assumes homogeneous,
  isotropic brane directions and does not solve the present inhomogeneous
  quadrupole. Its analytic comparison cases suggest useful solver benchmarks,
  not evidence that the HDBLAST cosmological hypothesis is correct.

The next model-specific step remains a quantitative derivative ledger for the
actual bulk source, affine shell source and feedback, including initial-corner
compatibility and constraints. These cited results do not supply that ledger.
