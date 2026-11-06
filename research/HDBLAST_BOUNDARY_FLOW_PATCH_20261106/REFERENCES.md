# Sources and method comparison

The model sources below are linked at immutable repository commit
`8fe1c4d99ac8a1c4a4a2ea4ec3cdf0e1dd12053c`. The cited bytes agree with the
retained source members. These are mathematical source checkpoints with
their stated hypotheses; citing them does not imply external peer review.

## Retained model sources

| Alias | Public source | Use here |
| --- | --- | --- |
| S1 | [Exact finite-regularity shell evolution](https://github.com/maldonado-research/HDblast/blob/8fe1c4d99ac8a1c4a4a2ea4ec3cdf0e1dd12053c/research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md), 15 October 2026 | Conformal PDE coordinate, relative gauge and conditional prescribed-slab theorem |
| S2 | [Exact nonlinear boundary straightening](https://github.com/maldonado-research/HDblast/blob/8fe1c4d99ac8a1c4a4a2ea4ec3cdf0e1dd12053c/research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/inputs/NONLINEAR_BOUNDARY_STRAIGHTENING.md), equations (1)–(5) | Proper-`y` lapse, full-collar background subtraction and the literal autonomous perturbation generator |
| S3 | [Initial harmonic completion and subsidiary system](https://github.com/maldonado-research/HDblast/blob/8fe1c4d99ac8a1c4a4a2ea4ec3cdf0e1dd12053c/research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md), §§1–4 | ADM time lapse, reference connection, shell-normal convention and matter Codazzi cancellation |
| S4 | [Registered exact-family binding](https://github.com/maldonado-research/HDblast/blob/8fe1c4d99ac8a1c4a4a2ea4ec3cdf0e1dd12053c/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/inputs/REGISTERED_EXACT_FAMILY_BINDING.md), §§1,3–4 | Actual `delta=1/1000` member, fixed background-labelled chart and coordinate-time qualification |
| S5 | [Actual potential and background collar bounds](https://github.com/maldonado-research/HDblast/blob/8fe1c4d99ac8a1c4a4a2ea4ec3cdf0e1dd12053c/research/HDBLAST_SECOND_INITIAL_NORMS_GAUGE_20261015/background/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json) | Conditional wider-collar enclosures, including `rho in [129,131]` and `phi0-1 in [-9/100000,-1/12500]` |

## Primary method references

1. Tatsien Li, Xing Lu and Bopeng Rao, [*Exact boundary controllability and
   exact boundary synchronization for a coupled system of wave equations
   with coupled Robin boundary controls*, arXiv:2202.00760v1](https://arxiv.org/html/2202.00760v1),
   1 February 2022; *ESAIM: COCV* 27 (2021), S7,
   [DOI:10.1051/cocv/2020047](https://doi.org/10.1051/cocv/2020047).
   Theorem 4.1, equations (4.3)–(4.7), uses an exponential normal
   conjugation to replace scalar Robin data by Neumann data and first-order
   bulk terms. Their coupled weak-solution/control results assume a boundary
   matrix similar to a real symmetric matrix. They do not supply an
   arbitrary nonsymmetric nonlinear twelve-row flow theorem or the present
   finite-order constants.

2. Grigorios Fournodavlos and Jacques Smulevici, [*The initial boundary value
   problem in General Relativity: the umbilic case*, arXiv:2104.08851v1](https://arxiv.org/html/2104.08851v1),
   18 April 2021. Theorem 1 and §2.2 provide the wave-gauge,
   boundary-orthogonal umbilic strategy; Appendix A gives scalar Robin energy
   bookkeeping. The theorem concerns smooth four-dimensional vacuum and
   constant umbilicity. Its boundary dimension, second-fundamental-form
   convention and absolute harmonic normalization differ from this
   five-dimensional Einstein–scalar model. Here the scalar junction cancels
   the variable-tension term in Codazzi; the vacuum identity cannot simply
   be copied.

3. Zhongshan An and Michael T. Anderson, [*Well-posed geometric boundary
   data in General Relativity, II: twisted Dirichlet boundary data*,
   arXiv:2507.15567v3](https://arxiv.org/html/2507.15567v3),
   revised 30 May 2026. Theorem 1.1 treats vacuum in arbitrary spatial
   dimension at least two through smooth tame/Nash–Moser theory. Its
   conformal-class and twisted-volume boundary data contain no metric
   derivatives. They differ from the Israel/scalar Robin rows here and do
   not provide this flow patch or a prescribed proper-time domain.

4. Dionysios Anninos, Raúl Arias and Damián A. Galante,
   [*Umbilic Observatories*, arXiv:2610.40783v1](https://arxiv.org/html/2610.40783v1),
   29 October 2026. Its linear analysis concerns a three-dimensional de
   Sitter boundary in four-dimensional Minkowski gravity with conformal or
   Dirichlet data. Section 3.2 discusses selected nonlinear boundary Weyl
   solutions. It does not establish the general nonlinear variable-tension
   Einstein–scalar shell problem.

These exact primary versions were read back on 6 November 2026. The useful
method comparison is Robin-to-Neumann conjugation followed by accounting
for all induced bulk terms, energy margins and subsidiary identities.
No cited external theorem literally supplies the present model-specific
construction. This scope distinction supports neither a novelty claim nor
an ill-posedness claim.
