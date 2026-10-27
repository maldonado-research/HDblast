# Focused comparison: smooth reconstruction and residual bounds

This is a bounded primary-source comparison, not an exhaustive literature
search or a novelty determination. The mathematical methods discussed below
are established; none establishes the HDBLAST second-order bulk and shell
solution or certifies its residual.

1. T. Chaumont-Frelet and A. Ern, *Damped energy-norm a posteriori error
   estimates using C²-reconstructions for the fully discrete wave equation
   with the leapfrog scheme*, ESAIM: Mathematical Modelling and Numerical
   Analysis 59 (2025), 1937–1972, DOI
   [10.1051/m2an/2025027](https://doi.org/10.1051/m2an/2025027),
   [arXiv full text](https://arxiv.org/html/2403.12954v2), especially §§1, 3,
   and Theorem 4.3. The authors use a globally C² piecewise-quartic time
   reconstruction of a leapfrog finite-element approximation and derive
   damped energy-norm a posteriori estimates with residual, reconstruction,
   and source-regularity terms. The scalar model assumes homogeneous Dirichlet
   boundary conditions and zero initial data; its hypotheses and norms differ
   from the mixed-boundary coupled HDBLAST formulation. Its useful methodological
   lesson is that reconstruction smoothness and its error/defect must both be
   accounted for; the theorem supplies no pointwise HDBLAST C² bound.

2. L. R. Scott and S. Zhang, *Finite element interpolation of nonsmooth
   functions satisfying boundary conditions*, Math. Comp. 54 (1990), 483–493,
   [original primary paper](https://www.ams.org/journals/mcom/1990-54-190/S0025-5718-1990-1011446-7/S0025-5718-1990-1011446-7.pdf),
   DOI [10.1090/S0025-5718-1990-1011446-7](https://doi.org/10.1090/S0025-5718-1990-1011446-7).
   The construction chooses local averaging appropriately to preserve
   discrete boundary values while approximating rough functions. Its finite
   element interpolant is a different construction; it does not prove that
   arbitrary smoothing preserves the current mixed shell conditions or gives
   the present approximation C² regularity. It illustrates why boundary and
   trace compatibility must be designed and checked.

3. T. Führer, G. Gantner, R. González and M. Karkulik, *Space-time finite
   element interpolation of nonsmooth functions satisfying boundary
   conditions*, [arXiv:2609.11496v1](https://arxiv.org/html/2609.11496v1)
   (2026), §§1–3. The authors develop Scott–Zhang-type space-time projections
   that preserve discrete Dirichlet boundary data in parabolic energy spaces.
   This is a recent methodological comparison; its parabolic setting and
   Dirichlet conditions do not cover the present hyperbolic mixed Robin problem.

References were read at the linked full primary texts. No theorem from these
papers is imported as a certificate for the 23-mode HDBLAST approximation.
