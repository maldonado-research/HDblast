# Independent review of the reflection-invariant formal continuation

14 October 2026. Internal AI-assisted review, not external peer review.

Reviewed: `/workspace/hdblast-research-work/referee-20261014/FORMAL_REFLECTION_INVARIANT_CONSTRAINTS.md`, as available before final support-wording correction.

## Verdict

The proposed finite-order initial-constraint result follows from the corrected localized linear solvability theorem and the stated smooth radial-background hypotheses. No substantive mathematical gap was found in the symmetry, covariance, combinatorics or finite-order scope. One support statement should be sharpened as below. This review does not certify newness in the mathematical literature.

## Checks

- The Einstein–scalar constraint map is a natural nonlinear differential operator. Treating `π` and `ϖ` as absolute densities makes orientation-reversing sign reflections legitimate; positive integration densities introduce no spurious determinant sign.
- The four boost lapses carry the characters `σ_a`; the six rotations carry `σ_aσ_b`. Each character is nontrivial under the full sixteen-element independent-reflection group. Thus its invariant KID subspace is zero. Equivariant higher-order scalar/covector-density sources are orthogonal to all ten KIDs by change of variables.
- The locality of the nonlinear differential operator ensures each higher source is supported in the finite union of the previously chosen supports. Constant/background terms do not generate an extra exterior source because the background solves the constraints.
- Averaging a localized solution is legitimate because the linear operator commutes with the finite group and a concentric annulus is group invariant. No averaging of a nonlinear solution is used.
- The ordered-tuple factorial formula has the correct second, third and fourth derivatives and involves only lower coefficients. Smooth potential and metric positivity suffice for every finite Taylor expansion; analyticity is not required.
- Simultaneously diagonalizing real commuting symmetric `A,B` supplies the necessary full independent-reflection invariance. Antipodal symmetry by itself fails to eliminate rotational KIDs.
- Opposite-copy first-order parity requires `u_j^−=(−1)^j u_j^+`; it does not impose odd copy parity on every coefficient or on a finite-`ε` nonlinear path.
- For each fixed finite order, a small enough parameter interval keeps the metric positive and gives the Taylor remainder in a chosen spatial norm after controlling the required spatial/potential derivatives. That interval and norm bounds may depend on the truncation order.

## Support clarification

At the induction step, apply the localized theorem **on** `A_n`, not merely on the original `A`, to obtain

`supp z_n ⊂ closure(A_n)`.

Choose

`closure(A_{n−1}) ⊂ A_n ⊂ closure(A_n) ⊂ A_infinity ⋐ A`.

The source is then compactly supported in `A_n`; averaging preserves the same support enclosure. This makes the common fixed collar away from shell and vertex explicit. Because each operator theorem gives infinite-order boundary vanishing, inclusion in the closed annulus is enough for smooth zero extension.

## What is not implied

A compatible sequence of smooth coefficients is a formal initial-data jet. Neither their growth nor a common convergent `ε`-series is controlled by this induction. Borel-type realization, if subsequently used, would at most make the residual flat in `ε`; flatness is not vanishing. Exact nonlinear initial data need a nonlinear inverse or gluing theorem with suitable function-space and support estimates. Shell-contact evolution additionally needs the hyperbolic and material-interface problem.

## Citation context and a directly relevant comparison

**Judith M. Arms, Jerrold E. Marsden and Vincent Moncrief, “Symmetry and bifurcations of momentum mappings,” Communications in Mathematical Physics 78 (1981), 455–478**, [DOI:10.1007/BF02046759](https://doi.org/10.1007/BF02046759). The [author's archived abstract](https://www.cds.caltech.edu/~marsden/bib/1981/01-ArMaMo1981/) was checked. It describes singularities of momentum-map zero sets at points with symmetry, with local quadratic models, and an application to Einstein solutions with spacelike Killing fields and compact Cauchy surfaces. This establishes the longstanding symmetry/constraint framework. It does **not** by itself prove this annular scalar-background, finite-reflection support theorem. The present mechanism should be described as standard equivariance and finite-group averaging applied to the model-specific local solver.

**Peter Hintz, “The linearized Einstein equations with sources,”** [arXiv:2306.07715](https://arxiv.org/abs/2306.07715), 13 June 2023. Canonical metadata and primary in-body passages were checked. On smooth connected globally hyperbolic vacuum spacetimes of dimension `n+1`, `n≥1`, Theorem 1.1 proves that a smooth spatially compact divergence-free source admits a spatially compact metric perturbation exactly when all its Killing charges vanish. The proof first solves support-controlled linearized constraints and then solves the gauge-fixed hyperbolic Cauchy problem. Remark 1.2 extends this to a cosmological constant. This is a particularly close precedent for combining charge conditions, localized constraints and bulk evolution. It does not cover nonconstant bulk scalar matter or a finite timelike material shell.

**De Vuyst–Eccles–Hoehn–Kirklin,** [arXiv:2411.19931](https://arxiv.org/abs/2411.19931), supplies current boundary-charge context and explicitly limits integration of all higher perturbative consistency conditions to cases where analyticity and convergence are assumed. Its context does not replace the model-specific estimates.

The focused contextual search did not establish whether an identical abstract finite-reflection proposition has already been published. No external novelty claim is warranted.

