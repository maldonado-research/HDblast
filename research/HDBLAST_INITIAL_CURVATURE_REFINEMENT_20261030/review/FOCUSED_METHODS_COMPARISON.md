# Focused methods and timelike-boundary literature check

Checked 30 October 2026. Two bounded searches were followed by primary-source
reading. This is not an exhaustive literature search or an originality review.
The useful numerical guidance is established methodology: preserve shared
dependence before taking ranges, and subdivide a low-dimensional domain where
coarse intervals dominate. No inspected reference supplies the missing physical
constraints or nonlinear shell solution for HDBLAST.

## Interval dependency and the current initial-data calculation

M. Neher, K. R. Jackson and N. S. Nedialkov, *On Taylor Model Based Integration
of ODEs*, author-hosted preprint dated 18 August 2006, 21 pages.
[Primary preprint](http://www.cs.toronto.edu/pub/reports/na/Taylor_Models_ODES_2006.pdf).
[Publisher DOI](https://doi.org/10.1137/050638448) was identified by search;
this check read the author preprint, not the publisher's version.

Sections 2.2 and 2.3 distinguish interval dependency from wrapping and explain
the use of a polynomial plus an interval remainder. Their definition encloses
the function over a specified domain; symbolic polynomial terms retain common
dependence while remainder terms are enclosed. Example 2.3 shows that blindly
composing Taylor models can fail when the outer model's domain does not contain
the inner model's entire range. The introduction also explains that subdivision
of initial or parameter intervals is practical primarily in low dimensions.
Later sections analyze naive and preconditioned ODE integration; they do not
state an Einstein-shell PDE existence theorem.

Implementable consequences for the current one-dimensional initial profiles:

- First split the few radial cells that dominate the existing second-derivative
  error majorant, keeping an exact complete cover and the same reconstruction.
- Preserve the common radial variable and correlated leading actual/nominal
  terms until after subtraction. The existing paired telescoping evaluation
  already addresses this problem; refinement should retain that cancellation.
- If a genuine interval-dependency floor remains, a shared local polynomial
  model with separately proved derivative remainders is a candidate next method.
  Every composition needs a domain-inclusion check. An enclosure for an unknown
  remainder is not permission to differentiate that remainder.
- A thinner upper enclosure is not a lower bound on true error. At a cutoff
  join, only an actual-minus-smoothed interval excluding zero proves a nonzero
  pointwise error. The exact even-kernel join identity must be proved directly
  from the actual piecewise regularity and retained smoothing width.

This check does not import software, assume a new Taylor-model implementation
is correct, or claim that subdivision must improve the final bound. Arithmetic,
coverage and remainder checks remain part of this checkpoint's own proof.

## Recent timelike-wall work and its limits

H. Helton, G. T. Horowitz and D. Marolf, *Dirichlet walls and the end of time*,
[arXiv:2606.05505v3](https://arxiv.org/html/2606.05505v3), revised 16 September
2026; the [abstract and revision record](https://arxiv.org/abs/2606.05505)
confirm the authors, title and version. The complete HTML body was retrieved;
the introduction, matter discussion and discussion of boundary well-posedness
were inspected for this limited comparison.

The paper considers Einstein-Hilbert gravity with a prescribed finite Dirichlet
wall and argues that open sets of initial data can reach finite-time singularities
at the wall. Its simplest cosmological constructions are vacuum with vanishing
or negative cosmological constant; it also studies BTZ walls in 2+1 dimensions
and trapped-surface constructions in higher dimensions. A scalar-potential
discussion is not a proof for arbitrary scalar-dependent material-shell laws.

Its cited An-Anderson boundary theorem is a different entry from the twisted
Dirichlet paper already reviewed in the preceding HDBLAST checkpoint:
Z. An and M. T. Anderson, *Well-posed geometric boundary data in General
Relativity, I: Dirichlet boundary data*,
[arXiv:2505.07128v2](https://arxiv.org/abs/2505.07128v2), revised 30 May 2026.
Only that paper's abstract and metadata were read in this check. The abstract
states local-in-time well-posedness for the vacuum Einstein IBVP with Dirichlet
data under a Lorentz-signature condition on the boundary Brown-York tensor,
and explicitly says the condition cannot generally be dropped. This should
not be conflated with part II, arXiv:2507.15567, on twisted Dirichlet data.

The practical requirement is to compare the actual boundary operators and
initial/corner constraints before importing a nonlinear theorem. A material
shell junction law has not been identified with prescribed induced-metric
Dirichlet data in this calculation. Local existence under one boundary law
also does not prove stable long-time cosmological evolution under another.
Neither reference establishes that the HDBLAST perturbative candidate is an
exact nonlinear family or that a higher-dimensional event caused our Big Bang.

## Coverage and access limits

The searches returned at most five and four web hits, respectively. Search
excerpts were used only to select full sources. An author-hosted PDF on the
umbilic Einstein IBVP returned 403 Forbidden; no theorem from it was used.
The earlier local HDBLAST literature notes were read before the new searches;
they already document Taylor-model provenance and the part-II An-Anderson
comparison. This note adds primary-text detail useful for error evaluation and
a bounded check of the recent Dirichlet-wall paper. It establishes no new
observational evidence, external specialist endorsement or novelty finding.

Raw retrievals remain in the working evidence directory. They need not be
copied into a public supplement; this original comparison and source links
are sufficient to describe the scope of the check.
