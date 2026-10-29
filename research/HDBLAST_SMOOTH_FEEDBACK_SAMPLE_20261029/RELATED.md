# Focused primary-source comparison for the feedback evaluation

This is a bounded source comparison, not an exhaustive literature or novelty
review. The local scalar quadrature and coupled boundary system retain their
own model-specific assumptions and error ledger.

* Fredrik Johansson, *Numerical integration in arbitrary-precision ball
  arithmetic*, arXiv:1802.07942v1, sections 2–2.1,
  https://arxiv.org/html/1802.07942. The primary paper implements rigorous
  quadrature in Arb using Petras-type adaptive bisection/Gaussian quadrature
  and validates analyticity on complex regions for the high-order error bound.
  Its output ball radius is a rigorous a posteriori enclosure even when the
  requested absolute/relative tolerance cannot be met. This supports the
  distinction between requesting precision and proving an integration error.
  The present real weighted-interval midpoint approach uses an explicit
  Lipschitz bound instead; it does not import this paper's analytic-domain
  estimates or claim to use Arb/Petras.

* Matthew Hassell and Francisco-Javier Sayas, *Convolution Quadrature for Wave
  Simulations*, arXiv:1407.0345v1, section 1.1 and the causal-operator framework,
  https://arxiv.org/html/1407.0345v1. These authored course notes explain causal
  convolutions, the Duhamel principle and operator-valued wave convolutions,
  then develop Lubich convolution quadrature for acoustic wave boundary
  integral formulations. Their capitalized Convolution Quadrature is a
  particular time-discretization method; ordinary integration of the local
  characteristic feedback kernel is not automatically that algorithm. The
  notes do not certify the present 23-mode mixed shell system.

* Théophile Chaumont-Frelet and Alexandre Ern, *Damped energy-norm a posteriori
  error estimates using C2-reconstructions for the fully discrete wave
  equation with the leapfrog scheme*, arXiv:2403.12954v2, sections 1 and 3–4,
  https://arxiv.org/html/2403.12954v2. The primary text constructs a globally
  C2 time reconstruction and uses the reconstructed error equation in a
  damped energy estimate. It assumes the stated scalar-wave/finite-element
  setting, homogeneous Dirichlet boundary, zero initial data and smooth
  forcing supported away from zero. A smooth nominal reconstruction and a
  rigorous numerical enclosure of one feedback sample are separate from a
  full PDE error theorem; this paper does not supply the present coupled
  mixed-boundary sup-norm second-derivative estimate.

All three bodies were read selectively in primary full text. No theorem is
transferred to the Einstein–scalar shell problem without checking its complete
hypotheses. The local integral identities in theory/LOCAL_FORMULAS.md are proved there
by elementary changes of variables and integration by parts, rather than
attributed as a new physical or mathematical discovery.
