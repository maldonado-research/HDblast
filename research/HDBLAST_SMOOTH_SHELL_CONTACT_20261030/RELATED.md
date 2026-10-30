# Focused primary-source comparison for shell contact

This is a bounded literature check, performed on 30 October 2026, rather than an exhaustive novelty review. It supplements the earlier smooth-feedback comparison. Primary full text was inspected selectively at the sections named below. None of the references supplies a theorem for the present coupled Einstein–scalar shell system without a separate verification of hypotheses.

## Robin reflection and the separate boundary contribution

J. D. Bondurant and S. A. Fulling, *The Dirichlet-to-Robin Transform*, arXiv:math-ph/0408054v1 (2004), published in *Journal of Physics A* **38** (2005), 1505–1532. [Primary text](https://arxiv.org/html/math-ph/0408054), especially §§II.1–II.2 and III.2, equations (1), (3), (21)–(23); [record](https://arxiv.org/abs/math-ph/0408054).

The paper uses the inward-normal convention `∂r u(0)=κu(0)` on a flat half-line. For constant scalar `κ>0`, its one-dimensional wave propagator contains the Neumann image term and a delayed reflected convolution. For zero initial velocity, equation (22) contains the correction

```
−κ 1_(t>r) ∫_0^(t−r) exp[−κ(t−r−y)] f(y) dy.
```

This is a directly useful comparison: a Robin boundary cannot generally be imposed by simply declaring the physical field even or odd. A reflection construction for the auxiliary homogeneous Dirichlet/Neumann propagator and a separate causal boundary term must be kept distinct. The sign convention also makes a useful elementary test of a lift: if `Lg(t,r)=−∫_0^(t−r)_+ g(s)ds`, then its inward derivative at `r=0` is `g(t)`.

The source does **not** justify commuting the present output parity projector through the coupled bulk matrices. Its exact transform relies on flat geometry and commuting constant-coefficient scalar operators; the principal elementary Green-function construction restricts to `κ>0`. The full 23-mode problem, with its mixed projectors, lower-order coupling, affine forcing and frozen versus variable coefficients, retains its own proof and error ledger. The delayed convolution is established classical mathematics, not a new HDBLAST discovery. This comparison motivates checking reflection and boundary lift separately; it does not certify a complete nonlinear solution.

## Recent reconstruction results and their limits

Marcus J. Grote, Omar Lakkis and Carina Santos, *A posteriori error estimates for the wave equation with mesh change in the leapfrog method*, arXiv:2411.16933v2 (revision 26 June 2025). [Primary text](https://arxiv.org/html/2411.16933v2), §§2.1, 3–4; [version record](https://arxiv.org/abs/2411.16933v2).

This paper constructs fully computable energy-norm estimates for a scalar wave equation discretized by conforming finite elements and leapfrog time stepping, including mesh transfer and local time stepping. Equations (2.1)–(2.4) specify homogeneous Dirichlet–Neumann boundary conditions, with a Dirichlet portion of positive measure. Its reconstruction yields a continuous error equation and accounts for the different residual contributions. This supports keeping approximation, reconstruction and PDE error bounds distinct. It does not provide the present matrix Robin boundary-defect estimate, a pointwise `C2` error certificate, or an Einstein-constraint estimate. Its recent mesh-change machinery does not change the elementary characteristic/interval method used for this checkpoint.

Marvin Knöller and Jörg Nick, *The temporal domain derivative in inverse acoustic obstacle scattering*, arXiv:2510.21471v1 (24 October 2025). [Primary text](https://arxiv.org/html/2510.21471v1), §1, Theorems 3.3 and 4.4.

The authors study the shape derivative of a scalar acoustic scattering map for a sound-soft obstacle in two or three dimensions. The derivative solves another wave problem with a derived boundary condition. Laplace-domain estimates give temporal regularity and Runge–Kutta convolution-quadrature convergence under sufficient regularity of the incident boundary trace. In Theorem 4.4 the pointwise constant depends inversely on distance from the boundary: the stated interior point estimate is not a boundary-trace certificate. The prime denotes a domain/shape derivative, not merely the time derivative used in a differentiated shell-error equation. This paper therefore does not fill the shell-contact or differentiated-boundary gap in the present work, and its named Convolution Quadrature method is not the ordinary interval integration used here.

## Consequence for this checkpoint

The useful lesson is structural: verify the physical-side source, its reflected **output**, and the causal Robin lift in the correct order, and report the residual boundary defect independently of the candidate value. Recent reconstruction and acoustic shape-derivative results do not remove the need for model-specific coefficient, source, initial-jet and boundary estimates. The search did not establish priority, novelty, observational support, or a transferable nonlinear existence theorem.
