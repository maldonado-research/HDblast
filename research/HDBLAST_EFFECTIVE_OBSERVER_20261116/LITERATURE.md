# Timelike Jacobi comparison: methods and source limits

16 November 2026 (America/Los_Angeles; 17 November UTC). This focused review concerns the quantitative normal-cut argument through q=25ell/16. The evaluated constants and specified target belong to the present proof; they are not imported from a generic comparison theorem.

| Primary source | Version and body inspected | Relevant method and limit |
| --- | --- | --- |
| Andersson–Howard, [Comparison and Rigidity Theorems in Semi-Riemannian Geometry](https://arxiv.org/abs/dg-ga/9707020v1) | v1, 25 July 1997; versioned HTML, §§1 and 3.1, Eq.(3.1), Theorem 3.2, Proposition 3.3 and adjacent examples | Matrix Riccati order comparison includes curvature and initial-shape hypotheses. Theorem 3.2 assumes both compared solutions already exist throughout the interval. Proposition 3.3 adds a two-sided continuation result. The first theorem alone cannot be used circularly to establish focal freedom. |
| Alexander–Bishop, [Lorentz and semi-Riemannian spaces with Alexandrov curvature bounds](https://arxiv.org/abs/0804.2493v1) | v1, 15 April 2008; versioned HTML, introduction and §4 | Their tensor inequality reverses the sectional-curvature direction on timelike planes. Their Weingarten sign is opposite to Andersson–Howard's. Curvature and shape conventions must be translated together. |
| Grant–Kunzinger–Ohanyan–Schinnerl–Steinbauer, [Timelike conjugate points in Lorentzian length spaces](https://arxiv.org/abs/2509.12855v4) | v4, 9 November 2026; first submitted September 2025; versioned HTML, Definition 2.7 and §§4.1–4.2 | The weak timelike Rauch theorem bounds locations of conjugate endpoints under its causal and curvature-comparison hypotheses. The authors explicitly leave quantitative Jacobi-growth comparison unestablished. Endpoint conjugacy does not replace a hypersurface-focal estimate with the actual initial cut shape. |

The direct equation in the actual normalized parallel screen is

```
S''-S=E S,       S(0)=Id,       S'(0)=0,
S(q)=cosh(q)Id+integral_0^q sinh(q-s)E(s)S(s) ds.
```

This Volterra identity permits a bound for the full matrix without diagonalizing E or assuming its values commute. A closed operator-norm estimate bounds the deviation from cosh(q)Id; a lower singular-value estimate then proves positive cuts. It avoids using the inverse Jacobi matrix before invertibility has been established. The actual phase retention, complete curvature/velocity/frame error bounds and initial time jet remain essential inputs.

Andersson–Howard's shape variable is B=-S'S^-1 for curvature coefficient R=-Id-E, wherever the inverse exists; its background is -tanh(q)Id. It is a different object from the Jacobi matrix S. The present background check S0=cosh(q)Id fixes the sign directly.

The clock is q=r/rho_b, where r is physical proper past duration. Velocities and second derivatives therefore carry rho_b and rho_b^2 factors. Equivalently, h(x_q,x_q)=-rho_b^2 while the normalized metric gives unit timelike q velocity.

A positive cut differential is one part of the global cylinder argument. Complete label spheres, the retained normal flow, Gauss orthogonality and the anchored degree argument supply additional injectivity information. None of this establishes an arrival theorem, a nonzero response, detector calibration or measured cosmological confirmation.

Firecrawl was used for focused searches, canonical version checks and primary-body retrieval. The earlier observer sources were reused read only. The named sections were inspected; returning full HTML is not a claim to have read every proof. The date-bounded 2025–2026 screen and the sources above do not constitute an exhaustive survey or a novelty determination.
