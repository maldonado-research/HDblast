# Primary-source comparison: generic warped-product TT equation

Prepared 2026-10-31. Bounded discovery: two web searches, followed by reading the actual bodies of two primary papers. This note is a synthesis and formula audit; it does not reproduce full third-party texts, establish novelty, or claim exhaustive literature coverage.

## Closest applicable primary result

Hideo Kodama, Akihiro Ishibashi, and Osamu Seto, *Brane World Cosmology — Gauge-Invariant Formalism for Perturbation*, Physical Review D **62**, 064022 (2000), arXiv:hep-th/0004160. Read version v3 (7 August 2000), primary body: <https://arxiv.org/pdf/hep-th/0004160>.

Relevant exact locations:

- Section II.A, equation (2): arbitrary local warped-product background `g_ab(y) dy^a dy^b + r(y)^2 gamma_ij(x) dx^i dx^j`, where the fiber has constant sectional curvature K. The symbol y here denotes *all m base coordinates*, so it permits both our radial coordinate y and our time coordinate t. No requirement that r depend on only one base coordinate, and no static-base hypothesis.
- Section II.A, equations (7)–(10): matter is allowed, with background `T_ai=0`, `T^i_j=P delta^i_j`. The background satisfies Einstein's equation, including a cosmological constant if present. The later vacuum specialization occurs in Section IV, not in the tensor result in Section II.
- Section II.C.1, equations (19)–(20): `Delta_hat T_ij = -k^2 T_ij`, trace zero, divergence zero.
- Equation (21): `h_ab=h_ai=0`, `h_ij=2 r^2 H_T T_ij`. H_T is gauge invariant because an infinitesimal diffeomorphism has no tensor-type component.
- Equation (22): `delta T^i_j = tau_T T^i_j`, where tau_T is the tensor anisotropic stress (a *mixed-index* matter coefficient).
- Equation (23):

\[
-\Box_B H_T-\frac{n}{r}D r\cdot D H_T
+\frac{k^2+2K}{r^2}H_T=\kappa^2\tau_T,
\qquad \Box_B=D^aD_a.
\]

The source uses a normalization factor 2 in h_ij. In the present notation `k_ij=R^2 U T_ij`, the binding is `H_T=U/2`.

## Substitution for the specified five-dimensional background

Use

\[
g_B=dy^2-\rho(y)^2dt^2,\quad R=\rho(y)\cosh t,
\quad a=\rho'/\rho,\quad n=3,\quad K=1,\quad k^2=6.
\]

Then

\[
\Box_B U=U_{yy}+aU_y-\rho^{-2}U_{tt},
\qquad R^{-1}D R\cdot D U=aU_y-\rho^{-2}\tanh t\,U_t.
\]

Consequently the left side of KIS (23), with `H_T=U/2`, is exactly

\[
\frac{1}{2\rho^2}\left[
U_{tt}-\rho^2 U_{yy}-4\rho^2aU_y
+3\tanh t\,U_t+8\operatorname{sech}^2t\,U
\right].
\]

Thus the stated `+8/cosh^2(t)` potential is `k^2+2K=6+2`; the first derivative coefficient 4a combines one a from the base d'Alembertian and three from the fiber volume factor.

This quantity is the coefficient of the *mixed* angular tensor Einstein residual. If the residual is written with both angular indices down, its coefficient multiplying `T_ij` acquires the factor `R^2`. A definition of delta E_TT should specify which normalization is used. On a background solution, the mixed TT part is unchanged when passing between the Einstein residual and its trace-reversed form; their difference is angular trace.

## Why the scalar field does not supply a tensor source

For a minimally coupled canonical scalar with background `phi_0=phi_0(y)`, `T^i_j=P delta^i_j`. At first order, the anisotropic product `nabla^i phi nabla_j phi` has zero variation because both angular background derivatives vanish. The remaining first variation is proportional to `delta^i_j`, even when the scalar perturbation depends on all coordinates. Hence the tensor anisotropic stress `tau_T` in KIS (22) is zero. In covariant indices, the variation includes the isotropic-pressure term `P h_ij`; it is removed by raising an index, not an extra TT source.

For an arbitrary metric perturbation, the scalar/vector/tensor harmonic splitting on the compact round S^3 and orthogonality isolate this tensor coefficient. KIS states after equations (19)–(20) that expansion coefficients with different eigenvalues decouple. To establish the full-mode statement rigorously in the project's setting, retain the project's explicit S^3 TT projection/integration-by-parts argument and identify its coefficient with this operator. The source is corroboration of that identity, not a substitute for specifying the project's projection and gauge/regularity conventions.

## Independent primary corroboration with explicit covariant matter correction

Rong-Gen Cai and Li-Ming Cao, *Generalized Formalism in Gauge-Invariant Gravitational Perturbations*, Physical Review D **88**, 084047 (2013), arXiv:1306.4927, primary body: <https://arxiv.org/pdf/1306.4927>; publisher DOI landing page discovered: <https://link.aps.org/doi/10.1103/PhysRevD.88.084047>.

Section II.A equation (2.1) allows an arbitrary Lorentzian m-dimensional base and Einstein n-dimensional fiber. Section III.B.1 equation (3.10) is the tensor equation *before spectral expansion*:

\[
-\Box_B(h^{(2)}_{Tij}/r^2)
-n(D^ar/r)D_a(h^{(2)}_{Tij}/r^2)
+r^{-2}[\widehat\Delta_L-2(n-1)K](h^{(2)}_{Tij}/r^2)
=2\kappa^2\left[\delta T^{(2)}_{Tij}/r^2-P h^{(2)}_{Tij}/r^2\right].
\]

The scrape's formula typography mislabels the superscript on the base d'Alembertian, so the above uses Box_B to avoid copying that parser artifact. The surrounding definition identifies the base operator.

Equations (3.11)–(3.12) define the Lichnerowicz operator. On a constant-curvature fiber, `Delta_L=-Delta_hat+2nK` on TT tensors, so `lambda_L-2(n-1)K=k^2+2K`.

Section IV.A equations (4.1)–(4.3) explicitly set

\[
h^{(2)}_{Tij}/r^2=2H_T\mathbb T_{ij},\qquad
\delta T^{(2)}_{Tij}/r^2-P h^{(2)}_{Tij}/r^2=\tau_T\mathbb T_{ij},
\]

and obtain

\[
\Box_B H_T+n(D^ar/r)D_aH_T
-r^{-2}[\lambda_L-2(n-1)K]H_T=-\kappa^2\tau_T.
\]

This is KIS (23) multiplied by -1. It makes the covariant stress correction particularly explicit and confirms applicability before imposing vacuum or static black-hole conditions.

## Limits

The audit assumes the Einstein equation with a minimally coupled canonical scalar, a smooth background solution, and the compact round unit S^3. It does not establish a global boundary/initial value theorem, nonlinear integrability of linear modes, regularity at zeros of rho, or classification of all perturbations. The project's differential operator applies on a regular region where rho>0. No claim is made that the exact-background obstruction or the project-specific TT projection is new to the literature.

## Retrieval provenance

The primary PDF bodies were returned as parsed Markdown on 2026-10-31, and the sections and equations listed above were inspected. The returned full-body representations were not retained after the reading session, so no full-body digest is asserted. The accompanying public receipt records this limit and pins this comparison note instead. This note and the mathematical verification do not depend on access to a private retrieval cache or authenticated retrieval service.
