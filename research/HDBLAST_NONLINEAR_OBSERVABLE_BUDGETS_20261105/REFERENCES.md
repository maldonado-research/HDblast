# Focused references for the conditional observable error budget

5 November 2026. Working methods comparison. This is a bounded
methods search, not an exhaustive novelty assessment. Earlier project
references and quotient-budget bookkeeping were consulted first.
The three primary papers below were inspected in the indicated **full-text
passages**; none received a complete proof audit.

1. **Catherine Bandle and Alfred Wagner**, *Second Domain Variation for
   Problems with Robin Boundary Conditions*, Journal of Optimization Theory
   and Applications **167** (2015), 430–463,
   [doi:10.1007/s10957-015-0801-1](https://doi.org/10.1007/s10957-015-0801-1).
   Read [arXiv:1403.2220v2](https://arxiv.org/pdf/1403.2220v2), especially
   §2.2.2, equations (2.17)–(2.20), and §3, equations (3.4), (3.14),
   (3.18)–(3.19). These give the determinant-based second surface-density
   expansion, fixed-domain pullbacks, and the distinction between material
   and untransported derivatives. They directly support the calculus used
   to budget actual area. Their elliptic Robin problem does not supply a
   dynamical shell family, observer-clock bounds, or the project's constants.

2. **Marc Dambrine and Jimmy Lamboley**, *Stability in shape optimization
   with second variation*, Journal of Differential Equations **267** (2019),
   3009–3045,
   [doi:10.1016/j.jde.2019.03.033](https://doi.org/10.1016/j.jde.2019.03.033).
   Read [arXiv:1410.2586v4](https://arxiv.org/pdf/1410.2586v4), §1.2–§1.3
   (in particular the improved-continuity assumption and integral Taylor
   formula on PDF pp. 5–6), and §4.1, Lemma 4.3 and its proof. The close
   methodological point is that a modulus controlling the second derivative
   gives a controlled Taylor remainder; smooth pullback Jacobians provide
   geometric expansions. Their shape-stability and elliptic estimates do
   not establish such a modulus for this Einstein–scalar observable.

3. **Leonardo Pires and Rodrigo A. Samprogna**, *A quasilinear problem with
   large diffusion and nonlinear dynamic boundary conditions*, Journal of
   Elliptic and Parabolic Equations **12** (2026), 1277–1295, published
   10 April 2026,
   [doi:10.1007/s41808-026-00460-w](https://doi.org/10.1007/s41808-026-00460-w).
   Read the publisher's full-text introduction, §2 operator formulation,
   Lemma 2.2 and Theorem 2.4,
   and §5 Theorem 5.1 with the opening energy comparison. This is a verified
   recent quasilinear dynamic-boundary parameter-limit source: it treats
   parabolic bulk/surface diffusion and convergence as diffusion becomes
   large. It does **not** establish twice differentiable amplitude dependence,
   a Hölder/Lipschitz second-parameter-derivative estimate, or a reversible
   Einstein–scalar junction problem. No result from it is imported.

The following are elementary specializations of the first two methods, not
new external theorems. For the positive-definite metric `h_e` on the **same
observer proper-time cut**, put `P=h_e^{-1}h'_e` and
`j_e=sqrt(det h_e/det h_0)`. Then

```
j'  = (j/2) tr P,
j'' = j[(1/2)tr(h^{-1}h'')+(1/4)(tr P)^2-(1/2)tr(P^2)].
```

The metric must already contain the full observer, embedding and clock
pullback. In a common coordinate chart, `A_e=Phi_e(X_e)-phi_b` gives

```
A'' = (partial_e^2 Phi)(X)
      +2 D(partial_e Phi)(X)[X']
      +D^2 Phi(X)[X',X']+D Phi(X)[X''].
```

For an odd `C^2` moment map into `C(I)`, with `S=M'(0)`, let
`omega(r)=sup_{|u|<=r} ||M''(u)-M''(0)||`. Oddness gives `M''(0)=0`, and
integral Taylor yields

```
||M(e)/e-S|| <= |e| integral_0^1 (1-t) omega(t|e|) dt.
```

Thus `C^2` alone gives `o(|e|)` for the normalized error, not a cubic
unnormalized remainder. If `omega(r)<=C r^alpha`, `0<alpha<=1`, the
normalized bound is `C|e|^(1+alpha)/[(1+alpha)(2+alpha)]`; Lipschitz
`M''` gives `C e^2/6`. These implications require the stated **uniform**
modulus on the admitted common domain.

The remaining project-specific work is still to certify that domain and
the pulled-back scalar, metric/area and modulus bounds. This comparison
evaluates no nonlinear constant or numerical amplitude radius. No focused
search result was found here that supplies the missing model-specific
`C^2`/modulus estimate; that is a report of this bounded search, not a
nonexistence or originality claim.
