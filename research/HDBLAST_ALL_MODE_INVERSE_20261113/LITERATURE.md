# Focused literature check for the all-mode weighted initial inverse

**Focused method note.** Scientific search and retrieval: **13 November 2026 UTC**. This supplements the frozen [continuation](CONTINUATION.md), using repository base `d26f530122f8bfa186291ce601e4876e3646d760`. Primary identifiers, passage coverage and access limits are recorded in the [source ledger](LITERATURE_SOURCES.json). Retrieved excerpts remain in private provenance.

The retained **norm shifts and weight powers match Delay’s primary theorem**: scalar order four gives Ḣ16→Ḣ12, momentum order two gives Ḣ14→Ḣ12, with the same weighted norm on each side. The theorem remains conditional on its operator, coercivity, kernel and scaling hypotheses; it supplies no evaluated inverse norm. A separate primary-source issue was found in the printed boundary Korn identity of Chruściel–Delay: its displayed powers need a corrected substitution before numerical use. This does not establish an error in the retained block mapping.

## Exact norm and operator dictionary

Delay defines the squared Sobolev norm by integrating the sum of `φ^(2j)|∇^j u|² ψ²`, including j=0, with the background geometric volume. The dotted space is the weighted completion of compactly supported sections. It is not an ordinary Dirichlet space obtained by imposing classical boundary traces on the potentials. [Delay v4, §2, (2.1)](https://arxiv.org/html/1003.0535v4#S2.E1).

For the selected four-dimensional member, the translation is exact:

```
x = d = (ell/100) r0,
a = 2, s = beta = ell/100,
φ = x² = qweight,
ψ = x^(2a−n) exp(−s/x) = exp(−1/r0) = w.

||u||²_(H^k_(qweight,w))
  = ∫_A Σ_(j=0)^k qweight^(2j) |∇^j u|² w² dμ_(gamma0).
```

Thus the unsquared differentiated quantity is `qweight^j ∇^j u`; the squared derivative weight is `qweight^(2j)=d^(4j)`. Replacing it by `d^(2j)` would change the space.

For an order-m underdetermined elliptic operator P, Delay’s normal operator and theorem are

```
L = ψ^-2 P ψ² φ^(2m) P*,
π_(K⊥) L : K⊥ ∩ Ḣ^(k+2m)_(φ,ψ) → K⊥ ∩ Ḣ^k_(φ,ψ),
K = ker(φ^m P*) in Ḣ^m_(φ,ψ).
```

Orthogonality and projection here use **L²_ψ**, with measure `ψ²dμ`. [Delay, Proposition 3.6 and Theorem 3.8, (3.6)–(3.9)](https://arxiv.org/html/1003.0535v4#S3.Thmtheorem8).

| Selected block | m | Correction multiplier | Energy multiplier | Domain → range at k=12 |
| --- | --- | --- | --- | --- |
| Scalar metric constraint | 2 | `w² qweight⁴ P_g*` | `qweight² P_g*` | `Ḣ16_(qweight,w) → Ḣ12_(qweight,w)` |
| Fixed-density momentum constraint | 1 | `w² qweight² P_m*` | `qweight P_m*` | `Ḣ14_(qweight,w) → Ḣ12_(qweight,w)` |

The normal-operator multiplier is the square of the energy multiplier; interchanging them would lose two powers of qweight in the scalar block or one in the momentum block.

The retained normalization `f=w^-2 f_physical` reconciles the theorem’s weighted projection with the usual physical obstruction pairing:

```
<f,k>_(L²_w) = ∫ (w^-2 f_physical) k w² dμ
             = ∫ f_physical k dμ.
```

On the invariant subspace, averaging over the finite reflection group makes every invariant source orthogonal to the full kernel provided the kernel is group-stable and has no invariant part. The inverse is then invariant by uniqueness. This argument must use the actual geometric action on one-forms, the radial weights and fixed positive pairings; averaging bare components under an unrelated scalar action would not prove it.

The scalar principal Hessian map is of the injective full-Hessian type used in [Delay Lemma 8.1](https://arxiv.org/html/1003.0535v4#S8.Thmtheorem1), after checking the retained lower-order terms. Symmetric gradient is not injective on an arbitrary first-derivative tensor: the momentum API must instead use the Killing/Korn route in [§9.3](https://arxiv.org/html/1003.0535v4#S9.SS3). The selected convention `P_m=-2 div` has `P_m*=L_(Z#)gamma0=2S(Z#)`, with indices raised as prescribed. Chruściel–Delay uses `S=(1/2)sym∇`; those factors must be converted before assigning a constant. The formal-adjoint calculation is also explicit in [Isett–Mao–Oh–Tao, Appendix A.4, (A.13)–(A.14)](https://arxiv.org/html/2509.04617v1#A1.E14).

## What is required beyond the theorem statement

Delay’s API is a boundary-collar estimate of the form `||u||_(L²_w)≤C||qweight^m P*u||_(L²_w)` outside a compact interior set. Surjective principal symbol and finite-dimensional kernel alone do not replace this assumption. Proposition 3.6 combines it with interior control and a compactness argument to obtain coercivity on the weighted kernel complement. That proof establishes existence of a finite constant, without evaluating its value for this model. [Delay, (1.1), Lemma 3.5 and Proposition 3.6](https://arxiv.org/html/1003.0535v4#S3.Thmtheorem6).

Appendix 10 requires uniformly bounded scaled derivatives of the weights, comparability of φ and ψ on φ-sized balls, balls retained inside the domain, and uniformly controlled scaled elliptic coefficients. Its Sobolev proof also uses bounded overlap and local chart estimates. In dimension four the local change-of-variables factor in the **unsquared** norm is `w_i q_i²`, agreeing with the retained scaled-cover mapping. Fixed tensor and density conversions preserve qualitative boundedness but their numerical multiplier and derivative constants still have to be charged. [Delay, (10.1)–(10.6), Lemma 10.2](https://arxiv.org/html/1003.0535v4#S10.Thmtheorem2).

A useful quantitative bookkeeping target follows directly. Suppose an all-mode coercivity constant and a high-order regularity constant have actually been evaluated, in these precise norms:

```
||u||_(Ḣ^m) ≤ C_gap ||qweight^m P*u||_(L²_w),
||u||_(Ḣ^(12+2m))
  ≤ C_reg [||Lu||_(Ḣ12)+||u||_(L²_w)].
```

The weak energy identity gives `||u||_(L²_w)≤C_gap²||f||_(L²_w)` for `Lu=f`, hence

```
||L^-1 f||_(Ḣ^(12+2m))
  ≤ C_reg (1+C_gap²) ||f||_(Ḣ12).
```

This is a derived interface, not an externally supplied number. Apply it separately to scalar and momentum blocks, then combine using the registered direct-sum norm. The constants must control every allowed angular sector and both annulus ends. A large angular eigenvalue at a fixed interior radius does not give an endpoint-uniform estimate when the weighted angular principal factors vanish there.

The dotted potential domain permits weighted boundary growth. Physical zero extension concerns `A u`, not a freely imposed zero trace of u. Delay’s smooth-extension corollary assumes its corresponding smooth-source regularity; it cannot turn arbitrary Ḣ12 data into a C∞ correction. The retained finite-regularity multiplier/zero-extension argument remains the appropriate separate join.

## A source identity that must be corrected before numerical use

The exact primary [Chruściel–Delay v2 PDF](https://arxiv.org/pdf/gr-qc/0301073v2) was inspected in memory; Appendix D, p.67 was also viewed as a small image. Corollary D.5 prints `V=∇x/x` on its left side but a leading `x^(2t−4)` weight on its right side. Their Proposition D.2 gives leading `x^(2t−3)` with that V and `u=−s/x+t log x`. Using `V=∇x/x²` produces the printed right-side power.

An independent flat one-dimensional check makes the mismatch explicit. Put `W=x^(2t)exp(−2s/x)`. The printed left side is exactly

```
(3/2) ∫ W Y'Y/x
  = −(3s/2) ∫ W Y²/x³
    −(3/4)(2t−1) ∫ W Y²/x².
```

It cannot have the displayed leading x^-4 factor. With `V=∇x/x²`, the valid weighted identity instead has leading

```
−∫ W x^-4 (s/2)
       [|dx|²|Y|² + 2<Y,∇x>²],
```

plus controlled lower-order terms. Choosing t=2 and estimating the corrected left side yields the required `||Y||_(L²_w)≤C||x²S(Y)||_(L²_w)` on a sufficiently small collar, for fixed s≠0. A numerical proof must evaluate that collar and the remainder. This is an apparent typographical mismatch in the cited version, with an independently derived correction; it is not an asserted published erratum.

The compact-boundary assumptions include a defining function with |dx| bounded away from zero near the boundary, a metric continuous to the boundary and locally W1,∞, and `|Hess x|=o(x^-1)`; the smooth retained background satisfies the qualitative version. Proposition C.4 has the explicit leading Hardy coefficient `s²−epsilon`, but the sufficiently small collar is still an input to quantify. [Chruściel–Delay, Appendix C, pp.63–64, and Appendix D, pp.66–67](https://arxiv.org/pdf/gr-qc/0301073v2).

Two further distinctions prevent a mistaken source transfer. Their full coupled constraint theorem has staggered momentum/scalar Sobolev indices, whereas Delay 3.8 applies directly to the separated fixed-order blocks above. Their nonlinear Theorem 5.9 chooses the growing source weight exp(s/x), s>0; its linear boundary argument works for either nonzero sign, but the nonlinear theorem’s sign convention must not be copied wholesale into the selected decaying-potential formulation.

## Related 2025–2026 methods and actual applicability

These are focused comparisons, including sources already identified in the retained weighted-regularity review.

| Source and verified primary version | Relevant inspected evidence | What can be used here |
| --- | --- | --- |
| Isett, Mao, Oh and Tao, [arXiv:2509.04617v1](https://arxiv.org/html/2509.04617v1), **4 September 2025**; retrieved primary record lists no journal reference. | Theorems 1.1–1.2, Theorem 1.15, Remark 1.18, quantitative averaging hypotheses, and Appendix A.4. | A constructive support-controlled inverse and dual coercivity strategy with finite-rank obstructions. The inspected bounded-domain estimates are in ordinary Sobolev spaces. Its weighted extension at infinity has additional asymptotic hypotheses. Neither gives this exponential-boundary normal inverse automatically. |
| Le Floch and LeFloch, [arXiv:2607.28599v1](https://arxiv.org/html/2607.28599v1), **30 July 2026**; retrieved primary record lists no journal reference. | Theorem 2.4, (2.19)–(2.21), its geometric-constant dependencies, and passages locating the weighted Korn/Poincaré/Hardy and shell estimates. Appendix A was not independently read in this update. | Related all-angular stability organization. Its sphere-at-infinity and conic radial regime differ from the compact two-end annulus. Its own geometric inequalities must be checked, and no numerical constant transfers to the registered operator. |
| Mengesha and Scott, [arXiv:2610.33325](https://arxiv.org/abs/2610.33325), primary v1 **23 October 2026**; no journal reference retrieved. | Indexed in-body introduction, (1.1)–(1.3), and Corollary 4.1. The indexed text carries an internal 27 October date without a version marker; that is not treated as a revision. | Fractional Euclidean Korn estimates and rigid-motion projection are related methodology. They do not supply the requested integer-order curved exponential-weight estimate. No constant is imported. |

A different support-preserving right inverse cannot simply replace the selected variational inverse. It could help bound the latter only after proving the relevant weighted estimates and transferring them through the actual adjoint pairing and projection.

## Primary three-sphere vector spectrum and the Killing gap

The retained [Lindblom–Taylor–Zhang primary source, arXiv:1709.08020v1](https://arxiv.org/html/1709.08020v1#S3), **23 September 2017**, was reread specifically for the complete one-form spectrum. Its §3 contains three vector families A=0,1,2, index ranges in Table 1, rough-Laplacian equations (15)–(17), divergence equations (18)–(20), orthonormality (21), and the full vector expansion (22)–(23). Both transverse families are needed.

On the **unit** round S3, write positive scalar eigenvalues as λ=k(k+2). The paper gives:

| One-form family | Allowed indices | Positive rough eigenvalue −∇² | Positive Hodge eigenvalue |
| --- | --- | --- | --- |
| Exact gradients, A=0 | k≥1; ell=0,...,k | λ−2 | λ |
| Coexact families, A=1,2 | k≥1; ell=1,...,k | k(k+2)−1 | μ=(k+1)² |

The Hodge column is a derivation from the paper’s rough eigenvalues and the Weitzenböck identity Δ_H=−∇²+Ric=−∇²+2 on unit-S3 one-forms. The two coexact families have total multiplicity 2Σ_(ell=1)^k(2ell+1)=2k(k+2). Thus μ=4 at k=1 has multiplicity six, and the next coexact eigenvalue is μ=9 at k=2.

To identify the first sector, integration by parts gives the explicit convention check

```
||L_V kappa||²
  = 2[||∇V||²+||div V||²−2||V||²].

coexact eigenmode:
||L_V kappa||² = 2(μ−4)||V||².

exact eigenmode df:
||L_(grad f) kappa||² = 4(λ−2)||df||².
```

Consequently the six coexact μ=4 modes are Killing and are precisely the ambient rotation space. Orthogonality to this space leaves the conservative coexact gap μ≥9. This does not assert that μ=9 occurs in the reflection-invariant subspace; further symmetry may enlarge its gap.

These are integrated angular identities, not a pointwise estimate of an arbitrary full derivative tensor by its symmetric part. For radius R3, restore the R3^-2 factors displayed in the primary equations. The actual warped-cylinder connection, radial derivatives, weights and canonical-density factors still require their separate joins. No new finite-mode restriction of the nonlinear family is inferred.

## Publication, access and stopping limits

Delay’s mathematical text is pinned to **arXiv v4, 17 November 2011**. Publisher metadata verifies Communications in Partial Differential Equations 37(10), 1689–1716 (October 2012), [DOI](https://doi.org/10.1080/03605302.2012.711794); its retrieved body was a navigation stub. Chruściel–Delay is pinned to **arXiv v2, 11 July 2003**, with the primary journal reference Mémoires de la Société Mathématique de France 94 (2003). The PDF’s internal 2018 date and the paper index’s later “updated” field do not establish a new version.

Corvino’s [primary publisher page](https://link.springer.com/article/10.1007/PL00005533) verifies Communications in Mathematical Physics 214, 137–189 (October 2000), received 8 November 1999 and accepted 27 March 2000. Its body was unavailable: the PDF URL returned the subscription preview, and the paper index returned no passages. Therefore no exact norm shift or theorem statement from Corvino’s unread body is claimed verified. The exact matching norm map is supported directly by Delay’s fully retrieved theorem.

The ledger records **27 Firecrawl calls**, **18 explicitly reported credits**, and **12 calls with unreported usage**; the total actual cost is unknown. No broad crawl or autonomous research job was launched. Primary-paper bytes and inspection images were held in memory; only these two private deliverables were written.

The model-specific all-mode C_gap and C_reg must be evaluated separately, with the fixed density, physical pairings, weights and zero-extension domain retained; the four companion proof notes supply that calculation for this checkpoint. Literature collection stops here because further retrieval would not calculate those model-specific numbers. No nonzero amplitude, nonlinear tube, observational confirmation or novelty theorem follows from this review.

## Public-summary-ready paraphrase

The selected all-mode weighted normal inverse uses the same source and target weight, with scalar gain four derivatives and momentum gain two. Its qualitative isomorphism matches the primary gluing theorem. A numerical certificate must additionally evaluate global weighted coercivity and all-mode scaled-cover regularity in the registered background; recent related inverse and Korn methods do not supply those constants automatically.
