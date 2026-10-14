# Independent continuation review — 14 October 2026

This is internal AI-assisted mathematical review, not external peer review.
Reviewed input: `HDBLAST_SECOND_ORDER_LOCALIZED_20261014` recovered from the
uploaded ZIP, plus the published shell–Noether manuscript. No claim of
worldwide novelty or cosmological confirmation follows.

## Main verdict

The localized **second-order initial-constraint** result is defensible after
repairing the operator-specific justification of Delay's hypotheses. The
restriction to an annulus away from the shell and vertex is substantive.
Vanishing corrections give initial junction compatibility, not an evolved
second-order solution through shell contact. The revised gauge scope of the
monopole identity is appropriate.

The formula tests are evidence for the stated explicit formulas. Finite
sample-point exact tests alone are not proofs of general rational identities;
symbolic reductions and hand arguments must bear any universal claim.

## Required correction: Delay hypotheses

Primary source read in full relevant passages:
[Delay, arXiv:1003.0535v4](https://arxiv.org/html/1003.0535v4), particularly
Definitions 1.1–1.2, Theorems 1.3 and 5.1, Lemmas 8.1–8.2, Corollary 8.4,
and Section 9.3.

Lemma 8.1 assumes injectivity of the linear map on the **entire highest
derivative tensor**, not just injectivity of the principal symbol on a pure
covector power. It proves the asymptotic Poincaré inequality (API), not the
kernel restriction condition (KRC). For the Killing operator, the full
first-derivative map is symmetrization and has a nonzero antisymmetric
kernel. Consequently the manuscript sentence invoking Lemma 8.1 for both
operators and both properties is unsupported.

This is a repairable proof gap, not a counterexample to solvability.

Replacement paragraph:

> The scalar Hessian map `A ↦ ½(A − γ₀ tr_{γ₀} A)` is invertible on
> symmetric two-tensors in dimension four (its trace multiplies by `−3/2`),
> so `P_g*` has the full-derivative injectivity required by Delay's
> Lemma 8.1 and satisfies (API). The momentum adjoint is the Killing
> operator. Its principal symbol is injective, but symmetrization of a full
> first-derivative tensor has an antisymmetric kernel; use Delay's
> Section 9.3, which records the established (API) for divergence on
> symmetric two-tensors. By Lemma 1 both adjoint kernels are spanned by
> smooth boosts or rotations on a neighbourhood of the closed annulus,
> so (KRC) holds (also Delay's Corollary 8.4 for the Killing operator).
> Theorems 1.3 and 5.1 applied on `𝒜′` then give smooth solutions flat at
> its boundary and hence smooth zero extensions with compact support in
> `𝒜`.

The theorem numbering and the infinite-order boundary vanishing in Theorem
5.1 match the text of v4.

## Checks of the local symmetry arguments

The metric-only lapse classification can be verified directly without
the supplied code. Write `u = Ω′/Ω`, `f = Φ′` and
`u′ = u² − f²/3`, the background relation. For radial `N=n(s)`, twice the
adjoint is

```
2 P_g* n = C_δ δ + C_x xxᵀ,
C_δ = −6 n′ − 4s u n′ − 4s n″
      + (12u + 4s u′ + 8s u²)n,
C_x = 4 n″ − 8u n′ + (8u′ − 8u² + 4f²)n.
```

The `xxᵀ` component gives `n″−2u n′+f²n/3=0`. Eliminating it from
the trace component gives
`(1+2su)n′ = 2u(1+su)n`. The compatibility condition excludes
`l=0`, including the pointwise exceptional set `1+2su=0`. For
`l≥2`, the tangential trace-free Hessian argument excludes all harmonics;
for `l=1`, the mixed component forces `N=cΩx_a`. These are local
classifications on a concentric annulus, not just statements about global
background symmetries.

The full boost scalar-adjoint identity follows from
`U′=Ω⁻²(10f+4sf′+12suf)` and
`Δ_{γ₀}Φ=Ω⁻²(8f+4sf′+8suf)`. For `N=n(s)x_a`,

```
P_ψ* N = 4s f x_a Ω⁻² (n′ − u n),
```

so `N=Ωx_a` lies in the full scalar adjoint kernel as required for
necessity. For radial `N`, the same relation gives
`P_ψ* N=Ω⁻²f(4s n′−(2+4su)n)`, with kernel `N=cΩ√s`.

The Killing vector classification via Euclidean conformal Killing fields
uses the strict sign `(1/Ω)″>0`, excluding dilations, translations and
special conformal combinations. The rotations preserve the background
scalar. The density pairing for the momentum source is consistent.

## Scope, remaining gaps and housekeeping

1. For two independent copies, cancellation of **sums** of local moments
   does not supply compactly supported corrections on each separate copy.
   The per-copy requirement and fixed-parity specialization are correctly
   stated. No interface adjoint classification follows from this method.
2. Initial shell compatibility follows because every initial bulk and shell
   perturbation jet used in the junction equations can be chosen identical
   to the background near the cut. This does not solve the nonzero
   quadratic junction sources generated at later shell contact.
3. The areal-radius monopole is invariant under a linearized gauge change
   of the second-order correction with the first-order packet fixed; it is
   not invariant under gauge changes of the first-order packet. The revised
   manuscript makes this distinction. A sign is not a physical energy.
4. Explicit correction choice in Section 4 has scalar and scalar-momentum
   parts; Theorem 2's existence of a separate correction with these parts
   zero is compatible with that choice. Do not conflate the constructions.
5. Section 4 still contains the nonexistent path
   `basis/results_snapshot/construct_cc_M.json`; use
   `basis/results/construct_cc_M.json`.

## Most useful next scientific step

For the selected commuting quadrupole, provide an actual finite-slab
second-order bulk evolution in a specified hyperbolic gauge and derive the
quadratic pullback, Israel and scalar-junction sources. Check their
constraint propagation and corner compatibility before claiming a solution
through shell contact. A formal time jet or a compact initial correction
is useful intermediate progress but does not supply that evolution theorem.

A distinct, valid mathematical continuation is formal initial-data
solvability through every finite order in the coordinate-reflection
invariant sector. A proof is supplied in
[the formal reflection proof](../continuation/REFLECTION_FORMAL_INITIAL_DATA.md); it makes no convergence or
nonlinear spacetime claim.

## Review of the additional precontact evolution note

`PRECONTACT_SECOND_ORDER_DYNAMICS.md` uses the standard normally hyperbolic
wave-map reduction, the finite-order constraint propagation argument and
finite propagation. This is valid with Cauchy gauge data chosen to agree
with the background outside the containing annulus. Arbitrary global
pure-gauge Cauchy data would invalidate the claimed support property;
the gauge-data support choice is an essential sentence to add.

The causal bound is for the union of first- and second-order supports.
The moving shell is at `R_b(T)=√(b²+T²)`, so the support upper edge
`R_++|T|` leaves the strict gap
`b²−R_+²−2R_+|T|>0`. A positive-`s` controlled collar also requires
`R_-−2|T|>0`. These bounds correctly support a finite precontact
spacetime-coefficient theorem; they supply no theorem through contact.

The Gaussian-time kinetic acceleration factors in equation (9) can be
checked independently by canonical inversion:

```
K₁=(p/Ω)(sδ−3xxᵀ),   tr_{γ₀}K₁=sp/Ω³,
K₁·γ₀⁻¹·K₁=(p²/Ω⁴)(s²δ+3sxxᵀ).
```

The second parameter derivative of the metric kinetic terms is therefore
`8K₁·γ₀⁻¹·K₁−4(tr K₁)K₁=4sp²(sδ+9xxᵀ)/Ω⁴`, as displayed.
Similarly the scalar kinetic forcing is `−2(tr K₁)ν₁`, yielding
`3s𝒦(2s𝒦′+7𝒦)H²/(Ω⁶Φ′)`. The remaining displayed terms are the
required linear responses to the second-order spatial correction, not
terms that may be dropped. This uses Gaussian proper time, not ambient
conformal time.

The Weyl-tensor angular witness is appropriately described as the
projection of one quadratic **scalar-matter** source term. It does not
prove a nonzero projection of the full gravitational forcing; a possible
cancellation must be calculated before claiming that conclusion.

## Independent proof of the exact Hamiltonian angular reduction

For `γ₀=Ω(s)²δ` and a covariant spatial tensor `h` with
`tr_δh=0` and `h_{ij}x_j=0`, the scalar-metric term in `P_H(h,0)`
vanishes, as does `Ric_{γ₀}:h`, since both background tensors are radial
combinations of `δ` and `xxᵀ`. The trace-Laplacian term vanishes too.
Writing `h^{ij}=Ω⁻⁴h_{ij}` and
`Γ^i_{jk}=2(Ω′/Ω)(δ^i_kx_j+δ^i_jx_k−δ_{jk}x^i)`, every
connection or first conformal-factor derivative in
`D_jh^{ij}` contracts `h` with `x` or its trace. Thus

```
D_jh^{ij}=Ω⁻⁴∂_jh_{ij}.
D_iD_jh^{ij}=Ω⁻⁴∂_i∂_jh_{ij}       (four spatial dimensions),
P_H(h,0)=(2Ω⁴)⁻¹∂_i∂_jh_{ij}.
```

For `P=sδ−xxᵀ` and
`h=sP(∇²χ)P−tr(P(∇²χ)P)P/3`, the tensor in the unit sphere
frame is `r⁴(Hess_{S³}χ)^TF`. Its two flat divergences reduce to
`s div_{S³}div_{S³}(Hess_{S³}χ)^TF`. Since `Ric_{S³}=2g`,

```
div(Hess χ − gΔχ/3)=(2/3)d(Δ+3)χ,
P_H(h,0)=s/(3Ω⁴) Δ_{S³}(Δ_{S³}+3)χ.
```

This proves the conformal reduction and angular sign/normalization without
the inherited operator code. The new angular exact audit then correctly
checks the polynomial identity for the regenerated source and stored
coefficients by reduction modulo `|x|²=s`; its source denominators have
no spatial-coordinate dependence, so clearing denominators is valid.
This verifies the correction against that source exactly. It does not
constitute a second independent derivation of the source Hessian itself.
The audit runner should use explicit runtime checks rather than Python
assertions before claiming reliable optimized-mode or failure detection.
