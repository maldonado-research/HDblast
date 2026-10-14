# Assessment of an exact nonlinear initial-data upgrade

Date: 14 October 2026. Internal AI-assisted mathematical assessment. This file is an isolated continuation note; no frozen repository or published checkpoint was changed.

## Conclusion

A **finite-Sobolev-regularity exact initial-data family** can be justified from the corrected annular linear solver, using an explicitly checked weighted Banach map. This is stronger than a formal Taylor jet and does not require convergence of that jet. It is still an initial-constraint result: it proves neither evolution through shell contact nor a nonlinear shell spacetime theorem. It does not claim that the family is spatially C-infinity.

The precise linear input is Delay's global weighted isomorphism, not merely his existence theorem for individually compactly supported sources. The mixed orders two and one fit the required weighted products in four spatial dimensions. The fixed scalar and scalar-momentum seed fields introduce lower-order metric dependence and compactly supported parameter dependence, which do not spoil those products.

This conclusion remains conditional on the checkpoint's verified local KID classification and the corrected API/KRC hypotheses. An external novelty claim is not justified.

## Primary-source facts checked

1. **Delay, arXiv:1003.0535v4**, [primary text](https://arxiv.org/html/1003.0535). Theorem 3.8 states a global isomorphism
   `K-perp ∩ H̊_{φ,ψ}^{k+2m} → K-perp ∩ H̊_{φ,ψ}^k`
   for the projected operator `L=ψ^-2 P ψ²φ^{2m}P*`, assuming API. This includes noncompact-in-the-open-domain weighted sources. Section 2 defines `H̊` by completion of compactly supported fields. Sections 4 and 5 control physical zero extensions. Its nonlinear illustration is **§9.6** (a Weyl-connection equation); §10 is the scaling appendix. It is not an Einstein–scalar nonlinear theorem.
2. **Chruściel–Delay, “On mapping properties of the general relativistic constraints operator in weighted function spaces, with applications,” Mémoires de la SMF 94 (2003)**, [DOI:10.24033/msmf.407](https://doi.org/10.24033/msmf.407); [primary PDF](https://www.numdam.org/item/10.24033/msmf.407.pdf). Theorems 3.7/3.9 are conditional on differentiability of the nonlinear weighted map; Remark 3.8 explicitly requires weight-specific checking. Theorem 5.9 proves the compact-boundary exponential-weight vacuum construction using weighted Moser/scaling estimates with `k>n/2`. Proposition 5.10 and Corollary 5.11 give the corresponding Hölder improvement and smooth zero extension for smooth input. Proposition 3.17 permits Hölder bootstrap; the earlier warning about insufficient bootstrap concerns Sobolev solutions. These vacuum results contextualize the argument below but are not cited as an automatic matter extension.

## Formulation and actual density convention

Choose one smooth concentric annulus `D ⋐ A` containing the compact support of the first-order datum `v`; its closure avoids both shell and vertex. Let `d` be a smooth, positive, radial boundary defining function of `D`. Fix `β>0`, a nonnegative integer `a`, and weights

```
φ = d²,             ψ = d^{2a−4} exp(−β/d) .
```

All derivatives of these weights measured at the scale `φ` have the usual uniform bounds and local comparability on balls of radius proportional to `φ`; use the scaling cover of Delay's appendix. Fixed smooth background coefficients and the reference volume have uniform rescaled bounds because the background is regular on the closed annulus.

Rescale the defining function if necessary so `B(p,φ(p))⊂D`, as required by Delay's equation (10.4); equations (10.1)–(10.6) give the multiplier and comparability statements used below.

Normalize the momentum output by the **fixed** background density `√det γ₀`. Write

```
Ĉ = (H, M/√det γ₀).
```

Use fixed background bundle identifications when converting a symmetric tensor to the contravariant momentum density. Let `P_g` be the metric-only Hamiltonian derivative and let `P_m` be the actual normalized momentum derivative on symmetric tensors. Its adjoint is the Killing operator with the normalization implied by `M=-2 div π+ϖ dφ`; this normalization must not be replaced by an arbitrary convention for “symmetrized gradient.”

For invariant potentials `u=(N,Y)`, define

```
h(u)       = ψ² φ⁴ P_g* N ,
δπ(u)      = √det γ₀ · ψ² φ² P_m* Y ,
A u        = (h(u), 0, δπ(u), 0) .
```

The tensor/index conversion in the second line is through `γ₀`. Keep matter fixed at the seed

```
φ_ε = Φ + ε v_φ,           ϖ_ε = ε v_ϖ ,
x_seed(ε) = x₀ + εv .
```

The nonlinear map is

```
F(ε,u) = ψ^-2 Ĉ(x_seed(ε)+A u).
```

At `(0,0)), its derivative in the potentials is exactly the block-diagonal weighted operator `L_g ⊕ L_m` from Delay's Theorem 3.8.

## The Banach spaces and absence of a residual projection

For a fixed integer `k≥10`, take

```
X = H̊_{φ,ψ}^{k+4}(scalar) ⊕ H̊_{φ,ψ}^{k+2}(vector),
Z = H̊_{φ,ψ}^{k}(scalar)   ⊕ H̊_{φ,ψ}^{k}(covector),
```

with their closed `Γ≅(Z₂)^4` invariant subspaces. The radial domain, weights, background volume and operators commute with `Γ`.

The scalar adjoint kernel is spanned by `Ωx_a`, carrying characters `σ_a`; the momentum adjoint kernel is spanned by rotations, carrying `σ_aσ_b`. No such character is trivial. Hence invariant elements are orthogonal to the full kernel in `L²_ψ`, and the projection in Theorem 3.8 acts as the identity on the invariant range. The unique inverse on the kernel complement commutes with the group. Thus

`L_g ⊕ L_m : X^Γ → Z^Γ`

is a bounded isomorphism. This is a statement about the derivative at the fixed background; it does not require classifying KIDs at every perturbed datum.

## Explicit mixed-order weighted estimates

Use a countable scaling cover `B_i`, with bounded overlap, radius comparable to `φ_i`, and locally comparable `ψ_i`. In four dimensions the weighted norms are equivalent to square sums of the rescaled local norms

```
V_H,i = ψ_i φ_i² ||Ñ_i||_{H^{k+4}},
V_M,i = ψ_i φ_i² ||Ỹ_i||_{H^{k+2}} .
```

A tilde denotes pullback to a fixed-size ball. The fixed adjoint correction formulas give

```
||h̃_i||_{H^{k+2}}   ≤ C ψ_i V_H,i ,
||δπ̃_i||_{H^{k+1}} ≤ C ψ_i φ_i^-1 V_M,i .
```

Derivatives of the weights do not lose these factors because the scale is `d²). Smooth density/index conversions do not change them.

Since `k>4/2`, the local Sobolev product and composition estimates apply uniformly. The highest derivative quadratic terms in the Hamiltonian and momentum remainders have the forms

```
Q_H : h D²h, (Dh)², (δπ)² ,
Q_M : h Dδπ, δπ Dh .
```

Here `D` is an unscaled spatial derivative. On each ball every one has local source norm bounded by

`C ψ_i² φ_i^-2 (V_H,i+V_M,i)²`.

Normalizing by `ψ_i²` and multiplying by the output norm factor `ψ_i φ_i²` gives

```
W_i(ψ^-2 Q) ≤ C ψ_i (V_H,i+V_M,i)² .
```

Consequently

```
||ψ^-2 Q(u)||_Z ≤ C ||u||_X² ,
||ψ^-2(Q(u)−Q(w))||_Z
    ≤ C (||u||_X+||w||_X) ||u−w||_X
```

on a sufficiently small ball. The sequence estimate uses bounded `sup ψ_i` and
`||V²||_{ell²} ≤ ||V||_{ell-infinity} ||V||_{ell²} ≤ ||V||_{ell²}²`.
Higher derivatives of the map are controlled by the corresponding multilinear product estimates.

The rational factors from `γ^-1` and `1/det γ` are smooth Nemytskii maps on the uniform positive-metric neighborhood. They alter the constants but not the derivative orders or weight factors. First-order seed terms give an `O(|ε|)||u||_X` contribution on their compact support, where all weights are bounded above and below.

Matter terms are specifically checked rather than regarded as arbitrary prescribed stress-energy:

- `|Dφ_ε|_γ²` depends algebraically on `γ^-1` and on a known smooth one-derivative scalar field. Its linear metric part is already in `P_g`; the remainder has no extra metric derivative.
- `U(φ_ε)` is independent of the unknown corrections and differs from the background only on the seed support.
- `ϖ_ε²/det γ` is algebraic in `γ` and compactly supported in the seed region.
- The normalized term `ϖ_ε dφ_ε` is a known compactly supported source; the momentum's gravitational part has only the mixed terms listed above.

All background and potential derivatives needed for smooth parameter dependence are available under the manuscript's smoothness assumption. Therefore `F` is a smooth map between the stated fixed Banach spaces near `(0,0)`. The estimates concern the actual Einstein–scalar constraint formulas, not just a vacuum citation.

For the dotted target space, first define the map on compactly supported smooth potentials. Its constraint output is compactly supported, because the seed also is. The uniform product/composition estimates extend this map to the weighted completion; the output therefore belongs to the dotted target closure, rather than merely having finite weighted norm.

## Proper zero extension and physical regularity

The physical correction map has bounded zero-extension mappings

```
A_H : X_H → H^{k+2}(R⁴),
A_M : X_M → H^{k+1}(R⁴).
```

Indeed, rescaling an ordinary derivative back to the physical ball multiplies the local norm by powers of `φ_i^-1`; in dimension four both resulting estimates are bounded by

`C sup_D(ψ φ^-k) ||u||_X`.

This supremum is finite because `β>0`; the exponential dominates every power at the boundary. Approximate the potentials in `H̊` by compactly supported smooth potentials, symmetrizing the approximants if necessary. Their physical corrections extend by zero; the above bounds make those extensions converge in the indicated unweighted spaces. Thus the limiting corrections have genuine Sobolev zero extensions, without distributional boundary jumps. Their support is contained in the closed annulus.

For `k=10`, these spaces imply substantially more than `C^4` spatial regularity by four-dimensional Sobolev embedding. A claim of spatial `C^∞` is unnecessary and is not made here.

## Implicit-function consequence

The Banach implicit-function theorem now applies to `F(ε,u)=0`: `F(0,0)=0` and its potential derivative is a bounded isomorphism on the invariant spaces. It supplies one interval `|ε|<ε₀` and a smooth map `ε↦u(ε)∈X^Γ`.

Because the seed solves the linearized constraints, `∂_ε F(0,0)=ψ^-2 Pv=0`; hence `u(0)=0` and `u′(0)=0`. The exact initial-data family

`x(ε)=x₀+εv+A u(ε)`

therefore has tangent `v`, solves the full nonlinear constraints, and differs from the background only on the fixed closed annulus. Its metric stays positive for a sufficiently small parameter interval. Its parameter dependence is smooth into the physical finite Sobolev spaces on that same interval.

For two copies one can use identical families or the parameter-reversed family `x^−(ε)=x^+(−ε)`. Keeping all shell initial variables at their background values gives exact initial junction compatibility because the bulk family equals the background on a fixed shell collar. This statement still does not supply an evolved material-interface solution.

## What should remain open

- Spatial `C^∞` regularity and smooth dependence in every spatial seminorm require a separate mixed-order Hölder bootstrap/parameter argument; the cited vacuum results suggest the route but do not automatically establish the scalar extension.
- Evolution through shell contact requires the coupled hyperbolic and junction initial-boundary problem.
- Any claim of global nonlinear spacetime integrability, stability, cosmological origin, external novelty or readiness for submission remains unsupported by this initial-data result.

No new numerical or symbolic tests were run for this assessment. It is a functional-analytic proof check and source audit, not an empirical PDE solution.
