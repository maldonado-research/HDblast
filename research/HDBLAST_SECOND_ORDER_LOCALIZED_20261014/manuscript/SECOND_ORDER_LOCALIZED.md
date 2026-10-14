---
title: Localized second-order constraint solvability for Einstein–scalar shell perturbations, and an areal-radius monopole identity
author: Ricardo Maldonado
affiliation: Independent Researcher
orcid: 0009-0009-3937-6527
email: ricardomaldonadoandrade@gmail.com
date: 14 October 2026 (America/Los_Angeles) — technical draft for author review
---

**Status.** Technical continuation of the 13 October shell-Noether checkpoint
[1], for author review. Not an APS submission and not external peer review.
All earlier capsules are unchanged. The [subsequent reflection theorem](../continuation/REFLECTION_FORMAL_INITIAL_DATA.md) proves formal initial-data solvability through every finite order in a symmetry sector; the [pre-contact note](../continuation/PRECONTACT_SECOND_ORDER_DYNAMICS.md) supplies finite-slab spacetime coefficients. This manuscript retains its localized second-order scope.

## Abstract

The shell-Noether checkpoint [1] proved that ten Lorentz (Taub-type) quadratic
conditions are *necessary* for a linearized Einstein–scalar shell perturbation
to be tangent to an exact solution family, and that within a fixed-profile
quadrupole family they reduce to `[A,B]=0`. It left open whether they are also
*sufficient* for the second-order initial constraint problem, with the
interface adjoint domain flagged as the obstacle. For every first-order
datum supported in an annulus that avoids the shell and the cone vertex, we
settle this initial-data question within corrections that vanish near the
shell. For fixed reflection parity the ten global tests are then necessary
and sufficient. On such an annulus of the time-symmetric slice, the
local Killing initial data (KIDs) are exactly the four boost lapses `Ω x_a`
and the six rotations. A localization theorem of Delay [2] then gives
compactly supported second-order data if and only if the ten per-copy Taub
integrals vanish. Because the correction vanishes near the shell, every
second-order junction condition holds with zero second-order shell
displacement, and the interface adjoint problem never arises. For the focused
commuting state `A=B=diag(1,−1,0,0)` the correction is written in closed form,
up to one explicit radial quadrature, and checked by exact symbolic and
exact-rational tests and by an independent numerical code. Finally, an exact
identity with the areal-radius weight `Ωr` shows that a monopole moment `m`
of every compactly supported correction is fixed by the first-order packet,
once the packet is written in the standard reconstruction gauge. For
velocity-only packets in that gauge `m` is strictly negative. The spherically
symmetric part of the correction can then be neither a compactly supported
gauge mode nor a scalar correction. `m` is not invariant under first-order
gauge changes of the packet, so this is a gauge-fixed statement, not a
physical observable. Second-order evolution through shell contact and nonlinear existence remain open; the subsequent note addresses finite-order formal initial data in the reflection-invariant sector.

## 1. Setting

The action, conventions and background are those of [1] and its inputs:
five-dimensional Einstein–scalar gravity with `κ₅ = 1`,
`S = ∫√−g [R/2 − (∇φ)²/2 − U(φ)]`, two regular interiors glued at one
timelike shell with scalar-dependent tension. The fixed regular background is
`g₀ = Ω(s)² η`, `φ₀ = Φ(s)`, `s = R² − T²`, with `Ω > 0` and `Φ_s ≠ 0`.

On the `T = 0` slice of each copy the data are time symmetric:
`γ₀ = Ω(r²)² δ` on the ball `r < b` (shell at `r = b`), `φ₀ = Φ(r²)`,
`π₀ = 0`, `ϖ₀ = 0`. The constraint map is
```
H   = R/2 + (2/det γ)(tr π²/3 − |π|²) − ϖ²/(2 det γ) − |Dφ|²/2 − U(φ),
M_i = −2 D_j π^j_i + ϖ ∂_i φ          (π, ϖ densitized).
```
At a time-symmetric background the linearization decouples:
`P_H(h, ψ) = DH[h, ψ]` does not involve momenta, and
`P_M(π′, ϖ′) = −2 div π′ + ϖ′ dφ₀` does not involve `(h, ψ)`.

Along `x(ε) = x₀ + εv + ε²w/2`, the second-order constraint equation is
```
P[w] = S := −D²C[v, v],      P = (P_H, P_M).
```
First-order data use the inherited ADM reconstruction (eq. (22) of [1]):
for master data `q`, `p`,
`h = 2q(sδ − 3xxᵀ)`, `δφ = −3(Dq + 4q)/(2Ω²Φ′)`,
`δπ = −3p xxᵀ/(2Ω)`, `δϖ = −3Ω(Dp + 5p)/(2Φ′)`, with `D = x·∇`.
For quadrupoles, `q = G(s) xᵀAx` and `p = K(s) xᵀBx` with arbitrary radial
profiles `G`, `K`. The frozen family of [1] is `K = ωG`.

Throughout, `𝒜 = {r₁ < r < r₂}` with `0 < r₁ < r₂ < b`, and first-order data
are supported in a compact subset of `𝒜` in each copy and solve the
linearized constraints.

## 2. Local KIDs on an annulus

**Lemma 1.** On `𝒜`:

(a) The Killing fields of `γ₀` are exactly the rotations `so(4)`.

(b) The kernel of the metric-only adjoint
`P_g* N = ½(Hess N − γ₀ ΔN − N Ric) + ½ N dφ₀⊗dφ₀`,
the adjoint of `P_g(h) = P_H(h, 0)`, is exactly `span{Ω x_a : a = 1..4}`.

(c) The ten are full KIDs: they are the `T = 0` data of the background
boosts and rotations, so `Ωx_a` also annihilates the scalar component
`P_ψ* N = div(N∇φ₀) − N U′(φ₀)`, and the rotations preserve `φ₀`. This
inclusion is used for the necessity half of Theorem 2. By (a) and (b) there
are no other KIDs.

*Proof of (a).* `γ₀` is conformally flat, so by Liouville's theorem a Killing
field on the connected set `𝒜` is the restriction of a conformal Killing field
of `R⁴`: `Y = a + λx + 2(b·x)x − s b + (rotation)`. It is Killing for
`Ω²δ` exactly when `Y(ln Ω²) = −½ div Y`. Separating angular degrees gives
`λ(2sΩ′/Ω + 1) = 0` and `(Ω′/Ω)(a + s b) + b = 0`. A nonzero `λ` forces
`Ω = c s^{−1/2}` with `c > 0`, so `(1/Ω)″ = −s^{−3/2}/(4c) < 0`. A nonzero
`(a, b)` forces `Ω = C/(a_i + s b_i)`, so `1/Ω` is linear in `s` and
`(1/Ω)″ = 0`. Both contradict the background equation
`(1/Ω)″ = (Φ′²/3)(1/Ω) > 0`, which holds because `Φ′ ≠ 0`. ∎

*Proof of (b).* `P_g*` commutes with `SO(4)`, so its kernel is a sum of
isotypic pieces `n(s) ⊗ V_l`. For `l ≥ 2`, take the part of `P_g* N` that is
tangential and trace-free. `Ric` and `dφ₀⊗dφ₀` are combinations of `δ` and
`xxᵀ`, and the Christoffel terms of a radial conformal factor are tangentially
pure trace. Hence that part is `½ n(s) r^{l−2} (∇²_{S³}Y_l)^{TF}`, whose
squared `L²(S³)` norm is `(2/3)λ_l(λ_l − 3)∫Y_l² > 0` for `λ_l = l(l+2)`,
`l ≥ 2`. So `n = 0`. For `l = 1`, `N = n(s)x₁`, the mixed component of
`P_g* N` is `−(nΩ′ − n′Ω)/Ω`, so `(n/Ω)′ = 0` and `N = cΩx₁`, a boost lapse.
For `l = 0`, the two radial components give the first-order relation
`n′ = λn`, `λ = 2Ω′(Ω + sΩ′)/(Ω(Ω + 2sΩ′))`. Compatibility with the
second-order component leaves the residual `−Φ′²Ω² n/(3(Ω + 2sΩ′)²)`, so
`n = 0` wherever `Ω + 2sΩ′ ≠ 0`. Where `Ω + 2sΩ′ = 0`, the eliminated
relation reads `nΩ′(Ω + sΩ′) = 0` with `Ω′ = −Ω/(2s) ≠ 0`, so again `n = 0`.
∎ The `l = 0, 1` reductions are checked symbolically for all `s`
(`new/kids_symbolic_l01.json`), and independently by the referee's own code.
The `l = 2, 3, 4` rank statements are checked at two sample points
(`new/kids_tensor_only.json`); the general `l ≥ 2` case rests on the
trace-free Hessian argument above.

## 3. Localized solvability

**Theorem 2.** Let `v` be smooth with support in a compact `K ⊂ 𝒜` in each
copy, solving the linearized constraints, and `S = −D²C[v, v]`. There is
smooth second-order data `w` with compact support in `𝒜` and `P[w] = S` if
and only if in each copy
```
∫_𝒜 Ω x_a S_H dμ_γ₀ = 0   (a = 1..4),      ∫_𝒜 Y_J · S_M d⁴x = 0   (six rotations J).
```
The solution can be chosen with `ψ = 0` and `ϖ′ = 0`.

*Proof.* Necessity: integrate by parts against the full KIDs of
Lemma 1(c); boundary terms vanish because `w` has compact support. This uses
arbitrary `w`, including `ψ ≠ 0` and `ϖ′ ≠ 0`. Sufficiency:
solve the two decoupled equations separately. (i) `P_M(π′, 0) = −2√γ₀ div τ`
for `π′ = √γ₀ τ` is first order, and its adjoint is the Killing operator.
(ii) `P_g` is second order. Choose a smooth annulus `𝒜′` with
`K ⊂ 𝒜′ ⋐ 𝒜`; Lemma 1 holds on every annulus. The scalar Hessian map
`A ↦ ½(A − γ₀ tr_γ₀ A)` is invertible on symmetric two-tensors in dimension
four (trace multiplication by `−3/2`), so `P_g*` has the full-derivative
injectivity required by Delay's Lemma 8.1 and satisfies his estimate (API).
The Killing operator has injective symbol, but symmetrization of a full
first-derivative tensor has an antisymmetric kernel, so Lemma 8.1 alone does
not establish (API) for the momentum operator. Instead use Delay's §9.3,
which records (API) for divergence on symmetric two-tensors, citing
Chruściel–Delay [5]. By
Lemma 1, both local adjoint kernels are spanned by smooth boosts or rotations
on a neighbourhood of the closed annulus, hence kernel regularity (KRC)
holds (also Delay's Corollary 8.4 for the Killing operator). His Theorems 1.3
and 5.1 then give smooth solutions on `𝒜′`, vanishing to infinite order at
`∂𝒜′`, whenever the source is `L²`-orthogonal to the adjoint kernel on `𝒜′`.
Extended by zero, they have compact support in `𝒜`. By Lemma 1 the
orthogonality conditions are exactly the ten displayed integrals. Theorem and
lemma numbers follow arXiv:1003.0535v4, the version published in Comm. PDE 37
(2012). ∎

**Corollary 2.1 (interface and vertex).** Assume in addition that the
first-order shell displacement and tilt, and their time derivatives, vanish
at `T = 0` (as for data vanishing near the shell). The solution of Theorem 2
vanishes in a neighbourhood of the shell and of the vertex. So do the
first-order data and every quadratic interface source built from them. Every
second-order junction and compatibility condition at `T = 0`, including
induced matching, the Israel and scalar junctions and the shell's normal
force, therefore holds with zero second-order shell displacement and tilt. The interface adjoint-domain
problem identified in [1] does not arise for this class.

**Corollary 2.2 (relation to the global tests).** For fixed reflection
parity the two copies carry equal quadratic sources, so each per-copy
integral is half the corresponding global one. The ten global tests of [1]
are then necessary *and* sufficient for the second-order initial-data problem
in this localized class. For two copies with different packets the localized
method needs ten conditions per copy, while only the ten sums are known to be
necessary. That gap is not resolved here.

**Corollary 2.3 (quadrupoles).** For `q = G H_A`, `p = K H_B` with arbitrary
profiles supported in `𝒜`, the boost moments vanish sphere by sphere by
parity. `O(4)` equivariance forces the per-sphere rotational moment to be
`c(s; G, K)·tr(M[A, B])`: the symmetric square of the degree-two harmonics
contains no `so(4)` component, and reflections tie the self-dual and
anti-self-dual parts together. Hence `[A, B] = 0` implies compactly supported
second-order initial data for every such profile pair. Combined with [1]:
for fixed-parity quadrupoles with `K = ωG`, `ω ≠ 0` and `B_rad(G, G) ≠ 0`
(for example thin support, as in [1]), `[A, B] = 0` is necessary and
sufficient for second-order initial data in this localized class.

## 4. Closed form for `A = B = diag(1, −1, 0, 0)`

Write `Q = (x₁² − x₂²)²`, `P = sδ − xxᵀ` (the tangential projector times `s`),
and
```
Q_rot(s) = 4Φ′²GKs² + 63GK + 18GK′s + 21G′Ks + 6G′K′s².
```
The same polynomial gives the rotational moment of the non-commuting witness.
The correction `w = (h_T + h₀, Ψ, π′, ϖ′)` is:

* **Momentum.** `π′ = β(s) Q P` with `β = 3Q_rot/(4Φ′²Ω³s)`, and
  `ϖ′ = ν(s) Q`. Both are explicit rational functions of the background and
  profile jets (`basis/results/construct_cc_M.json`). Because
  `Px = 0`, no radial derivative of `β` enters `P_M`. The equation is solved
  sphere by sphere, with no ODE.
* **Hamiltonian, `l ≥ 1` part.** `h_T = [s P (∇∇χ) P − tr(P(∇∇χ)P) P/3]`
  with `χ = C₁(s)(x₁² + x₂²) + C₂(s) Q`. This tensor is tangential and
  trace-free, so no radial derivative of `C₁` or `C₂` survives. It solves
  `P_H[h_T] = S_H − f₀`, where `f₀` is the `S³` average of `S_H`.
* **Hamiltonian, `l = 0` part.** `h₀ = a(s) xxᵀ`, `ψ = Ψ(s)`. Exact
  identities, valid for every background:
  ```
  s^{3/2}Ω⁵ P_H(a xxᵀ, 0) = d/ds[3s^{5/2}(Ω + 2sΩ′)a] − (3/2)s^{3/2}Ω a,
  s^{3/2}Ω⁵ P_H(0, Ψ)     = −4 d/ds[s^{5/2}Ω³Φ′Ψ].
  ```
  Take `a = κ b(s)` with `b ∈ C_c^∞((s₁, s₂))`, `b ≥ 0`, `b ≢ 0`, where
  `[s₁, s₂] ⊂ (r₁², r₂²)` contains `supp f₀`, and
  `κ = −(2/3)∫s^{3/2}Ω⁵f₀ ds / ∫s^{3/2}Ω b ds`. Then
  `4s^{5/2}Ω³Φ′Ψ = 3s^{5/2}(Ω + 2sΩ′)a − ∫_{s₁}^{s}[σ^{3/2}Ω⁵f₀ + (3/2)σ^{3/2}Ωa]dσ`
  is smooth and vanishes outside `(s₁, s₂)`.

Every coefficient is a product of profile jets, so `w` has the packet's
support, up to the chosen bump.

## 5. An areal-radius monopole identity

**Theorem 3.** Let `w` be any compactly supported solution of `P[w] = S` on
`𝒜`, and let `h_{l=0} = α(s)δ + a(s)xxᵀ` be the spherically symmetric part of
its metric. Then
```
m[w] := ∫ Ω√s (3α + s a) ds = −(2/3) ∫ s^{3/2} Ω⁵ f₀ ds,
```
where `f₀` is the `S³` average of `S_H`. The functional `m` does not involve
`ψ`. It is invariant under compactly supported linearized gauge
transformations of `w` at fixed `v`, including radial ones `ξ = ζ(s)x`. It
does depend on the gauge in which `v` is written: it is fixed once `v` is
given in the reconstruction gauge (eq. (22) of [1]).

*Proof.* Let `N = Ωr = Ω√s`. A direct computation gives
`P_ψ* N = 0` and `P_g* N = −Ω(2sδ + xxᵀ)/(2s^{3/2})`, so
`⟨P_g* N, αδ + a xxᵀ⟩_γ₀ = −3(3α + sa)/(2Ω³√s)`. For compactly supported `w`,
`∫ N P_H[w] dμ = ∫⟨P_g* N, h⟩ dμ + ∫ (P_ψ* N) ψ dμ`, with no boundary terms.
Only the `l = 0` parts survive because `N` is radial. With
`dμ = (s/2)Ω⁴ ds dΩ₃`, the two sides give
`(1/2)∫s^{3/2}Ω⁵f₀ ds = −(3/4) m[w]`. The same reduction, done by integrating
the exact `l = 0` operator by parts, is checked symbolically
(`new/l0_full_sector.json`) and numerically (`independent/indep_l0.json`).
Invariance under linearized gauge changes of `w` follows from
`P_H(L_ξ γ₀, ξφ₀) = 0`, also checked directly. ∎

*Gauge dependence on `v`.* Under a first-order gauge change
`x̃(ε) = exp(ελ L_ξ) x(ε)`, the correction becomes
`w̃ = w + 2λ L_ξ v + λ² L_ξ² x₀`. For radial `ξ = ζ(s)x` the referee
computed
`m[L_ξ²γ₀] = ∫[8s^{5/2}Ω³ζ′² + 8s^{3/2}Ω(sΩΩ″ − sΩ′² + ΩΩ′)ζ²] ds`, which is
nonzero in general (`+6.34×10⁻³` on the test background with a bump on
`(0.5, 1.5)`). For a velocity packet `L_ξ v` has no metric part, so a
suitable `λ` makes `m[w̃] = 0` for the same physical packet in another gauge.
Theorem 3 is therefore a statement about the gauge-fixed packet.

*Remark.* `Ωr` is the areal radius of the background spheres. Up to scale it
is the only radial lapse annihilating `P_ψ*`, which is why `m` contains no
`ψ`. It is also the lapse, at `T = 0`, of the radial boost `R∂_T + T∂_R`.
That field preserves `s`, hence `Φ` and `Ω`, and its deformation
`L_K g = 2RTΩ² dΩ²` vanishes at `T = 0`; it fails to be a KID only through
the metric component of the adjoint.

**Corollary 3.1 (velocity packets, reconstruction gauge).** For `G = 0` the
weighted density is, exactly,
```
s^{3/2}Ω⁵ f₀ = s^{7/2}[8Φ′²s²K² + 3(2sK′ + 7K)²] / (8Φ′²Ω),
```
which is pointwise non-negative and strictly positive where `K ≠ 0`. (For a
general velocity matrix `B` it scales with `tr B²`.) Hence, for velocity
packets in the gauge (22) of [1], `m[w] < 0` for every compactly supported
correction. In that gauge the spherically symmetric part of the second-order
metric correction cannot be a compactly supported linearized gauge mode, and
cannot be replaced by a scalar correction. A first-order gauge change of `v`
shifts `m`, so this is not a gauge-invariant physical statement.

**Corollary 3.2 (displacement and the frozen family).** The displacement part
reduces, modulo total derivatives, to `c₂G″² + c₁G′² + c₀G²` with
`c₂ = 6s^{13/2}/(Φ′²Ω) > 0`. The cross term between `G` and `K` vanishes,
because `H` is even in the momenta. For `K = ωG` the moment is
`m = −(2/3)(I_G + ω²I_K)` with `I_K > 0`. It is nonzero except for at most one
value of `ω²`. For `G` supported in a sufficiently thin interval (by the
Wirtinger inequalities `∫G″² ≥ (π/w)²∫G′²` and `∫G′² ≥ (π/w)²∫G²`), or
sufficiently rapidly oscillating, `c₂` dominates and `m < 0` for every `ω`.
The sign of `I_G` for broad profiles depends on the background and is not
determined here. All statements in this corollary are in the reconstruction
gauge.

## 6. Verification

Counts and outcomes are in `RESULTS.md`, filled only from completed
receipts. In brief:

* the original construction's symbolic checks (background, first-order
  constraints, sources, KIDs, `l = 0` identities, sphere-by-sphere solves);
* an exact symbolic proof that no radial derivative of any closed-form
  coefficient survives, and that `P_M[w] = S_M` (`new/fast_verify.py`);
* an exact-rational identity test of `P_H[h_T] = S_H − f₀` and
  `P_M[w] = S_M` at random points, with three wrong-formula controls
  (`new/exact_random_point.py`);
* the all-radius kernel analysis of Lemma 1(b);
* an exact symbolic verification of `P_H[h_T] = S_H − f₀`, reducing the
  independently derived angular operator modulo `|x|² = s` before
  substituting radial coefficients (`new/angular_exact.py` and the
  [companion derivation](ANGULAR_EXACT_VERIFICATION.md)); the source is
  regenerated by the original engine, so this does not constitute an
  independent symbolic derivation of `−D²C[v,v]`;
* an independent numerical check sharing no code with the original engine.
  It re-derives the background relations from the 5D Ricci tensor, uses an
  exact closed-form test background, differentiates every field numerically
  at 70 digits, carries exact `ε`-series and evaluates in 40-digit
  arithmetic. It also checks the `l = 0` weighted identity and gauge
  invariance (`independent/`);
* an independent AI referee report, with every point accepted
  (`reviews/`).

The companion session's original exact verifier (`basis/code/verify_cc.py`)
terminated without output and is not counted.

## 7. Prior work and scope

Solvability of an underdetermined elliptic equation with compact support,
exactly when the source is orthogonal to the local adjoint kernel, is Delay's
theorem [2]. It is part of the Corvino–Schoen localized gluing programme
[3–5], with recent support-preserving solution operators [6, 7]. Taub
conditions and linearization stability go back to Fischer–Marsden, Moncrief
and Arms–Marsden–Moncrief [8–10]. The contribution here is the application
to this shell model: the annulus KID classification, the explicit commuting
solution, and the monopole identity with its areal-radius weight. No
external novelty is certified.

**Not established:** second-order evolution through shell contact (the
second-order junction initial-boundary problem with quadratic interface
sources); higher orders outside the reflection-invariant formal construction; convergence of the formal series; packets that already touch the shell at
`T = 0`; nonlinear existence; stability; positivity of a physical energy; any
cosmological or Big Bang inference.

**Author and disclosure.** Ricardo Maldonado, Independent Researcher, ORCID
0009-0009-3937-6527. No external funding; no known competing interests.
Substantial AI assistance (Anthropic Claude) was used for derivations,
symbolic and numerical code, literature location and drafting. Internal
automated checks are not external peer review. The author retains
responsibility. The reviewed publication and separate follow-on proofs used OpenAI Codex with collaborating internal agents.

## References

1. R. Maldonado, *A rotational obstruction to nonlinear continuation of
   Einstein–scalar shell perturbations*, HDBLAST_SHELL_NOETHER_20261013
   (GitHub checkpoint, 13 Oct 2026).
2. E. Delay, *Smooth compactly supported solutions of some underdetermined
   elliptic PDE, with gluing applications*, Comm. PDE 37 (2012) 1689–1716,
   [arXiv:1003.0535](https://arxiv.org/abs/1003.0535).
3. J. Corvino, Commun. Math. Phys. 214 (2000) 137–189.
4. J. Corvino and R. Schoen, J. Differential Geom. 73 (2006) 185–217.
5. P. T. Chruściel and E. Delay, Mém. Soc. Math. Fr. 94 (2003).
6. Y. Mao and Z. Tao, *Localized initial data for Einstein equations*,
   [arXiv:2210.09437](https://arxiv.org/abs/2210.09437).
7. Y. Mao, S.-J. Oh and Z. Tao, Bogovskii-type solution operators for the
   constraint equations (2023).
8. A. Fischer and J. Marsden, Gen. Rel. Grav. 7 (1976) 915.
9. V. Moncrief, J. Math. Phys. 16 (1975) 493; 17 (1976) 1893.
10. J. Arms, J. Marsden and V. Moncrief, Ann. Phys. 144 (1982) 81.
