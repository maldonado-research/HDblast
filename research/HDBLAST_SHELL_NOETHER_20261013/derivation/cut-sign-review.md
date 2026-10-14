# Independent check of the GHY corner and tangent Noether cut signs

This calculation uses a fixed timelike hypersurface, an outward spacelike unit normal `n`, `K_ab = h_a^A h_b^B ∇_A n_B`, and the Einstein normalization `L_EH = (1/2) R ε`. The metric variation below is covariant: `q_AB = δg_AB`, while `q^AB = g^AC g^BD q_CD` has **no** extra minus sign. Let `ε_Σ = i_n ε`. All identities are exact at the nonlinear background; no Gaussian gauge is imposed on the variation.

## 1. Fixed-boundary variation

Gaussian normal coordinates for the background at the boundary suffice to check the signs. Write

```
g = dr² + h_ab(r,x) dx^a dx^b,
n = ∂_r,
λ = q_rr,       v_a = q_ra,       q = h^ab q_ab.
```

The varied unit normal is `δn^r = −λ/2`, `δn^a = −v^a`. Direct calculation gives

```
δK = ½ ∂_r q − ½ λ K − D_a v^a,
θ_grav|_Σ = ½ [D_a v^a − K^ab q_ab + K λ − ∂_r q] ε_Σ,
δ(K ε_Σ) = [½ ∂_r q − ½ λ K − D_a v^a + ½ K q] ε_Σ.
```

Consequently, with `p = n^A ∇_A φ`,

```
θ|_Σ + δ(K ε_Σ)
  = −½ (K^ab − K h^ab) q_ab ε_Σ − p δφ ε_Σ + d_Σ β,
β(δ) = −½ i_v ε_Σ,
v^a = h^ab n^A q_Ab.
```

In particular, the proposed minus sign and factor `1/2` in `β` are correct. Neither `q_ra` nor `q_rr` may be set to zero before checking this identity.

## 2. The extension-independent tangent charge

Let `ξ` be tangent to the fixed boundary, and define the tangent vector

```
A^a = h^a_B n^A ∇_A ξ^B,
(K·ξ)^a = K^a_b ξ^b.
```

Tangency implies

```
h^a_B n_A ∇^B ξ^A = −(K·ξ)^a.
```

The normal-tangent part of `dξ♭` is therefore `A + K·ξ`, whereas the normal-tangent part of `L_ξ g` is `A − K·ξ`. The Hodge identity in the stated orientation is

```
(*dξ♭)|_Σ = i_(A + K·ξ) ε_Σ.
```

For `Q_ξ = −½ *dξ♭`, this yields

```
Q_ξ|_Σ = −½ i_(A + K·ξ) ε_Σ,
β(L_ξ g) = −½ i_(A − K·ξ) ε_Σ,

Q_ξ|_Σ + i_ξ(K ε_Σ) − β(L_ξ g)
  = i_[(K δ^a_b − K^a_b) ξ^b] ε_Σ
  = −i_(T·ξ) ε_Σ,
T^a_b = K^a_b − K δ^a_b.
```

All dependence on the normal extension of `ξ` cancels. A Killing assumption is unnecessary. Minimal scalar coupling adds no separate bulk Noether charge of this form.

## 3. Two sides, one tension, and orientation bookkeeping

Use the same induced metric and scalar on the two copies of the shell. If all boundary forms have already been represented in one common outward-oriented boundary convention, put

```
ℓ = (K_+ + K_− − σ) ε_Σ,
β = β_+ + β_−,
I_ab = (K_+ab − K_+ h_ab) + (K_−ab − K_− h_ab) + σ h_ab,
P = p_+ + p_−.
```

There is one shell action `−σ ε_Σ`, not one tension term per bulk region. Then

```
Σ_i θ_i|_Σ + δℓ = F + d_Σ β,
F(δ) = −[½ I^ab δh_ab + P δφ] ε_Σ,

B_ξ := Σ_i Q_iξ|_Σ + i_ξ ℓ − β(L_ξ)
     = −i_(I·ξ) ε_Σ.
```

For a shell potential `σ(φ)`, the same formulas hold with `P = p_+ + p_− + σ'(φ)` and the nonlinear `I_ab` shown above.

The qualification about common orientation is substantive. If both bulk regions inherit a single glued-manifold orientation, their raw outward-induced forms pull back with opposite signs: `j_i^* ε_Σi = s_i ε_Σ`, where `s_i = ±1`. The physical sum written on the common oriented shell uses `s_i j_i^*` on **every** boundary form (`θ_i`, `Q_i`, `β_i`, and `K_i ε_Σi`). Thus `s_i j_i^*(K_i ε_Σi) = K_i ε_Σ`. One must not add raw pullbacks with opposite volume forms and then insert an unsigned Israel sum. Equivalently, use side-specific bulk orientation conventions that identify both outward boundary forms with the same abstract `ε_Σ` from the outset.

The cut orientation also matters. For a future timelike vector `u` tangent to the shell and normal to a spacelike cut,

```
i_u ε_Σ = i_u i_n ε = −i_n i_u ε.
```

Thus the future-end orientation of the timelike shell is opposite to the outward-boundary orientation of a future-oriented spacelike slice. The form identities above remain valid; integration must retain this orientation reversal.

## 4. Exact variation and the corner correction

Fix `ξ` in field space and keep it tangent to the fixed shell. Define

```
k_ξ(δ) = Σ_i [δQ_iξ − i_ξ θ_i(δ)]|_Σ,
C_ξ(δ) = (δ_field β)(δ,L_ξ)
       = δ[β(L_ξ)] − L_ξ[β(δ)].
```

The second equality uses covariance and fixed `ξ`. Cartan's identity gives the exact algebraic relation

```
k_ξ − C_ξ − d_Σ(i_ξ β(δ)) = δB_ξ − i_ξ F(δ).
```

This statement is stronger and more precise than setting the bulk `δQ_ξ − i_ξ θ` to zero by itself: the field-space variation of the GHY corner is required. Expanding the right side gives

```
δB_ξ − i_ξ F
 = −i_[(δI^a_b) ξ^b] ε_Σ
   − ½(h^cd δh_cd) i_(I·ξ) ε_Σ
   + [½ I^ab δh_ab + P δφ] i_ξ ε_Σ.
```

At a nonlinear shell solution `I = P = 0`,

```
k_ξ − C_ξ − d_Σ(i_ξ β) = −i_[(δI^a_b) ξ^b] ε_Σ.
```

For a tangent variation satisfying the full linearized junction conditions, `δI = δP = 0`, the corrected local cut form vanishes exactly. In particular, the index-raising and volume-form variations vanish because they multiply the *background* residual `I`; they must not be discarded before imposing the nonlinear shell equation. On a closed cut, or with compact support away from its edge, integration removes the remaining `d_Σ(i_ξ β)` term. With an edge, that term must be retained as an edge contribution.

To turn this cut identity into presymplectic degeneracy one also uses the bulk equations, their linearization, and the corresponding Noether constraint identity. It does not by itself establish degeneracy for a vector that acts at an exterior boundary or changes prescribed exterior data. If `ξ` is field dependent, the fixed-parameter formulas require the usual subtraction of the charge associated with `δξ`; the derivation above deliberately avoids that complication.

## 5. Regularity and scope

A sufficient classical setting is a smooth timelike shell, one-sided `C²` metrics and scalars up to it, matching smooth induced fields, `C²` variations obeying the linearized matching and junction conditions, and a `C²` tangent vector field with suitable cut support. The basic boundary identities only use first traces of the fields and variations, but the bulk equations and Noether-current argument require the corresponding higher derivatives or a justified weak formulation. A bulk finite-energy estimate alone does not supply all normal traces, corners, or products appearing here. At Sobolev regularity, prove a trace/approximation extension before invoking these identities. The calculation introduces no horizon, infinity, or finite-radius endpoint condition.

Verdict: the proposed `β = −½ i_v ε_Σ` is correct, the completed tangent charge is `−i_(I·ξ) ε_Σ`, and exact cut cancellation requires its variation, the field-space corner term, and explicit treatment of cut edges and orientations.
