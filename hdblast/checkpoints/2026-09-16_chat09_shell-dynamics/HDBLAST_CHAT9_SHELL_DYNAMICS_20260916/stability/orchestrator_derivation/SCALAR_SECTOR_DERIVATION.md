# Scalar-sector perturbations of the registered de Sitter shell (orchestrator derivation)

16 September 2026. Status: analytic derivation + floating-point shooting. Not an interval certificate.
Independent of the three workflow derivations (A, B, C); used as a fourth cross-check.

## 1. Background

Five dimensions, signature (−++++), κ₅² = 1, `R_AB = ∂_Aφ∂_Bφ + (2/3)U g_AB`, `□φ = U_φ`,
`W = 1 − φ + φ³/3`, `U = W_φ²/2 − (2/3)W²`.

    ds² = dy² + ρ(y)² γ_μν dx^μ dx^ν ,   γ = unit dS₄ (R_μν[γ] = 3γ_μν),   φ = φ(y)
    ρ'² = 1 + ρ²(φ'²/12 − U/6),   φ'' + 4(ρ'/ρ)φ' = U_φ,   (ρ'/ρ)' = −1/ρ² − φ'²/3 .

`y` is proper distance from the regular first cone (ρ ≈ y). The Z₂ shell sits at `y_b` with
`ρ'/ρ = σ_t/6`, `φ' = −σ_t'/2`, `σ_t = 2W + t(1 + cφ)`. Brane Hubble rate `H = 1/ρ_b`.

Conformal coordinate `dz = dy/ρ`, `a = ρ`, `ℋ = a_z/a = ρ'`:

    6(ℋ² − 1) = φ_z²/2 − a²U,    ℋ² − ℋ_z − 1 = φ_z²/3,    φ_zz + 3ℋφ_z = a²U_φ .

## 2. Perturbations in generalized longitudinal gauge

    ds² = a²[(1 + 2ξ)dz² + (1 + 2ψ)γ_μν dx^μ dx^ν],   δφ = χ,   □_γ Y = μ² Y.

A 4D mode has physical mass `m² = μ²H²`.

**Traceless μν equation** (for harmonics with `(∇_μ∇_ν − ¼γ_μν□)Y ≠ 0`, i.e. `μ² ≠ −4`): `ξ = −2ψ`.

**Codazzi (zμ) equation.** The slices z = const have unit normal `n = N⁻¹∂_z`, `N = a(1+ξ)`, and umbilic
extrinsic curvature `K^μ_ν = k δ^μ_ν`, `k = (ℋ + ψ_z − ℋξ)/a`. `R_{nμ} = D_νK^ν_μ − D_μK = −3∂_μk` and
`R_{nμ} = (φ_z/a)∂_μχ` give

    ψ_z − ℋξ = −φ_z χ/3      ⇒      χ = −3(ψ_z + 2ℋψ)/φ_z .                         (C)

**Gauss (nn) equation.** For a spacelike normal `2G_nn = −R̃ + K² − K_μνK^μν` (checked on the background:
it returns `6(ℋ²−1)/a²`). With `R̃ = (12/a²)(1 − 2ψ) − (6/a²)□ψ` (4D conformal rescaling) and
`K² − K_μνK^μν = 12k²`:

    3□ψ + 12ψ + 12ℋ(ψ_z − ℋξ) = φ_zχ_z − ξφ_z² − a²U_φ χ .                           (H)

Eliminating χ with (C), ξ = −2ψ, and the background identities gives the master equation

    ψ_zz + (3ℋ − 2φ_zz/φ_z)ψ_z + [4ℋ_z − 4ℋφ_zz/φ_z + 6 + μ²]ψ = 0 .                 (M)

The flat-wall limit (drop the "+6", which is 6K with K = 1 the slice curvature) is the standard
DeWolfe–Freedman–Gubser–Karch / Giovannini equation. In proper distance:

    ψ'' + 2(ρ'/ρ − φ''/φ')ψ' + [ −(4/3)φ'² − 4(ρ'/ρ)(φ''/φ') + (2 + μ²)/ρ² ]ψ = 0 .   (M_y)

**Schrödinger form.** With `ψ = (φ_z/a^{3/2})u`, `g = φ_zz/φ_z`:

    −u_zz + V(z)u = μ²u,   V = g² + ℋg + (9/4)ℋ² − g_z − (5/2)ℋ_z − 6 .

At the cone (`z → −∞`, `ℋ → 1`, `g → 2`) `V → 9/4`: the continuum starts at `μ² = 9/4`, the same gap as the
tensor sector. Bound states have `μ² < 9/4` and the normalizable cone behaviour
`ψ ∝ y^α`, `α = 1/2 + sqrt(9/4 − μ²)`.

## 3. Shell boundary condition

Transform to Gaussian-normal gauge (`ξ̄ = 0`, no shift, shell at fixed `y_b`) with `ε^y`, `ε^μ = ∇^μ ε_∥`:
`∂_yε^y = ξ`, `∂_yε_∥ = −ε^y/ρ²`, `Ē = −ε_∥`, and the shell displacement in longitudinal gauge is `ζ = ε^y(y_b)`.

* Traceless Israel condition `∂_yĒ|_b = 0` ⇒ `ζ = 0` for every harmonic with `μ² ≠ −4`
  (no brane bending in longitudinal gauge; this is also Frolov–Kofman's statement).
* Trace Israel condition `∂_yψ̄ = σ_t'χ̄/6` becomes `ψ' − (ρ'/ρ)ξ = σ_t'χ/6`, which is identically the
  bulk constraint (C) once the background junction `φ' = −σ_t'/2` is used. It carries no new information.
* Scalar junction `∂_yχ̄ = −σ_t''χ̄/2` becomes the one non-trivial condition

      χ' + 2φ'ψ + (σ_t''/2)χ = 0     at y = y_b⁻ .                                  (BC)

Because χ contains ψ', (BC) involves ψ'' and therefore depends on μ² through (M_y).

## 4. Result for the registered shell

`scalar_spectrum_orchestrator.py` integrates (M_y) together with the background from `y₀ = 10⁻³` and
finds the roots of (BC). For `t = 10⁻³`, `c = c★`:

    exactly one bound state below the continuum,   μ² = m²/H² = −7.71787  (step-halving stable to 5×10⁻⁹).

The mismatch function is almost exactly linear in μ², `B ∝ (μ² + 7.718)`.

Growth rate on the shell: modes behave as `e^{sHτ}` with `s² + 3s + μ² = 0`, so `s = 1.657`:
the static shell is unstable with an e-folding time of 0.60 Hubble times.

See `../../effective_theory_orchestrator/` for the closed-form 4D explanation,
`μ²(t→0) = −4(3c² − 4c + 8)/(c(3c + 4)) = −7.719796`.

## 5. Caveats

* Floating point; RK4 with step halving. The background reproduces the certified M462 root
  (δ to 1×10⁻⁹, η to 4×10⁻⁸) but is not itself certified.
* `μ² = −4` harmonics are excluded from the no-bending argument; for a single Z₂ brane these are the
  usual pure-gauge "radion" configurations. The bound state found here is far from −4.
* Positivity of the mode norm (tachyon rather than ghost) is supported by the 4D effective action
  (`Z_E > 0`) and by the Schrödinger form with standard L² norm; a full second-order-action proof is not included.
