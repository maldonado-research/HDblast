# An exact entropy identity for the registered shell instanton

16 September 2026. Status: exact identity (short proof below) + floating-point confirmation.
The physical interpretation (creation probability) is conditional and is not a claim about our Universe.

## Statement

Continue the registered SO(4,1) solution to Euclidean signature: two mirrored 5-balls
`ds² = dy² + ρ(y)² dΩ₄²`, `0 ≤ y ≤ y_b`, glued along the shell `y = y_b` with the registered junction
conditions. This compact O(5)-symmetric geometry is the "creation of the shell universe" instanton
(the bulk-scalar analogue of the Garriga–Sasaki brane-world instanton, hep-th/9912118).

Its on-shell Euclidean action is **exactly**

    S_E = −16π² ∫₀^{y_b} ρ(y)² dy = −8π² M₄²/H² ,        M₄² ≡ 2∫₀^{y_b} (ρ/ρ_b)² dy ,   H = 1/ρ_b ,

i.e. minus the four-dimensional de Sitter entropy evaluated with the graviton zero-mode Planck mass, and
also minus the Bekenstein–Hawking entropy `A₅/(4G₅)` of the five-dimensional horizon
(`A₅ = 2·4π∫ρ² dy`, `1/(4G₅) = 2π` for κ₅² = 8πG₅ = 1).

## Proof

On shell `R = (∂φ)² + (10/3)U`, so the bulk Lagrangian `R/2 − (∂φ)²/2 − U` equals `(2/3)U`. The curvature
delta-function at the shell plus the Gibbons–Hawking terms contribute `[K] = −4σ_t/3`, so the shell term is
`σ_t − 4σ_t/3 = −σ_t/3`. With `Ω₄ = 8π²/3` and the Israel condition `σ_t = 6ρ'_b/ρ_b`:

    S_E = −(4/3)Ω₄ ∫₀^{y_b} ρ⁴U dy − 2Ω₄ ρ_b³ ρ'_b .

The background equations `ρ'² = 1 + ρ²(φ'²/12 − U/6)` and `(ρ'/ρ)' = −1/ρ² − φ'²/3` give
`ρ³ρ'' = ρ⁴(−φ'²/4 − U/6)` and therefore the exact total derivative

    d/dy [2ρ³ρ'] = 6ρ²ρ'² + 2ρ³ρ'' = 6ρ² − (4/3)ρ⁴U .

Since `ρ³ρ' → 0` at the regular centre, `2ρ_b³ρ'_b = ∫₀^{y_b}[6ρ² − (4/3)ρ⁴U]dy`, and the `U` terms cancel:
`S_E = −6Ω₄∫ρ² dy = −16π²∫ρ² dy`. ∎

The identity holds for every regular O(5) solution of this Einstein–scalar system with a Z₂ shell obeying
the Israel condition; it does not use the scalar junction condition or the specific tension function.

## Numbers (`shell_instanton_action.py`)

| t | S_E from the 5D action | −8π²M₄²/H² | ratio |
|---|---:|---:|---:|
| 1e-2 | −1.00186621e5 | −1.00186571e5 | 1.0000005 |
| 3e-3 | −3.37088814e5 | −3.37088224e5 | 1.0000018 |
| 1e-3 (registered) | −1.01448911e6 | −1.01448846e6 | 1.0000006 |

(For smaller t the direct 5D evaluation loses digits to a large bulk/shell cancellation; the identity is exact.)
The leading effective-theory value `−24π² f²/V` with `f = 2I₊`, `V = t` agrees to 0.2% at the registered t,
the difference being the O(t) corrections to `h` and `M₄²`.

## Reading

* The registered shell universe has de Sitter entropy `S_dS = 1.0145×10⁶` in registered units.
* In a no-boundary/tunnelling-wavefunction reading, `exp(±S_dS)` weights its creation; which sign applies is
  the long-standing Hartle–Hawking versus Linde–Vilenkin question and is not resolved here.
* Because the shell's modulus mass is `μ² = −7.72 < −4`, the 4D fluctuation operator on S⁴,
  `ℓ(ℓ+3) + μ²`, is negative for ℓ = 0 and ℓ = 1 (1 + 5 modes) in the effective-theory counting. The instanton is
  therefore not a standard single-negative-mode bounce. The precise 5D negative-mode count (where the ℓ = 1
  harmonics are the special conformal-Killing ones) has not been computed.

Prior art: brane-world creation instantons (Garriga–Sasaki), de Sitter brane entropy as bulk horizon area
(Hawking–Maldacena–Strominger hep-th/0002145). The identity with a bulk scalar is elementary; no novelty is claimed.
