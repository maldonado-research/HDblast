# A closed-form 4D effective theory for the HDBLAST shell from the radial Hamilton–Jacobi equation

16 September 2026. Status: analytic derivation (two-derivative order, t → 0), validated against the full 5D
linear spectrum. Method has prior art (holographic RG / Hamilton–Jacobi, de Boer–Verlinde–Verlinde
hep-th/9912012; low-energy braneworld expansions, Kanno–Soda; BPS moduli actions, Brax–van de Bruck–Davis–Rhodes).
What is new for this project is the explicit result for the registered model and its consequences.

## 1. Radial Hamiltonian constraint

Bulk action `S = ∫√−G [R/2 − (∂φ)²/2 − U] + GHY`, unit-lapse radial slicing `ds² = dy² + g_μν dx^μdx^ν`,
`K_μν = ½∂_y g_μν`. For a spacelike normal, `R⁽⁵⁾ = R̃ + K² − K_μνK^μν + (total derivative)`, so

    L = √−g [ ½(R̃ + K² − K_μνK^μν) − ½φ_y² − ½(∂φ)² − U ],
    π^μν = ½√−g (K g^μν − K^μν),   π_φ = −√−g φ_y ,

and the radial Hamiltonian constraint is

    ½[(4/3)π² − 4π_μνπ^μν]/|g| − ½π_φ²/|g| − ½R̃ + ½(∂φ)² + U = 0 .                 (HJ)

(On the registered background this is `6(ρ'² − 1)/ρ² = φ'²/2 − U`.)

## 2. Derivative expansion of the on-shell action

The on-shell action of one bulk copy, as a functional of the data on the shell, is expanded as

    S₁[g, φ] = ∫d⁴x √−g [ 𝒲(φ) + Φ(φ) R̃ + ½ M(φ)(∂φ)² ] + (non-local remainder),

with `π^μν = δS₁/δg_μν`, `π_φ = δS₁/δφ`. Order by order (HJ) gives:

* **Order 0:** `U = ½𝒲'² − (2/3)𝒲²`. The solution selected by regularity of the flow into the
  φ = −1 AdS₅ throat is the registered superpotential, `𝒲 = W`.
* **Order 2, coefficient of R̃:**  `W_φ Φ' − (2/3) W Φ = −1/2`.
* **Order 2, coefficient of □φ:**  `M = 2WΦ'/W_φ`.
* **Order 2, coefficient of (∂φ)²:**  `−2WΦ'' + (1/3)WM + ½W_φM' + ½ = 0`, which is *identically satisfied*
  once the previous two hold (checked algebraically) – the expansion is consistent.

The IR-regular solution of the transport equation is

    Φ(φ_b) = I(φ_b)/2,   I(φ_b) = ∫_{y_b}^{∞} e^{2[A(y) − A(y_b)]} dy   along the BPS flow φ = −tanh y, A' = −W/3,

because `dI/dy_b = −1 + (2W/3)I`, i.e. `W_φ I' − (2/3)W I = −1`. At `φ_b = 0`, `I = I₊ = 1.0357712571567`.

## 3. The shell effective action

Two mirror copies plus the shell action `−∫√−g σ_t(φ)` give

    S₄ = ∫d⁴x √−g [ f(φ)R/2 − ½Z(φ)(∂φ)² − V(φ) ] + (CFT-like non-local part),

    f = 2I,     Z = −W f'/W_φ = 2W(1 − 2WI/3)/W_φ²,     V = σ_t − 2W .             (EFT)

Remarks.
* `f = 2I` is the 4D Planck mass² (κ₅ = 1); `V` is exactly the detuning of the tension from its BPS value.
* For pure AdS (`I = 3/(2W)`) `Z = 0`: a single Randall–Sundrum brane has no radion, as it must.
* `Z > 0` iff `I < 3/(2W)`; this holds along the whole registered wall (no ghost).
* The non-local remainder encodes the Kaluza–Klein continuum above `9H²/4` (the holographic CFT).

## 4. de Sitter shells and their stability

Einstein frame: `V_E = V/f²`, `Z_E = Z/f + (3/2)(f'/f)²`. A static de Sitter shell is a stationary point,

    V'/V = 2f'/f ,                                                                  (S)

and its modulus mass in Hubble units is

    μ² ≡ m²/H² = 3 (ln V_E)''/Z_E ,    (ln V_E)'' = V''/V − (V'/V)² − 2[f''/f − (f'/f)²] .   (μ²)

**Registered tension** `V = t(1 + cφ)` at `φ_b = 0` (`W = 1, W_φ = −1, W_φφ = 0`):

* (S) gives `c = 2/I₊ − 4/3` – *exactly the project's registered `c★`*, whose origin is thereby explained:
  it is the condition that `φ_b = 0` is an equilibrium of the modulus.
* `I' = 1 − 2I₊/3 = cI₊/2`, `I''/I = (2 − c)/3`, `Z = f' = cI₊`, `Z_E = c/2 + 3c²/8`, and

      μ²(t → 0) = −4(3c² − 4c + 8) / (c(3c + 4)) = −7.719795918 .

  The discriminant of `3c² − 4c + 8` is negative: **μ² < 0 for every c > 0**.
* Full 5D linear spectrum at `t = 10⁻³` (independent code): `μ² = −7.71787`. Difference `1.9×10⁻³` = O(t).

**General linear detuning.** With `V'' = 0`, (S) turns the stability condition `(ln V_E)'' > 0` into
`(f²)'' < 0`. Along the registered wall `(I²)'' > 0` everywhere (table in `HJ_EFT_TABLE.json`), and
`μ²` ranges from −6.6 to −16: no equilibrium position of a linearly detuned shell is stable or slow-roll.

**Quadratic detuning** `V = t(1 + cφ + dφ²/2)` at `φ_b = 0`: `μ²(t→0) = 3(d − d₀)/Z_E` with
`d₀ = c² + 2(ln f)'' = 1.1134966`.

* `d > d₀` gives a **stable** de Sitter shell. Checked in the full 5D linear problem at `t = 10⁻³`
  (referee track, re-run independently by the orchestrator, `../stability/orchestrator_derivation/stable_shell_recheck.py`):
  `d = 1.3 → μ² = +1.2949033`, `d = 1.4 → +1.987133`, and for `d = 1.6, 2.0` there is no scalar bound state below
  `9/4` at all. The instability is a property of the registered *linear* detuning, not of de Sitter shells in this model.
* **Correction (post-review).** An earlier version of this note said that `d = 1.105924` gives `n_s ≈ 0.965`.
  That was wrong twice. (i) The canonical potential has a large cubic term at the hilltop (`V'''/V = +7.7`), so
  slow roll gives `n_s − 1 = −2|η| coth(|η|N/2) ≤ −4/N`, i.e. `n_s ≤ 0.93` for every `d` in this family
  (tensor track: 0.910/0.918/0.924 at N = 50/55/60) – excluded by Planck. (ii) At finite `t` the hilltop mass obeys
  `μ² ≤ −2.31√t` (slope track); at `t = 10⁻³`, `d = 1.105924` gives `μ² = −0.0899`, not −0.0525. An inflection
  variant (`c = c★ − 9.44×10⁻⁵`, `d = d₀`) reaches `n_s = 0.9649` only with a slope tuned to ~2×10⁻⁵, and
  inflation still never ends (max ε_H ≈ 0.5). No viable inflationary model is obtained from this modulus.

## 5. Post-review additions (16 Sept 2026, after the agent tracks and referee)

* **O(t) coefficient derived analytically** (`../slope_analytic/`): `μ² = μ₀ + s₁t + O(t²)` with
  `s₁ = 1.9243896400`, a closed form in `c`, `I₊`, `ℓ₋³/8` and two explicit wall quadratures `K₁`, `K₂`;
  5D shooting gives 1.9243896(5). The exact shell mass formula is `μ² = −4 − E_b/(3X_b s_b)`.
* **Kinetic normalisation confirmed in 5D** (`../verification/math/`): the pole residue of the brane response to a
  source is −6932.17 against the effective-theory value `−ρ_b²/(fZ_E) = −6932.10`; positive norm (tachyon, not ghost).
* **Landscape corrections.** The modulus space does *not* end at `φ_b = −1`: `Z(−1) = 9/23` is finite and the effective
  theory continues onto the `φ = −coth y` branch, where `V_E` turns negative beyond `φ_b = −1/c = −1.673`; the 5D
  fate there is unknown. Toward `φ_b → +1` the space ends at finite distance `Θ = −3.4653` in a quadratic minimum with
  `m² = 2H²` exactly and `V_E → (1+c)t/81 = 0.019723 t` (the earlier "0.01976" was a finite-distance table entry).
* **The "−4".** Its identification with the Frolov–Kofman / Gen–Sasaki `−4H²` is by value and structure only;
  FK's quantitative analysis assumes a rigid boundary condition that forbids `μ² < −4`, whereas the registered shell
  is in the opposite, soft regime (`σ_t'' ≈ 10⁻⁴`).

## 6. What this does and does not establish

Established (at the stated level): a consistent two-derivative effective action for the shell modulus, the
explanation of `c★`, the closed-form tachyon mass, and agreement with the 5D spectrum through O(t) (~10⁻⁶).

Not established: control of higher-derivative and non-local (CFT) corrections away from `t → 0`; the
non-linear fate of the instability in the full 5D theory; any reheating, perturbation spectrum or
observational statement; external novelty of the method.
