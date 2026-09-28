# Exact analytic structure of the static "+1" branch

HDBLAST research checkpoint, 27–28 September 2026. This is a local calculation. It has not been published or externally reviewed. It deals only with the registered static branch. It does not test whether a five-dimensional event produced our Big Bang.

**Summary.**
- The small-detuning expansions of φ_b, H² and ρ_b are now known exactly in c through O(δ⁸), four orders beyond the previous checkpoint. Solutions of the full nonlinear boundary-value problem at about 31–41 significant digits (solver-setting agreement) confirm the new coefficients. Neither fit nor tuning was used.
- The degree-14 Gegenbauer polynomial comes from the superpotential: Δ = 3W''/W = 18 at φ = 1 is an integer. It has a closed elementary form. The same construction is elementary for every bulk mass.
- At leading order, the tensor spectrum around the φ = 1 shell has an exact closed-form quantization condition: one massless graviton, then a continuum starting at m = (3/2)H. This agrees with Garriga–Sasaki.
- The decoupled bulk scalar has no discrete mode when H/k < 16.6.
- Scalar–metric mixing is not included in this workstream. The stability of the branch therefore remains **open**.

Status labels used below: **exact** means symbolic or rational-arithmetic verification. **numerical** means a floating-point or multiprecision computation with stated convergence evidence. **conditional** means the result depends on a stated approximation.

## Model and conventions

These are the registered model and the conventions of the 22 September package (`frozen/registered_solver.py`, `static_branch/solve_plus_branch.py`). They were not changed.

- Superpotential and potential: W = 1 − φ + φ³/3 and U = ½W_φ² − ⅔W².
- Shell tension: σ = 2W + δ(1 + cφ), with c = 2/1.0357712571566784 − 4/3. This is evaluated at 60 digits from that decimal string. It differs from the float64 value by about 2×10⁻¹⁷.
- Static metric: ds² = dy² + ρ(y)² ds²(dS₄, unit), with a regular cone at y = 0 and H = 1/ρ_b.
- Bulk equations:
  - ρ_y² = 1 + ρ²(φ_y²/12 − U/6)
  - φ_yy + 4(ρ_y/ρ)φ_y = U_φ
- Junctions: ρ_y/ρ = σ/6 and φ_y = −σ'/2.
- Working variables: k = 1/9 (the AdS scale at φ = 1), u = ky, R = kρ, η = φ − 1.
  - V ≡ −U/(6k²) = 1 − 21η² − 25η³ + (9/4)η⁴ + 6η⁵ + η⁶
  - U_η/k² = 252η + 450η² − 54η³ − 180η⁴ − 36η⁵

## 1. Small-δ expansions, exact in c: a new result for this project

### Method (`series_expansion.py`)

1. **Series variables.** Cone regularity selects the growing bulk mode. Write η and R as a formal double series in X = αe^{14u} (the regular growing mode) and T = e^{−2u} (the dS curvature):
   - η = N(X,T), with N(X,T) = X(1 + …) + O(X²)
   - R = (e^u/2)P(X,T), with P = 1 − T + O(X²)
2. **Operator.** d/du acts as D(X^mT^j) = (14m − 2j)X^mT^j. The bulk equations become algebraic recursions with the operator (p−14)(p+18), where p = 14m − 2j.
3. **Truncation.** At the shell, X_b and T_b are both O(δ). The two junction conditions are solved as power series in δ.
4. **Why cone regularity drops out through δ⁸.** It enters only through homogeneous pieces proportional to α = X_bT_b⁷ ∝ δ⁸. The first resonance, (m,j) = (2,7), sits at total degree 9 (Section 1c). So through O(δ⁸) the expansion is fixed entirely by the local analysis near the shell.

The computation uses exact rational arithmetic in c. It checks that every lower-order residual vanishes. It also checks that the linear part of N equals C₁₄⁽²⁾(cosh u)/15 term by term.

### Results (exact; all coefficients through δ⁸ are in `SERIES_COEFFICIENTS.json`)

- **Scalar at the shell:**

  φ_b = 1 − (9c/64)δ − [9c(2503c+1288)/942080]δ² − [27c(60469c²+49568c+4048)/331612160]δ³ − [81c(2674271029069c³+2637538943248c²+105518716992c−201399593472)/146654708419788800]δ⁴ + O(δ⁵)

- **Expansion rate:**

  H² = (1+c)δ/27 + [(29c²+64c+32)/1152]δ² − [c²(11767c+11362)/2826240]δ³ − [3c²(161863c²+330944c+113344)/1326448640]δ⁴ + O(δ⁵)

- **Shell radius:**

  ρ_b = √(27/((1+c)δ)) × {1 − [3(29c²+64c+32)/(256(1+c))]δ + [9(384281c⁴+1465672c³+2144336c²+1413120c+353280)/(15073280(1+c)²)]δ² − [27(236579153c⁶+…)/(42446356480(1+c)³)]δ³ + O(δ⁴)}

- **Cone displacement:** η_h = (136/3)X_bT_b⁷, up to relative O(δ⁸ log δ). Therefore

  η_h = −(111537/131072)c(1+c)⁷δ⁸ [1 + O(δ)]

  The previously "tiny" cone value −1.3367×10⁻²³ therefore has an exact leading closed form.

- **Shell position:** y_b = (9/2)ln[4/(3(1+c)δ)] + [27(35c²+64c+32)/(256(1+c))]δ + O(δ²).

**Known terms reproduced exactly.** The O(δ) and O(δ²) coefficients of the 22 September checkpoint are reproduced exactly: −9c/64, (1+c)/27, and (1+c)²/36 − c²/384 = (29c²+64c+32)/1152.

**Known-limit check.** At c = 0 the H² series stops at δ², as it must. In that case φ ≡ 1 solves the problem exactly and H² = σ²/36 − k².

Coefficients at the registered c:

| Order | φ_b − 1 | H² | ρ_b factor |
|---|---:|---:|---:|
| δ | −0.08403678774 | 0.05917018278 | −0.5912394192 |
| δ² | −0.01589265892 | 0.06996748900 | 0.5439862661 |
| δ³ (new) | −0.002688961364 | −0.002324227250 | −0.5490083550 |
| δ⁴ (new) | −0.0004536042830 | −0.0002979720758 | 0.5824790232 |

The previous package measured (η_b + 9δc/64)/δ² ≈ −0.015895 and (H² − 2nd order)/δ³ ≈ −0.0023245. These agree with the exact a₂ and h₃ above to the expected O(δ) accuracy.

### Numerical confirmation (`mp_bvp.py`, `run_bvp_scan.py`, `analyze_series_vs_bvp.py`)

**The solver is independent of the series code.** It uses mpmath at 50–65 digits and a Taylor-series integrator of order 50–64. It starts from a Frobenius cone series of order 90 at u = 0.25. The shell is located by an event on J1, and η_h is found by secant iteration on J2. The first integral is not enforced; its residual is recorded as a check.

**Coverage.** Four values of c: 0.5976 (registered), −0.4, 1.3 and 0. Detunings from 10⁻⁴ to 0.1: 31 solves in total.

| Check | Result | Status |
|---|---|---|
| Junction residual J2; first-integral residual | ≤ 3.6×10⁻⁵⁰; ≤ 2.8×10⁻³² absolute (≤ 4.9×10⁻³³ relative to its natural scale; the largest values occur at δ = 0.1) | numerical |
| Solver-setting convergence (dps 50/order 50/h 0.125 vs dps 65/order 64/h 0.1) | Maximum relative differences over φ_b, H², ρ_b, η_h, y_b: 3.6×10⁻⁴² (δ = 10⁻⁴), 2.4×10⁻⁴¹ (8×10⁻⁴), 1.0×10⁻³⁷ (0.0064), 1.8×10⁻³¹ (0.1) | numerical |
| Agreement with the 22 Sept float64 values at δ = 0.001 | φ_b = 0.99991594731691335388 (package 0.9999159473169134); ρ_b = 129.92476284968243 (package 129.9247628496781, relative 3×10⁻¹⁴) | numerical |
| **Blind extraction.** Polynomial extrapolation of (Q − S₂)/δ³ over δ = 10⁻⁴·2ᵏ, k = 0…6. Uses only the previously known two orders. | Recovers the exact new coefficients. Relative errors for q₃ / q₄ / q₅: registered c: 1.1×10⁻²⁶ / 1.7×10⁻²¹ / 1.0×10⁻¹⁶ (H²) and 1.4×10⁻²⁵ / 1.7×10⁻²⁰ / 6.8×10⁻¹⁶ (φ_b); similar for c = −0.4 and 1.3 | numerical confirmation of exact coefficients |
| Scaled remainders (Q − S_n)/δ^{n+1} → q_{n+1} | Deviation ∝ δ for every order n = 0…7, for all five quantities and all three nonzero c | numerical |
| Full series through δ⁸ vs BVP at δ = 0.001 | abs. diff 2.8×10⁻³⁴ (H²), 3.3×10⁻³³ (φ_b) | numerical |
| Full series through δ⁸ vs BVP at δ = 0.1 | 2.9×10⁻¹⁶ (H²), 3.4×10⁻¹⁵ (φ_b) | numerical |
| c = 0 (known limit) | η ≡ 0; H² − (δ/27 + δ²/36) ≤ 7×10⁻⁵³ | numerical, exact limit |
| **Control:** wrong coefficient, h₃ → h₃(1+10⁻⁴) | Scaled δ⁴ remainder is off by factors of 7.8, 152 and 5.1 (three values of c), vs 10⁻⁵ when correct | control fails as it should |
| **Control:** perturbed parameter, c → c(1+10⁻⁶) in the series | Scaled δ⁴ remainder off by about 10⁷–10⁹ | control fails as it should |

### 1c. Where the local expansion ends (`resonance_probe.py`)

**The resonance (exact).** At (m,j) = (2,7) the recursion operator vanishes. The source there is S₂,₇ = −24576/17 in exact rational arithmetic, which is nonzero. The regular solution therefore contains a secular term κ·u·α²e^{14u} with κ = 768/17.

**Consequence 1: no logarithm in the shell observables (numerical).** For φ_b and H², the log can be absorbed into the growing-mode amplitude that the shell equations determine. We fitted the numerical δ⁹ remainders with an added ln δ term. For the registered c and c = 1.3, the log coefficient is below 3×10⁻⁷ of the constant term. For c = −0.4, where the δ⁹ coefficient is about 10⁻¹¹ and near the precision floor, the ratio is ≤ 4×10⁻³. The δ⁹ coefficient itself depends on cone regularity; it is not given by the local recursion alone.
- Registered c: 3.249147×10⁻⁶ for φ_b and −2.759695×10⁻⁷ for H² (numerical).
- Computing it exactly would need global second-order matching. That was not attempted.

**Consequence 2: a logarithm in the cone value (prediction confirmed numerically).** η_h must acquire a term

η_h ⊃ δ⁸ · (81/4)c²[3(1+c)/4]¹⁴ · δ⁸ ln δ

Checked against the numerics:
- The fitted log coefficient agrees with this prediction to 2.2×10⁻⁵, 1.9×10⁻⁶ and 4.8×10⁻⁵ relative for c = 0.5976, −0.4 and 1.3.
- A pure-polynomial model fails an out-of-sample test by a factor of more than 10⁶.

## 2. Why a degree-14 Gegenbauer polynomial (`gegenbauer_structure.py`): exact, 28 checks

1. **Linear equation.** Around φ = 1 in pure AdS (R = sinh u), the linear radial equation is η'' + 4coth(u)η' = (m²/k²)η with m²/k² = U''(1)/k² = (28/9)·81 = 252. With x = cosh u it becomes

   (x² − 1)F'' + 5xF' − n(n+4)F = 0

   This is the Gegenbauer equation with λ = 2 and n(n+4) = 252, so n = 14.

2. **Why the number is an integer.** For a potential built from a superpotential, at a critical point of W: U'' = W''² − (4/3)WW'' and k = W/3. Hence m²/k² = s(s−4) with s = 3W''/W. This is the conformal dimension Δ of the AdS/CFT dictionary, and n = Δ − 4.
   - At φ = +1: s = 3·2/(1/3) = 18, so Δ = 18 and **n = 14 exactly**.
   - At φ = −1: W = 5/3 and W'' = −2, so s = −18/5, Δ = 38/5 and n = 18/5, which is not an integer. There is no polynomial there.
   - The superpotential-to-dimension relation is standard in this "fake supergravity" setting (DeWolfe–Freedman–Gubser–Karch, https://arxiv.org/pdf/hep-th/9909134). No novelty is claimed for it.

3. **Quantization condition.** The solution regular at the cone is 2F1(−n, n+4; 5/2; (1−x)/2). The coefficient ratio contains (j − n), so the series terminates, giving a polynomial, iff n is a non-negative integer.

4. **Closed forms (verified exactly):**
   - C₁₄⁽²⁾(cosh u) = Σ_{j=0}^{14} (j+1)(15−j) e^{(14−2j)u}, with every coefficient positive.
   - C₁₄⁽²⁾(cosh u) = ½ U′₁₅(cosh u) = ½ (1/sinh u) d/du[sinh 16u / sinh u].
   - Explicit polynomial: C₁₄⁽²⁾(x) = 245760x¹⁴ − 745472x¹² + 878592x¹⁰ − 506880x⁸ + 147840x⁶ − 20160x⁴ + 1008x² − 8.
   - C₁₄⁽²⁾(1) = 680 = C(17,3).

5. **Every bulk mass is elementary, not just this one (5D is odd).** f_ν = (1/sinh u) d/du[sinh(νu)/sinh u] with ν = √(4 + m²/k²) solves the radial equation for any ν. The cone-singular partner is (1/sinh u) d/du[cosh(νu)/sinh u] ≈ −u⁻³. The polynomial property is the special case of integer ν − 2.

6. **Consequences:**
   - η_h/α = 680/15 = 136/3 (used in Section 1).
   - The positive exponential form shows d ln C/du > 0 for u > 0. So the finite-curvature linear static response η_b = −(9/2)δc/(d ln C/du + 18) never has a pole. In particular, **no scalar zero mode exists at any shell radius** at linear order.
   - The producer's log-derivative formula, which uses dC₁₄⁽²⁾/dx = 4C₁₃⁽³⁾, is confirmed.
   - Wrong-index (C₁₄⁽¹⁾) and wrong-degree (C₁₃⁽²⁾) controls fail as they should.

## 3. Leading-order fluctuation spectrum around φ = 1 (`spectrum_leading_order.py`)

**Setting (conditional).** The background is the leading order in δ: pure AdS₅ in dS slicing, R = sinh u, with the shell at u_b and H = k/sinh u_b. Modes are f(u)Y(x) with M² = m²/H², governed by

f'' + 4coth u f' + (M²/sinh²u − m₅²/k²)f = 0

- **Tensor modes:** m₅ = 0, with the Neumann condition f'(u_b) = 0.
- **Bulk scalar:** m₅²/k² = 252, with the Robin condition f' = −(σ''/2k)f = −18f. The coefficient 18 = 3W''/W = Δ exactly, because the tension is σ = 2W at leading order.
- **Mixing omitted:** scalar–metric mixing through the O(δ) background gradient (brane bending/radion) is **not included**.

### Exact reduction (checked symbolically)

Substituting f = sinh^{−3/2}(u) w(cosh u) gives the associated Legendre equation with degree ν(ν+1) = m₅²/k² + 15/4 and order μ² = 9/4 − M². This gives ν = 3/2 for the tensor and ν = 31/2 for the scalar. Solutions normalizable at the horizon are w = P^{−μ}_ν. The exact quantization condition, for any bulk mass and Robin coefficient b, is

(ν − 3/2)cosh u_b P^{−μ}_ν(cosh u_b) + b sinh u_b P^{−μ}_ν(cosh u_b) = (ν − μ)P^{−μ}_{ν−1}(cosh u_b)

Here the physical mass is m² = (9/4 − μ²)H². The continuum threshold m² = 9H²/4 does not depend on the bulk mass. The Legendre derivative identity was checked with mpmath to 4×10⁻²⁹.

| Sector | Exact / numerical result | Status |
|---|---|---|
| Tensor | The condition reduces to **(μ − 3/2) P^{−μ}_{1/2}(cosh u_b) = 0**. Apply Pfaff's transformation, P^{−μ}_{1/2} ∝ (1−z)^{1/2} 2F1(μ−½, −½; 1+μ; z/(z−1)). For μ > ½ the series is strictly decreasing toward a positive Gauss value, Γ(1+μ)/(Γ(3/2)Γ(μ+3/2)); for μ ≤ ½ every term is ≥ 0. So P^{−μ}_{1/2} > 0 for x > 1 and μ > 0. **The only discrete mode is μ = 3/2 (m = 0, f = const), for every u_b. The continuum is m ≥ (3/2)H.** There is no tensor tachyon. | exact (numerical sweep: min normalized P = 1.0000027 > 0; no sign change of f' for 0 < M² < 9/4 at 8 shell radii) |
| Tensor, full nonlinear branch | f = const is an exact zero mode for any warp factor. It is normalizable because the doubled bulk has finite warped volume. Float shooting on the full branch at δ = 0.001, 0.01, 0.1 finds no sign change of f'(u_b) for 0.02 ≤ M² ≤ 2.24. | exact (zero mode) / numerical (no gap mode) |
| Bulk scalar (decoupled) | **No zero mode and no tachyon for any u_b.** Proof: for M² ≤ 0 the regular solution stays positive and increasing (maximum principle), so f'/f > 0 > −18. For large u_b, f'/f → Δ − 4 = 14 for every M², while the shell requires −Δ = −18. Mismatch table (F + b > 0 throughout 0 < M² < 9/4): 31.3 at u_b = 1; 31.99 at u_b = 3.36 (the registered shell). **Result: no discrete scalar mode below 9H²/4 for u_b > 0.0600, i.e. H/k < 16.65.** One bound state appears only for H/k ≳ 16.6, far outside the small-δ regime; for example, M² = 1.004 at H/k = 50. | exact (M² ≤ 0) / numerical (scan, 90 points × 16 radii) |
| Calibration | Direct ODE shooting in pure AdS vs the exact Legendre formula: tensor abs. error ≤ 2.1×10⁻¹² and scalar relative error ≤ 8.3×10⁻¹⁴. A wrong-degree Legendre control (ν = 1/2) misses by ≥ 0.30. | numerical control |

### Literature (established; not new here)

- **Garriga & Sasaki,** https://arxiv.org/abs/hep-th/9912118. An inflating brane bounding an AdS₅ ball, with a massless graviton and a Kaluza–Klein gap m = (3/2)H. In Euclidean signature, our leading-order configuration (S⁴ shell bounding an H⁵ ball, doubled) is their "brane-world creation from nothing" instanton. That is an established idea related to the hypothesis. It is **not evidence for it**.
- **Karch & Randall,** https://arxiv.org/pdf/hep-th/0011156 (the AdS₄-brane counterpart).
- **Hawking, Hertog & Reall,** https://arxiv.org/abs/hep-th/0003052.
- **Himemoto & Sasaki,** https://arxiv.org/abs/gr-qc/0010035 (bulk scalar with a dS brane).
- **Frolov & Kofman,** https://arxiv.org/abs/hep-th/0309002. They find that the radion of inflating, stabilized braneworlds is typically tachyonic. This is exactly the scalar–metric sector omitted here.
- **Garriga & Vilenkin,** https://journals.aps.org/prd/abstract/10.1103/PhysRevD.44.1007 (wall fluctuations in dS; tachyonic wall mode not always an instability).
- **Langlois, Maartens, Sasaki & Wands,** https://arxiv.org/pdf/hep-th/0012044.

**New to this project, not claimed as new mathematics externally:**
- the δ-expansion through δ⁸ and the δ⁸ ln δ structure of η_h;
- the compact tensor condition and its positivity proof, written out here. We did not check whether this exact form appears in the literature;
- the observation b = Δ for the superpotential shell;
- the decoupled-scalar statements for this registered model.

## Limitations

- **Scope of the expansion.** It is a formal asymptotic expansion, verified numerically to δ = 0.1. No rigorous remainder bound or radius of convergence was proved. The δ⁹ coefficient depends on cone regularity and is only measured, not derived.
- **Numbers are not interval-certified.** Multiprecision agreement between two solver settings is convergence evidence, not a certificate.
- **Stability is not decided.** The spectrum analysis is leading order and omits the brane-bending/radion mixing. That mixing is O(δ) and is where Frolov–Kofman-type tachyons appear. Stability of the +1 branch is therefore **not** established. The absence of light decoupled scalar and tensor modes is necessary, not sufficient.
- **Radion work elsewhere.** Gauge-invariant radion analysis belongs to the separate stability workstreams (`../stability_gauge_invariant`, `../stability_time_domain`). It was not duplicated here.
- **Nothing here bears on the Big Bang claim.** Nothing addresses the connection from the φ ≈ −1 shell to this branch, a radiation era, or observations.

## Reproduction

Requirements: python3 with mpmath 1.3, sympy 1.14, numpy, scipy.

```sh
cd research/HDBLAST_CHECKPOINT_20260927/analytic_structure
./run_all.sh          # everything, about 20-25 min single core
# or individually:
python3 series_expansion.py 8          # -> SERIES_COEFFICIENTS.json (~50 s)
python3 gegenbauer_structure.py        # -> GEGENBAUER_STRUCTURE.json (28 exact checks)
python3 run_bvp_scan.py reg reg 50 50 0.125 0.25 0.0001,0.0002,0.0004,0.0008,0.0016,0.0032,0.0064,0.0003,0.001,0.003,0.01,0.03,0.1
python3 run_bvp_scan.py cm04 -0.4 50 50 0.125 0.25 0.0001,0.0002,0.0004,0.0008,0.0016,0.0032,0.0064
python3 run_bvp_scan.py cp13 1.3 50 50 0.125 0.25 0.0001,0.0002,0.0004,0.0008,0.0016,0.0032,0.0064
python3 run_bvp_scan.py c0 0 50 50 0.125 0.25 0.0001,0.001,0.01,0.1
python3 run_bvp_scan.py reg_conv reg 65 64 0.1 0.2 0.0001,0.0008,0.0064,0.1
python3 analyze_series_vs_bvp.py       # -> SERIES_VS_BVP.json
python3 resonance_probe.py             # -> RESONANCE_PROBE.json
python3 spectrum_leading_order.py      # -> SPECTRUM_LEADING_ORDER.json (~30 s; needs runs/BVP_SCAN_reg.json)
python3 make_manifest.py               # -> MANIFEST.sha256.json
```

Machine-readable outputs are `SERIES_COEFFICIENTS.json`, `SERIES_VS_BVP.json`, `RESONANCE_PROBE.json`, `GEGENBAUER_STRUCTURE.json`, `SPECTRUM_LEADING_ORDER.json` and `runs/BVP_SCAN_*.json`. The last store all BVP values as 50–65 digit strings. Logs are in `logs/`. No file outside this folder was modified.
