# Linear stability of the static "+1 branch": gauge-invariant perturbation theory (Method A)

HDBLAST research checkpoint, 27 September 2026. This is a new local calculation for this project. It has not been published or externally reviewed, and nothing in this folder claims observational support for the hypothesis.

## Question

The 22 September 2026 checkpoint (Zenodo record 22922928) found a static de Sitter configuration, the "+1 branch". It satisfies both shell junction conditions of the registered model, and its regular cone sits near φ=+1. That checkpoint left open whether the configuration is **linearly stable**. This folder answers that question for the scalar sector (including brane bending and the shell scalar condition) and the tensor sector. It uses a gauge-invariant formulation and first calibrates it on the original shell, which is known to be unstable.

## Short answer

| Sector | Result on the +1 branch, δ ∈ {0.0003, 0.001, 0.003, 0.01, 0.03, 0.1} | Status |
|---|---|---|
| Scalar (bulk scalar, metric and brane bending, both junctions) | **No normalizable mode with μ² < 9/4 at any δ**, so there are no unstable and no marginal scalar modes. The mismatch stays ≥ 0.9927 on its normalized [−1, 1] scale, far from its zero. | numerical, with analytic support |
| Scalar, complex μ² | No zeros inside two contours, [−60,2]×[−30,30] and [−400,2]×[−200,200]. Complex eigenvalues are excluded analytically whenever B>0, and B≈3.5 here. | numerical and conditional-analytic |
| Scalar, μ² < −4 | Excluded analytically, since B>0 at every δ | conditional (see §4) |
| Scalar ℓ=1 harmonic (μ²=−4 translations) | Pure gauge in the bulk. A relative brane displacement would need Bφ′_b=0, but Bφ′_b ≠ 0, so there is no physical mode. | exact-verified (symbolic) + numerical value |
| Scalar ℓ=0 static zero mode | None: the static junction Jacobian is non-singular (condition number 15–17) | numerical |
| Tensor | Only the massless graviton, μ²=0, which is normalizable because ∫ρ²dy is finite. No Kaluza–Klein bound states in (0, 9/4), so the gap is m ≥ 3H/2. No tachyons: μ² ≥ 0 is proven by a Sturm–Liouville argument. | exact-verified (argument) + numerical |
| Calibration: original unstable shell, δ=0.001 | μ² = −7.717871625260, growth p = 1.657193631259, exactly one bound state | numerical; passes |

**Conclusion.** Within the sectors computed (scalar, including brane bending and both junction conditions, and tensor), the static +1 branch is **linearly stable** at all six detunings tested. The massless graviton is the only marginal mode, and it is not an instability. Vector perturbations, nonlinear stability, quantum tunnelling and dynamical attraction were not studied (see Limitations).

## 1. Model and conventions

We use the registered model without changes, with κ₅²=1:

S = ∫√−g [R/2 − ½(∂φ)² − U] over two Z₂ copies, minus ∫_shell √−h σ(φ).

- W = 1−φ+φ³/3 and U = ½W_φ² − (2/3)W².
- σ = 2W + δ(1+cφ), with c = 0.5975949350280132.

The background metric is ds² = dy² + ρ(y)² γ, where γ is a unit dS₄. The regular cone is at y=0 (ρ≈y), the shell is at y_b, and the bulk is y<y_b. The junction conditions are ρ′/ρ = σ/6 and φ′ = −σ′/2. The shell Hubble rate is H = 1/ρ_b.

Perturbations are decomposed into dS₄ harmonics with □_γY = μ²Y, so the physical 4D mass is m² = μ²H².

## 2. Derivation (`derive_linearized.py` → `DERIVATION_RESULTS.json`; 55/55 symbolic checks pass)

The script computes the full 5D linearized Einstein–scalar equations with sympy, starting from the Christoffel symbols. It uses the most general scalar perturbation (N, B, ψ, E, χ) and a concrete harmonic Y=e^{pt}, with μ² = −(p²+3p). The key steps are:

- The background equations and both background junctions are reproduced.
- Every perturbation equation depends on p only through μ², which was checked by polynomial reduction.
- The gauge shifts δN=T′, δB=L′+T/ρ², δψ=HT, δE=L and δχ=φ′T equal the Lie derivative. All five linear equations are gauge invariant on-shell.
- The gauge-invariant variables are X = χ − φ′ρ²(B−E′) and Ψ = ψ − Hρ²(B−E′). The trace-free equation gives Ñ = −2Ψ. The momentum and Gauss constraints give a closed first-order system. The two second-order equations (trace, scalar) and the yy equation then follow from it identically, a Bianchi-type consistency check.
- The system is written in the variable **Z = Ψ/φ′**, which stays finite even though φ′ ~ 10⁻²³ near the +1 cone:

  Z′ = −(2H+g) Z − X/3,  X′ = g X + (3λ/ρ² − 2φ′²) Z,  g = φ″/φ′,  λ = μ²+4.

- **Junctions on the displaced shell** y = y_b + ζY. The script computes the unit normal, extrinsic curvature and induced metric to first order in a general gauge, and imposes K_mn = (σ/6)h_mn and n·∂φ = −σ′/2. All three conditions are invariant under the shift ζ → ζ − T.
  - In longitudinal gauge the trace-free condition forces ζ=0 (no brane bending for μ² ≠ 0, −4).
  - The trace condition is identical to the bulk momentum constraint.
  - The scalar condition becomes
    **M(μ²) = X′ + (σ″/2)X + 2φ′Ψ = B X + 3λZ/ρ² = 0 at y_b, with B = φ″/φ′ + σ″/2.**
- **ℓ=1 harmonic (p=1, μ²=−4).** The bulk solution space is exactly the two-parameter gauge family. A relative brane displacement Δ leaves the metric junction satisfied identically and changes the scalar junction by BΔφ′. So a physical ℓ=1 mode exists only if Bφ′_b=0.
- **Tensor sector** (h_{x1x2}, TT, homogeneous). The equation is h″ + 4Hh′ + μ²h/ρ² = 0 with the junction h′(y_b)=0 (Neumann).
- **Wrong-formula controls:** four deliberate errors were all detected:
  - wrong gauge shift δψ=−HT;
  - wrong sign of the ζ shift;
  - scalar junction without the 2φ′Ψ term;
  - wrong sign in the Z equation.

These equations match the independent Gaussian-normal and longitudinal derivations of Chat 9 (`D-Blast 3/untitled folder 146/.../stability/`). That agreement is a consistency check with this project's earlier work, not an external validation.

## 3. Stability criterion

A mode with shell dependence Q(x) obeys (□_dS₄ − m²)Q = 0. At late times it behaves as e^{pHτ} with p² + 3p + μ² = 0.

- The mode grows (Re p > 0) exactly when Re μ² < (Im μ²)²/9. For real μ², that means μ² < 0.
- The regular cone exponent of the bulk solution is s₊ = −3/2 + √(9/4−μ²). This is the same p, so normalizability at the cone and the growth rate are the same condition.
- Bound states require μ² < 9/4. Above that is the continuum, which is a stable (decaying) band.

The linear criterion used here is therefore: **the configuration is stable if there is no normalizable mode with Re μ² < (Im μ²)²/9, apart from pure-gauge ℓ=1 translations.**

This is the standard framework for dS-sliced branes and bubble walls:
- continuum gap m = 3H/2 above the massless graviton for a dS brane in AdS ([Garriga & Sasaki](https://arxiv.org/abs/hep-th/9912118));
- tachyonic radion of inflating stabilized branes with a bulk scalar ([Frolov & Kofman](https://arxiv.org/abs/hep-th/0309002));
- m² = −4H² radion / wall-fluctuation modes ([Gen & Sasaki](https://arxiv.org/pdf/gr-qc/0011078); [Garriga & Vilenkin, PRD 44, 1007](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.44.1007));
- the wall-fluctuation mode with p² = −4, interpreted as the Goldstone mode of O(4,1)→O(3,1) that is removed from the physical spectrum once gravity is included ([wall fluctuation modes in open inflation](https://arxiv.org/pdf/astro-ph/9702174));
- flat-wall scalar fluctuation equations ([DeWolfe, Freedman, Gubser & Karch](https://arxiv.org/pdf/hep-th/9909134)).

The growth-rate identification is confirmed on the calibration branch. The recorded rate 1.657193631 follows from the recorded μ² through this same relation.

## 4. Analytic statements (conditional on the derived equations; identities checked symbolically)

The following identities hold for real solutions of the (X,Z) system; the Wronskian one also holds with complex conjugates:

- [ρ²ZX]′ = 3λZ² − ρ²X²/3 − 2ρ²φ′²Z²
- [ρ²(Z₂X₁ − X₂Z₁)]′ = 3(λ₁−λ₂)Z₁Z₂

Cone terms vanish for normalizable modes. Using the shell condition X_b = −3λZ_b/(Bρ_b²):

1. **Energy bound:** λ·[3∫Z² + 3Z_b²/B] = ∫ρ²(X²/3 + 2φ′²Z²) > 0. If **B>0**, then λ>0, i.e. **μ² > −4**.
2. **Reality:** (λ−λ̄)·[3∫|Z|² + 3|Z_b|²/B] = 0. If **B>0**, every eigenvalue is real.
3. **Tensor:** −(ρ⁴h′)′ = μ²ρ²h with a Neumann condition. Integrating gives ∫ρ⁴h′² = μ²∫ρ²h², so μ² ≥ 0. The zero mode h=const is normalizable because the bulk is compact in y.

On the +1 branch, B = 3.478 to 3.555, dominated by σ″/2 = 2φ_b ≈ 2 plus g_b ≈ 14/9. The value 14/9 is the growth rate of the regular Δ=18 AdS scalar, and B → 32/9 as δ → 0 (a leading-order observation, not proved). The theorems therefore rule out μ² < −4 and complex modes. Only the window −4 < μ² < 0 needs numerics for instability; the full range up to 9/4 was scanned. On the original shell, B = −2.68×10⁻⁴ < 0, which is consistent with its tachyon at μ² < −4.

## 5. Calibration first (`spectrum.py` → `SPECTRUM_RESULTS.json`)

These are the same operator and code applied to the original unstable shell, whose background comes from `frozen/registered_solver.py`:

| Quantity | This folder | Recorded | Difference |
|---|---:|---:|---:|
| μ², δ=0.001 (rtol 1e-13, y₀ 1e-5) | −7.717871625260 | −7.717871625176 (Chat 9 GN) | −8.4×10⁻¹¹ |
| | | −7.717871618737 (Chat 9 longitudinal) | −6.5×10⁻⁹ |
| growth p | 1.657193631259 | 1.657193631245 (prior); 1.657193663767 (folder 152 spectral) | 1.4×10⁻¹¹; −3.3×10⁻⁸ |
| μ², δ=0.003 | −7.714023603 | −7.7140236 | < 10⁻⁸ |
| μ², δ=0.01 | −7.700561459 | −7.700561458 | ~6×10⁻¹⁰ |
| bound states below 9/4 | 1 | 1 | — |
| winding number (both contours) | 1.0000000000 | 1 | — |
| unreduced longitudinal cross-check, δ=0.001 | −7.717871624 | — | 1.4×10⁻⁹ |

For the δ=0.001 root, the spread over rtol ∈ {1e-10, 1e-12, 1e-13} × y₀ ∈ {1e-3, 1e-4, 1e-5} is 9.7×10⁻⁹ and dominated by rtol=1e-10. At rtol=1e-13 the spread is 8.6×10⁻¹¹. The recorded value was known in advance, so this is a validation, not a blind prediction.

## 6. +1 branch results (δ from the checkpoint scan)

The backgrounds were reproduced by calling the package's `solve_plus_branch.solve`. At δ=0.001 this gives φ_b = 0.9999159473169134, ρ_b = 129.92476284968, H² = 5.924014794329×10⁻⁵ and H/H₀ = 0.606721732, all matching the stated digits (`BACKGROUNDS.json`). Each background was then re-integrated and root-polished by independent code. Junction residuals are ≤ 4×10⁻¹⁷ after polishing, and ρ_b agrees with the package to ≤ 6×10⁻¹⁵ relative.

| δ | H_brane | B | min M̂ over 360 μ² in [−400, 2.2499] | M̂(0) | scalar roots | longitudinal roots | tensor roots | max spread over 9 precision variants |
|---:|---:|---:|---:|---:|---|---|---|---:|
| 0.0003 | 0.0042139 | 3.55532 | 0.999981 | 0.999988 | none | none | μ²=0 only | 3.3×10⁻¹⁶ |
| 0.001 | 0.0076968 | 3.55478 | 0.999937 | 0.999960 | none | none | μ²=0 only | 2.2×10⁻¹⁶ |
| 0.003 | 0.0133469 | 3.55323 | 0.999812 | 0.999880 | none | none | μ²=0 only | 4.4×10⁻¹⁶ |
| 0.01 | 0.0244683 | 3.54779 | 0.999366 | 0.999594 | none | none | μ²=0 only | 3.3×10⁻¹⁶ |
| 0.03 | 0.0428721 | 3.53231 | 0.998038 | 0.998744 | none | none | μ²=0 only | 4.4×10⁻¹⁶ |
| 0.1 | 0.0813286 | 3.47856 | 0.992741 | 0.995355 | none | none | μ²=0 only | 3.3×10⁻¹⁶ |

M̂ = M/(|BX_b| + |3λZ_b/ρ_b²|) ∈ [−1, 1]. An eigenvalue requires M̂ = 0. On the +1 branch M̂ stays within 0.73% of +1 at every tested δ. The gravitational mixing term 3λZ/ρ² is 6×10⁻⁶ (δ=0.0003) to 2.3×10⁻³ (δ=0.1) of BX at μ²=0.

Physically, the scalar sector is close to a heavy bulk scalar (U″(1)=28/9, i.e. Δ=18 in AdS₅ of radius 9) with a stabilizing shell term σ″/2≈2. It has no bound state below the continuum. Other diagnostics:

- Gauss-constraint residuals of the unreduced longitudinal integration are ≤ 3×10⁻¹⁵.
- The tensor zero mode's norm ∫ρ²dy is finite (e.g. 7.49×10⁴ at δ=0.001).
- Tensor precision spreads are ≤ 7.5×10⁻¹⁰, with the largest at μ²=2.2499 near threshold.
- The ℓ=1 condition value Bφ′_b is −1.4×10⁻⁴ (δ=0.0003) to −4.5×10⁻² (δ=0.1), which is nonzero.

## 7. Deliberate controls (`controls.py` → `CONTROLS_RESULTS.json`; `winding.py` → `WINDING_RESULTS.json`)

| Control | Result |
|---|---|
| C1 known limit: μ²=0 solution vs the Gegenbauer polynomial C₁₄⁽²⁾(cosh y/9) (decoupled AdS scalar) | Relative difference −3.52×10⁻⁵, −1.175×10⁻⁴, −3.53×10⁻⁴ at δ=0.0003, 0.001, 0.003. Difference/δ is constant to 0.16%, so the deviation vanishes linearly as δ→0, as expected from O(η_b) nonlinearity. Passes. |
| C2 detection power: wrong sign of σ″/2 on the +1 background | The wrong operator yields tachyons: μ² = −28555 (δ=0.001) and −283.5 (δ=0.1). The scan can see unstable modes on this background. |
| C3 wrong formula on the calibration: drop 3λZ/ρ² | The −7.718 root disappears; the gravitational term is essential |
| C4 small-δ limit of the calibration | Linear extrapolation of δ=0.001 and 0.003 roots gives −7.7197956 vs the closed-form 4D value −7.7197959 (diff 2.8×10⁻⁷) |
| C5 perturbed parameter c→c(1±10⁻³) | Original root shifts smoothly (dμ²/dc = −2.89). +1 branch re-solved: still no roots, min M̂ = 0.99994 |
| Winding calibration | Original shell: winding 1 on both contours. +1 branch: 0 (|w| < 3×10⁻¹⁵) at all δ |
| Renormalized-integrator consistency | Calibration root −7.7178716314 vs −7.7178716262 (direct) |
| Symbolic wrong-formula controls | 4/4 detected (§2) |

## Limitations

- **Floating point only.** Nothing here is interval-certified. The +1 conclusion rests on M̂ staying ≥ 0.9927 over a 360-point grid in [−400, 2.2499], on winding counts, and on analytic bounds that hold *given* the derived equations and B>0.
- **Grid coverage.** μ² < −400 and complex μ² outside the contours are covered only by the conditional theorems (B>0 ⇒ μ² > −4 and real spectrum). μ² → 9/4 was approached to 2.2499. The continuum itself is stable under the criterion but was not resolved spectrally.
- **Sectors not computed.** Vector perturbations were not computed. Only single-harmonic linear stability was tested: there is no nonlinear stability, basin of attraction, or proof that evolution from the original shell reaches this branch. The two configurations have cones near opposite vacua, and connecting them involves the outgoing wall and global structure.
- **Vacuum decay not studied.** U(+1) = −2/27 lies above U(−1) = −50/27, so the +1 vacuum is the higher of the two AdS critical points. Quantum decay, i.e. bubble nucleation near the shell, was not studied. For potentials of this superpotential form, positive-energy arguments protect both critical points in the bulk without branes ([Boucher, Nucl. Phys. B 242, 282](https://www.sciencedirect.com/science/article/abs/pii/0550321384903948); [Townsend, Phys. Lett. B 148, 55](https://www.osti.gov/etdeweb/biblio/5746769)). Whether this protection survives with the detuned shell has not been checked.
- **Interpretation of the criterion.** It uses the dS-invariant harmonic decomposition with regularity at the cone (horizon). This is standard in the cited literature but is an assumption about which states count.
- **Scope of claims.** The +1 branch being linearly stable makes it a consistent candidate late-time state of the matter-free model: an empty inflating brane. That strengthens the earlier "relaxation toward a de Sitter brane" reading. It is **not** evidence for a radiation era or for the higher-dimensional-origin hypothesis, and it is not a discovery claim.
- **Novelty.** The methods are standard. The gauge-invariant equations reproduce Chat 9's. The B>0 reality and lower bound are simple consequences of the Chat 9 identities; external novelty has not been assessed.

## Reproduction

Run the scripts with python3, numpy ≥ 2, scipy ≥ 1.16 and sympy ≥ 1.14, single-threaded. Everything runs in this folder and never writes outside it; bytecode writing is disabled so the package folders stay untouched.

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant
export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 derive_linearized.py      # ~12 s   -> DERIVATION_RESULTS.json (55 checks + 4 wrong-formula controls)
python3 reproduce_backgrounds.py  # ~6 s    -> BACKGROUNDS.json (imports package solvers read-only)
python3 spectrum.py               # ~4 min, 2 processes -> SPECTRUM_RESULTS.json (calibration + +1 branch)
python3 winding.py                # ~4 min, 2 processes -> WINDING_RESULTS.json
python3 controls.py               # ~5.5 min -> CONTROLS_RESULTS.json
```

`spectrum.py`, `winding.py` and `controls.py` read `BACKGROUNDS.json`, so run `reproduce_backgrounds.py` first. Logs from the recorded runs are in `*_log.txt`.

## Files

- `derive_linearized.py`: sympy derivation and all symbolic checks. Output `DERIVATION_RESULTS.json` holds the equations as strings.
- `gi_core.py`: backgrounds in shifted variables (exact rational potential coefficients), plus the gauge-invariant scalar solver, tensor solver, unreduced longitudinal cross-check and renormalized integrator.
- `reproduce_backgrounds.py`, `spectrum.py`, `winding.py`, `controls.py` and their JSON outputs, which record script SHA-256 hashes.

Other literature used for context: [Hawking, Hertog & Reall, Brane new world](https://arxiv.org/abs/hep-th/0003052) (dS brane bounding AdS balls, the same topology as the doubled bulk here); [BraneCode](https://arxiv.org/pdf/hep-th/0309001); a recent thick-brane stability study ([arXiv:2609.21421](https://arxiv.org/abs/2609.21421)), seen only as a search result and not read in full.
