# Linear stability of the static "+1 branch": time-domain and discrete-operator study (Method B)

HDBLAST research program (Ricardo Maldonado). Work package dated 27–28 September 2026. It uses floating-point numerics only. Nothing here is an interval certificate, and the work has not been externally reviewed.

## Question

The 22 September checkpoint (Zenodo record 22922928) found a static regular-cone solution of the registered model, the "+1 branch": φ_b = 0.9999159473169134, ρ_b = 129.9247628497, H² = 5.92401479433×10⁻⁵ and H/H₀ = 0.606721732. It left the stability of that solution open. Is the solution linearly stable, marginal or unstable under the registered 5D Einstein–scalar dynamics? What are its leading physical growth or decay rates?

The method has to be calibrated first, by recovering the known growth rate of the original unstable shell.

## Model and scope

This study uses the unchanged registered model: W = 1 − φ + φ³/3, U = ½W_φ² − (2/3)W², σ = 2W + δ(1 + cφ), δ = 0.001 and c = 0.5975949350280132. The bulk is doubled (Z2) with one shell. The metric is e^{2B}(−dt² + dz²) + e^{2A}dx₃², with the shell at z = 0 and the bulk at z < 0.

The PDE system, gauge, stretched grid, quintic-Hermite shell closure and nonlinear right-hand side are those of the frozen `registered_solver.py`. A copy is in `frozen_input/`; its SHA-256 `4147239c…cd6b6` is checked on import and matches the checkpoint's `frozen/registered_solver.py`. Only the static background is new. `td_core.plus_background` re-solves the +1 branch with a second-order (ρ, ρ_y, η, η_y, z) system and both junctions, without the first-integral square root. It reproduces the archived ρ_b to a relative difference of 3×10⁻¹⁵ and φ_b to all printed digits. Junction residuals are 6×10⁻¹⁷ and 1×10⁻¹⁹; on the evolution grid they are ≤ 6.5×10⁻¹⁵.

**Sector tested:** only perturbations that keep the brane homogeneous and isotropic (the 1+1 reduction that the registered solver evolves). Inhomogeneous (k ≠ 0) brane perturbations and tensor modes are **not** tested.

**Time unit:** the coordinate time t is measured in the Hubble time of the background being perturbed (1/H = ρ_b). On the +1 branch, a rate s per unit t equals s·H₊ = 0.6067·s·H₀.

## Method

1. **Discrete operator.** The linear operator of the registered PDE system is `Solver.linear_matrix()`, inherited unchanged. It was checked against a complex-step Jacobian of the unchanged nonlinear RHS on every grid, with relative error ≤ 3×10⁻¹⁶. The frozen outer nodes were removed.
   - Dense eigendecompositions were computed at h_min = 4, 3 and 2×10⁻⁴ (up to 7,188 unknowns). Further dense runs varied L = 6, 8 and 10 and the stretch parameter (0.05 or 0.1).
   - Shift-invert Arnoldi was run at h_min = 2×10⁻⁴, 1×10⁻⁴ and 5×10⁻⁵, with 12 shifts on the real axis from 1.66 to −4. It was also run at h_min = 1×10⁻⁴ with L = 8 and at δ = 0.003 and 0.01.
2. **Mode classification.** Each eigenpair was classified with four diagnostics:
   - (a) The linearised Hamiltonian and momentum constraints, relative to the sum of the absolute values of their terms.
   - (b) The two shell observables that are invariant under the residual shell-preserving conformal gauge: the shell scalar f_b and the shell Hubble perturbation (ȧ − b)_b.
   - (c) The least-squares residual against the analytic residual gauge mode ξ^t = e^{λt}cosh λz, ξ^z = e^{λt}sinh λz. This mode generates δa = cosh λz + H_c sinh λz, δb = λ cosh λz + H_c sinh λz and δf = φ_z sinh λz.
   - (d) Dependence on the outer boundary through the far-region weight and the L scan.
3. **Time domain.** The evolutions start from time-symmetric, constraint-satisfying initial data (`td_evolve.two_bump_data`).
   - The scalar disturbance is two smooth compact bumps in the bulk. a = b is obtained by solving the discrete linearised Hamiltonian constraint with the solver's own stencils and shell closure.
   - The second bump's weight is tuned so that a = 0 at the shell. With f = 0 near the shell, this removes the junction-compatibility (corner) obstruction identified in the 22 September controls. The momentum constraint holds identically.
   - Evolution uses RK4 with dt = 0.4 h_min, with either the verified linear matrix or the **unchanged nonlinear RHS** at small amplitude ε.
   - Rates were extracted in three ways: least-squares fits of ln|f_b|; matrix-pencil (Hua–Sarkar) poles of f_b(t); and, on the +1 branch, the decay of a positive weighted scalar energy E = Σ dz ρ³(ḟ² + f_z² + ρ²U″f²). The amplitude rate is half the slope of ln E.
   - Resolutions used: +1 branch at 4, 2, 1 and 0.5×10⁻⁴; original shell at 2, 1 and 0.5×10⁻⁴. There were two initial profiles, two amplitudes, L = 6 and 10, and constraint monitoring.
4. **Fully-discrete check.** The effective RK4 rate ln|R(λ dt)|/dt was computed for every dense eigenvalue.
5. **Controls.**
   - Calibration on the original shell.
   - Wrong-formula control: the sign of σ″(φ_b) flipped in the scalar junction.
   - A σ″ = 0 control.
   - Perturbed-parameter runs at δ = 0.003 and 0.01.
   - L and stretch variation.
   - A potential-precision control (see Limitations).
   - Nonlinear-versus-linear amplitude scaling.

## Results

| # | Result | Status | Evidence (JSON) |
|---|---|---|---|
| 1 | **Calibration.** On the original shell, the time-domain matrix-pencil growth rate is **1.657193368** (linear, h_min = 5×10⁻⁵) and **1.657193408** (nonlinear, ε = 10⁻⁸). The reference is 1.6571936312 from the independent Gaussian-normal shooting (Chat 9), a difference of −2.6×10⁻⁷ and −2.2×10⁻⁷. The sequence at h_min = 2, 1 and 0.5×10⁻⁴ is 1.657043, 1.657193 and 1.657193. Late-window log fits differ by −8.6×10⁻⁷. The semi-discrete eigenvalue is 1.6572110, 1.6571944 and 1.6571937 on the same grids. | numerical | `results/KEY_RESULTS.json` → `calibration_*` |
| 2 | **No growing physical mode on the +1 branch.** At every resolution, L and stretch tested, the only eigenvalues with Re λ > −3/2 are of three kinds. (i) λ = 0.9999966 (L = 6), which approaches 1 as L grows (0.99999985 at L = 8, 0.99999991 at L = 10); it matches the analytic gauge template to 1.5×10⁻⁸ and has f_b ≤ 4×10⁻¹² (0.26 and 0.016×10⁻¹² at finer h), so it is a **pure gauge mode**. (ii) Far-boundary modes with λ = 0.1040, 0.0753 and 0.0594 for L = 6, 8 and 10 (≈ 0.6/L), their complex ladder, and their mirror partners at −3 − λ̄; they have f_b ≤ 10⁻¹¹ and far-region weight ≈ 1, so they are **artifacts of the truncated domain (gauge or constraint sector)**. (iii) Grid-scale oscillations with \|Im λ\| ≥ 3,200 that scale as 1/h; they are **numerical**, and RK4 damps them. Under RK4 the fully-discrete maximum growth rate is exactly the gauge value 1.000 on every +1 grid. | numerical | `plus_dense`, `plus_shift_invert_modes_Re_gt_-3`, `far_boundary_mode_vs_L` |
| 3 | **Every shell-supported mode lies on Re λ = −3/2** to within 1×10⁻⁸ (599 to 1,198 upper-half-plane eigenvalues per grid). This is the de Sitter relation λ(λ+3) = −μ² with μ² > 9/4, a discretised continuum. The shift-invert scans at h_min = 5×10⁻⁵ find no real eigenvalue in (−3, 1) with shell support. So **no scalar bound state with μ² < 9/4 was found**, whether tachyonic or light. | numerical | same |
| 4 | **The time-domain rate converges to −3/2.** The amplitude rate from the scalar energy is −1.5370, −1.5393, −1.5062 and −1.5011 at h_min = 4, 2, 1 and 0.5×10⁻⁴ (observed order ≈ 2.7; Richardson estimate **−1.50012**). A second, shell-localised initial profile gives −1.5062, −1.5012 and −1.5003 (Richardson estimate **−1.50017**). L = 10 gives the same as L = 6 to 1×10⁻⁴. The nonlinear RHS at ε = 10⁻⁶ and 2×10⁻⁶ agrees with the linear evolution; the f_b/ε difference is 1.8×10⁻⁷ and 3.6×10⁻⁷, scaling linearly with ε, which is the expected quadratic effect. The observed departure from −1.5 is numerical dissipation. The compensated-energy check E_X loses 62%, 34%, 11% and 6.5% over t = 6 on the four grids, converging slowly. | numerical | `plus_td_*`, `nonlinear_vs_linear` |
| 5 | **The shell scalar decays promptly.** The pulse's impact on the gauge-invariant shell scalar peaks at an RMS of 3.6×10⁻⁹ (relative to the unit bulk bump), identically on all grids. It falls to about 2×10⁻¹⁴ in the window t ∈ [0.4, 0.5] (2.7 and 1.5×10⁻¹⁴ on the two finest grids). After that the signal sits at a floor that shrinks under refinement: 2.2×10⁻¹³, 2.1×10⁻¹³, 5.1×10⁻¹⁴ and 3.9×10⁻¹⁵ at t ≈ 6. No resolution-independent persistent or growing component was seen at the shell. | numerical | `plus_td_late_shell_scalar_rms_last_window`, `f_b_rms_windows_0p1` in ANALYSIS.json |
| 6 | **Wrong-formula control** (σ″ → −σ″ in the scalar junction on the +1 background). Both methods detect a strong instability: λ = **167.48925** from shift-invert at h_min = 1 and 0.5×10⁻⁴ and from the time-domain pencil (dense h4: 167.48899). Its partner is −170.48899, so the pair sums to −3 exactly. Setting σ″ = 0 leaves the branch stable (time-domain rate −1.5062). The method therefore can see an instability on this background when one is present, and the stability depends on the sign of the brane mass term σ″(φ_b) = 4φ_b > 0. | numerical (control) | `control_*` |
| 7 | **Perturbed parameter.** At δ = 0.003 and δ = 0.01 the +1 branch shows the same structure: only gauge and far-boundary modes above the line, with f_b < 10⁻¹¹. The time-domain rate at δ = 0.003 is −1.5048 at h_min = 10⁻⁴. | numerical | `perturbed_parameter_*` |
| 8 | **Constraint monitoring.** Near the shell (z > −0.2), the linearised constraints of the final state converge under refinement: the relative Hamiltonian residual is 1.4×10⁻², 2.6×10⁻³, 1.6×10⁻⁴ and 3.2×10⁻⁵. **In the far stretched region (z < −1), whole-domain residuals stay at O(0.1–0.6) and do not converge.** Outgoing high-frequency radiation (ω ≈ 200–1,500) is under-resolved where the grid spacing grows. This is the same limitation the earlier checkpoints found with outgoing fronts. The nonlinear runs start with an O(ε²) Hamiltonian mismatch, because the data solve only the linear constraint. That mismatch decays to 3×10⁻⁷ relative, independently of ε. | numerical / limitation | `plus_td_near_shell_constraints_final`, `plus_td_far_region_constraints_final` |

All 24 automated checks in `results/KEY_RESULTS.json` pass. One check was first written with a tolerance of 10⁻⁹, which is below the dense-eigensolver accuracy (eigen-residuals ~2×10⁻⁸). It failed; this is recorded under Limitations.

## Conclusion

For homogeneous (FRW-symmetric) perturbations in the registered model, the +1 static branch is **linearly stable** (numerical result). No growing or slowly decaying physical mode was found. The leading physical rate is the continuum edge, **Re λ = −3/2 in units of the branch's own Hubble rate** (≈ −0.910 H₀). All physical scalar-sector modes found have Re λ = −3/2 exactly, corresponding to μ² > 9/4. The time-domain rate extrapolates to −1.5001 to −1.5002 from two initial profiles at four and three resolutions.

This is ordinary damping of heavy modes on a de Sitter brane, not marginality. There is no discrete bound state below the μ² = 9/4 threshold. In our evolutions the brane-localised scalar signal leaves the shell within about 0.4 Hubble times.

A plausible mechanism, supported by the controls but not proved here: the bulk scalar near φ = 1 is heavy (U″(1) = 28/9 against an AdS curvature scale of 1/9), and the brane mass term σ″(φ_b) = 4φ_b > 0 is positive. Flipping its sign produces a tachyon at λ ≈ 167.

The only positive rates in the discrete system are the gauge mode (λ → 1 as L → ∞) and far-boundary artifacts that shrink as ≈ 0.6/L. Neither affects the gauge-invariant shell observables.

**What this does not establish:**
- stability against inhomogeneous or tensor perturbations;
- nonlinear stability, or a basin of attraction;
- that the Chat 14 roll-off actually reaches this branch;
- any rigorous bound.

The late-time dark-radiation-like mode at λ ≈ −4.000 (localised at the far boundary) and the far-region constraint non-convergence are noted but not resolved.

## Established literature versus new results

A spectral gap of (3/2)H between the zero mode and the Kaluza–Klein continuum on an inflating (de Sitter) brane is established for RS-type brane worlds ([Garriga & Sasaki, hep-th/9912118](https://arxiv.org/abs/hep-th/9912118)). The decay e^{−3Ht/2} of fields with μ² > 9/4H² in de Sitter is textbook. The first-order "superpotential" construction behind W and U follows [DeWolfe, Freedman, Gubser & Karch, hep-th/9909134](https://arxiv.org/abs/hep-th/9909134).

New to this project, and nothing beyond it:
- the numerical stability classification of this specific +1 branch of the registered model;
- the time-domain calibration against the Chat 9 rate;
- the constraint-satisfying two-bump linear data on the evolution grid;
- the explicit gauge-mode template check.

No novelty beyond this project is claimed.

## Limitations

- **Sector:** homogeneous perturbations only. The spectrum is truncated at z = −L with a Dirichlet outer boundary. The true cone is at z → −∞, and the continuum is represented by discrete box modes.
- **Resolution:** the stretched grid resolves the shell region well and the far region poorly for high frequencies. Hence the non-converging far-region constraints, the slow E_X convergence, and the unconverged classification of far-region box modes. Constraint labels ("constraint-violating"/"physical-candidate") in the raw spectra JSON are threshold-based and **not** reliable for far-region box modes. Items 2–3 rely on shell observables, gauge-template fits and refinement, not on those labels.
- **Classification threshold:** "shell-supported" means \|f_b\|/max\|a, b, f\| > 10⁻⁶ (or 10⁻³ where stated).
- **Initial data:** the data satisfy the linearised constraints exactly on the grid, and the nonlinear constraint only to O(ε²). A radiated O(ε²) Hamiltonian mismatch is unavoidable for time-symmetric data on this branch, because the scalar sources the metric only through the tiny background gradient.
- **Precision:** the registered φ-polynomial evaluation of U derivatives loses accuracy where |φ − 1| < 10⁻¹² (the deep bulk of the +1 branch). This study evaluates them as exact polynomials in η = φ − 1 (`td_core.eta_derivatives`, identical to 1.4×10⁻¹⁴ for moderate η). The registered evaluation (`--phi-poly`) gives the same eigenvalues above the line to < 10⁻⁷ but shifts a far-boundary partner mode (−3.104) by 2.3×10⁻⁶. The first version of this check, with tolerance 10⁻⁹, failed and was replaced by a 10⁻⁷ check on eigenvalues above the line, justified by the dense eigen-residuals.
- **ARPACK:** the first control shift-invert search (k = 6, real shifts 50–200) did not converge. It was stopped and rerun with k = 2 at shift 167.5, which converged. The failing configuration was replaced in `run_shift.sh`.
- **Runs discarded:** a preliminary set of evolutions, from before energy recording was added and with t_f = 3.5 for the original shell, was superseded by the batches below and deleted. A stray run from that set, finished by a process that was still running, was likewise removed.

## Files

- `td_core.py`: frozen-solver import with hash check, +1 background, general solver, linear and nonlinear constraints, gauge template, classification.
- `spectra.py`: dense and shift-invert spectra.
- `td_evolve.py`: constraint-satisfying data and RK4 evolutions.
- `analyze.py` writes `results/ANALYSIS.json`; `summarize.py` writes `results/KEY_RESULTS.json` (24 checks).
- `results/spectra/*.json|npz`, `results/td/*.json|npz`: raw outputs. Every run records its arguments, grid, background meta and runtime.
- `logs/`: batch logs.
- `frozen_input/`: frozen solver, archived +1 branch results, archived spectra and the Chat 9 reference rate. These are copies; the originals were not modified.

## Reproduction (Linux, python3, numpy 2.3.5, scipy 1.16.3; single-threaded)

```bash
cd research/HDBLAST_CHECKPOINT_20260927/stability_time_domain
./run_spectra.sh     # dense spectra, ~15 min total (h2 runs ~4 min each)
./run_td.sh &        # batch A: calibration + controls, ~22 min
./run_td2.sh         # batch B: +1 branch evolutions, ~23 min
wait
./run_shift.sh       # shift-invert scans, ~2 min
python3 analyze.py   # -> results/ANALYSIS.json
python3 summarize.py # -> results/KEY_RESULTS.json, prints pass/fail
```

Single run examples:
`python3 spectra.py --bg plus --hmin 4e-4 --mode dense --out results/spectra/plus_h4_L6`
`python3 td_evolve.py --bg plus --hmin 1e-4 --tf 6 --linear --out results/td/plus_lin_h1e-4`
