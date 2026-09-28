# Constraint control in the balanced-disturbance evolutions

HDBLAST research program (R. Maldonado). Workstream checkpoint, 27–28 September 2026.
Every number below comes from a script in this folder and has a machine-readable source: `ANALYSIS.json`, `CONTROLS.json`, or `runs/**/*.json`. These are floating-point numerical experiments. None of them is an interval certificate.

## Question

The 22 September checkpoint (Zenodo 22922928) evolved the compensated two-bump seed (ε=0.01, "wide" map with stretch 2 and L=3). It reported an outgoing Hamiltonian-constraint error front whose maximum barely improved on the finest grid: at t=0.5 it went from 0.00754 at h=1e-4 to 0.00655 at h=5e-5. This workstream asked four things:

1. Can that behaviour be reproduced?
2. Where does it come from?
3. Is there a change that makes the Hamiltonian residual converge at t=0.5 and t=1?
4. At what convergence order?

## Short answer

- **Reproduced (numerical).** The copied checkpoint wrapper gives H_max(t=1, h=2e-4) = 0.300596, against 0.300474 in the archive (0.04 % difference). The new laboratory code agrees with the copied wrapper to a relative 7.7e-10. The baseline stall reappears on four grids: at t=0.5 the observed order between h=1e-4 and h=5e-5 is 0.24, and at t=1 it is 0.27.
- **Diagnosis (numerical, with an exact transport law behind it).** The front is not generated in the outer or coarse region. Almost all of it is made in the first ~0.05 time units, within about 0.02 of the shell, where the steep wall (φ_z≈79, A_z≈17) meets the edges of the compact bumps.
  - Once made, the outgoing density e^{3A}(H+2M) is conserved along z≈−t, which is the exact continuum identity from folder 152.
  - The unweighted, background-normalised H therefore grows by the geometric gain e^{-3t}[ρ(z0)/ρ(z0−t)]³. From z0=−0.003 this gain is about 2.4×10³ at t=0.5 and 5.8×10³ at t=1.
  - A relative constraint defect of order 10⁻⁸ near the shell is enough to produce the reported percent-level errors in the bulk.
- **Why the stall (numerical).** Near the shell the seeded outgoing density has two parts:
  - a truncation part, which converges at roughly order 3;
  - a floor from the initial data, which does *not* decrease with h. The discrete initial H of the sampled profiles near the bump edge rises from 2.9e-4 at h=1e-4 to 5.5e-4 at h=5e-5, in weighted units.

  At h=5e-5 the floor dominates.
  - The floor is not caused by the coordinate-inversion tolerance: xtol=1e-300 gives bit-identical diagnostics near the shell.
  - It is not caused by the time step: with CFL 0.2 the outgoing density at t=0.1 is the same.
  - It is removed by solving the discrete constraint for the initial data, which is Remedy A below.
- **Remedies that worked (numerical):**
  - **A. Discrete-constraint projection of the initial data.** The Newton solve takes 2–3 iterations. It changes the warp by at most about 5e-8 and leaves the scalar data untouched. It removes the stall: at t=0.5 the observed orders are 2.23, 2.78, 2.91, and at t=1 they are 2.22, 3.15, 3.13. On the finest grid H_max(t=0.5) drops from 6.38e-3 to 9.35e-4.
  - **B. Outgoing-characteristic constraint damping.** A source −κ(H+2M)/(6(A_t+A_z)) is added to B_tt, switched off smoothly within 0.02 of the shell.
    - SymPy verifies the damped law exactly: (∂_t−∂_z)[e^{3A}C₊] = −κ e^{3A}C₊.
    - With κ=10 this lowers H by factors of 13 to 100.
    - Without A it still stalls at the finest pair at t≤0.5.
  - **A+B together** gives the smallest residuals: H_max(t=0.5) of 5.4e-3, 1.2e-3, 1.8e-4 and 2.3e-5, with orders 2.2, 2.7, 3.0.
  - **C. Order-6 interior with projection** gives the highest observed orders at the finest pair: 3.85 at t=0.5 and 3.67 at t=1 for the maximum; 3.97 and 3.91 for the causal-window L2 norm.
- **Correct order (numerical, partly conditional).**
  - A smooth, corner-free calibration pulse gives exactly order 4.00 and 3.97 for H in the bulk until it reaches a floor of about 1e-9 to 1e-10. So the scheme is fourth order.
  - The seed runs show orders rising towards 4 (2.2 → 2.8 → 3.0 for A; up to 3.85 for C), but four grids do not reach the asymptotic regime. The scalar velocity φ_t self-converges at only about 2.1–2.4 in every variant.
  - **An asymptotic order of 4 is therefore *not* demonstrated for the seed runs.** What is demonstrated is a convergent residual with order about 3 or better, and 3.7–3.9 for Remedy C.
- **Remedies that failed or did not help (negative):**
  - Tighter coordinate inversion near the shell: no change.
  - Kreiss–Oliger dissipation 0.2: the early outgoing density fell only from 2.2e-4 to 1.8e-4, and with projection it stayed at 2.2e-5.
  - Order 6 without projection: worse initial noise.
  - Uniform damping without the shell switch-off: a boundary instability, with H of order 56 at t=1.
  - The naive "H-only" damping: unstable, as predicted analytically.
  - A 6th-order Hermite shell closure: amplifies roundoff below h≈0.01.
  - Projection alone at h≥2e-4: no gain, because truncation dominates there.

None of this bears on the physical fate of the shell, reheating, or the Big Bang hypothesis. It is a numerical-method result about controlling constraint error in this particular 1+1 model.

## Setup and definitions

- **Model.** The registered 5D Einstein–scalar model and junctions exactly as in the checkpoint: W=1−φ+φ³/3, δ=0.001, c=2/1.0357712571566784−4/3.
  - Metric: ds² = e^{2B}(−dt²+dz²) + e^{2A}dx₃².
  - Perturbation variables (a,b,f,pa,pb,pf), measured about the zero-bump reference.
  - Seed: ε=0.01, L=3, stretch=2. The outer C² taper sits on [−0.99L, −0.85L].
- **Constraints.** H and M are the checkpoint's residuals.
  - Background normalisation: S0 = 1+6H_c²+φ_{s,z}².
  - "H_max" is max|H|/S0 over z>−0.8L, excluding the first six nodes. This is the checkpoint mask.
  - The "causal" window additionally excludes the future of the taper region, z < −0.85L + t + 0.05. It differs from the original mask only at t=0.25 on the finest grids, where the incoming taper violation enters the mask.
  - Weighted characteristic densities: C₊^w = (ρ/ρ_b)³e^{3(t+a)}(H+2M), and C₋^w likewise with H−2M.
- **Time and grids.** RK4 with Δt ≤ 0.4h_min, where 1/Δt is a multiple of 4 so that t = 0.25, 0.5, 0.75 and 1 are hit exactly. The grids h = 4e-4, 2e-4, 1e-4 and 5e-5 are exactly nested from the shell, so self-differences need no interpolation.

## Files

| File | Role |
|---|---|
| `src/` | Verbatim copies of the checkpoint's `registered_solver.py`, `balanced_constraint_seed.py`, `constraint_seed.py`, `evolve_balanced.py` (SHA-256 in each run JSON). Never edited. |
| `lab.py` | Rebuilt operators (orders 2/4/6 with a degree p+1 Hermite shell ghost); constraint damping `kappa` with a C² shell switch-off (`damp_off`); KO dissipation (frozen-solver form, applied to the velocity equations only); inverse-map tolerance; `project_initial_hamiltonian` (Newton, analytic sparse Jacobian); low-memory construction of the frozen operators. |
| `run_case.py` | Runs one evolution. Output: JSON records every 0.025 time units, plus NPZ fields and constraints at t = 0.25, 0.5, 0.75, 1. |
| `calib.py` | Smooth, corner-free order calibration: a Gaussian bulk pulse whose warp comes from integrating the continuum Hamiltonian constraint (DOP853, rtol 1e-12). |
| `controls.py` → `CONTROLS.json` | 17 pass/fail controls (all PASS) plus 2 recorded findings. |
| `analyze.py` → `ANALYSIS.json` | Every table below: orders, self-differences, transport gains, early densities, damping controls, calibrations, reproduction. |
| `reproduce_all.sh` | Regenerates everything sequentially, in roughly 4–5 h on one core. |
| `runs/final` | Main matrix: 5 variants × 4 resolutions, t ∈ [0,1]. |
| `runs/short`, `runs/damp`, `runs/ctrl`, `runs/calib`, `runs/calib2`, `runs/calib3`, `runs/baseline_original`, `runs/recheck` | Short diagnostics, damping controls, calibrations, the copied-wrapper reproduction, and reproducibility checks. |

## Results

### 1. Reproduction of the baseline (numerical)

| Quantity | Value |
|---|---|
| Archived checkpoint H_max(t=1), h=2e-4 wide | 0.300474 |
| Copied wrapper rerun here | 0.300596 (0.04 % from archive; floating-point/library differences) |
| This folder's `lab.py` (order 4, κ=0) | 0.300596; relative difference from the copied wrapper 7.7e-10 |
| Baseline H_max(t=0.5), h = 4e-4, 2e-4, 1e-4, 5e-5 | 0.245, 0.0449, 0.00755, 0.00638 (archive: —, 0.0447, 0.00754, 0.00655) |

Default-path runs are bit-identical across the development versions of `lab.py`: records and final states match exactly at h=4e-4 (`runs/recheck`). The finest baseline run (h=5e-5, t=1, 17.8 min) is new.

### 2. Hamiltonian residual H_max at t=0.5 and t=1, observed orders log2(E_h/E_{h/2}) (numerical)

| Variant | t | h=4e-4 | 2e-4 | 1e-4 | 5e-5 | Orders (successive pairs) |
|---|---|---:|---:|---:|---:|---|
| Baseline (checkpoint scheme) | 0.5 | 2.45e-1 | 4.49e-2 | 7.55e-3 | 6.38e-3 | 2.45, 2.57, **0.24** |
| | 1.0 | 1.33 | 3.01e-1 | 3.36e-2 | 2.78e-2 | 2.14, 3.16, **0.27** |
| A: discrete projection | 0.5 | 2.27e-1 | 4.83e-2 | 7.01e-3 | 9.35e-4 | 2.23, 2.78, 2.91 |
| | 1.0 | 1.25 | 2.70e-1 | 3.03e-2 | 3.46e-3 | 2.22, 3.15, 3.13 |
| B: C₊ damping κ=10 (off near shell) | 0.5 | 5.96e-3 | 1.16e-3 | 1.78e-4 | 6.50e-5 | 2.36, 2.70, **1.45** |
| | 1.0 | 4.34e-2 | 1.15e-2 | 1.19e-3 | 1.43e-4 | 1.92, 3.27, 3.06 |
| A+B | 0.5 | 5.39e-3 | 1.17e-3 | 1.81e-4 | 2.28e-5 | 2.20, 2.69, 2.99 |
| | 1.0 | 4.34e-2 | 1.15e-2 | 1.19e-3 | 1.44e-4 | 1.92, 3.27, 3.05 |
| C: order 6 + projection | 0.5 | 6.88e-2 | 9.73e-3 | 1.21e-3 | 8.41e-5 | 2.82, 3.01, **3.85** |
| | 1.0 | 2.42e-1 | 3.59e-2 | 4.54e-3 | 3.56e-4 | 2.75, 2.99, **3.67** |

The momentum residual behaves the same way, with M_max ≈ H_max/2 at the front. For example, the baseline M order is 0.24 at t=0.5 for the finest pair, and A+B gives 2.94.

Causal-window L2 orders for the finest pair:

| Variant | t=0.5 | t=1 |
|---|---:|---:|
| Baseline | 0.58 | 0.63 |
| A | 3.48 | 3.38 |
| B | 1.71 | 3.10 |
| A+B | 3.40 | 3.10 |
| C | 3.97 | 3.91 |

The pointwise self-difference of H on nested nodes gives orders 2.94–3.74 for A, A+B and C at the last pair.

Runs take 22–27 s at h=4e-4, 76–96 s at 2e-4, 283–355 s at 1e-4, and 998–1347 s at 5e-5. The nodes are 4583, 9164, 18327 and 36653.

### 3. Mechanism: where the front is seeded (numerical)

Weighted outgoing density max|C₊^w| near t=0.05, just after the near-shell passage. These values are from `early_outgoing_density_t0.05`.

| h | Baseline initial | Baseline t≈0.05 | Projected (A) t≈0.05 | Order 6 + A t≈0.05 |
|---:|---:|---:|---:|---:|
| 4e-4 | 8.6e-3 | 1.6e-2 | 1.45e-2 | 1.14e-2 |
| 2e-4 | 1.2e-3 | 2.6e-3 | 2.3e-3 | 9.2e-4 |
| 1e-4 | 2.9e-4 | 3.0e-4 | 2.6e-4 | 5.3e-5 |
| 5e-5 | **5.5e-4** | **2.2e-4** | **2.3e-5** | **4.0e-6** |

- In the baseline, the initial discrete constraint defect is largest at the bump edge nearest the shell, around z≈−0.0024. It *grows* between h=1e-4 and h=5e-5, which is a sampled-profile noise floor.
- Afterwards C₊^w stays of the same size as it moves out; for example it is 1.2e-4 at t=0.5 on the finest baseline grid.
- The unweighted H front then follows the background gain, tabulated in `transport_gain_unweighted`: 3016 at t=0.5 and 7345 at t=1 for a ray leaving z0=0.
- Projection removes the floor. What remains converges at order 3.2–3.5 at t≈0.05.

Short diagnostics at h=5e-5 (t≤0.1), giving C₊^w at t=0.05:

| Change from baseline | C₊^w(0.05) | Result |
|---|---:|---|
| None | 2.22e-4 | — |
| xtol=1e-300 | 2.22e-4 | identical at printed precision, so **negative** |
| CFL 0.2 | 3.95e-4; at t=0.1 1.81e-4 vs 1.82e-4 | not limiting by t=0.1 |
| KO 0.2 | 1.80e-4 | **negative** (marginal) |
| KO 0.2 + projection | 2.17e-5 | same as projection alone |
| Order 6 without projection | 1.78e-4 (initial 1.29e-3) | **negative** |
| Projection | 2.29e-5 | improvement |
| Order 6 + projection | 4.04e-6 | best |

### 4. Damping laws and controls (exact-verified and numerical)

- **Exact (SymPy, `CONTROLS.json`).** With the source s = −κ(H+2M)/(6(A_t+A_z)) added to B_tt:

  (∂_t−∂_z)[e^{3A}C₊] = −κ e^{3A}C₊
  (∂_t+∂_z)[e^{3A}C₋] = −κ e^{3A}C₊ (A_t−A_z)/(A_t+A_z)

  Wrong sign and the wrong denominator (A_t−A_z) are rejected. κ may depend on z, because no derivative of it appears. This builds on the standard constraint-propagation identity in folder 152. The damping term is a project-specific application of standard constraint-damping ideas, not a new method.
- **Numerical, at h=4e-4, uniform κ.** The front ratio at t=0.25 compared with baseline is 0.090 for κ=+10, against e^{-2.5}=0.082. For κ=−5 it is 3.31, against e^{1.25}=3.49, and at t=0.5 it is 11.2, against e^{2.5}=12.2.
- **Negative results.**
  - Uniform κ=10 becomes unstable at the shell boundary: a grid-scale pb sawtooth grows at rate ≈2κ, giving H≈56 at t=1.
  - H-only damping, s = −κH/(6A_t), stops on the amplitude guard. A constant-coefficient analysis predicts growth ∝ κ(A_z/A_t−1)/2 where A_z>A_t, which holds near this shell.
  - The damping term contains A_zz and A_tz, so it modifies the principal part. The shell switch-off (κ=0 for z>−0.02, full κ for z<−0.05) is an empirical stabilisation. It is not a proof of well-posedness.
- **κ=30 with switch-off.** At t≥0.5 the damped runs keep a residual of about equal parts C₊ and C₋. It is generated at the pulse in the bulk and is not reduced by larger κ: κ=30 gives H_max(t=1)=2.1e-2, against 4.3e-2 for κ=10 at h=4e-4. That part converges with h (Table 2).

### 5. Order calibration (numerical)

`calib.py` uses the ε=0 reference plus a Gaussian scalar pulse (amplitude 0.01, width 0.02, centre z=−1.4). The warp comes from the continuum constraint, so the shell data are exactly static and there is no corner problem. The window is kept outside the taper's future and away from the shell.

| t | H_max h=4e-4 | 2e-4 | 1e-4 | 5e-5 | Orders |
|---|---:|---:|---:|---:|---|
| 0.25 | 1.68e-6 | 1.05e-7 | 6.70e-9 | 4.93e-9 | **4.00, 3.97**, 0.44 (floor) |
| 0.5 | 4.22e-6 | 2.65e-7 | 1.67e-8 | 1.65e-9 | **3.99, 3.99**, 3.33 |

Field self-difference orders at t=0.5 on the first three grids are 3.96–4.01 for all six fields. At t≥0.75 the pulse left in the window is below a floor of about 1e-10. With width 0.05 the floor is reached already at h=2e-4.

**Negative control:** a *single* uncompensated Gaussian at z=−0.3 with projected warp gives H_max(t=0.5) = 0.098, 0.222 and 0.782 at h = 4e-4, 2e-4 and 1e-4. It *grows* under refinement, and the peak sits exactly on the corner characteristic z=−t. This agrees with the corner obstruction documented in folder 152 and shows why only compensated seeds can be tested for convergence.

### 6. Controls (`CONTROLS.json`, status PASS)

- Low-memory operators, order-4 operators, RHS, constraint arrays and KO term all equal the frozen solver to roundoff (relative ≤ 1e-13).
- Order-2 and order-4 stencils, including the shell ghost, converge at orders 2.0 and 4.0. A 1 % wrong shell slope destroys boundary D2 convergence, as it should.
- The order-6 first derivative converges at 6.0.
- Projection brings the initial H down to about 1e-14, keeps M exactly 0, and leaves ε=0 data exactly zero. The correction is 4.7e-8, 6.9e-9, 1.5e-9 and 1.5e-9 at the four resolutions.
- **Finding:** the degree-7 shell closure for order 6 loses D2 accuracy to roundoff below h≈0.01. At the evolution spacings, Remedy C's shell closure is therefore not truly sixth order. Its observed order near 4 is consistent with that.

## Result labels

| Claim | Label |
|---|---|
| Baseline stall reproduced (orders 0.24 at t=0.5 and 0.27 at t=1 for the finest pair) | numerical |
| Front seeded near the shell at t≲0.05, then carried out by e^{3A}(H+2M) with geometric gain | numerical (the transport law is exact-verified in folder 152) |
| Stall caused by a non-decreasing initial-data defect floor near the shell (not xtol, not Δt) | numerical |
| Damped transport law for the C₊ source | exact-verified (SymPy, with negative controls) |
| Remedy A restores a convergent residual (orders ≈2.9–3.1 at the finest pair) | numerical |
| A+B: smallest residuals, 2.3e-5 at t=0.5 and 1.4e-4 at t=1 on the finest grid | numerical |
| Remedy C: orders 3.85/3.67 (max) and 3.97/3.91 (L2) at the finest pair | numerical |
| Scheme is fourth order on smooth, corner-free data | numerical (calibration, before the floor) |
| Asymptotic fourth order for the ε=0.01 seed | **inconclusive**; φ_t self-converges at only ≈2.1–2.4 |
| Uniform damping, H-only damping, KO, tighter inversion, order 6 without projection | negative |
| Stability of the switched-off damping beyond t=1, or for other seeds | not tested |

## Limitations

- **Evolved time.** All runs stop at t=1, well inside the small-amplitude regime: the maximum |f| is about 0.011 and the shell scalar deviation is about −0.003. Nothing here addresses the nonlinear roll-off, the fate of the shell, particle production or reheating.
- **Order.** The seed's compact bumps are C^∞ but not analytic, with strong gradients at the edges, so four resolutions sit in a pre-asymptotic regime. An order of about 3–3.9 is measured, not proven. The φ_t field converges more slowly than H. Its origin, whether near the shell or in the data smoothness, was not isolated.
- **Projection.** Remedy A changes the initial warp by a discrete-constraint correction (≤4.7e-8, not decreasing below h=1e-4) and leaves the scalar data unchanged. At h≥1e-4 this is a truncation-size change; at h=5e-5 it removes the sampling-noise defect of about 1.5e-9. It is not a new solution of the continuum problem. It is the unique discrete-constraint solution closest to the sampled seed under the stated boundary conditions: far-left values fixed and the junction slope at the shell.
  - Solving the discrete problem on a *local* window, rather than up to the taper, created a kink at the window edge and was discarded.
  - Using the projection for uncompensated data does not fix corner incompatibility (§5).
- **Damping.** Remedy B changes the principal part. Stability is empirical: t≤1, κ ∈ {10, 30}, one switch-off profile, four grids.
- **Runtime environment.** Two finest-grid jobs were killed by the memory cgroup (OOM) while running concurrently. The cause was a dense n×n temporary in the frozen `R.matrices` (LIL assignment of `sparse.eye`), about 10.7 GB at h=5e-5. `lab.py` now builds the same operators in COO form, checked equal to roundoff, and all affected runs were repeated. The sparse projection solve now uses natural ordering; this changed projected runs by about 1e-10 relative (`runs/recheck`). The logged failures are in `runs/queue_exit_codes.txt`.
- **Precision.** Double precision throughout, with ODE profiles at rtol 8e-14. The calibration floor of about 1e-9 to 1e-10 in normalised H, and the h=5e-5 initial-data floor, show that roundoff and profile noise are within a factor of about 10–100 of the best residuals reported here.

## Reproduction

Environment: Python 3.11.15, NumPy 2.3.5, SciPy 1.16.3, SymPy 1.14.0, run single-threaded (`OMP_NUM_THREADS=1`). Run from this folder:

```sh
# full regeneration (roughly 4-5 h on one core)
./reproduce_all.sh
# or individual pieces, e.g.
python3 -B run_case.py --hmin 1e-4 --out runs/final/base_h1e-4                                  # baseline
python3 -B run_case.py --hmin 1e-4 --project --out runs/final/proj_h1e-4                        # remedy A
python3 -B run_case.py --hmin 1e-4 --kappa 10 --damp-off 0.02 0.05 --out runs/final/k10_h1e-4   # remedy B
python3 -B run_case.py --hmin 1e-4 --project --kappa 10 --damp-off 0.02 0.05 --out runs/final/k10p_h1e-4
python3 -B run_case.py --hmin 1e-4 --project --order 6 --out runs/final/o6p_h1e-4              # remedy C
python3 -B calib.py --w 0.02 --hmin 1e-4 --out runs/calib3/g_h1e-4                            # calibration
python3 -B controls.py      # -> CONTROLS.json
python3 -B analyze.py       # -> ANALYSIS.json (reads runs/)
```

Copied-wrapper reproduction:

```sh
python3 -B src/evolve_balanced.py --epsilon .01 --hmin .0002 --stretch 2 --L 3 --tf 1 --output runs/baseline_original/balanced_epsp01_wide_h2
```

## Literature context (standard methods only)

- Constraint damping by adding lower-order constraint terms: C. Gundlach, G. Calabrese, I. Hinder, J. M. Martín-García, "Constraint damping in the Z4 formulation and harmonic gauge", Class. Quantum Grav. 22 (2005) 3767, https://arxiv.org/pdf/gr-qc/0504114. Constraint-damped generalized harmonic evolution: F. Pretorius, https://arxiv.org/pdf/gr-qc/0407110.
- Kreiss–Oliger dissipation as a standard high-frequency filter in numerical relativity, for example GRChombo, https://arxiv.org/pdf/1503.03436.
- The damping source used here is a project-specific application of these standard ideas to the checkpoint's 1+1 system. No methodological novelty beyond this project is claimed.
