# Independent residual audit of the compensated-seed evolutions

22 September 2026. Analysis of seven completed runs using saved full fields and velocities. The residual reconstruction itself is read-only. A separately documented ten-step replay reconstructs an unrecorded coarse state at the exact comparison time; no additional complete PDE run was performed. An initialization-only sampling-tolerance experiment is described below.

**Verdict:** the full-domain nonlinear constraint gate remains open. The large outgoing-front errors in the narrow-refinement runs are not merely a consequence of normalizing by the old background. Refining a much wider region substantially reduces the errors, but the largest constraint errors nearly stall in the third refinement. The finest run retains actual-term front residuals of about 0.12% in H and 0.67% in M at t=0.5. Clean near-brane evolution alone does not establish a resolved five-dimensional solution throughout the retained causal domain.

## 1. What was recomputed

The saved perturbations `(a,b,f,pa,pb,pf)` and reference arrays give

\[
A_t=1+p_a,\quad B_t=p_b,\quad\phi_t=p_f,
\quad A_z=H_c+a_z,\quad B_z=H_c+b_z,\quad\phi_z=\phi_{s,z}+f_z.
\]

I independently constructed polynomial ghosts, evaluated the spatial differences explicitly, and reconstructed all terms of the nonlinear Hamiltonian and momentum constraints:

\[
\begin{aligned}
H={}&-2e^{2B}U+6A_t^2+6A_tB_t-12A_z^2+6A_zB_z
      -6A_{zz}-\phi_t^2-\phi_z^2,\\
M={}&-3A_{tz}-3A_tA_z+3A_tB_z+3A_zB_t-\phi_t\phi_z.
\end{aligned}
\]

The producer's background denominator remains `S0=1+6Hc^2+phi_s,z^2`. I also report

\[
r_H=\frac{|H|}{\sum_i|H_i|},\qquad
r_M=\frac{|M|}{\sum_i|M_i|},
\]

where the sums contain the eight and five separate terms displayed above. These are transparent measures of the failure of the actual nonlinear terms to balance. They are **additional diagnostics**, not replacement normalizations, unique physical norms, or certified solution-error estimates. Changing the grouping of terms would change these denominators; the exact grouping is recorded in the JSON.

The analysis uses the producer's interior mask `z>-0.8L` with the first six points excluded, and reports the near-brane core `z>-0.2` separately. It retains the independent one-sided velocity derivative at the brane. Recomputed constraints agree with the saved deviation-form constraints to about `2.2e-8` absolute Hamiltonian error or better and `4e-12` momentum error or better in these runs, far below the problematic front residuals. Tiny differences are expected from explicit full-term cancellation and separately assembled ghost arithmetic.

## 2. Residual sizes with both denominators

The percentages below are separate maxima of the actual-term ratios, so the Hamiltonian and momentum maxima need not occur at exactly the same node.

| Run | End time | Original max `abs(H)/S0` | Actual-term max `r_H` | Actual-term max `r_M` |
|---|---:|---:|---:|---:|
| epsilon=0.01, narrow refinement, h=0.0002 | 2.5 | 1.1378 | 25.68% | 41.01% |
| epsilon=0.01, narrow refinement, h=0.0001 | 2.5 | 1.7380 | 34.79% | 51.56% |
| epsilon=0.01, wide refinement, h=0.0002 | 1.0 | 0.30047 | 5.124% | 8.757% |
| epsilon=0.01, wide refinement, h=0.0001 | 1.0 | 0.034033 | 0.7766% | 2.357% |
| epsilon=0.01, wide refinement, h=0.00005 | 0.5 | 0.0065507 | 0.1231% | 0.6729% |
| epsilon=+0.001, wide refinement, h=0.0002 | 1.0 | 0.025596 | 0.5784% | 1.782% |
| epsilon=-0.001, wide refinement, h=0.0002 | 1.0 | 0.024588 | 0.5622% | 1.680% |

The narrow map uses stretch `0.05`, `L=6`; the wide map uses stretch `2`, `L=3`. These maps have different spacing away from the brane even at the same `hmin`. The comparison therefore tests the importance of resolving the outgoing region, not just the boundary spacing.

At the original Hamiltonian-error peak in the narrow fine run, `z≈-2.465`, the actual-term residuals are approximately `34.79%` in H and `51.56%` in M. Thus the principal negative conclusion survives the change of denominator. In the wide fine run the original H peak is near `z=-0.99276`, with `0.7766%` H and `2.357%` M using actual terms. The wide coarse run has a slightly different node maximizing its actual-term H ratio; all peak locations and term values are retained in the JSON.

The core improves much more strongly. At the end of the wide fine run its maximum actual-term ratios are about `2.53e-8` in H and `3.51e-7` in M. In the narrow fine run at t=2.5 they are about `1.60e-5` and `9.37e-5`. These support a better-resolved core while leaving the outgoing bulk region unresolved at the tested accuracy.

Reducing epsilon by ten on the coarse wide map reduces the absolute/background-normalized front defect by approximately an order of magnitude, with similar results for both signs. This is consistent with a substantial linear numerical-error component. It does not establish convergence to the continuum solution: smaller absolute perturbations naturally produce smaller absolute residuals, and the actual-term ratios still require resolution and amplitude interpretation.

## 3. Transport versus continued numerical generation

For the exact nondissipative continuum equations,

\[
(\partial_t-\partial_z)C_+=0,\qquad
(\partial_t+\partial_z)C_-=0,
\quad C_\pm=e^{3A}(H\pm2M).
\]

The script uses the same densities divided by the constant `rho_ref,b^3`, which does not alter their transport. All runs here use zero KO coefficient, so there is no intentional KO source term to add. Spatial/time truncation, boundary error, and profile approximation can still generate a discrete constraint residual.

In the narrow maps, the maximum weighted density grows far above its initial value: approximately **1,270 times** on the coarse grid and **668 times** on the fine grid by t=2.5. Moreover the late weighted maximum is near the brane (`z≈-0.00267` coarse and `-0.04445` fine), rather than at the far outgoing front. Pure continuum transport and geometric amplification of the initially measured constraint residual cannot explain these increases. There is continuing finite-resolution constraint generation or boundary inconsistency in addition to transport; this analysis does not uniquely isolate the offending operator or distinguish diagnostic truncation from evolved constraint error at every node.

For the wider runs, the end-time weighted maxima are **below** their initial maxima: ratios about `0.78` coarse, `0.45` fine, `0.66` for epsilon=+0.001, and `0.32` for epsilon=-0.001. Their front is overwhelmingly outgoing: at the weighted front peak, the incoming density is only about `0.3–0.6%` of the outgoing density. This is consistent with strong transport of early numerical error toward the cone.

However a weighted maximum that decreases does not prove the absence of later numerical generation. Transport comparisons between saved snapshots show shape/phase changes, dispersion, and nonzero differences from simply advecting the previous density. These comparisons also involve spatial interpolation of a narrow wave, and are diagnostic rather than exact source integrals. At several peaks the backward initial characteristic foot lies outside the initial domain, so the peak could involve early boundary-generated information or numerical dispersion; a direct initial-value comparison there is deliberately omitted from the JSON.

Thus: **the narrow runs show more than geometric amplification of initial error; the wide runs substantially improve that behavior but do not yet pass the global accuracy gate.**

## 4. Three-grid comparison at exactly t=0.5

The archived coarse state is at t=0.4992. To avoid comparing different physical times, `matched_time_refinement.py` reconstructs its solver from the saved background arrays and advances **exactly ten original RK4 steps of dt=0.00008** to t=0.5. No new background solve, seed, discretization, or full run is introduced. Its output `MATCHED_TIME_COARSE_STATE.npz` and receipt record the original data and solver hashes. The two finer states are read directly at t=0.5.

The three maps are exactly nested: strides 1, 2, and 4 select identical coarse-node coordinates. There is no spatial interpolation in the state comparisons. At this matched time:

| Spacing | Maximum original Hamiltonian ratio | Maximum original momentum ratio |
|---:|---:|---:|
| 0.0002 | 0.0446611 | 0.0221036 |
| 0.0001 | 0.00754067 | 0.00374193 |
| 0.00005 | 0.00655073 | 0.00326694 |

The largest constraint errors improve strongly in the first refinement and then nearly stall. The evolved-field coordinate differences nevertheless decrease. Using separate coordinate L2(dz) norms on the same coarse nodes gives:

| Component | Coarse–middle difference | Middle–fine difference | Reduction factor |
|---|---:|---:|---:|
| A-t | 7.2875e-8 | 4.7711e-9 | 15.27 |
| B | 4.7579e-8 | 4.3041e-9 | 11.05 |
| phi | 5.8786e-7 | 7.5486e-8 | 7.79 |
| A_t | 5.3094e-6 | 4.1532e-7 | 12.78 |
| B_t | 2.0266e-5 | 4.7345e-6 | 4.28 |
| phi_t | 0.00151585 | 0.000342577 | 4.42 |

These are useful self-convergence diagnostics. They are not a physical energy norm, a certified error bound, proof of an asymptotic order, or proof that the limiting fields satisfy all constraints. The derivative-based constraint stall remains material even though the sampled fields become closer.

## 5. Initial inverse-map tolerance experiment

`check_inverse_tolerance.py` overrides the sampling function only within its own process. All six initialization cases share the same two dense radial solutions and keep the ODE tolerance at `8e-14`. It compares the original inverse-map `brentq` absolute tolerance `2e-14` with `1e-300`, retaining the legal relative tolerance `4*machine_epsilon`. No producer file is modified and no time evolution is executed.

The original initial global H residual worsens by approximately four as the grid spacing halves, with its peak near z=-2.37 on the first two grids and z=-2.30 on the finest. Tightening the inverse tolerance reduces this outer error while leaving the analytic brane junction data unchanged:

| Spacing | Original initial global H ratio | Tighter-inversion global H ratio | Core H ratio, unchanged |
|---:|---:|---:|---:|
| 0.0002 | 3.0668e-7 | 1.4800e-7 | 1.4800e-7 |
| 0.0001 | 1.3147e-6 | 3.5264e-7 | 1.0064e-7 |
| 0.00005 | 5.1897e-6 | 1.4706e-6 | 5.6996e-7 |

The core residuals and initial weighted maxima are exactly unchanged in all three comparisons. Other outer sampling/differentiation error remains and still worsens with refinement. Full values are in `INVERSE_TOLERANCE_CHECKS.json`; the test therefore confirms one contributing inverse-sampling floor without establishing it as the entire source of the evolved front error.

The mechanism is numerically plausible and directly tested: changing an initial warp field by only about `3.5e-14` can appreciably change a second derivative on a fine grid. Consequently refinement alone cannot be assumed to improve constraints once it differentiates a fixed sampling-error floor. A small inverse residual rounded to zero is not an exact inversion certificate.

## 6. Next gate

The next useful improvement is initial-profile sampling accuracy or a directly evolved profile-difference representation, followed by a matched-time constraint check. The tolerance experiment does not establish that improved initialization alone will resolve the outgoing front, and no improved-tolerance PDE run has been performed here. A grid-wide stability analysis or a discretization with a constraint/energy estimate remains appropriate if the front errors persist.

No result here establishes the late roll-off fate, a singularity, visible particle production, or reheating. The useful result is a more precise numerical diagnosis and an explicit record of why the global gate remains open.

## Reproduction

Run `analyze_evolution_residuals.py` with NumPy. It reads completed `runs/balanced_*.json` and corresponding NPZ files, reconstructs the constraints and term sums, and writes `EVOLUTION_RESIDUAL_REVIEW.json`. It executes no ODE background shooting and no PDE evolution. Source hashes, actual snapshot times, peak locations, complete term decompositions, initial/previous characteristic comparisons, and reconstruction discrepancies are recorded for each included run.

`matched_time_refinement.py` additionally uses SciPy to reconstruct the original coarse spatial operators and perform only the ten specified RK4 steps; it writes `MATCHED_TIME_COARSE_STATE.npz` and `MATCHED_TIME_REFINEMENT.json`. `check_inverse_tolerance.py` uses SciPy for the two shared radial shooting solutions and six initialization checks; it performs no PDE evolution and writes `INVERSE_TOLERANCE_CHECKS.json`. The three analyses have deliberately distinct scopes.
