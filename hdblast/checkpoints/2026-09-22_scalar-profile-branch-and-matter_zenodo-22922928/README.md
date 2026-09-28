# HDBLAST: a corrected static branch and a consistent matter extension

Research checkpoint, 22 September 2026, America/Los_Angeles. New calculations, independent agent checks, and primary-source review. No publication or external peer review is claimed.

**The strongest new result is a correction to the proposed late-time endpoint.** The constant bulk field φ=1 used for the earlier Randall–Sundrum comparison fails the scalar boundary condition of the registered model. We calculated a nearby nonconstant scalar profile that satisfies both boundary conditions, independently reproduced it with a different form of Einstein's equations, and derived its small-detuning expansion.

We also evolved the balanced initial disturbances, found a remaining numerical limitation in the outgoing disturbance, and derived a matter-production extension with consistent gravitational backreaction and energy exchange. This is substantive progress in testing the hypothesis. A radiation-filled universe has not been demonstrated.

## 1. New static solution candidates in the unchanged model

The registered five-dimensional Einstein–scalar model has

\[
W=1-\phi+\phi^3/3,\quad U=\tfrac12W_\phi^2-\tfrac23W^2,
\quad \sigma=2W+\delta(1+c\phi),
\]

with δ=0.001 and c=0.5975949350280132. The scalar junction condition is nφ=−σ′/2. At φ=1, σ′=δc≠0. A bulk scalar identically equal to one has nφ=0, so it fails this condition. Satisfying the gravitational junction alone does not produce a solution of the full model.

Using η=φ−1, we solved the full nonlinear static boundary-value problem from a regular cone to the shell. The tiny cone displacement must be represented directly: at the registered setting it is about −1.34×10⁻²³, which would disappear if stored as 1+η in ordinary floating-point arithmetic.

| Registered candidate | Result in the model's dimensionless units |
|---|---:|
| Shell scalar φ_b | 0.9999159473169134 |
| Shell radius ρ_b | 129.9247628497 |
| Expansion rate squared H² | 0.0000592401479433 |
| H/H₀, relative to the original unstable shell | 0.606721732 |
| Cone displacement η_h | −1.336692365×10⁻²³ |

The producer found candidates at six detunings from 0.0003 to 0.1, with a tolerance comparison for each. A separate checker made 18 integrations using the second-order Einstein equation rather than enforcing the producer's radial first integral. Its 28 checks passed. Maximum junction residuals were approximately 3.61×10⁻¹⁶ for the metric and 1.57×10⁻¹⁴ for the scalar. These are measured floating-point residuals, not certified error bounds.

The small-detuning expansion is

\[
\boxed{\phi_b=1-\frac{9c}{64}\delta+O(\delta^2),}
\]

\[
\boxed{H^2=\frac{\delta(1+c)}{27}
+\delta^2\left[\frac{(1+c)^2}{36}-\frac{c^2}{384}\right]+O(\delta^3).}
\]

The coefficient of the scalar-profile correction was checked algebraically and against the numerical scan. The regular linear radial solution is a degree-14 Gegenbauer polynomial in cosh(y/9); it provides an accurate starting point for the nonlinear solve.

This correction changes the expansion-rate benchmark by only about **7.87 parts per million** at δ=0.001. It does not explain the much larger finite-time gap in the older roll-off runs. Its significance is that the candidate now satisfies the missing boundary condition.

**What remains open:** stability, uniqueness, rigorous existence bounds, and dynamical attraction. The new candidate's cone is near φ=+1, while the original unstable shell's cone is near φ=−1. Connecting the two configurations requires demonstrating what happens to the departing wall and the global spacetime; a static calculation does not establish that connection.

Read the [independent branch review](static_branch/INDEPENDENT_BRANCH_REVIEW.md), [numerical results](static_branch/PLUS_BRANCH_RESULTS.json), and [matter report's analytic derivation](matter/MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md).

## 2. Reconciliation with the newly completed Chat 14

Your new reference correctly points to additional work completed after our previous source audit. The supplied earlier package is unchanged: its SHA-256 matches our delivered archive. The external Chat 14 archive now contains the completed registered-detuning runs; all 84 of its manifest payloads pass hash verification.

Its finite-time reversal is supported across three grids: H₀τ≈5.8944–5.8947, at φ_b≈−1.9412 to −1.9417. On the positive side, H/H₀≈0.638642 at H₀τ=6.9 is a useful late comparison, still above the new static candidate's rate.

Some source conclusions remain stronger than the stored evidence. The final positive value near 0.628690 occurs after signals from the seed taper and the far boundary could have reached the shell. The terminal collapse is not converged, and the archive lacks full velocity-field snapshots and Hamiltonian histories needed for an independent full-constraint audit. Agreement among several runs does not resolve those omissions. The archive supports finite-time tendencies and reversal, not a proven final singularity or a settled asymptotic endpoint.

The copied source archive and the [independent Chat 14 audit](source_audit/COMPLETED_CHAT14_INDEPENDENT_AUDIT.md) preserve the chronology, exact versions, event calculations, and limits. Earlier documents' instructions were treated as reference material; no publication recommendation in them was executed.

## 3. What happened when the balanced disturbances were evolved

We moved beyond initial-data construction and ran finite-amplitude evolutions. The initial profiles are sampled at the same shell-anchored conformal coordinate, and the reference background and all cached coefficients use the same zero-bump construction. This avoids comparing different physical radial locations or silently mixing two background integrations.

The disturbance amplitude ε is separate from the fixed model detuning δ. Tested amplitudes include ε=0.01 and ε=±0.001. The sign of ε is not itself a physical prediction of which final branch is selected.

The main finding is a limitation worth preserving: **satisfying the tested initial boundary conditions does not by itself fix the evolution solver.** On the grid concentrated near the shell, an outgoing constraint-error front becomes large in the coarse region. Expanding the region with fine spacing reduces the error substantially. At coordinate t=1, refining that wider grid reduces the reported background-normalized Hamiltonian maximum from about 0.3005 to 0.0340. A near-shell value alone would conceal the larger outgoing error.

The independent review also normalizes residuals by the magnitudes of the actual nonlinear terms. It confirms that the problem is not merely an unfortunate background normalization. A three-grid comparison at exactly coordinate t=0.5 shows self-convergence of the fields and velocities, while the maximum Hamiltonian residual improves only from about 0.00754 to 0.00655 on the last refinement. This is encouraging early-time field convergence together with an unresolved constraint-accuracy limit. It does not establish globally controlled late-time evolution, and no final fate is inferred from these short runs.

An initialization-only experiment also isolated one source of error: a tiny uncertainty in converting radial coordinates becomes larger when the fields are differentiated twice on a fine grid. Tightening that conversion reduces an outer initial residual, but leaves the core and weighted maxima unchanged. No improved-tolerance time evolution is claimed.

Read [initialization checks](evolution_review/INITIALIZATION_AND_EVOLUTION_REVIEW.md) and [the evolution residual review](evolution_review/EVOLUTION_RESIDUAL_REVIEW.md). Raw outputs retain all six fields and sampled Hamiltonian and momentum constraints, so this limitation can be investigated rather than hidden.

## 4. A matter extension that can actually be tested

We specified a shell scalar χ with a scalar-dependent mass, and an optional decay channel. Varying the same doubled-bulk, single-shell action gives

\[
nA=\frac{\sigma+\kappa_5^2\rho}{6},\qquad
nB=\frac{\sigma-\kappa_5^2(2\rho+3p)}{6},\qquad
n\phi=-\frac{\sigma'+\kappa_5^2j}{2},
\]

\[
\boxed{\dot\rho+3H(\rho+p)=j\dot\phi.}
\]

Here ρ and p are actual shell-matter density and pressure, j is the scalar source obtained from that matter action, and dots denote shell proper time. Matter changes both gravitational boundary conditions and the scalar boundary condition. Its energy gain must be balanced by the same coupling in the bulk/shell system.

The derivation passed 47 exact checks and four deliberately incorrect-formula controls, followed by a separate action and geometry review. The resulting four-dimensional identities still contain bulk scalar and Weyl information; they do not close the coupled problem by themselves. A radiation equation appended to the old vacuum trajectory would not be a self-consistent solution.

Established instant-preheating physics supplies a possible nonoscillatory production mechanism. We screened five hypothetical crossing locations using three archived grids and two derivative windows. Among the tested locations, crossings near φ=0.5 and 0.75 require the least demanding coupling ratios for the specified local approximation checks. These are conditional screening results, not a chosen physical coupling, a particle-production simulation, or proof of thermalization. [Felder, Kofman and Linde](https://arxiv.org/abs/hep-ph/9812289).

A positive remaining vacuum energy creates another requirement: radiation must dominate it for a sufficient duration, not merely be nonzero. The report derives conditional radiation thresholds with the high-energy density-squared correction retained. It also explains why a simple fixed-Planck-mass four-dimensional energy budget cannot be imposed blindly across these five-dimensional configurations. The formal particle-energy proxy used in such estimates is explicitly distinguished from a renormalized stress tensor during production.

Read the [full matter derivation and primary-source review](matter/MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md), [independent review](matter/INDEPENDENT_MATTER_REVIEW.md), and [crossing-screen results](matter/CROSSING_ELIGIBILITY_AND_TOY_BUDGET.json).

## 5. The next discriminating tests

1. Test the perturbation spectrum of the new static branch and its relation to the original cone and outgoing wall. Numerical existence alone does not make it the endpoint.
2. Control both constraints along the moving disturbance, with further spatial and initial-profile convergence and a benchmarked late-time coordinate system. Shell observables alone are insufficient.
3. Select explicit physical scales and matter couplings, evolve quantum mode functions with a consistent renormalized stress/source, and feed all three modified junctions back into the bulk.
4. Demonstrate energy transfer, thermalization and a sustained radiation era before computing observational predictions.

This checkpoint supplies new derivations and reproducible numerical candidates for this project. External mathematical novelty and a higher-dimensional origin of our universe remain unestablished. All calculations retain Einstein's equations within the stated model; matter is an explicitly separate extension.

## Package scope

The source review was targeted to the relevant new archives, nearby research summaries and references. The web review used primary papers and checked recent applicable work; it was not an exhaustive search of the web or all computer files. Original research folders, prior delivered packages and synced project sources were not edited.

The [reproduction guide](REPRODUCE.md) distinguishes saved-data checks, independent integrations, initial-data reconstruction, and full PDE reruns. The delivery archive includes the source snapshots, scripts, numerical results, independent reviews and a payload manifest. Passing package checks establishes precisely those checks, not a hot Big Bang mechanism.
