# Independent review of the compensated-seed evolution setup

22 September 2026. Reviewed frozen `registered_solver.py`, `balanced_constraint_seed.py`, and the new `evolution/evolve_balanced.py` wrapper. This review executed initialization and algebra controls, not a PDE time evolution. The original checkpoint was left unchanged.

**Finding:** the wrapper correctly samples the compensated profile and matched zero-bump reference at the same shell-anchored conformal coordinate, and resets the frozen solver's background coefficients consistently. Initial acceleration checks converge under refinement. Important remaining qualifications concern interpolation noise, finite resolution of compact bumps, initial-mode interpretation, and the distinction between short evolution tests and a resolved nonlinear fate.

## 1. Coordinate matching is essential

For each shooting solution separately define

\[
z(y)=\int_{y_0}^{y}\frac{d\tilde y}{\rho(\tilde y)}-z_b,
\]

where `z_b` is the integral at that solution's brane. At each target evolution node `z_i`, invert that solution's own monotone map to obtain its own `y_i`. Comparing profiles at the same proper `y` would introduce a coordinate displacement of the same order as the desired physical perturbation.

For the earlier archived `epsilon=1e-4` candidate, `delta y_b≈-6.2310e-7` and `phi_y≈1`. The scalar correction from this coordinate shift is approximately `-6.2292e-7`, comparable to the actual brane scalar displacement `+8.9932e-7`. Ignoring coordinate matching would therefore significantly change the seed, rather than merely adding a tiny interpolation error.

The wrapper's use of independent dense-output inversions is correct. It also uses the shifted scalar `psi=phi+1` for scalar differences and `log1p((rho_seed-rho_ref)/rho_ref)` for the warp perturbation. These avoid avoidable subtraction errors. The equally-spaced-proper-distance archived profiles are suitable for scale estimates, not for linearly interpolating the finite-difference PDE's initial data.

## 2. Matched reference and gauge

The frozen solver originally built its background through a different static ODE integration and cone-series truncation. Those small differences could act as a false seed. The wrapper instead computes its zero-bump reference through the same mass-transport formulation as the compensated seed and replaces all relevant arrays and cached quantities: `rho`, `Hc`, `phi`, `phiz`, `rho^2`, potential derivatives, endpoint `rho_b/phi_b`, and tension derivatives.

The independent zero-seed control returns exactly zero perturbation arrays and exactly zero right-hand side. This confirms consistency of that subtraction path; it does not establish stability of nonzero perturbations.

Initial `A_t=1`, `B_t=phi_t=0` is retained. The reported `H/H0=(1+pa_b) exp(-b_b)` and `H0 tau=integral exp(b_b) dt` use the **reference** brane Hubble normalization `H0=1/rho_ref,b`. A compensated seed can have a different initial physical brane Hubble rate because its `rho_b` differs. That is consistent with the chosen initial data and should not be mistaken for an initialization bug.

## 3. Local causal scales and resolution

In the earlier small-amplitude transport profiles, the two compact bump supports map approximately to

\[
(-0.0179034,-0.00777863),\qquad
(-0.00676364,-0.00196988)
\]

in conformal `z`. The brane has an exactly static initial neighborhood only out to the nearest support edge. Conditional on the continuum initial-boundary problem being well posed, no nonstatic disturbance from the compact modification reaches the brane before approximately `t=0.00196988`. This is a short interval, not a Hubble-scale period.

For the newly chosen `epsilon=0.01` profile, the closest edge is approximately `z=-0.001963`. The earliest-arrival value obtained here is a grid-interpolation diagnostic; it is not an interval bound. The analytic expectation before arrival is constant `phi_b` and `B_b`, with `A_b` increasing at unit coordinate-time slope. Finite differences have truncation error and a wider numerical stencil domain, so this is an accuracy test rather than a requirement of literal zero numerical change at every step.

The upper bump occupies only about 22 near-shell cells at `hmin=2e-4`, and approximately twice that at `1e-4`. Its smooth, compact edges can require more resolution than a simple wall-width estimate suggests. A fixed fine region also does not automatically follow the later outgoing wall.

The wrapper tapers the seed in `[-0.99L,-0.85L]` and enforces zero deviation at the outer boundary. This modifies the exact constrained initial data in that outer region. The earliest physical influence of the taper on the shell is about `0.85L`; the wrapper correctly rejects requests extending to that time. Core diagnostics and an eventual domain-size comparison remain useful because finite-difference causality is approximate. The quintic taper is C2 when extended by constants, not C-infinity.

## 4. Independent initial-acceleration control

The mass-transport seed provides an analytic check distinct from merely evaluating its radial constraint ODE. Let `I=rho^2 D` be the integrated defect. At the initial slice the unmodified evolution equations predict

\[
(A_{tt},B_{tt},\phi_{tt})
=\left(-2I/\rho^2,\ 4I/\rho^2,\ \epsilon\rho^2 b\right).
\]

I independently regenerated dense profiles and compared these predictions with the finite-difference right-hand side on `z>-0.15`, where no taper is active.

| Initial control at epsilon=0.01 | hmin=0.0002 | hmin=0.0001 |
|---|---:|---:|
| Maximum absolute A acceleration error | 0.0001825 | 0.00003740 |
| Maximum absolute B acceleration error | 0.001285 | 0.0001160 |
| Maximum absolute scalar acceleration error | 0.07585 | 0.01440 |
| Maximum normalized initial Hamiltonian residual in the core | 1.8004e-7 | 2.3080e-8 |
| Initial momentum residual | 0 | 0 |

The corresponding maximum analytic acceleration magnitudes are approximately `(1.56, 3.12, 31.85)`. These observed reductions support the sampling and derivative implementation, but do not establish a precise asymptotic convergence order from only two grids.

Analytic endpoint slopes agree with target junction data to about `3.2e-13`. Exact endpoint accelerations vanish for precisely compensated data; the sampled direct scalar acceleration is about `5.35e-5` on the coarse grid and `3.60e-5` on the fine grid. This slower endpoint improvement is consistent with a dense-profile/differentiation error floor and deserves separate reporting. Over the approximately `0.002` pre-arrival interval its implied scalar displacement is of order `1e-10`. Much finer grids can amplify fixed interpolation or shooting error through second derivatives; tightening the spatial grid alone need not reduce that floor.

## 5. Mode interpretation and fit readiness

- Evolve the compensated nonlinear seed itself. Replacing it with its projection onto a linear eigenvector would discard the finite-amplitude constraint and corner construction.
- A Euclidean dot product with a right eigenvector is not generally the correct projection for this non-self-adjoint system containing physical, gauge, and constraint degrees of freedom. A justified projection would require a suitable left eigenvector or a separately derived physical-mode pairing.
- Vanishing integrated defect does not guarantee nonzero overlap with the unstable mode, and `sign(epsilon)` need not select the same future branch as the previous shifted-tension seed. Indeed the new positive `epsilon=0.01` seed has a negative initial brane scalar displacement. Measure the actual trajectory.
- Prefer reported scalar-growth windows after initial transients, signed-amplitude controls, and comparison across grids. Agreement with an expected growth rate is useful only alongside constraint convergence and small enough perturbations.
- Retain complete fields and velocities, absolute and normalized constraints, weighted characteristic residuals, and analytic-versus-independent boundary derivative defects. Small absolute errors in a tiny seed alone do not establish accurate linear dynamics.

## 6. Candidate difference equations if sampling cancellation becomes limiting

The current wrapper is adequate for the displayed initialization controls. If much smaller amplitudes are required, one can integrate differences directly rather than repeatedly subtract two independently approximated profiles. The following equations are an algebraically verified candidate, **not an implemented or validated replacement integrator**.

At common conformal `z`, let the reference quantities be `(rho0,v0,s0,psi0)` with `I0=0`; define

\[
a=\ln(\rho/\rho_0),\quad f=\psi-\psi_0,\quad k=s-s_0,
\quad E=e^a,\quad \Delta U=U(\psi_0+f)-U(\psi_0).
\]

In this subsection `U(psi)` denotes the model potential evaluated at `phi=psi-1`. Write `Delta U'` similarly. Then

\[
\begin{aligned}
\Delta X={}&\rho_0^2\left[(e^{2a}-1)
 (s_0^2/12-U_0/6)
 +e^{2a}\{(2s_0k+k^2)/12-\Delta U/6\}\right]
 -\frac{I}{\rho_0^2 e^{2a}},\\
\Delta v={}&\frac{\Delta X}{\sqrt{v_0^2+\Delta X}+v_0},\\
a_z={}&\Delta v,\\
f_z={}&\rho_0[\operatorname{expm1}(a)s_0+e^a k],\\
k_z={}&\rho_0[\operatorname{expm1}(a)U'_0+e^a\Delta U']
 -4[v_0k+\Delta v(s_0+k)]+\epsilon\rho_0e^a b,\\
I_z={}&\frac{\epsilon}{6}\rho_0^5e^{5a}(s_0+k)b.
\end{aligned}
\]

The rationalized square-root difference uses the positive radial branch and requires its radicand to remain positive. Polynomial `Delta U` and `Delta U'` can be evaluated by finite Taylor/Horner sums, as already done in the time-evolution solver. A proper-coordinate difference, if needed, satisfies `delta y_z=rho0 expm1(a)`.

The integration direction, regular-cone starting data, shell anchoring, and shooting parameters still need a stable treatment. Integrating backward from the brane can magnify unwanted cone solutions. The algebra alone does not settle these numerical issues. Standard cancellation-aware arithmetic is being applied here; no new mathematics is claimed.

## 7. Focused primary-source guidance

The relevant lesson from [Calabrese, Lehner and Tiglio's Einstein–scalar boundary study](https://arxiv.org/abs/gr-qc/0111003) is to distinguish incoming physical/gauge information from incoming constraint violations using characteristic propagation. Their formulation and geometry differ from this brane model; their successful evolutions do not validate the present boundary scheme.

[Frittelli and Gómez](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.69.124020) connect normal Einstein projections at boundaries to constraint preservation. This supports investigating a constraint-preserving discretization, while keeping the physical Israel/scalar junction conditions intact.

[Wang and Petersson's fourth-order wave method](https://arxiv.org/abs/1809.04310) provides a concrete summation-by-parts and ghost-point route with a discrete energy estimate. Its significance here is methodological: polynomial accuracy of a Hermite closure is not itself an energy-stability proof. Adapting such a construction to the coupled nonlinear brane boundary would be further work. These are established methods found in a bounded search, not new cosmological evidence.

## Verification files and scope

- `check_wrapper_initialization.py` regenerates the zero-seed and two-grid initial-acceleration controls in `WRAPPER_INITIALIZATION_CHECKS.json`; no time evolution.
- `verify_sampling_and_difference.py` verifies six exact paired-profile identities, rejects three incorrect variants, and records causal-scale diagnostics from five archived profiles in `SAMPLING_DIFFERENCE_CHECKS.json`.
- The small-amplitude source profiles and metadata are copied, unchanged, under `inputs/`, with hashes in the verification output. They are reference inputs distinct from the new epsilon=0.01 initialization calculations.

None of these controls establishes matter production, reheating, a resolved nonlinear end state, or an origin of the Big Bang.
