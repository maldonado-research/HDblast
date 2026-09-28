# Independent review of the constrained compact initial seeds

22 September 2026. Reviewed `constraint_seed.py`, `verify_constraint_seed_algebra.py`, `balanced_constraint_seed.py`, and their saved numerical results. This review did not edit those files, repeat their shooting solves, or evolve the initial data.

**Verdict:** the defect transport, initial accelerations, and second-time corner obstruction are algebraically correct. A precisely compensated bump gives a valid conditional route to all-order *local formal* compatibility at the brane. The supplied numerical candidates approximate the compensation and junction conditions to floating-point tolerances; they do not establish an exact compensated root or a subsequent nonlinear evolution.

## 1. Initial geometry and constraints

On the initial slice, take proper-distance coordinate `y`, set `dz=dy/rho`, and prescribe

\[
A=B=\ln\rho(y),\quad A_t=1,\quad B_t=\phi_t=0.
\]

Here `A=B` refers to the initial slice, not their proposed time dependence. Let `v=rho_y` and `s=phi_y`. The initial momentum constraint vanishes because `A_z=B_z=v`. The Hamiltonian constraint reduces to

\[
\rho_{yy}=\frac{1-v^2}{\rho}-\frac{\rho s^2}{6}-\frac{\rho U}{3}.
\]

The scalar profile is independently selected through

\[
\phi_{yy}=U'-4\frac{v}{\rho}s+\epsilon b(\phi).
\]

The additional term chooses an initial spatial profile. It is not an additional term in the subsequently intended Einstein–scalar evolution equations. No evolution has yet been demonstrated for these profiles.

Define

\[
D=1-v^2+\rho^2\left(\frac{s^2}{12}-\frac{U}{6}\right).
\]

`D` is a static radial first-integral defect, **not** the Hamiltonian constraint itself: the prescribed radial metric ODE already sets the initial Hamiltonian constraint to zero even when `D` is nonzero. Direct differentiation gives

\[
D_y=-2\frac{v}{\rho}D+\frac{\epsilon\rho^2s b}{6},\qquad
(\rho^2D)_y=\frac{\epsilon\rho^4s b}{6}.
\]

In the exact regular-cone limit `rho^2 D -> 0`, the compensated condition is therefore the vanishing of the weighted bump integral. At a finite numerical start `y0`, the identity includes the initial integration constant `rho(y0)^2 D(y0)`. Truncated cone-series data and floating-point integration do not make that constant identically zero, although it is small for the starts used here.

## 2. Accelerations and corner compatibility

Evaluating the **unmodified** evolution equations at the initial slice yields

\[
A_{tt}=-2D,\qquad B_{tt}=4D,\qquad\phi_{tt}=\epsilon\rho^2b.
\]

The initial target junction conditions are

\[
v=\frac{\rho\sigma}{6},\qquad s=-\frac{\sigma'}2.
\]

The first time derivatives of the junction residuals vanish because the initial `B_t` and `phi_t` are zero and `A_t=1` is spatially constant. This first-order result alone is insufficient for smooth evolution through the initial corner.

For completeness, differentiating the junction conditions twice gives the following general residuals, after imposing the initial junction conditions:

\[
\begin{aligned}
R_A^{(2)}&=0,\\
R_B^{(2)}&=-12vD+\epsilon\rho^3s b,\\
R_\phi^{(2)}&=-4\rho sD
 +\epsilon\rho^2b\left(2v+\frac{\rho\sigma''}{2}\right)
 +\epsilon\rho^3s b_\phi.
\end{aligned}
\]

Thus, when the bump is identically zero on an **open neighborhood** of the brane,

\[
\boxed{(R_A^{(2)},R_B^{(2)},R_\phi^{(2)})
=(0,-12\rho_yD_b,-4\rho\phi_yD_b).}
\]

Merely setting `b(phi_b)=0` is insufficient: a nonzero `b_phi(phi_b)` can still enter the scalar second-time corner residual. The selected compact profiles have a genuine open gap from their support to the brane, so the stronger condition applies.

## 3. What exact compensation would establish

If `b=0` on a connected open neighborhood of the brane and `D_b=0` **exactly**, then `D_y=-2vD/rho` implies `D=0` everywhere in that neighborhood. The metric equation becomes

\[
\rho_{yy}=-\rho\left(\frac{s^2}{4}+\frac{U}{6}\right),
\]

while the scalar equation becomes the ordinary static radial equation. These local profiles therefore furnish an explicit local static solution

\[
A(t,z)=t+\ln\rho(z),\quad B(t,z)=\ln\rho(z),\quad\phi(t,z)=\phi(z),
\]

with the same initial fields, velocities, and target junction conditions. All time derivatives of the junction residuals vanish for this local solution. This supplies all-order formal corner compatibility without separately computing an unbounded list of time derivatives.

It does **not** prove a smooth global initial-boundary solution, its numerical convergence, its long-time behavior, or heating. It also does not say the perturbed interior remains static. Disturbances from the compactly modified interior may subsequently reach the brane; the exact local static solution supplies compatibility and, conditional on a well-posed hyperbolic initial-boundary problem, a finite local domain of dependence.

## 4. Reading the saved compensated candidates

The two-bump profile has support at shifted scalar `psi=phi+1` in `(0.20,0.50)` and `(0.55,0.85)`. The saved branes are near `psi=1.00002`, leaving an open scalar gap of about `0.15002`. The solved compensation coefficient is approximately `0.40933` and changes smoothly over the displayed signed amplitudes.

Both numerical resolutions solve the two target junction conditions to roughly `10^-14`. The integrated defect contribution `I_b/rho_b^2` is around `10^-14`, but the independently evaluated direct defect is about `1.8e-11` to `5.7e-11`. Their difference exposes the integration/cancellation floor; the integral value must not be substituted for the direct defect without reporting it. The direct second-time corner residuals are about `6e-9` to `2e-8`, rather than exactly zero.

The numerical solve tunes `I_b=0`, not a certified exact `D_b=0`. It is reasonable to call these **floating-point compensated initial-data candidates**, with separately reported junction, direct-defect, and corner residuals. It would overstate the evidence to call them an exact all-order compatible family or a proven evolved solution. The conditional analytic implication in Section 3 is exact; existence of the required compensated shooting root has only been supported numerically here.

The positive-only bump is useful as a negative control: its nonzero weighted integral produces a nonzero `D_b` and the predicted second-time obstruction. For example, the saved refined `epsilon=1e-6` single-bump case has `D_b≈1.51958e-4` and `R_B^(2)≈-0.047923`. The compensated candidates reduce this defect to the numerical floor. That is a concrete improvement in the initial-data design, not a cosmological outcome.

## 5. Independent verification

`verify_seed_compatibility_independent.py` re-derives the radial identities and differentiates the complete junction twice, including general nonzero bump terms. It does not import either producer. It passes **13 positive symbolic identities**, rejects **four incorrect variants**, and checks the brane gap and junction tolerances in **ten saved compensated candidates** (four signed amplitudes and the zero reference, each at two tolerances). Source hashes and results are in `SEED_COMPATIBILITY_CHECKS.json`.

The checks of saved numerical data are read-only consistency checks. They do not independently rerun the shooting algorithm, prove interval bounds, or evolve the PDE.

Release verification note: the final saved-data check was rerun after the zero-reference row was added. All ten checks passed. The mass-transport formulation is assessed separately by the sampled-profile checks; the ten checks above cover the nominal and refined second-order formulations.
