# Constraint-compatible initial data: a corner obstruction and a compensated construction

22 September 2026. New calculations in the unchanged registered bulk and brane model, δ=10⁻³, d=0, c=2/1.0357712571566784−4/3. This note concerns initial data only. No nonlinear spacetime evolution is performed here, and no particle production or cosmological claim follows.

## 1. Initial constraint equations

Use coordinate time T to distinguish it from the detuning. At T=0 choose

\[
A=B=\log\rho(y),\quad A_T=1,\quad B_T=\phi_T=0,\quad dz=dy/\rho.
\]

Then A_z=B_z=ρ_y, A_zz=ρρ_yy, and φ_z=ρφ_y. The momentum constraint vanishes. The Hamiltonian constraint reduces to

\[
6-6\rho_y^2-6\rho\rho_{yy}-\rho^2\phi_y^2-2\rho^2U=0,
\]

or

\[
\boxed{\rho_{yy}=\frac{1-\rho_y^2}{\rho}-\frac{\rho\phi_y^2}{6}-\frac{\rho U}{3}.}
\]

To parameterize initial profiles, choose

\[
\phi_{yy}+4\frac{\rho_y}{\rho}\phi_y-U_\phi=\epsilon b(\phi).
\]

This is an initial-data selection equation, **not an alteration of the scalar's subsequent equation of motion**. The profile b is smooth and compactly supported away from the regular cone and brane. Integrate from the original regular cone series and adjust φ_h and y_b to satisfy the original brane junctions

\[
\rho_y/\rho=\sigma/6,\qquad \phi_y=-\sigma_\phi/2.
\]

The zero time derivatives of B and φ, and spatially constant A_T, satisfy the first time derivative of all three junction conditions. Thus the construction removes the old initial jump obtained by seeding with a different tension slope.

The direct Hamiltonian residual reported by the producer is an algebraic substitution check, not an independent discretized spatial-derivative test. Independent crosschecks are the mass-defect integral, tolerance changes, and the second equivalent integration formulation below.

## 2. A one-sign bump still fails at second time order

Define

\[
D=1-\rho_y^2+\rho^2(\phi_y^2/12-U/6).
\]

The initial profile equations give the exact identities

\[
D_y=-2\frac{\rho_y}{\rho}D+\frac{\epsilon\rho^2\phi_yb}{6},
\qquad
\boxed{(\rho^2D)_y=\frac{\epsilon\rho^4\phi_yb}{6}.}
\]

The **unmodified** Einstein–scalar evolution evaluated on this initial slice gives

\[
A_{TT}=-2D,\quad B_{TT}=4D,\quad \phi_{TT}=\epsilon\rho^2b.
\]

Because b vanishes in a neighborhood of the brane, differentiating the boundary conditions twice in time gives residuals

\[
\mathcal C_A^{(2)}=0,\qquad
\boxed{\mathcal C_B^{(2)}=-12\rho_yD_b},\qquad
\boxed{\mathcal C_\phi^{(2)}=-4\rho_b\phi_yD_b}.
\]

Therefore initial constraints and the first time corner conditions are not sufficient for smooth evolution. For ε≠0, a nonnegative bump and monotone φ_y>0 produce a nonzero integral defect of the sign of ε; that specific one-sign construction cannot remove the second-order obstruction. This is a restriction of this initial-data ansatz, not a general obstruction to perturbing the model.

Eight exact symbolic checks in `verify_constraint_seed_algebra.py` verify the Hamiltonian reduction, defect identity, three accelerations and three second-order corner formulas. All pass.

The one-bump example uses b=exp[1−1/(1−x²)] for |x|<1 and zero otherwise, x=(φ+.5)/.3. Seven numerical roots were found at ε=0, ±10⁻⁸, ±10⁻⁶, ±10⁻⁴, each refined in integration tolerance and cone start. At ε=+10⁻⁶, D_b≈1.5195815×10⁻⁴ and the B/φ corner residuals are approximately −.04792/−.04790. This is a resolved incompatibility, not numerical round-off. These one-bump files must not be presented as fully compatible evolution seeds.

## 3. Signed compensation repairs the continuum obstruction

Use two separated C∞ bumps:

\[
b(\phi)=b_1(\phi)-k b_2(\phi),
\]

where b₁ is supported on −.8<φ<−.5 and b₂ on −.45<φ<−.15, each of the same smooth exponential form and width .15 around centers −.65 and −.30. Solve for the three quantities (φ_h,y_b,k) with the two original junction conditions and

\[
\boxed{\int_0^{y_b}\rho^4\phi_y[b_1(\phi)-k b_2(\phi)]\,dy=0.}
\]

Regularity makes the cone contribution ρ²D vanish. The integral condition therefore enforces D_b=0. On the open neighborhood between the second bump and the brane, b=0, so the homogeneous defect equation also gives D=0 throughout that neighborhood. There the data obey the original static Einstein–scalar equations. Consequently A=logρ+T, B=logρ, φ=φ(y) is a local exact static-brane evolution: all formal local time-corner conditions are satisfied in the exact construction.

This is a continuum compatibility argument conditional on an exact root. It is **not** an existence theorem for that root, a proof of global smooth evolution, or a proof that finite differences preserve constraints. The finite floating-point approximants retain small nonzero residuals.

`balanced_constraint_seed.py` finds compensated roots for ε=±10⁻⁶ and also ±10⁻⁴, plus a zero reference. For the zero reference the unused compensation coefficient is fixed conventionally at .40932878; it is not uniquely measured when ε=0. For ε=+10⁻⁶, k≈.4093288451; for ε=−10⁻⁶, k≈.4093287185. Junction residuals are about 10⁻¹⁵ to 10⁻¹⁴. All roots remain well outside the bump support at the brane.

The calculation uses two formulations:

1. Integrate (ρ,ρ_y,φ,φ_y) with the derived second-order Hamiltonian equation and an independent defect integral. Refine relative tolerance from 3×10⁻¹³ to 8×10⁻¹⁴.
2. Independently evolve the integral I=ρ²D and use the positive branch
   \(\rho_y=\sqrt{1+\rho^2(\phi_y^2/12-U/6)-I/\rho^2}\).
   This avoids subtractive loss when evaluating a near-zero D from large terms.

The second formulation reaches defect-integral values around 10⁻¹⁴ or smaller at the brane. Direct subtraction of the first formulation leaves a residual of order 10⁻¹¹, and the different formulations agree on φ_b to approximately 10⁻¹⁰. This is useful numerical evidence, not outward-rounded error certification. All exact machine values and tolerance comparisons are preserved in `BALANCED_CONSTRAINT_SEED_RESULTS.json`.

As a separate check, `check_sampled_constraints.py` differentiates the saved profile values with fourth-order finite differences instead of substituting the ODE right-hand side. Across all five compensated/reference profiles, away from the excluded cone and endpoint stencils, the normalized Hamiltonian residual is below 2.11×10⁻¹⁰. The two independent first-derivative comparisons also pass. This is a sampled control with 15 explicit tolerances, not a global residual enclosure.

## 4. Reproduction and use

The scripts are portable Python with NumPy/SciPy for the ODE roots and SymPy for the eight algebra identities. Run them in their folder or pass the one-bump script an output directory. No original user files are loaded, no installation or network access is required, and the new results reside only in this checkpoint.

`CONSTRAINT_SEED_PROFILES.npz` contains the one-bump counterexamples; `BALANCED_CONSTRAINT_SEED_PROFILES.npz` contains the best mass-transport compensated profiles. Each profile is sampled at 4,097 proper-distance points. Array rows are y, ρ, ρ_y, ψ=φ+1, φ_y, z with brane z=0, and I=ρ²D. The samples are an exchange format, not a high-order interpolator. For subsequent evolution, resample from the ODE dense solution or construct and validate an interpolator on the target grid, keep boundary derivatives consistent, and test both discrete constraints and the propagated boundary residuals before interpreting growth or nonlinear fates.

The next physical computation is an evolution of these compensated data with demonstrated numerical constraint control. This checkpoint has not performed that evolution.
