# Pre-execution implementation choices

This addendum records remaining numerical choices before executing any new archive audit or control. Registration commit: 34a2be98b5e972ce6fa0263cda201bae56488720.

Reconstruction matrix: scipy CubicHermiteSpline using stored values and derivatives, CubicSpline using values, and make_interp_spline of degrees 5 and 7 using values and their default not-a-knot boundary conditions. The higher-degree reconstructions do not impose the archived derivatives. Measure that discrepancy rather than treating it as zero.

Evaluate the union of saved times within [0,6.9] and 2001 uniform points. Evaluate all actual interior spline breakpoints by exact one-sided polynomial coefficients, including the B-spline piecewise-polynomial conversion. Require relative continuity residual <=1e-9 wherever continuity is guaranteed by the stated C class; derivative jumps beyond that class are reported rather than failed. Conformal time is integrated per Hermite interval with scipy quad epsabs=epsrel=1e-11.

Report U, R, Rdot and Rdd from each reconstruction and local metric variations of a unit R² action. These are sensitivity diagnostics, not a fitted quantum stress or measured R² coefficient. Report absolute and range-scaled reconstruction differences without a post hoc physical closeness threshold.

The exact abrupt-step control uses U_minus=1,U_plus=4 and k in {8,16,32,64,128,256}; require the final relative deviation of k^4 n from DeltaU²/16 below 1e-3. Distinct-knot interference may be checked by integrating the leading asymptotic beta expression analytically with cosine-integral functions. That check is a leading-tail algebraic diagnostic, not a full high-k mode integration. If used, high bands are set from inverse minimum conformal-knot spacing, with logarithmic-slope relative target 1e-3.

The independently derived generalized adiabatic subtraction on a smooth FRW background may be checked symbolically or with analytic jets. It will be reported as a construction/control only; no archived absolute renormalized source is inferred from the readiness audit.
