# Independent review of the scalar-profile de Sitter branch

22 September 2026. This review did not edit or import `solve_plus_branch.py`. It checked its equations and performed separate four-variable integrations using an independently expanded potential and the second-order Einstein equation. It did not rerun nonlinear time evolution.

**Verdict:** the calculated profiles are consistent floating-point candidates for regular-cone static de Sitter configurations satisfying both registered junction conditions. The analytical small-detuning expansion is correct. No code error was found. The calculation does not prove stability, dynamical attraction, uniqueness, or interval-certified existence.

## Independent equations and regular-cone start

The producer evolves `(ρ,η,η_y)` with η=φ−1 and the positive square root of the first integral. The independent checker evolves `(ρ,ρ_y,η,η_y)` using

\[
\rho_{yy}=-\rho\left(\frac{\eta_y^2}{4}+\frac{U}{6}\right),
\qquad
\eta_{yy}=U_\eta-4\frac{\rho_y}{\rho}\eta_y.
\]

It does not enforce the first integral during evolution. Its independent constraint test is

\[
\mathcal C=\rho_y^2-1-\rho^2\left(\eta_y^2/12-U/6\right).
\]

The independently expanded potential used by the checker is

\[
U=-\frac2{27}+\frac{14}{9}\eta^2+\frac{50}{27}\eta^3
-\frac16\eta^4-\frac49\eta^5-\frac2{27}\eta^6.
\]

Both this expression and the producer's factored potential follow from the registered W. The producer's first and second derivatives are correct.

The cone series coefficients also check. With u=U(η_h), v=U′(η_h), w=U″(η_h),

\[
\rho=y-\frac{u}{36}y^3+
\left(\frac{u^2}{4320}-\frac{v^2}{750}\right)y^5+O(y^7),
\]

\[
\eta=\eta_h+\frac{v}{10}y^2+
v\left(\frac{w}{280}+\frac{u}{630}\right)y^4+O(y^6).
\]

The numerical start is necessarily at finite y. The independent integrations varied that start over **2×10⁻⁴, 5×10⁻⁵, 10⁻⁵**, using the same archived shooting parameters rather than refitting each case. Their consistency tests sensitivity to the truncated cone start; it is not a rigorous bound on the omitted series tail.

The Gegenbauer initialization is also consistent. At η=0, k₊=1/9 and U″=28/9. Setting x=cosh(y/9), the regular linear solution is proportional to C₁₄²(x), because 14(14+4)=252=U″/k₊². Its logarithmic derivative in the producer is correct. This supplies only an initial guess; the reported profiles solve the nonlinear ODE.

## Numerical checks

`independent_branch_checks.py` made **18 independent integrations**: three cone starts for each of δ=0.0003, 0.001, 0.003, 0.01, 0.03, 0.1. It compared 4,097 archived profile samples per detuning and sampled its independently evolved constraint at 8,001 positions per integration. All **28 assertions pass**.

Maximum values across those integrations:

| Diagnostic | Maximum magnitude |
|---|---:|
| Metric junction residual ρ_y/ρ−σ/6 | 3.61×10⁻¹⁶ |
| Scalar junction residual φ_y+σ′/2 | 1.57×10⁻¹⁴ |
| Independent first-integral residual, absolute | 1.60×10⁻¹² |
| Independent first-integral residual, normalized | 2.02×10⁻¹⁵ |
| Relative H² difference from producer | 3.35×10⁻¹⁵ |
| Maximum relative ρ profile difference | 2.19×10⁻¹⁵ |
| Maximum η profile difference / |η_b| | 4.04×10⁻¹³ |

The constraint normalization is explicitly `1+ρ_y²+ρ²(η_y²/12+|U|/6)`. Its small ratio should not be reported as an interval error bound. The independently checked H² identity is

\[
\frac1{\rho_b^2}=\frac{\sigma_b^2}{36}
-\frac{(\sigma'_b)^2}{48}+\frac{U_b}{6}.
\]

This identity follows by substituting **both** junctions into the radial first integral; retaining the scalar-gradient term is essential. It is a consistency test, not an independent theorem of branch existence.

At the registered detuning δ=0.001 the refined candidate has

\[
\phi_b=0.9999159473169134,\qquad
\phi_{y,b}=-0.0001306991661942445,
\]

\[
\rho_b=129.9247628496781,\quad
H^2=5.924014794329276\times10^{-5},\quad
\eta_h=-1.336692365158688\times10^{-23}.
\]

Using the original unstable-shell normalization ρ₀=78.82817714224423 from the external archive gives **H/H₀≈0.606721732**. The old constant-φ metric benchmark was 0.606726506. Their fractional difference is only about **7.87×10⁻⁶**. Therefore the scalar-junction correction is logically important but does not explain the roughly 5% gap between the archived finite-time roll-off rate near 0.639 and a possible asymptotic rate near 0.607.

## Analytical expansion

The positive regular AdS exponent solves λ²+4k₊λ−U″=0 and equals 14/9. The linearized scalar junction implies

\[
\eta_b=-\frac9{64}\delta c+O(\delta^2).
\]

Substituting η_b=aδ+O(δ²), a=−9c/64, into the full H² identity gives

\[
H^2=\frac{\delta(1+c)}{27}
+\delta^2\left[\frac{(1+c)^2}{36}-\frac{c^2}{384}\right]
+O(\delta^3).
\]

The coefficient of c² was independently checked using exact rational arithmetic, not a numerical fit. At second order the a² terms cancel; the remaining scalar correction is `ca/27−ca/6−c²/48 = −c²/384`.

The data support the stated remainder orders over the small-detuning end of the scan:

| δ | (η_b+9δc/64)/δ² | (H²−H²_second_order)/δ³ |
|---:|---:|---:|
| 0.0003 | -0.01589347 | -0.002324315 |
| 0.001 | -0.01589535 | -0.002324521 |
| 0.003 | -0.01590073 | -0.002325122 |

This is numerical consistency with the expansion, not a uniform remainder bound in δ.

## Interpretation and remaining work

This branch addresses the exact scalar-junction inconsistency in a constant φ=1 endpoint. Its regular cone has φ_h very close to **+1**, whereas the original unstable shell's cone is near **−1**. Both satisfy the same local field equations and junction law, but they are distinct configurations. Connecting them dynamically may involve the original domain wall leaving the observable patch; that connection and the global extension have not been demonstrated here.

The branch need not be an attractor merely because it exists. A perturbation-spectrum analysis, a controlled evolution from compatible original-shell data, and a careful comparison across a nonsingular chart are needed to assess that question. Root convergence and agreement between floating-point IVPs also do not give an interval existence certificate. The proper status is **a reproducible static scalar-profile candidate**, with both junctions satisfied to the displayed numerical accuracy.

## Reproduction

Run `python independent_branch_checks.py` with NumPy and SciPy available. Inputs are the producer's script, `PLUS_BRANCH_RESULTS.json`, and `PLUS_BRANCH_PROFILES.npz`; output is `INDEPENDENT_BRANCH_CHECKS.json`. The script verifies the producer hash, preserves all producer files, and recomputes its own second-order IVPs. It records input hashes and explicit tolerances, so the test's scope is reviewable.
