# E1: matched quantum source through a smooth mass-zero crossing

Registered 2 October 2026 before numerical mode integration for this experiment. Earlier HDBLAST finite-band results and exact pulse formulas are known. Analytical design and counterterm derivations precede this registration.

## Purpose and scope

Test a consistently matched absolute scalar stress/current through a mass-zero crossing, including the pressure counterterm invisible to an energy test on flat spacetime. This is a Minkowski-space prerequisite benchmark, not a new archived-shell evolution, a radiation era or a solution of the original delta=0.001 model. The external field supplies the quantum-sector work.

Use units T=1, m_infinity=2; x(t)=M²(t)=4 tanh²(t), with primary interval [-10,10]. The physical scalar is real and minimally coupled. Define S=1/2<chi²> so rho_dot=xdot*S. If phi=tanh(t), J_phi=8phi*S; this is a declared test profile, not the archived phi solution.

## Regulator and matching

Covariant Pauli–Villars masses y_j=x+j Lambda² and coefficients (1,-3,3,-1), j=0,1,2,3. All species use the same minimal coupling and derivative partial_x y_j=1. These auxiliary weights are not physical species; Lambda is not the HDBLAST five-dimensional gravity cutoff.

For one common effective action with local terms -V+F R, use
V_PV=Sum c_j y_j² ln(y_j/mu²)/(64pi²),
F_PV=Sum c_j y_j[ln(y_j/mu²)-1]/(192pi²).
The common arbitrary scale cancels in these PV sums.

At positive r=x_ref=4 impose V_R(r)=V_R'(r)=V_R''(r)=0 by subtracting the quadratic Taylor polynomial C_V. Impose F_R(r)=F_R'(r)=0 by subtracting the linear Taylor polynomial C_F. These are explicit finite matching choices, not observational determinations.

Compute stable dynamic integrals after subtracting each instantaneous zero-point integrand, then restore the analytically integrated matched potential V_m=V_PV-C_V:
rho=Sum c_j rho_dyn,j+V_m;
S=1/2 Sum c_j Q_dyn,j+V_m';
p=Sum c_j p_dyn,j-V_m-2F'_PV(r) xddot.
The last term is required even at R=0 because varying F(x)R contributes to pressure. Curvature-squared counterterms have zero first variation on this exact flat background; their finite matching remains necessary for a future curved extension.

Do not use a time-local heavy-mass expansion at x=0, a field-dependent subtraction scale, or an uncorrected hard three-momentum vacuum cutoff. Any high-k tail expansion is confined to positive high momentum and must be documented.

## State and numerical formulation

Use the exact asymptotic in-vacuum of the tanh-squared frequency, initialized at the finite start by a convergent small-z Jost hypergeometric series, z=(1+tanh(t))/2. Check series stability and normalization. The intended state is the asymptotic scattering state, not a freshly chosen finite-time instantaneous vacuum.

The hypergeometric parameters obey a+b=1, ab=4, c=1-i omega_infinity,j. The exact asymptotic occupation is
n_exact,j=cosh²(pi sqrt(4-1/4))/sinh²(pi sqrt(k²+4+j Lambda²)).
Use stable logarithmic evaluation where needed. Finite output time, phase and truncation errors are retained.

Evolve instantaneous Bogoliubov variables and phase with DOP853, primary rtol=1e-10, atol=1e-13, and a phase-resolving maximum step. Compute dynamic energy, pressure, Q and the scalar work using common fixed k weights. Report the actual integration method and work normalization.

## Prespecified matrix and tests

Primary: Lambda=16, K=128, 24 Gauss–Legendre nodes per dyadic panel [0,2,4,8,16,32,64,128]. Comparison cases:
1. Lambda=4,K=32, same panel rule.
2. Lambda=8,K=64, same panel rule.
3. Primary 48 nodes per panel.
4. Primary K=64, retaining the identical lower panels.
5. Primary K=192, adding [128,192].
6. Primary tighter rtol=1e-12, atol=1e-15.
7. Primary exact-in preparation evaluated from start -8 to end +8.
8. Constant x=4 matched-vacuum null control, whose numerical/algebraic implementation must be stated.

Compare source curves on the common [-6,6] interval with 2001 samples or a documented refinement. Retain full-interval endpoints for work and spectrum tests. All cases and failures are reported.

Control thresholds:
- maximum normalization drift <=1e-7;
- maximum absolute late-time occupation discrepancy from the exact formula <=1e-6, with finite-time basis corrections reported;
- normalized integrated rho/work residual <=1e-5, using absolute energy and absolute accumulated work terms plus a stated scale floor;
- quadrature, momentum-cutoff and integration-tolerance curve differences target <=1e-5 after scaling rho,p by m_infinity^4 and S by m_infinity²;
- regulator-sequence source differences target <=1e-3 on those same scales, reported separately from numerical accuracy. Failure does not become a successful regulator-removal claim.

Pressure receives an independent trace check:
Q_R=2S;
-rho+3p=Q_R''/2-x Q_R+A_Lambda,
A_Lambda=-4V_m+2x V_m'-Sum_{j>=1} c_j j Lambda² Q_dyn,j.
Compare resolved finite-difference Q_R'' with the direct pressure, account for the finite-K vacuum tail, and report stencil/refinement sensitivity. Target residual <=1e-4*m_infinity^4. The continuum-regulator limit predicts
A_infinity=(x-r)²/(32pi²)+xddot/(96pi²)
in this matching scheme; test convergence instead of assuming it.

A fixed K/Lambda sequence alone cannot establish regulator removal. Use the fixed-Lambda K controls to bound the pressure tail before interpreting Lambda sensitivity. Late-time pressure can retain coherent oscillations; occupation is not a thermalization test.

## Interpretation

The specified finite conditions remove scheme ambiguity for this benchmark only. They are not a unique HDBLAST vacuum prediction. A successful experiment would validate a matched source calculation and expose the limitations of energy-only checks. It would not complete FRW renormalization, shell matching, decay, thermalization, energy affordability or five-dimensional backreaction. Changes after results are seen will be labeled exploratory.
