# Pre-execution numerical specification

Prepared after the first protocol commit (32f1fc617b0691bbc10783e7e12b337afb4a3af9), before any mode integration in this experiment. The original REGISTRATION.md is retained unchanged. These specifications address analytical/numerical design, not observed outcomes.

- Add infrared panel edges 0.125,0.25,0.5,1 to every original momentum grid. This resolves the integrable low-k behavior when x=0. The primary 24 and comparison 48 nodes per panel remain unchanged; no zero-momentum endpoint is sampled.
- The integrated work normalizer is |rho(t)|+|rho(initial)|+|work(t)|+m_infinity^4, with m_infinity=2. Also report absolute residual/m_infinity^4. The stationary asymptotic vacuum makes a vanishing relative denominator unsuitable.
- Use an analytic phase primitive and a stable constant-reference vacuum subtraction. This is algebraically equivalent to the registered instantaneous subtraction plus analytic continuum static restoration. Verify the static functions independently at high precision.
- Pressure trace uses both five- and nine-point finite differences on the 2001-point common grid, plus a coarsened nine-point comparison. Include the exact finite-K vacuum-tail correction D_K; the 1e-4*m_infinity^4 threshold is unchanged.
- Add a negative control computed from the same modes: omit only the curvature counterterm and report the resulting pressure drift with Lambda. Its leading predicted difference is -xddot*ln(Lambda2/Lambda1)/(48pi²). Energy and scalar work cannot detect this pressure error at H=0.
- Add an independent representation of the regulator limit using the physical modes alone. The matched pressure integrand adds xddot*(2k²/3+r)/[16*(k²+r)^(5/2)], r=4, to physical p minus its instantaneous zero point, and subtracts V_R. The potential is V_R=[x²ln(x/r)-1.5x²+2rx-.5r²]/(64pi²). The positive reference avoids an infrared subtraction singularity. Compare with PV results and retain finite-cutoff sensitivity; no extra physical model is being fitted.
- The constant-mass null is an algebraic prescription control. It is not reported as a nontrivial dynamical integration.
- All numerical changes following execution will be identified, including failed attempts. Numerical convergence evidence is distinct from a mathematical proof of regulator removal.

No result is known at this addendum's preparation.
