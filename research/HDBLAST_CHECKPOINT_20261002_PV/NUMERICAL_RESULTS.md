# Executed matched-source benchmark

One real scalar evolves through x(t)=M²(t)=4tanh²(t) on flat spacetime, with an exact asymptotic incoming state. Both potential and curvature matching are fixed at r=4. All numerical values below use the declared T=1 units; they are not HDBLAST H0-based cosmological predictions.

## Eight cases and independent checks

All eight dynamic cases and the algebraic constant-mass null pass. The largest exact scattering-occupation error is 5.74e-10. The largest Wronskian error is 1.37e-9. The paired energy/work residual is at most 5.59e-14 with the stated non-cancelling normalization and natural floor. The independent nine-point pressure-trace check has residual at most 2.34e-11*m_infinity^4.

Seventy-digit analytic checks verify the matching functions, static vacuum tails, Jost initialization and phase quadrature. Full metadata and tolerances are in [the audit](AUDIT_AND_EXECUTION.md), [registration](REGISTRATION.md) and [pre-execution addendum](PROTOCOL_ADDENDUM.md).

## Pressure changes qualitatively when the missing matching term is included

At t=0, x=0 and xddot=8:

| Auxiliary regulator Lambda | Matched rho | Matched p | Matched S | Potential-only pressure |
|---|---:|---:|---:|---:|
| 4 | -0.0065899634 | +0.0090611554 | 0.0107389505 | +0.0022929902 |
| 8 | -0.0087665234 | +0.0068873167 | 0.0128505404 | -0.0091715645 |
| 16 | -0.0094394883 | +0.0065477090 | 0.0134528976 | -0.0205254361 |

The potential-only column deliberately omits -2F'_PV(r)xddot. The scalar work and energy are unchanged by that omission, so the flat-space energy ledger cannot detect the pressure error. The wrong pressure continues to drift with Lambda. The corrected sequence approaches the independent direct-limit calculation.

The negative matched energy near the crossing is a quantum expectation in a stated finite matching convention with the asymptotic static vacuum set to zero. It is not an assertion of a negative-energy particle gas, a cosmological instability or free energy. The external mass profile supplies or extracts the work.

## Numerical refinement and remaining regulator effects

Differences below are maximum absolute curve changes on [-6,6], divided by m_infinity^4=16 for rho,p and m_infinity²=4 for S:

| Comparison with primary Lambda=16,K=128 | rho | p | S |
|---|---:|---:|---:|
| 48 instead of 24 nodes/panel | 4.26e-14 | 1.81e-12 | 4.62e-14 |
| K=64 | 9.32e-11 | 3.88e-8 | 6.30e-10 |
| K=192 | 5.00e-13 | 2.00e-10 | 3.14e-12 |
| Tighter ODE tolerance | 1.11e-16 | 2.41e-13 | 7.88e-17 |
| Exact Jost preparation at -8 instead of -10 | 2.22e-16 | 4.10e-13 | 1.21e-16 |
| Lambda=8, with its prescribed momentum band | 4.21e-5 | 3.78e-5 | 1.51e-4 |

The first five comparisons pass the numerical target. The regulator comparison passes the separate 1e-3 natural-scale target, but is not negligible relative to the signal: at the crossing, Lambda=8 to 16 changes rho by 7.13%, p by 5.19%, and S by 4.48% relative to the Lambda=16 values. The report's generic numeric_target_pass=false for regulator comparisons is not a failed numerical-refinement case; those rows have a different registered target, shown alongside it.

The direct physical-mode prescription uses convergent positive-reference integrands in the same formal matching convention. In the primary band its crossing values are rho=-0.0096746050, p=+0.0064457340 and S=0.0136595508. These are finite-K approximations to that representation, not a proof of its continuum limit. The direct cutoff/refinement comparisons and matched-trace residuals are recorded in outputs/analysis_summary.json.

| Direct representation at the crossing | rho_K | p_K | S_K |
|---|---:|---:|---:|
| K=64 | -0.00966533956 | +0.00643957650 | 0.01365144080 |
| K=128 | -0.00967460501 | +0.00644573396 | 0.01365955079 |
| K=192 | -0.00967632225 | +0.00644687775 | 0.01366105352 |

For K=192 versus 128, maximum curve differences divided by the natural scales are 1.07e-7, 1.09e-7 and 3.76e-7 for rho,p,S. Direct quadrature/tolerance/prehistory differences are far smaller. The leading cutoff behavior is 1/K².

Post hoc derivation gives an exact finite-K trace factor v_r³, v_r=K/sqrt(K²+r), multiplying the local continuum offset. A separately executed nine-point check of saved primary columns gives maximum absolute residual 1.673e-9 (1.046e-10*m_infinity^4). Omitting this factor leaves 2.164e-5 in absolute units; that is the known finite-cutoff contribution. See [the appendix](DIRECT_CUTOFF_APPENDIX.md).

Leading UV-tail addbacks yield crossing estimates at K=64,128,192 that agree within about 3e-8 in all three sources. These are exploratory asymptotic estimates, not rigorous continuum error bars. The mathematical identity and addbacks were developed after the registered outcomes, and no original pass criterion was changed.

## Precision qualifications

The extraordinarily small internal refinement differences apply to the unit-Wronskian estimator. Conservative differences from raw bilinear reconstruction reach 1.69e-7 in primary energy and 1.61e-6 at K=192—approximately 8.5 and 81 parts per million when compared with energy 0.02. These are diagnostics, not full numerical error bars. Remaining finite-Lambda corrections are substantially larger.

Pressure trace is checked independently by finite-differencing Q_R=2S and including the analytically derived finite-K tail. The formal regulator-limit trace offset is scheme dependent. Its primary mismatch from the finite-Lambda offset is 5.65e-5*m_infinity^4 and decreases across the tested regulator sequence.

## Late-time energy and interpretation

At t=10 the primary energy is 0.01996432657 and pressure is -0.00207606169. The instantaneous pressure can retain coherent oscillations; particle counting does not by itself supply it. The exact asymptotic spectrum gives number density 0.008831115308, energy 0.019964326635 and gas pressure 0.001405514039. Its gas p/rho is 0.0704013 and mean particle energy is 2.2606801. The primary endpoint energy differs from this spectrum's energy by -3.36e-9 fractionally. Full pressure includes coherent terms and is not the gas-pressure column.

S is the source conjugate to x=M². The corresponding scalar force for x=4(phi-phi_*)² is 8(phi-phi_*)S, zero exactly at the crossing.

This is a matched absolute one-loop source on a specified external flat background. It is not a new five-dimensional evolution, decay calculation, thermal radiation bath, energy-affordability result or hot Big Bang. The main advance is a reproducible source and pressure benchmark that the future curved-shell calculation must pass.
