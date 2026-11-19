# Two-time state-independent localization on the stationary hilltop

Input: `prior/hilltop/HILLTOP.json`. All times are x = Htop t. n=2 means eigenvalue 8 on the unit S³.

Equal interval halfwidths are epsilon; probabilities concern one initial density operator, evolved without intermediate measurement.

| branch | detuning | kappa | homogeneous x_match (reference) | absolute kappa B |
|---|---:|---:|---:|---:|
| phi_m=0.5 | 0.001 | 3.93515010e-06 | 3.38174878 | 3.24996185e+12 |
| phi_m=0.5 | 1e-10 | 3.93515010e-13 | 8.75476389 | 6.12938749e+35 |
| phi_m=0.5 | 1e-20 | 3.93515010e-23 | 16.43004758 | 1.08400001e+69 |
| phi_m=-0.5 | 0.001 | 3.93515010e-06 | 3.02840670 | 1.90513220e+07 |
| phi_m=-0.5 | 1e-10 | 3.93515010e-13 | 8.40136219 | 5.07497155e+24 |
| phi_m=-0.5 | 1e-20 | 3.93515010e-23 | 16.07664586 | 3.94755514e+49 |

## Last time admitting both one-mode marginal probabilities >= 0.95

This uses the sharp interval/prolate bound, evaluated numerically. The exact operator identity and rigorous simpler bounds are in THEORY.md.

| branch | detuning | epsilon=0.01 | epsilon=0.001 | epsilon=0.0001 |
|---|---:|---:|---:|---:|
| phi_m=0.5 | 0.001 | 0.43755574 | 0.107226 | 0.0015218653 |
| phi_m=0.5 | 1e-10 | 1.65644044 | 1.30391304 | 0.953748091 |
| phi_m=0.5 | 1e-20 | 3.42885966 | 3.07398487 | 2.71918359 |
| phi_m=-0.5 | 0.001 | 0.521654762 | 0.115723731 | 0.00152190424 |
| phi_m=-0.5 | 1e-10 | 2.03891866 | 1.60056467 | 1.16445302 |
| phi_m=-0.5 | 1e-20 | 4.23758915 | 3.79764977 | 3.35773823 |

## Last time admitting both one-mode marginal probabilities >= 0.99

This uses the sharp interval/prolate bound, evaluated numerically. The exact operator identity and rigorous simpler bounds are in THEORY.md.

| branch | detuning | epsilon=0.01 | epsilon=0.001 | epsilon=0.0001 |
|---|---:|---:|---:|---:|
| phi_m=0.5 | 0.001 | 0.402090869 | 0.0776740533 | 0.000939045911 |
| phi_m=0.5 | 1e-10 | 1.6194049 | 1.26705907 | 0.917240564 |
| phi_m=0.5 | 1e-20 | 3.39164716 | 3.03677787 | 2.68198776 |
| phi_m=-0.5 | 0.001 | 0.477647108 | 0.0816282731 | 0.000939055059 |
| phi_m=-0.5 | 1e-10 | 1.99289687 | 1.55469757 | 1.11895227 |
| phi_m=-0.5 | 1e-20 | 4.19145864 | 3.75152117 | 3.31161423 |

## At the homogeneous departure reference time

The rigorous HS upper bound on min(p0,pT) equals 1/2 plus the values below. The exact sharp prolate ceiling is never larger. This is not a bound on one specified marginal alone.

| branch | detuning | epsilon=0.01 | epsilon=0.001 | epsilon=0.0001 |
|---|---:|---:|---:|---:|
| phi_m=0.5 | 0.001 | 2.212947e-09 | 2.212947e-10 | 2.212947e-11 |
| phi_m=0.5 | 1e-10 | 5.095673e-21 | 5.095673e-22 | 5.095673e-23 |
| phi_m=0.5 | 1e-20 | 1.211701e-37 | 1.211701e-38 | 1.211701e-39 |
| phi_m=-0.5 | 0.001 | 9.140028e-07 | 9.140028e-08 | 9.140028e-09 |
| phi_m=-0.5 | 1e-10 | 1.770897e-15 | 1.770897e-16 | 1.770897e-17 |
| phi_m=-0.5 | 1e-20 | 6.349594e-28 | 6.349594e-29 | 6.349594e-30 |

If the initial unmeasured marginal probability is already at least 0.99, the stronger Landau-Pollak angle relation bounds the final unmeasured marginal by 0.01 plus at most 3.637718e-07 across these reference cases. The JSON preserves the small corrections separately.

## Validation

Direct DOP853 transfer vs independent hypergeometric solution: maximum relative discrepancy 3.406e-11.
100-digit analytic Wronskian check: maximum relative discrepancy 2.011e-66.
Prolate 96- vs 160-node critical-c difference: maximum 7.905e-14.

The prolate quadrature convergence is numerical evidence, not a rigorous interval eigenvalue certificate. The Hilbert-Schmidt probability inequality itself is rigorous.
