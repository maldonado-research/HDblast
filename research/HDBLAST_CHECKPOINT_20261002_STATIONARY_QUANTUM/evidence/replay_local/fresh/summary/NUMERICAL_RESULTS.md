# Registered stationary shell results

Numerical validator status: **PASS**, 1584 gates, 0 failed. The producer completed 48 registered integration/root cases and 3 wrong-model controls, with 0 solver failures.

These are stationary regular-cone roots of the declared one-field mathematical model. The sampled hierarchy screen is reported separately; it does not certify an EFT cutoff or physical viability.

| b | gamma | Delta phi | Delta H²/H0² | phi/H² shift resolved | phi/H² nonlinear resolved | hierarchy screen |
|---:|---:|---:|---:|:---:|:---:|:---:|
| -1 | 0.01 | 4.167864242e-14 | 2.553434782e-11 | yes/yes | no/no | pass |
| -1 | 1 | 4.16787114e-12 | 2.553372899e-09 | yes/yes | no/no | fail |
| -1 | 100 | 4.167873931e-10 | 2.553370533e-07 | yes/yes | no/no | fail |
| -1 | 10000 | 4.168112976e-08 | 2.553564055e-05 | yes/yes | no/no | fail |
| -1 | 1e+06 | 4.192302211e-06 | 0.002573097526 | yes/yes | yes/yes | fail |
| 0 | 0.01 | 3.930232875e-19 | 2.553400466e-11 | no/yes | no/no | pass |
| 0 | 1 | -2.90430657e-17 | 2.553373699e-09 | yes/yes | no/no | fail |
| 0 | 100 | -2.937388288e-15 | 2.553370682e-07 | yes/yes | no/no | fail |
| 0 | 10000 | -2.996311288e-13 | 2.553563672e-05 | yes/yes | no/no | fail |
| 0 | 1e+06 | -3.01975675e-11 | 0.002573047459 | yes/yes | yes/yes | fail |
| 1 | 0.01 | -4.16788186e-14 | 2.553400466e-11 | yes/yes | no/no | pass |
| 1 | 1 | -4.167930364e-12 | 2.553369009e-09 | yes/yes | no/no | fail |
| 1 | 100 | -4.167932687e-10 | 2.553370779e-07 | yes/yes | no/no | fail |
| 1 | 10000 | -4.168173107e-08 | 2.553564272e-05 | yes/yes | no/no | fail |
| 1 | 1e+06 | -4.192383249e-06 | 0.002573097748 | yes/yes | yes/yes | fail |

Shifts subtract each integration method’s own solved gamma-zero control. The registered nonlinear detection requires departure from the archived linear prediction of at least 0.1% of the measured shift and more than ten empirical refinement spreads plus the fixed floating-point floor. Detection classification is separate from root acceptance.

## Largest registered amplitude

- b=-1: Delta phi=4.19230221095e-06; Delta H²/H0²=0.00257309752638; nonlinear fractions phi=0.0058276551, H²=0.0076673306; Ehat/M5=11.137806.
- b=+0: Delta phi=-3.01975674961e-11; Delta H²/H0²=0.00257304745948; nonlinear fractions phi=0.0076430919, H²=0.0076479792; Ehat/M5=11.137806.
- b=+1: Delta phi=-4.192383249e-06; Delta H²/H0²=0.00257309774833; nonlinear fractions phi=0.0058325765, H²=0.0076673315; Ehat/M5=11.137806.

## Controls, failures and provenance

Maximum producer residuals: |C|=1.0495e-18, |B|=5.4503e-15. Rejected Newton trial steps: 1.

- flip_current: correct-equation (C,B) at the wrong root = [-2.779326858178173e-22, 2.9804527594655147e-05].
- omit_metric_source: correct-equation (C,B) at the wrong root = [-1.5126304639714704e-07, -3.7608262858090935e-19].
- freeze_sources: correct-equation (C,B) at the wrong root = [-1.1570167388155482e-09, 8.621284182929228e-08].

The original frozen validator exited 1 with 48 identical diagnostic-keyword TypeErrors. Its code, CHECKS.json, log and exit are retained. A separately published one-keyword metadata amendment was run on unchanged producer results; no gate, model or grid was changed.

Public freeze: `b4f77f5826e0a603400513832aa53b0673df7533`. Producer results SHA-256: `427d3270a5b1faaf9f3fc87624c54cfc4b3ed13690268d7a21a8ee2ada6568f5`. Amended checks SHA-256: `9ad0bb2c859db98124b6c45cdf8f31e628d386f7b2a9ad21f8d0d46f8d0c3dc6`.

No uniqueness, dynamical or quantum stability, heating, particle production, history of exit, observational fit or physical parameter inference is established.
