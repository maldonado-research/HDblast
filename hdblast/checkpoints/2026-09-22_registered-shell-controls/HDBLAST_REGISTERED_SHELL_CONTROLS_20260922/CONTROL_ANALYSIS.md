# Completed registered-shell controls

Analysis of six spectra and six short nonlinear-PDE evolutions. All numbers are floating-point results. The expected rate was supplied before the new eigensolves; comparison is a code-validation test, not a blind physical prediction.

| Detuning | Shell spacing | L | Nodes | Selected growth | Difference from prior rate |
|---:|---:|---:|---:|---:|---:|
| 0.001 | 0.0004 | 6.0 | 601 | 1.658701798553 | 1.508e-03 |
| 0.001 | 0.0002 | 6.0 | 1200 | 1.657211039351 | 1.741e-05 |
| 0.001 | 0.0001 | 6.0 | 2399 | 1.657194372000 | 7.408e-07 |
| 0.001 | 5e-05 | 6.0 | 4797 | 1.657193663767 | 3.252e-08 |
| 0.001 | 0.0001 | 8.0 | 2542 | 1.657194372340 | 7.411e-07 |
| 0.1 | 0.0004 | 6.0 | 601 | 1.627021369743 | 2.729e-08 |

The finest selected eigenpair has absolute residual 6.723e-08; the observed rate difference is not a rigorous error bound. Other near-shift eigenvalues are not independently validated.

| Evolution | Spacing | Fitted growth | Fit time window | Maximum normalized H | Maximum normalized M |
|---|---:|---:|---|---:|---:|
| reg_mode_plus_h2 | 0.0002 | 1.657210815291 | [1, 2.5] | 3.401e-09 | 1.701e-09 |
| reg_mode_minus_h2 | 0.0002 | 1.657211264420 | [1, 2.5] | 3.401e-09 | 1.701e-09 |
| reg_mode_ko_h2 | 0.0002 | 1.657210814572 | [1, 2.5] | 3.401e-09 | 1.701e-09 |
| reg_mode_h1 | 0.0001 | 1.657194158956 | [1, 2.5] | 2.922e-10 | 1.422e-10 |
| reg_plus_h2 | 0.0002 | 1.657169774085 | [3, 4.5] | 2.771e-06 | 1.065e-06 |
| reg_plus_h1 | 0.0001 | 1.657155492040 | [3, 4.5] | 1.502e-05 | 6.791e-06 |

Mode seeds use shell scalar amplitude ±10⁻⁸. They test the linear limit of the nonlinear equations, ending at coordinate time 2.5. Their full nonlinear junction mismatch starts at second order. Shifted-tension seeds use dc=10⁻⁷ and end at 4.5; their constraints worsen under refinement, so they do not pass nonlinear initial-data acceptance.

The L6→L8 selected-rate change is 3.404e-10. The KO 0→0.02 fitted-rate change is -7.193e-10; the opposite-sign fitted-rate difference is 4.491e-07. These are bounded controls, not guarantees of late-time behavior.

H,M are divided by 1+6Hc(z)^2+phi_z(z)^2 on z>−0.8L excluding first 6 nodes. This is a background normalization, not a relative perturbation-error bound. Raw arrays are in NPZ snapshots.

Time and rates are dimensionless: one unit of the static coordinate time equals one initial-shell Hubble time. No observed cosmological time or energy scale has been calibrated.

All 27 analysis checks passed, including the check that preserves the negative shifted-seed finding. This script did not rerun the spectra or evolutions. Full fits, inputs, hashes and final rows are retained in CONTROL_ANALYSIS.json.
