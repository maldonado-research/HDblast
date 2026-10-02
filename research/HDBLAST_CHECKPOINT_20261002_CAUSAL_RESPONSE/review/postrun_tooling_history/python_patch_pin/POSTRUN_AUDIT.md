# Independent audit of the preserved forced-mode results

**Passed, normally and with Python optimization enabled.** This is a post-run
archive audit of the experiment frozen publicly in commit
`a01d17e015f070be2ad2018745e877148fe2f4e4`. It adds no new source or mode evolution
and changes no frozen implementation or scientific acceptance gate.

The audit verifies all 32 frozen input hashes, the validator's bindings to the
primary and mode result files, and the four archived mode-file hashes. It
reconstructs all 72 coarse/fine finite-cutoff responses directly from individual
archived momentum nodes, mode amplitudes, counterterms, and weights, with a
different summation grouping from the producer. All 24 stored observations
satisfy the amplitude and reconstructed canonical-mode Wronskian checks.
Every panel satisfies its Gauss-Legendre polynomial moments through degree 31:
98,304 panel/moment checks in total.

| Check | Largest observed discrepancy | Acceptance |
| --- | ---: | --- |
| Archived integral versus producer | 1.735e-18 | Post-run arithmetic allowance 2.843e-14 |
| Fine modes versus primary finite-K formula | 6.939e-17 | Frozen 2e-9 plus primary quadrature estimate |
| Coarse versus fine modes | 1.284e-16 | Frozen 2e-10 |
| Observation amplitude Wronskian / epsilon | 9.758e-14 | Frozen 1e-10 |
| Observation canonical Wronskian / epsilon | 9.692e-14 | Frozen 1e-10 |
| Gauss-Legendre panel moment | 5.554e-14 | Post-run arithmetic allowance 1e-10 |

The archived integrands reproduce the pointwise sum of the raw mode response
and the fixed-reference subtraction exactly. For the positive pulse at
eta=-4, over 128<=k<256, their separate maximum magnitudes are approximately
7.813e-7, while the combined integrand's maximum magnitude is 3.333e-11. This
exposes the ultraviolet cancellation directly in the preserved arrays.
The archived high-k coefficients also satisfy the separate analytic coefficient
bound with the declared diagnostic arithmetic allowance.

All 36 fine-mode comparisons with the removed-cutoff memory response lie inside
the actual analytic UV-tail bound plus the frozen numerical allowances. The
largest observed continuum difference divided by its nonzero tail bound is
0.09643. The audit independently reconstructs the bound from each recorded
source-derivative budget. It binds the original validator's 17 rejected
wrong-formula controls to these exact result files; it does not rerun those
controls.

Finite-step and finite-momentum integration accuracy remains supported by
refinement and independent-route agreement. The tiny discrepancies are
empirical evidence, not certified integration error bars. Only the omitted
UV tail has an analytic enclosure, using the separate interval derivative
budget. Stored observations allow an independent Wronskian audit at those
times; maxima between observations remain the frozen producer's recorded
checks, rather than fully archived trajectories.

The audit used Python 3.12.14 and NumPy 2.2.6. Its portable command interface is:

```sh
python audit_archived_modes.py \
  --checkpoint /path/to/checkpoint \
  --primary /path/to/primary/results.json \
  --modes /path/to/modes/results.json \
  --checks /path/to/CHECKS.json \
  --output /path/to/new/ARCHIVE_AUDIT.json
```

`--output` must not already exist. No producer module is imported. The script
pins the original public freeze and full registration; it can audit replay
outputs when supplied their matching validator record. All checks use explicit
exceptions and remain active under `python -O`.

The scope remains one fixed-geometry linear scalar-variance susceptibility.
This package does not compute a stress response, stability, particle yield,
or heating.
