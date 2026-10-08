# Conditional physical-model progress

This directory defines a new flat 4+1-dimensional scalar half-space testbed,
separate from the current HDBLAST de Sitter pulse calculation.

- `HALF_SPACE_RESPONSE_THEOREM.md`: full action, exact retarded reduction,
  stability proof, spectral density, Gaussian state/noise rules, energy ledger,
  restricted-rival theorem and exact four-dimensional continuum replica.
- `MODEL_TESTBED_CONTRACT.json`: conditional model domain and assessment against
  the seventeen inherited physical-model fields. It is not an observational
  test registration.
- `verify_half_space_controls.py`: manufactured algebra, sign, bound-mode,
  energy, spectral and covariance controls. No physical input or retained array
  is opened. Numerical quadrature is explicitly diagnostic, not certified.
- `CONTROLS_NORMAL.json` and `CONTROLS_OPTIMIZED.json`: replay receipts.
- `SOURCE_MANIFEST.json`: hashes the review candidate's source files, receipts
  and the inherited documents read to delimit scope. This is a review freeze,
  not a public physical-execution GO.

Replay from any directory with SymPy and mpmath installed:

```sh
python /workspace/hdblast-research-work/continuation-trajectory-20261008/physical-model/verify_half_space_controls.py --output /tmp/half-space-controls.json
```

The precise advance is a joint response/noise restriction from an explicit
stable action, together with an explicit microscopic four-dimensional rival
that reproduces it. No geometric dimensional origin is identified. A prescribed
bulk-to-brane blast, a coupled cosmology, gravity, backreaction, reheating, an
observational likelihood and a finite-resolution discrimination bound remain
undeveloped. The continuum bath's equal-time raw noise variance diverges;
smoothed measurements are required. KMS equilibrium is a state premise.

The scalar numerical track is unchanged: metric calibration **FAIL**; full
continuous certificate **UNRESOLVED**; higher-dimensional Big Bang cause
**NOT_ESTABLISHED**; external novelty **NOT_ASSESSED**.
