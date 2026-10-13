# Portable exact-bound calculation

This directory is a standard-library arithmetic capsule. Its only retained data
dependency is `public_baseline/EXACT_UNIFORM_CONSTANTS.json`, a public accepted
baseline fixture used for regression. No archive checkout, numerical solver,
network access, or third-party Python library is required.

Run from any directory:

```sh
python -B /path/to/evaluate_rescaled_bounds.py --out /path/to/new-output/RESULT.json
```

To authenticate a frozen companion proof, supply both options:

```sh
python -B /path/to/evaluate_rescaled_bounds.py \
  --proof /path/to/PROOF.md --proof-sha256 THE_PROOF_SHA256 \
  --out /path/to/new-output/RESULT.json
```

A wrong proof pin fails before writing a certificate. Hash checking verifies
document identity. Separate internal proof and arithmetic review receipts are
required to assess the mathematical claims; neither this execution nor a source
hash establishes external peer review.

The primary interval is every `0 < delta <= 1/500` for the declared exact-decimal
polynomial model. Readable constants are `K_eta < 3.196`, `K_H < 0.254`, and
`q - K_H/500 > 0.000423`. Additional historical binary64 coefficient models are
declared as distinct exact-rational models. No historical ODE trajectory error
is certified.

The primary theorem requires the eleven construction checks. The stored explicit
bulk and shell convexity checks are additional sufficient checks for the chosen
models. Bulk convexity on the actual tube already follows from scalar
contraction; separate shell convexity is not required by the monotonicity proof.
The quantitative factored shell-slope check is also redundant under that proof.
