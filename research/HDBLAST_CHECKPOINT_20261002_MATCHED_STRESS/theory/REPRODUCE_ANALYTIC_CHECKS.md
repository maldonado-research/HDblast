# Analytic preparation checks

These commands verify symbolic identities and rational derivative envelopes.
They evaluate no physical source sample, response integral, mode, stress, or
deferred interval budget. Use fresh output filenames; evidence is not
overwritten. The validated interpreter is
`/workspace/hdblast-cloud-setup/venv-frw/bin/python` (Python 3.12, SymPy 1.14).

The coefficient/finite-cutoff proof checks the inherited source hashes listed
in `SOURCE_PINS.json`; provide a full checkout as `--repo-root`:

```bash
python verify_stress_tail_algebra.py --repo-root /path/to/HDblast --output /path/to/fresh/TAIL.json
python -O verify_stress_tail_algebra.py --repo-root /path/to/HDblast --output /path/to/fresh/TAIL_OPTIMIZED.json
```

The saved normal and optimized results both pass 25 exact identities and six
wrong-formula controls. The pulse-norm directory contains its independent
polynomial/root verifier and rational envelope proof, with their normal and
optimized results: 18 identities, six complete root certificates and seven
mutations, plus the exact rational global envelopes. Consult
`pulse-norms/PULSE_DERIVATIVE_NORMS.md` for those proof interfaces.

The global stress ceiling script performs rational arithmetic only and uses
the independently proved uniform derivative ceilings. Its archived output
and log show the K=256 feasibility result. It writes a fresh
`GLOBAL_STRESS_CEILINGS.json` beside itself; replay it in a copied temporary
directory with the pulse certificate and without an existing output file.

The response-facing libraries `stress_tail_bounds.py` and
`pulse-norms/pulse_derivative_budgets.py` are intentionally not exercised by
these analytic checks. Their source-jet and interval evaluations belong to
the publicly frozen numerical run. The main library accepts its pulse
dependency either beside itself or in the `pulse-norms/` subdirectory.

Only the omitted momentum band is bounded by these libraries. Time and
momentum quadrature, mode evolution, derivative evaluation, subtraction
cancellation and any Ward-ledger discretization need independent evidence.
