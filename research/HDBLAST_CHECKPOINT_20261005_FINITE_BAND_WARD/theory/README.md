# Analytic finite-band Ward/stress preparation

The current derivation is `FINITE_BAND_ERROR_TRANSPORT.md`. Use the release
receipts `SYMBOLIC_RELEASE_PHASE_NORMAL.json` and
`SYMBOLIC_RELEASE_PHASE_OPTIMIZED.json`, which bind its final bytes. The verifier
has 55 exact symbolic/rational checks and 12 rejected symbolic mutations.
It evaluates no source callback, mode, physical source integral or saved array.

The result consists of the complete off-shell conservation-error identity,
exact finite-band retarded kernels for independently defined density and
pressure, initial-state/residual error transport and conditional analytic
source-Taylor-tail bounds. At the unchanged K=256 and throughout
[-9/2,-7/2], the degree-24 **analytic truncation contribution alone** is bounded
by 6e-28 for scaled density, 2.2e-25 for scaled pressure and 6e-28 for the
integrated source mismatch work. Exact rational bounds for K=64/128/256 are
stored in the receipts. These are not measured errors or new acceptance gates.

The full twelve-case numerical pressure/contact certificate remains
unresolved. Incoming-state, coefficient/arithmetic, residual, contact,
momentum, time and geometry errors remain separate. Conservation and nine
momentum probes cannot exclude a normalized incoming state perturbation
between probes; section 7 proves this explicitly.

`INPUT_PINS.json` hashes seven inherited analytic documents read as plain
text. The full subtraction inventory and Cauchy/calculus premises are
conditional inputs, not proof-assistant formalizations by this script.

To repeat the pure proof, with SymPy 1.14.0 and fresh output paths:

```sh
python -B verify_finite_band_theory.py --output NEW_NORMAL.json
python -B -O verify_finite_band_theory.py --output NEW_OPTIMIZED.json
```

The program refuses to overwrite an existing receipt. It has no scientific
data input option and imports no project module. `development/` and the
earlier top-level receipts preserve preparation history; see
`PROOF_HISTORY.md`. The independent root implementation supplies a separate
standard-library exact-polynomial/rational check, when completed.

This directory has no remote publication authority or credential material.
The checkpoint can be published only with its full conditional scope and
unchanged historical numerical failures visible.
