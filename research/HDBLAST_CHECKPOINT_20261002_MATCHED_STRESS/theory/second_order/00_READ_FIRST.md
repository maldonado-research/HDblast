# Proposed next task: quadratic canonical energy

This folder is an **unregistered, analytically developed follow-up**, prepared
after the matched-stress input freeze. No new source, response, Fourier
spectrum, mode, work integral, or energy yield has been numerically evaluated
for it. The concrete future settings in the proposal are candidates that
require a new public code/input freeze before use.

Read `PROSPECTIVE_SECOND_ORDER_CANONICAL_ENERGY.md` for the proposed three-way
comparison: positive Fourier energy, independently evolved first-order
Bogoliubov amplitude, and matched finite-cutoff source work with its exact
background drift removed. It includes finite grids, diagnostic allowances,
failure gates, and constructive analytic energy tails.

The main scientific distinction is perturbative and physical: the preceding
first-order minimal stress can be signed coherent polarization. Positive
canonical excitation energy occurs at second order and does not by itself
give the second-order physical stress. For the positive pulse, an
occupation-only variance is even infrared divergent; the full fixed-time
coherent response cancels that divergence. The independent note in
`independent/IR_COHERENCE_AUDIT.md` checks this point analytically.

The proposal does not infer heating, thermalization, an equation of state,
actual endpoint or metric kernels, a coupled solution, or stability. It does
not change the inherited finite action or the prior experiment's inputs.

## Exact checks actually executed

`verify_second_order_energy_algebra.py` checks exact algebraic normalizations,
Bogoliubov phase conventions, the finite-K contact derivative, the work drift,
physical/minimal post-pulse Ward and trace reductions, the infrared
cancellation, and endpoint-bound coefficients. It uses explicit exceptions
and checks the nine faithful inherited input copies in `reference_inputs`.
It does not prove every analytic theorem or validate future numerical
quadratures. Proofs of Fourier positivity and smooth compact-source decay
are given in the accompanying prose with their assumptions.

The new script passed 21 identities and rejected 8 algebraic mutations in
normal and optimized Python. Full JSON, stdout,
and exit records are `EXACT_ALGEBRA{,_OPTIMIZED}.{json,log,exit}`. The local
FRW interpreter used was
`/workspace/hdblast-cloud-setup/venv-frw/bin/python`; its Python and SymPy
versions are recorded in the reports. There was no failed execution of this
new proof. Earlier inherited development failures remain preserved in their
original checkpoints; this folder does not replace their record.

The separate independent bounded formula review is recorded in
`independent/FORMULA_REVIEW.md`; it found no mathematical blocker in the
energy, drift, post-pulse stress, and endpoint-bound formulas. It does not
certify the future numerical implementation or its candidate allowances.

Portable replay with Python and SymPy installed:

```sh
python verify_second_order_energy_algebra.py --output local-check.json
python -O verify_second_order_energy_algebra.py --output local-check-optimized.json
```

Add `--repo /path/to/HDblast` to check the original inherited paths as well.
The archived runs used that option; it is optional for portable replay. The
proof does not depend on the old checkout when checking its faithful local
copies. `SOURCE_PINS.json` identifies every original path and SHA-256.
`ARTIFACT_SHA256.json` records this staged bundle, excluding itself.

No environment configuration change or new dependency installation was
needed; the existing FRW environment supported the analytic work.
