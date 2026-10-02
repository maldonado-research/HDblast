# Reproduce the fixed-geometry causal-response calibration

Use Linux and Python 3.12 with the exact versions in `requirements.txt`:
NumPy 2.2.6, SciPy 1.15.3, mpmath 1.3.0 and SymPy 1.14.0. Plotting uses
Matplotlib 3.10.1; the complete replay environment is pinned separately in
`requirements-replay.txt`. Linux is specified because the independent resource
gate measures `ru_maxrss` in KiB.

The scientific checkpoint contains its numerical code and faithful copies of
the required inherited protocol/reference inputs. Its ZIP can be extracted
and replayed without the author's private archive or older raw arrays.
Repository paths in `PREREGISTRATION.json` record provenance; local copies
are pinned in `FULL_REGISTRATION.json` and the outer manifest.

From a checkout containing this completed package, or its extracted ZIP:

```bash
python3.12 -m venv /tmp/hdblast-causal-venv
/tmp/hdblast-causal-venv/bin/python -m pip install -r /path/to/HDBLAST_CHECKPOINT_20261002_CAUSAL_RESPONSE/requirements-replay.txt
/tmp/hdblast-causal-venv/bin/python /path/to/HDBLAST_CHECKPOINT_20261002_CAUSAL_RESPONSE/code/replay_checkpoint.py \
  --output /fresh/external/output
```

Choose an unused output path outside the source checkpoint. The replay verifies
exact file membership, all outer payload hashes, the prospective hashes and
inherited local-copy hashes before running fresh copies. It reruns the exact
symbolic checks, primary formulas, independent forced modes, validation normally
and with optimization, post-run archive audit, summary and registered-point
figures. Read the generated `VALIDATION.json`, execution records and logs.
Matching archived timestamps, absolute paths or ZIP timestamps inside fresh
NPZ files is not required. Required scientific gates and coverage must pass.

## Original experiment and immutable choices

Public freeze: `a01d17e015f070be2ad2018745e877148fe2f4e4`.
Full registration SHA256:
`f59abda221a6a0296e839bc3022298bdf81fac913feb3877ac8f1485f4d71da6`.
Independent manifest SHA256:
`5bf042c7fe686e23dc469031eef81f0e682d56b4db4ff4c0381cf0d8a2ecee22`.
The model, two pulses, six times, three cutoffs, integration settings, tolerances,
state diagnostics and stop rules are in `EXPERIMENT.json` and the registrations.
Do not rewrite those files or substitute a different reference/state.

The original primary run and all four independent runs completed. The
validator records 36 comparisons and 17 wrong-formula/control detections.
Normal and optimized results agree scientifically. Exact theory has 46
identities and 14 wrong-formula detections. The failed symbolic-development
attempt is separate: its original output capture survives, and the old source
is explicitly a reconstructed version rather than an asserted original file.

## What numerical agreement does and does not mean

The independent route evolves forced canonical modes with an exact free
propagator and local-step Gauss quadrature. It combines mode and subtraction
terms before summing the finite momentum integral. Its two settings refine
time and momentum quadrature together, so their spread is an empirical
combined refinement check. It does not isolate every error contribution.

The primary finite-cutoff expression exchanges two finite integrals; it does
not use the removed-cutoff result to manufacture its answer. The logarithmic
and regular finite-part memory formulas provide a separate internal identity
check. At each point, the actual combined omitted-tail bound is computed from
the pulse derivatives with directed intervals. That bound covers only omitted
momentum. QUADPACK estimates, floating-point allowances and mode refinements
remain distinct, and a passing replay is not a certified total error enclosure.

The 17 controls include algebraic/reference mutations and selected wrong
response/state comparisons. They are not 17 complete reruns of mutated mode
solvers. The spectral energy/work companion is analytic and has no executed
energy-yield test. No output of this replay supplies the missing stress and
metric kernels, a coupled shell solution, thermalization or stability.

## Optional analytic follow-up

The [matched stress proof](theory/MATCHED_STRESS_PROOF_REPRODUCE.md) is separate from the nine-command numerical replay. It passes 21 exact identities, eight mutations and seven inherited-source hash checks. It requires a full HDblast checkout supplied with `--repo-root`; it neither adds dependencies to the standalone causal replay nor performs a numerical stress calculation.
