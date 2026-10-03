# Reproduce the separate post-run analytic stress proof

This proof checks the newly derived stress formulas by exact symbolic algebra.
It is separate from the completed, preregistered numerical variance experiment
and its nine-command replay. It performs no source sampling, mode evolution,
stress quadrature, or numerical spectral-energy calculation.

Use Python with SymPy installed. The verified cloud environment uses
`/workspace/hdblast-cloud-setup/venv-frw/bin/python`. A full HDblast checkout is
required for this optional proof. The script embeds the seven inherited
source hashes also recorded in `STRESS_PROOF_INPUT_SHA256.json`, using portable
repository-relative paths. Pass the checkout explicitly using `--repo-root`;
all seven source contents must match. The standalone causal numerical replay
does not acquire this prerequisite.

From the directory containing this document and the proof script, run:

```bash
python verify_matched_stress_proposal.py --repo-root /path/to/HDblast --output /path/to/fresh-proof/normal.json
python -O verify_matched_stress_proposal.py --repo-root /path/to/HDblast --output /path/to/fresh-proof/optimized.json
```

The implementation uses explicit exceptions, so optimization leaves every
check active. The output records actual source hashes, interpreter metadata,
the checked identities, and the distinction between this analytic follow-up
and the completed numerical experiment. Normal and optimized runs should
agree on all algebraic results; their optimization metadata differs.
Choose fresh output filenames: the script refuses to overwrite existing
results. Archived outputs from this run are named `MATCHED_STRESS_PROOF.json`
and `MATCHED_STRESS_PROOF_OPTIMIZED.json`.

The earlier unsuccessful general-r symbolic pressure comparison is preserved
separately in `development/` with its corrected result and explicit source
reconstruction provenance. Those records are not relabeled as this later
portable proof or as a failed physical experiment.
