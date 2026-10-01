# Reproduce the 1 October checkpoint

The analytical bounds and pulse controls were actually evaluated in JavaScript/V8 in this session. The scripts are saved as standalone JavaScript without third-party dependencies. Use Node.js 20 or later in this checkpoint folder:

```sh
node code/cohort_verify.js
node code/duration_verify.js
node code/mode_controls.mjs
```

Each script prints JSON to standard output. Compare it with the corresponding saved output under outputs/. Small last-digit differences across engines are acceptable; the scientific tolerances are recorded in REGISTRATION.md. The code defines its own numeric inputs, also supplied with provenance in inputs.json.

C3's exact source is the source executed for all registered cases. Registration commit ceaeb64adf76b901ded8fa7d1bf1c913a29d22e4 preceded its first execution. C1 and C2 were exploratory before registration and are explicitly disclosed as such.

## Full trajectory audit

code/trajectory_audit.py reads the original B3 .npz files, not just saved summary points. It requires Python 3 and NumPy, with pickle disabled. The public checkpoint does not redistribute the older large raw arrays. Their filenames and required columns are documented by the script; the author's private archive supplies them to the read-only runner. Its saved report and precise limitations are in TRAJECTORY_AUDIT.md.

This audit measures intervals in sampled trajectories; it does not rerun the field equations or establish continuum convergence. Read its sampling and reliability policy before interpreting any duration.

No new full five-dimensional simulation, renormalized coupled quantum calculation, or live literature search is claimed.
