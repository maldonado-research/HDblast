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

code/trajectory_audit.py reads the original B3 .npz files, not just saved summary points. It requires Python 3 and NumPy, with pickle disabled. The package includes the eleven numeric timeseries files and matching summaries needed for the nine selected histories, including restart parents. Their original SHA-256 values are preserved; data/PROVENANCE.md records the source. From this checkpoint folder, install NumPy 2.2.6 in your Python environment and run:

```sh
python3 code/trajectory_audit.py --scan-root data/B3_Y_scan --source-ref bundled-september-baseline --output replay_trajectory.json
```

Output path/source-ref metadata differs from the original private-run report; compare numerical values by each run's tag. Its saved report and precise limitations are in TRAJECTORY_AUDIT.md.

This audit measures intervals in sampled trajectories; it does not rerun the field equations or establish continuum convergence. Read its sampling and reliability policy before interpreting any duration.

No new full five-dimensional simulation or self-consistent renormalized quantum calculation is claimed.

## Download integrity

The curated ZIP contains this folder's public scientific files and MANIFEST.sha256.json. Verify its CRCs, exact inventory and per-file hashes with Python's standard library:

```sh
python3 code/verify_package.py HDBLAST_CHECKPOINT_20261001.zip
```

The public package redacts one private archive-path metadata string in the pulse source/output and corrects its command comment. Integration, inputs, acceptance logic and numeric outputs are unchanged. The original source and its root replay remain in the private checkpoint.
