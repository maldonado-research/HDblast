# Reproduce the curved-source readiness audit

Use Python 3.12, NumPy 2.2.6, SciPy 1.15.3 Matplotlib 3.10.1 and mpmath 1.3.0. No private Drive or private GitHub access is needed for the curated package.

```sh
python -m pip install numpy==2.2.6 scipy==1.15.3 matplotlib==3.10.1 mpmath==1.3.0
python code/verify_package.py HDBLAST_CHECKPOINT_20261002_FRW.zip
python -m compileall -q code
python code/verify_protocol.py
python code/verify_inputs.py
python code/verify_frw_counterterms.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python code/smooth_controls.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python code/independent_smooth.py
python code/geometry_audit.py --registration REGISTRATION.md --registration-sha256 939e7f6ad178c32cdcf638d01b8694d2194e805f0aad38964416518e3b72a2b5 --source-ref local-replay
python code/verify_archive_interference.py
python code/check_reconstruction_fidelity.py
python code/compare_results.py
python code/smoothing_limit.py
python code/plot_results.py
```

The input NPZ files are loaded with allow_pickle=False. Original summary files retain their original bytes, including NaN values describing the eventual failed late-time run; the loader excludes the unreliable interval and the audit uses only [0,6.9]. Outputs use strict JSON with nonfinite archival metadata represented as null. Nothing extrapolates the archived physical solution.

The exact-rational verifier uses the Python standard library. It checks selected local jets rather than proving a universal identity. The two oscillator integrators are independent implementations. Their toy frequency step tests the scattering controls, not full archived high-momentum evolution.

geometry_audit.py writes all eight reconstruction reports, per-knot jumps, full CSV/NPZ curves and a selected set of permanent report samples. plot_results.py renders those executed arrays. Actions artifacts retain full regenerated arrays temporarily; all are reproducible.

The constant-reference counterterm construction is checked algebraically. No absolute renormalized FRW stress, admissible global state, quantum backreaction, thermal radiation bath or new 5D solution is produced by these scripts.

The publication ZIP can be checked with code/verify_package.py after extraction; it contains the original inputs as well as the complete scripts and curated outputs. The manifest hashes every payload file.

The independent JavaScript leading-jump calculation can also be reproduced with Node.js: node code/archive_interference.js. Its frozen output is cross-checked with scipy.special.sici by verify_archive_interference.py. Neither is exact archived high-k mode evolution. The frozen primary report and figures came from the repaired source; generated output files are overwritten by a replay, while the checked-in ZIP preserves the original record. The archive manifest intentionally excludes the manifest itself and the outer ZIP.

Run the commands from the checkpoint directory. In a repository checkout the ZIP is in that directory; after extracting a separately downloaded ZIP, pass its actual path (often ../HDBLAST_CHECKPOINT_20261002_FRW.zip) to verify_package.py. Verification also checks that the extracted/checked-out payload matches the frozen ZIP; run it before regenerating outputs.
