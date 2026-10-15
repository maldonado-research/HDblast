# Reproduce the inertial first-variation checkpoint

Use the tested Python 3.12.14 environment and SymPy 1.14.0. The scientific replay uses integer,
rational and symbolic calculations. Matplotlib is needed only for the optional
illustrative plot, whose pixels are not mathematical validation inputs.

Obtain the capsule ZIP and its external SHA256 and manifest SHA256 from the
reviewed publication receipt or the dedicated GitHub Actions workflow. Verify
the ZIP hash before extraction. The capsule contains the scripts, full generated
data, reviewed notes and minimal pinned prior mathematical inputs. No private
millennium-project audit or raw retrieved paper body is needed.

In the extracted capsule's parent directory:

```sh
python -I -B -m pip install --only-binary=:all: -r HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015/requirements.txt
python -I -S -B HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015/verify_manifest.py \
  --expected EXTERNAL_MANIFEST_SHA256
python -I -B HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015/reproduce.py \
  --manifest-sha256 EXTERNAL_MANIFEST_SHA256 --output-dir fresh-replay
```

Supply the complete external pin in place of `EXTERNAL_MANIFEST_SHA256`. The
output directory must be new or empty and outside the frozen capsule. Each
process receives its own authenticated copied script/data tree; no archived
producer is imported or executed. Source authentication runs before any local
capsule execution and again after all checks.

The eight positive processes in each execution mode are the primary certificate
checker, independent complete-cell reconstruction, full JSON/CSV generator,
Sol free-comparator checker, Astra bound checker, independent referee residual
checker, archive derivative producer and archive symbolic identity checker.
The runner compares all scientific JSON fields against frozen goldens and
compares the complete generated JSON/CSV bytes. It verifies execution-mode
metadata before normalizing those fields for the normal/optimized comparison.
The Sol wrapper is retained as provenance; the runner directly performs its
normal/optimized comparator processes and named controls to avoid counting the
same nested runs twice.

The accepted roster is recorded in `REPLAY_ROSTER.json`: 16 positive processes,
22 external corrupted-formula/input/certificate rejections, 18 embedded controls,
and 112 explicitly enumerated positive groups per mode. Generator grid cells
and archive subdivision panels are reported separately. Each external control
must exit with code 1 and its exact named failure; unrelated exceptions or stale
passing output files fail the replay. Full data comprise 4,353 continuous cells:
one causal prefront cell, 4,096 support cells and 256 postpulse cells.

The output receipt lists every result/log/CSV artifact by relative path, length
and SHA256. Its copied executable workspaces are scratch, not delivery evidence.
These executions validate finite calculations and conditional enclosures. They
do not formally prove the retained analytic PDE assumptions, establish a
finite-amplitude nonlinear waveform or effective nonlinear epsilon interval,
calibrate an observed physical time, or infer a higher-dimensional cosmological
origin. The wider state is newly specified; no narrow-state response coefficient
is transferred to it.
