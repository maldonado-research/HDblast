# Reproduce the action, junction and stationary bridge checks

Use Python 3.12 and the exact packages in `requirements.txt`. The cloud setup
already supplies them in `/workspace/hdblast-cloud-setup/venv-frw`; activation
is `source /workspace/hdblast-cloud-setup/activate-frw.sh`. Elsewhere, install
the requirements in a virtual environment. No quantum or gravitational
simulation is needed for these algebra checks.

From the repository root, choose an output path that does not exist and lies
outside the checkout:

```sh
python research/HDBLAST_CHECKPOINT_20261002_JUNCTIONS/code/replay_junctions.py \
  --repo "$PWD" --output /tmp/hdblast-junction-replay
```

The runner verifies the checkpoint manifest and inherited source hashes before
copying the package into that fresh directory. It executes the primary action
derivation normally and with `python -O`, the separate action review, the
stationary bridge, and its separate review. Each command's exit code and log
are retained; a failing check stops the replay. `VALIDATION.json` records a
completed replay. The runner requires nonempty expected check counts and
unchanged action reports under optimization. Historical outputs are preserved.

The action checks give 21 identities and 14 detected negative controls. The
independent action derivation gives 21 identities and nine detected controls.
The bridge gives 33 algebra/archive checks and five detected controls. The
separate static review's counts and residuals are in its retained JSON report.
These sets overlap intentionally and must not be added up as independent
physical experiments. The recorded runtime used Python 3.12.14, NumPy 2.2.6,
SymPy 1.14.0 and mpmath 1.3.0.

`evidence/SOURCE_PINS.json` and `static_bridge/SOURCE_PINS.json` identify inputs
at `e17a01b428bb8049e919c42376ab0359e41d613c`. A later checkout is acceptable
only if those files retain their pinned bytes. Archive extraction without Git
metadata is supported by the static verifier; observed HEAD is provenance,
not a substitute for the individual file checks.

The next quantum expectation calculation, constrained initialization and
dynamic shell/bulk experiment remain UNRUN. This replay neither changes that
status nor supersedes the original failed smooth-FRW cutoff result.
