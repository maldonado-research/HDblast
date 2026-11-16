# Reproduce the evolved-parameter capsule

Use an externally expected manifest digest from the trusted workflow
at the tested commit. Run every command with fresh outputs outside the
capsule. Finite replay and written continuum proofs have separate roles.

Use Python 3.12.14 and the exact binary-wheel hashes in
`requirements-checks.txt` (SymPy 1.14.0 and mpmath 1.3.0):

```
python3.12 -m venv /tmp/hdblast-parameter-env
/tmp/hdblast-parameter-env/bin/python -m pip install --require-hashes --only-binary=:all: -r research/HDBLAST_EVOLVED_PARAMETER_20261116/requirements-checks.txt
```

Obtain the final manifest SHA256 from
`.github/workflows/evolved-parameter-20261116.yml` at the trusted tested
commit. An externally expected digest authenticates the capsule; a digest
computed from an untrusted local capsule does not supply that expectation.
Paste that trusted 64-character digest when this command reads a line:

```
read -r PIN
```

Use fresh outputs outside the capsule:

```
/tmp/hdblast-parameter-env/bin/python -I -B research/HDBLAST_EVOLVED_PARAMETER_20261116/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/hdblast-parameter-normal.json
/tmp/hdblast-parameter-env/bin/python -O -I -B research/HDBLAST_EVOLVED_PARAMETER_20261116/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/hdblast-parameter-optimized.json
cmp /tmp/hdblast-parameter-normal.json /tmp/hdblast-parameter-optimized.json
/tmp/hdblast-parameter-env/bin/python -I -B research/HDBLAST_EVOLVED_PARAMETER_20261116/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-parameter-auth-normal-fresh
/tmp/hdblast-parameter-env/bin/python -O -I -B research/HDBLAST_EVOLVED_PARAMETER_20261116/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-parameter-auth-optimized-fresh
/tmp/hdblast-parameter-env/bin/python -I -B research/HDBLAST_EVOLVED_PARAMETER_20261116/checks/run_math_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-parameter-math-normal-fresh
/tmp/hdblast-parameter-env/bin/python -O -I -B research/HDBLAST_EVOLVED_PARAMETER_20261116/checks/run_math_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-parameter-math-optimized-fresh
```

Every output path must not yet exist. Use a new private evidence root for
later campaigns rather than modifying an immutable completed evidence
index. The optimized harness must retain optimization in its child
checker processes. Compare normal/optimized main reports byte for byte,
and compare the mathematical case names, changed claims and intended
rejection diagnostics across both campaigns. Their receipts explicitly
record different child optimization levels, so those whole receipts need
not have identical bytes.

## Authentication and evaluator boundary

`checks/check_checkpoint.py` authenticates the external manifest,
exact closed file roster, regular nonsymlink members, safe paths, strict
JSON and every size/hash before importing the pinned local
`checks/evaluate_parameter.py`. That evaluator checks the immutable
source/member joins and invokes `checks/parameter_oracles.py`. The source inventory contains 50 snapshots and 141 qualified joins;
the executed report verifies those counts. No source producer
is imported or executed to supply an antecedent result.

The own finite oracle defines 194 baseline named checks and 81
changed-claim mathematical/interface controls. Thirteen separate negative
parser/type guard cases are recorded by the private translation review.
The full public evaluator has 556 named checks; its executed roster
includes source, provenance and regularity contracts. The own oracle's
standalone claim evaluation does not replace closed-capsule authentication
or fresh subprocess controls.

Five protocol faults change source bytes, add an unlisted member, remove
a member, replace a source with a symlink, or stale the external manifest
pin. Each rejects before evaluator import and reaches its intended
authentication diagnostic.

For each of the 81 mathematical faults, the separate harness copies a
complete capsule, changes exactly one genuine claim, rebuilds its
manifest and passes that new expected digest to a fresh full-checker
process. Authentication must succeed, followed by rejection at the
intended named mathematical/interface check. A parser, syntax,
schema/engine guard or authentication failure is not a substitute for
that target. Fresh fixtures and child receipts must preserve this
distinction in both Python modes. CI repeats the externally pinned
campaign with the same exact wheels and final manifest digest.

## What this reproduces

Hashes establish byte identity and qualified provenance. The exact
finite two-field/five-jet models test normalized-source derivatives,
principal cross terms, residual subtraction and declared derivative
interfaces. They do not evaluate all actual sixteen-field arrays,
prove the full mixed energy estimates, prove compactness or construct
the actual parameter derivatives. Those conclusions are supplied by
the four written proofs and their recorded internal review.

The replay evaluates no new numerical PDE-energy ceiling, parameter
interval, observable response, noise model or cosmological origin. The
actual analytic conclusion is lower-norm `C2` parameter dependence and
its fixed-base physical chain rule, with a qualitative field remainder.
The existing Zenodo edition remains unchanged.
