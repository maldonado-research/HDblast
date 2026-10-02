# Reproduce the stationary quantum checkpoint

Use Python 3.12 with NumPy 2.2.6, SciPy 1.15.3, SymPy 1.14.0 and mpmath 1.3.0,
as pinned in `requirements-replay.txt`; Matplotlib 3.10.1 supports the
optional result plots. The prospectively frozen
`requirements.txt` retains its original lower-bound dependency declarations.
The separate replay pin file records the tested environment without changing
those registered bytes.

Reproduction requires a full HDBLAST checkout containing this checkpoint and
its pinned inherited DESITTER, JUNCTIONS and classical background inputs. A
checkpoint-only ZIP is **not self-contained**. Use the stationary checkpoint
branch or a later commit containing it. The scientific baseline
`e24643d579f7560db719775b2998c0ffa22898b0` predates this checkpoint; checking
out that baseline alone is insufficient. The replay verifies the inherited
input bytes rather than requiring checkout HEAD to equal that baseline.

From the full repository root:

```bash
python3.12 -m venv /tmp/hdblast-stationary-venv
/tmp/hdblast-stationary-venv/bin/python -m pip install -r research/HDBLAST_CHECKPOINT_20261002_STATIONARY_QUANTUM/requirements-replay.txt
/tmp/hdblast-stationary-venv/bin/python research/HDBLAST_CHECKPOINT_20261002_STATIONARY_QUANTUM/code/replay_checkpoint.py \
  --repo /path/to/HDblast \
  --output /fresh/output
```

Replace the example repository and output paths. The output location must be
unused so that prior evidence is preserved. The replay performs the registered
producer run, the historical validator's expected failure, corrected primary
validation normally and with Python optimization, exact symbolic checks, independent
acceleration roots, proper-time endpoint checks and the post-run comparison
audit. Source/root runs and high-precision validation can take several minutes.
Read the generated status records and logs, not only the process exit code.

## Expected components

| Component | Coverage and interpretation |
|---|---|
| Primary roots | 48 solves covering the registered 18 model points and prescribed integration checks |
| Wrong-model solves | 3 converged wrong-model controls tested against the correct junctions |
| Historical validator | Expected exit 1 with the original 48 diagnostic-keyword TypeErrors, retained as failed implementation evidence |
| Corrected primary validator | `validate_stationary_v2.py` passes 1,584 gates normally and with `python -O`, and separately classifies shifts, nonlinear departure and hierarchy |
| Exact symbolic checks | 51 identities and 20 wrong-formula controls, also checked with `python -O` |
| Independent acceleration roots | 12 solves: 6 selected model points at 2 resolutions |
| Proper-time source checks | 8 evaluations at 4 independent endpoints; 48 gates |
| Independent comparison audit | 564 gates covering source/Hessian reconstruction, selected roots and wrong-model residuals |

An accepted root and a detected observable shift are distinct outcomes.
Inspect `shift_diagnostics` for resolution limits and nonlinear classification.
The positive-gamma hierarchy screen also has its own outcomes; a mathematical
root is retained when this physical-scale diagnostic fails. Do not replace a
reported unresolved shift or failed hierarchy screen with a blanket pass.

Timestamps, absolute paths and platform-sensitive floating-point values need
not match the archive byte for byte. Compare the specified numerical gates
and status records. A passing replay supplies empirical numerical evidence;
it does not create a certified error enclosure or a proof of existence.

## Original evidence and later corrections

`PREREGISTRATION.json` pins the producer, original validator, experiment,
model, dependencies and inherited scientific inputs. `FULL_REGISTRATION.json`
also pins the initial theory and independent files. Its 17 entries plus the
manifest itself form the 18 initial checkpoint files. Public commit
`b4f77f5826e0a603400513832aa53b0673df7533` preceded all new source evaluations
and radial/root calculations in this round. Symbolic algebra preceded that
commit and is identified accordingly. Do not run `stationary_shell.py --freeze`
when reproducing the preserved registration.

The initial frozen validator produced
[`outputs/primary/CHECKS.json`](outputs/primary/CHECKS.json) with `FAIL`: the
diagnostic metadata keyword `condition=` collided with the `gate` helper's
positional argument, yielding 48 TypeErrors. It could not finish the intended
per-root validation. The original source, failed checks and run evidence are
retained. `validate_stationary_v2.py` changes only that diagnostic keyword to
`condition_number=`. Repair commit
`d8a0b28fa7440215108ddfff51355816c932b541` was published after the original
run and before corrected validation. The replay retains the historical
validator's expected failure and uses v2 for completed scientific validation.
Scientific gates and model choices
are those of the original registration.

The independent comparison audit was developed after the root/source runs.
Its first version failed to serialize a NumPy boolean to JSON; the failed
script, log and exit status remain under `independent/`.
`audit_primary_v2.py` corrects JSON scalar serialization and adds portable
arguments while retaining the same scientific comparisons. This audit is
post-run evidence, not a prospectively registered calculation. The original
independent root and proper-time executables needed no scientific correction.

The original producer used Python 3.12.14, NumPy 2.2.6, SciPy 1.15.3 and
mpmath 1.3.0. Original independent roots record NumPy 2.3.5, SciPy 1.16.3
and mpmath 1.3.0. Preserve those version facts when interpreting the archived
run; the pinned replay uses the single environment specified above.

Literature search metadata and curated notes are archival context. Scientific
replay does not redownload copyrighted primary papers, reproduce prospective
timing, or constitute external peer review. Citation and Zenodo metadata are
prepared for reuse; no checkpoint DOI or completed deposition is implied.
