# Reproduce the registered smooth FRW control

This benchmark uses exact incoming static vacuum data on a prescribed smooth
background. It does not evolve HDBLAST shell backreaction or establish thermal
radiation. `REGISTRATION.md` fixes the numerical matrix and acceptance gates;
`COMMON_ACTION_SMOOTH_FRW.md` fixes the finite-action convention and its limits.
The original protocol, source record and failed first attempt remain preserved.
`REPAIRS.json` maps the registered primary implementation to the mechanical
NumPy-dispatch repair and the original source in `code/history/`.
The completed original matrix is **FAIL**: the flat control passed, while the
curved pressure's final K96-to-K192 cutoff change exceeded its fixed tolerance.
`FOLLOWUP_REGISTRATION.md` separately registers K384 controls without altering
that original verdict. `FOLLOWUP_SOURCES.json` pins the follow-up source and
portable historical baseline before those new physical modes are executed.

Use Linux on x86-64, Python 3.12, and the exact requirements in `requirements.txt`:
NumPy 2.2.6, SciPy 1.15.3, Matplotlib 3.10.1, SymPy 1.14.0 and mpmath 1.3.0.
The numerical implementation uses extended-precision `numpy.longdouble`
subtraction; the replay verifies that its epsilon is below 1e-18. Use one BLAS/OpenMP
thread per process and leave Python assertions enabled. No runtime credentials or
network service is required after package installation.

## Verify frozen bytes before replay

Run from this checkpoint directory. In the repository checkout its ZIP is adjacent
to `MANIFEST.sha256.json`; when using a separately downloaded ZIP, pass that ZIP's
actual location. The verifier reads the archive without extraction, rejects unsafe
paths, duplicate entries and symlinks, checks exact inventory and every payload's
SHA-256 against both ZIP and checkout, and validates the original source/repair
chain, follow-up source pins and portable baseline hashes against the follow-up
wrapper's literal pins. The manifest itself and the outer ZIP/checksum sidecar are excluded from the
payload hashes. Hash verification establishes consistency with the frozen package,
not a signature or an independent scientific result.

```bash
python3.12 -m venv /tmp/hdblast-smooth-frw-venv
source /tmp/hdblast-smooth-frw-venv/bin/activate
export PIP_CACHE_DIR=/tmp/hdblast-smooth-frw-pip-cache
python -m pip install --only-binary=:all: --disable-pip-version-check -r requirements.txt
python -m pip check
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONOPTIMIZE=0 MPLBACKEND=Agg
python code/verify_package.py HDBLAST_CHECKPOINT_20261002_SMOOTH_FRW.zip
```

The inherited PV and curved-readiness packages can also be checked from the
repository root before regenerating any of their outputs:

```bash
python research/HDBLAST_CHECKPOINT_20261002_PV/code/verify_package.py research/HDBLAST_CHECKPOINT_20261002_PV/HDBLAST_CHECKPOINT_20261002_PV.zip
python research/HDBLAST_CHECKPOINT_20261002_FRW/code/verify_package.py research/HDBLAST_CHECKPOINT_20261002_FRW/HDBLAST_CHECKPOINT_20261002_FRW.zip
```

## Complete replay in a fresh directory

The existing cloud checkout is already isolated; no Git worktree is required.
Keep the frozen checkpoint unchanged. The helper verifies it before copying it,
creates a fresh source copy and writes new logs, diagnostics, arrays and figures
to a separate runtime directory. It refuses an existing destination.

```bash
checkpoint="$PWD"
replay_parent=$(mktemp -d /tmp/hdblast-smooth-frw.XXXXXX)
HDBLAST_PYTHON=python bash code/replay_smooth_frw.sh "$checkpoint" "$replay_parent/replay"
```

The helper executes, in order:

1. Frozen package/source-history verification, copied package verification and
   exact runtime/dependency checks.
2. All three symbolic scripts: `derive_local_trace.py`,
   `verify_local_identities.py` and `verify_trace_integrand.py`.
3. `verify_numeric_counterterms.py`, which checks NumPy dispatch, local subtraction
   exchange and static/flat matching without physical mode evolution.
4. The complete original DOP853 matrix, its registered conditional K=192 extension only
   when triggered, and the exact-static negative control.
5. The independently implemented Radau solver at both amplitudes and all six
   momenta, followed by 600 comparisons at the ten registered time nodes.
6. Original figures from the generated finite arrays and the explicit original
   outcome validator: it requires actual command exit 1, original scientific
   `FAIL`, flat `PASS`, curved `FAIL` only for the final pressure cutoff gate,
   and all other original acceptance gates passing. Coarse exchange and sampled
   trace diagnostics remain report-only, as the original protocol specifies.
7. The separately registered K384 follow-up: curved tight/quadrature/half-step
   solves plus a new K384 static floor control, then follow-up figures and
   nonempty terminal checks requiring its 13 acceptance gates to pass.

The fixed protocol SHA-256 passed to every physical solver is
`37fda5b30332d0079ec44db708e266788fdf5859303707bb709750fdcf59d898`.
The primary solver records its actual current source digest in `provenance.json`.
This is separate from the original pre-execution digest in `PROTOCOL_SOURCES.json`.
The distinct follow-up protocol hash is
`b067afda194895255081553e7bd9a3833b781f5277d3011b619829b4b3033d51`.
Its wrapper imports the exact repaired core without editing its bytes. The new
experiment consumes `followup_baseline/` from the frozen package, rather than
the freshly regenerated original directory: the registered baseline's exact
summary and minimal NPZ hashes are intentionally fixed, while a fresh run's
timestamps and environment metadata can differ. The baseline contains only the
losslessly exported K192 observables consumed by the registered follow-up.

Each command's actual exit code is retained in `exit-codes.tsv` and its log in
`logs/`. A failed symbolic check or preflight stops physical execution. A failed
original matrix still permits the independent comparison and plotting to produce
reviewable evidence. Its expected historical failure is recorded as
`scientific_status: FAIL` and `original_command_exit_code: 1` in
`runtime/original-outcome.json`; recognizing that outcome does not turn the
original research verdict into `PASS`. Any additional original gate failure,
execution exception or incomplete matrix prevents the new follow-up from
executing and retains a nonzero exit. The follow-up must independently complete
with `PASS`; its `FAIL` or `IN_PROGRESS` remains a failing replay. All unexpected
failures retain the underlying command's exit code. An import-only or empty
comparison does not satisfy the terminal checks.

Both matrices are CPU work rather than a distributed/GPU workload. Each uses one
numerical thread and preserves the registered tolerances, grids and one possible
original UV extension, followed by the separately registered K384 experiment.
Budget up to 40 minutes on a standard CI Linux runner; actual runtime
depends on its CPU and the contingency trigger. Full physical mode arrays can use
hundreds of MB. Keep those arrays in local replay directories or temporary Actions
artifacts. The compact frozen publication retains numerical observables,
diagnostics, provenance and array size/hash inventory rather than embedding every
large complex-mode array.

The repository workflow `.github/workflows/hdblast_smooth_frw.yml` performs this
complete replay under Python 3.12, with a 40-minute job timeout and temporary
artifact upload even if an executed gate fails. A passing CI replay means the
original expected scientific failure was reproduced, the separately registered
follow-up passed, and all required independent/implementation checks passed.
Those distinct outcomes remain finite-regulator evidence about this prescribed
benchmark under its declared finite-action convention. They do not certify the
continuum limit or validate a coupled cosmological solution.

## Optional read-only audit of full replay artifacts

After the complete replay, independently audit the saved raw modes and their
integrated observables. These commands read existing runtime outputs and do not
evolve any new modes. They require the full replay arrays; the compact publication
directories intentionally omit the large complex-mode archives.

Using `checkpoint` and `replay_parent` from the replay commands above:

```bash
python "$checkpoint/code/audit_matrix.py" "$replay_parent/replay/runtime/primary" --expected-cases 13 --output "$replay_parent/replay/runtime/original-artifact-audit.json"
python "$checkpoint/code/audit_matrix.py" "$replay_parent/replay/runtime/followup" --expected-cases 4 --output "$replay_parent/replay/runtime/followup-artifact-audit.json"
```

The helper verifies finite numeric archives, agreement of their time/momentum
grids, raw Wronskians, reintegrated mode quantities, and future occupation energy.
The required case count rejects incomplete or empty matrices. An artifact audit
can pass while `matrix_status` remains `FAIL`: the original negative scientific
outcome is preserved. This read-only audit does not certify the uncomputed tail
or replace any registered acceptance gate.

The public helper reports the run directory's basename, without serializing an
absolute local path. Its numerical audit expressions and tolerances are unchanged
from the original executed helper retained at
`code/history/audit_matrix_executed.py` (SHA256
`1741011a2478a5d83ae715ae95382a3034839c9f3d26bf0182ef87fa862b155f`).
The public CLI adds the expected-case-count guard and requires a terminal matrix
status; it was checked against the same 13+4 saved cases and synthetic empty and
incomplete inputs. Original read-only audit outputs remain in the execution
record. Public audit reports record the actual helper's SHA256.
