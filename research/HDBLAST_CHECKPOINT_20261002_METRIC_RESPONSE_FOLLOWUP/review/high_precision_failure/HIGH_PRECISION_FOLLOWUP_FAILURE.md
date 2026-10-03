# High-precision metric followup: independent resource failure

The registered `HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_FOLLOWUP` has
status **FAIL**. Its high-precision primary completed all twelve observation
rows. The independent producer completed the positive_B coarse and fine
evolutions and saved their raw archives, then stopped because its recorded
peak memory exceeded the unchanged 256-MiB limit. This resource failure
prevented complete independent validation; it is not evidence of physical
instability or rejection of the HDBLAST hypothesis.

The followup public freeze is
`5ff571cdb7729398834ff3c025b50191c4ee29f2`, with registration SHA-256
`e9eee4aa5d0db84e84039d177a1ae12e8bf7ac1d41bb9ae032287daa8a667a3b`.
Its independent manifest remained
`4760f2a7e864dbe81ad8911797499dcd9fcc9e37a45421a607a101d35e6ce507`.
The saved replay directory is
`/workspace/hdblast-research-work/metric-followup-complete-replay-001`.

## Exact completed coverage

The replay started at 19:26:53.564014 UTC on 2 October 2026 and stopped at
19:31:02.686896 UTC. Of 32 planned commands, 28 were attempted and 27
passed: the same 26 pure/preflight commands, followed by the physical
primary producer. The independent_modes command then exited 1. No timeout
or operating-system out-of-memory kill is recorded; the producer deliberately
raised `RuntimeError` at its resource gate.

The primary completed both pulse profiles at all six observation times,
with continuum and K=64,128,256 outputs: **12 rows and 36 finite-K points**.
Its internal elapsed time was 37.526871408 seconds and its process elapsed
time was 38.187905336 seconds. The saved continuum algorithm used separate
50- and 70-digit mpmath evaluations with analytic integer source-polynomial
arithmetic, tanh-sinh quadrature, endpoint subtraction and the squared
endpoint transformation. Its output records empirical precision/quadrature
diagnostics. Producer completion is not a PASS of the unexecuted cross-route
validators or a certified total error bound.

The independent producer saved two complete positive_B run records:

| Setting | Time steps | Momentum nodes | Observation rows | Finite-K points |
|---|---:|---:|---:|---:|
| coarse | 576 | 8192 | 6 | 18 |
| fine | 1152 | 16384 | 6 | 18 |

Their preserved archives are:

- `metric_modes_positive_B_coarse.npz`: 33,643,203 bytes; SHA-256
  `0b366d473a672fa110a499759aa61fb8e879b36a7104f9049fee1c5147337102`.
- `metric_modes_positive_B_fine.npz`: 65,892,627 bytes; SHA-256
  `8bd8daac59d12cf0af5acfc7ff33185e0a18bd6910a27c7e30963fea1c89ac35`.

These run records retain their direct responses, histories, raw mode and
subtraction data, and per-step Wronskian diagnostics. They are useful partial
evidence but do not constitute the scheduled four-run independent response
calibration.

## Resource stop and unexecuted work

The independent failure capture records:

```
peak_rss_kib=306088,
frozen_budget_kib=262144,
excess_kib=43944,
producer_elapsed_seconds=153.807740940,
process_elapsed_seconds=153.936909779,
exception='Frozen 256-MiB peak-memory budget exceeded'.
```

The excess is approximately 42.91 MiB. The active-run marker was cleared
after the saved positive_B fine run, and the failure occurred at the main
loop's resource check. No signed_uB independent physical evolution began.
The independent combination of all four runs, numerical refinement and Ward
acceptance checks did not execute. There are no combined independent rows
or completed internal-gate status.

Normal and optimized cross-route validators, summary and figure renderers
were never attempted. Therefore no claim of completed cross-route accuracy,
accepted refinement/ledger gates, continuum-tail comparisons or physical
negative-control validation follows from this round.

Static source inspection identifies a plausible memory contributor:
`sha(path)` uses `Path.read_bytes()` and loads the whole archive while the
evolve function still retains its archive dictionary. Archive hashing occurs
after evolve's post-write resource check and before the main-loop check that
failed. The 65,892,627-byte fine archive makes this a substantial transient
allocation. The capture does not profile individual allocations, so the
contribution is a source-based explanation rather than a measured causal
allocation diagnosis. A storage-only followup must demonstrate that its
actual process stays within the same budget.

## Preserve both failures

The replay verifies sources before and after this attempt. The read-only
failure audit additionally rechecks all **331 followup frozen files**,
all **30 captured fresh-output hashes**, both raw archive hashes, and all
**265 original frozen files** in the earlier replay. The original
`57668b8...` calibration remains FAIL from its SciPy continuum roundoff
warning; the `5ff571c...` followup remains FAIL from its independent resource
gate. Their outputs, gates and histories must be published and retained
separately. A subsequent successful registration would not relabel either.

The evidence establishes neither a completed metric response matrix nor
actual shifted-root response, coupled bulk/shell dynamics, quantum stability,
heating, particle yield, observational support or a cosmological discovery.
