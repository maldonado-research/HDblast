# Portable metric checkpoint replay

Use Python 3.12 with the checkpoint's `requirements-replay.txt`. This driver
uses only checkpoint-local registered sources and preserved reference copies.
It does not require a parent repository, Git, credentials, or network access.

```
python code/replay_metric.py --output /an/unused/external/directory
```

`--plan-only` verifies the source payload and writes the complete command plan
without launching any subprocess. The optional `--freeze-commit` must match
the receipt; it cannot substitute a different prospective experiment.
The root coordinator independently verifies the public remote before creating
the post-freeze `FREEZE_RECEIPT.json`, with these required fields:

```
{
  "public_freeze_commit": "full 40-hex prospective commit",
  "registration_sha256": "FULL_REGISTRATION.json SHA256",
  "independent_manifest_sha256": "independent/MANIFEST.json SHA256"
}
```

The receipt is added after the public prospective freeze. A final package's
`MANIFEST.json` must include and hash the receipt. The replay authenticates the
receipt hash, full registration bytes, every registered file, independent
manifest and its sources, whole experiment/gates, and any final package
manifest before execution, after each command, and at completion. Every path
is confined to the checkpoint; missing files, symlinks, changed membership and
path escapes are fatal. All scripts are resolved uniquely from registered
filenames, so the standalone package preserves the theory directories and
their sibling reference inputs.

The fresh output contains a complete checkpoint copy, command plan, source
manifest, separate stdout/stderr logs and exit receipts, execution ledger,
fresh analytic and numerical results, and final validation status. No output
is overwritten, and failures remain available. Physical producers run normally
once each. The scientific wall budget remains 900 seconds; the replay wrapper
timeout is 1000 seconds. Every expected process exit is zero, and any failure
halts the replay without changing cutoffs, grids, tolerances or gates.

The 32 commands comprise normal and optimized primary algebra/preflights,
primary action/contact and Ward/trace proofs, independent generic-WKB,
primary-inventory and CSE proofs, actual-root prerequisite algebra, tail
algebra, synthetic validator guards, stable raw baseline proof and faithful
inherited-source verifier; an independent input preflight and derivative
certificate regeneration; two unoptimized physical producers; normal and
optimized saved-mode validation; a saved-data summary and figure renderer.
An optional separately authored `audit_saved_metric.py`, if registered, adds
normal and optimized saved-mode audit commands without a new physical run.

The driver verifies exact Python 3.12 major/minor and NumPy 2.2.6,
SciPy 1.15.3, SymPy 1.14.0, mpmath 1.3.0 and Matplotlib 3.10.1. Numeric threads
are one, bytecode writes are disabled for children, and Matplotlib uses Agg.
No new application service or environment configuration is needed.

The summary retains all 36 finite-K comparisons, three separately integrated
physical baselines, scalar-current normalization, direct-mode refinement,
independently reconstructed conservation evidence, and all 324 omitted-band
comparisons. It distinguishes the conservative directed UV envelopes from
empirical finite integration/refinement/arithmetic estimates. The three
figures each export PNG, SVG and PDF, for exactly nine artifacts, with input
and renderer hashes.

The implementation was checked before public freeze using synthetic source
trees and fabricated saved-input fixtures only. Fourteen replay guards passed
under normal and optimized Python with zero child calculations. Synthetic
postrun checks validate summary/figure plumbing; their values and figures are
not physical response evidence. These development artifacts stay separate
from the actual post-freeze numerical results.

This replay calibrates one homogeneous conformal metric channel at fixed
H=1, x=r=2, xi=0, physical scalar and incoming BD state. The actual-root proof
is algebra only. It executes no actual shifted-root modes, bulk/boundary/state
response, coupled evolution, stability or heating experiment, and it cannot
mark the complete automation prerequisite PASS.
