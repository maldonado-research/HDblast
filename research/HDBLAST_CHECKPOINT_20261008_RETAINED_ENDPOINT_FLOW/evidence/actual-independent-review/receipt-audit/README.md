This audit passed both actual frozen-005 runs and found the nine scientific exports byte-identical. It verifies exported source accounting, attempted-operation coverage, and authorization/custody bindings. It imports no production code and opens no retained archive.

`audit_actual_receipts.py` independently checks the 128 source rows, all 3,200 chosen dyadic coefficients and their vector hashes, the analytic tail identity, exact upward rounding of each complete source error, and all four transported source-error, polynomial-L1, and factorial-tail budgets. It checks the ordered 128-source and 36-array journals, derives 688,152 real slots from the frozen shapes, and recomputes all four canonical InputSpec digests. The 76 authenticated-member receipts match those digests.

The static worker and decoder review establishes the sequencing relied on here: `authenticate_inputs` returns only after all four `verify_inputs` calls succeed; `verify_inputs` hashes all members and parses every header before returning its snapshot object. The worker constructs its source models only after those snapshots exist, then requests the nine selected iterators per capsule. The successful, authenticated execution and journals support that sequence. This audit does not independently repeat archive authentication.

The auditor rehashes the exact registered file and directory rosters, binds the fixed public-GO and registration hashes, checks the actual entry and all ten output pins, reconstructs each complete custodian output-tree hash, and checks log custody, exact command arguments, nonroot execution, wait4 completion, and resource caps. Public-readback evidence is checked from the pinned GO receipt without making a new network request. Historical execution paths are retained and checked for internal command/custodian consistency; those recorded paths are not used to open new files when copied evidence is audited elsewhere.

The original represented values cannot be independently checked without the prohibited duplicate retained-array decode. The physical source coefficients and their coefficient-discrepancy premise are not reconstructed. This audit verifies their exported exact accounting; the underlying physical premise relies on the frozen implementation and reviewed analytic proof. It does not perform a second target evolution. All-node arithmetic and exact represented-state consistency are reviewed separately. Normal/optimized identity is a control diagnostic, not an independent mathematical proof.

Portable invocation, with an existing writable review output directory:

```bash
python -I -B audit_actual_receipts.py \
  --core /path/to/published/core \
  --frozen-source /path/to/frozen/core \
  --registration /path/to/FULL_REGISTRATION.json \
  --go /path/to/PUBLIC_GO.json \
  --normal /path/to/saved/normal \
  --optimized /path/to/saved/optimized \
  --output-dir /path/to/review-output
```

The same authenticated core can be supplied for both source locations when only one copy is available; each location is compared with the complete externally pinned registration. The registration SHA, GO SHA, and public freeze commit remain fixed inside the auditor and are not command-line overrides. No supplied retained-repository path is opened. Omitting `--optimized` produces a normal-only receipt and explicitly leaves paired identity unchecked.

`negative_controls.py` executes only the SHA-pinned independent auditor. Three positive baselines pass, followed by eight expected rejections: omitted rounding, doubled rounding, erased kernel tail, missing source attempt, wrong authenticated-member count, wrong InputSpec digest, wrong entry GO digest, and wrong entry registration digest. Each case uses copied small artifacts and never opens the node stream. All eight pass under normal and optimized Python.

For portable controls, supply `--actual /path/to/saved/normal --core /path/to/published/core --output /path/to/fresh-control-output`. The output directory must not already exist. Run with `python -I -B` and separately with `python -I -B -O`, using distinct output directories.

`RECEIPT_AUDIT_SEAL.json` binds the final auditor, controls, receipts, and this note by relative paths. The control receipts also bind every copied control artifact. Initial script versions and earlier controls remain preserved as development evidence; only the final auditor SHA is accepted by the final control suite. Historical metric calibration remains FAIL, and the continuous time/momentum/contact/UV certificate remains UNRESOLVED.
