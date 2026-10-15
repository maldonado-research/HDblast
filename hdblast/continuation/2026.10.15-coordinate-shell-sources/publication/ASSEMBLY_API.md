# Approved source and evidence assembly

All paths below are relative to the new task root unless an explicit worktree
or output is supplied. Do not edit a frozen capsule to repair a failed check.
Accept revised source/configuration, assemble a new directory and freeze again.

## Source roster

`assemble_source.py --config publication/SOURCE_ASSEMBLY_CONFIG.json` accepts:

```json
{
  "schema": "hdblast-approved-source-assembly-v1",
  "status": "APPROVED_FOR_FROZEN_SOURCE",
  "capsule": "HDBLAST_COORDINATE_SHELL_SOURCES_20261015",
  "files": [
    {"path": "README.md", "bytes": 123, "sha256": "complete actual SHA256"},
    {"path": "reproduce.py", "source": "publication/replay-tools/reproduce.py", "bytes": 123, "sha256": "complete actual SHA256"}
  ]
}
```

The displayed sizes/hashes are schema examples, not accepted pins. `source`
defaults to `path`. Every member needs its actual accepted byte count and hash.
Names are canonical relative paths; internal/private/cache/Git paths, symlinks,
reserved manifest names and duplicate output paths are rejected. The approved
roster must include `requirements.txt`, `reproduce.py`, `REPLAY_ROSTER.json` and
`verify_manifest.py`. Map the latter three from `publication/replay-tools/`.
The scientific roster, notes, prior inputs, goldens, audits, citations and any
figures must be explicit. No directory is copied recursively.

The default fresh output is `assembled-source/`; the separate receipt is
`publication/SOURCE_ASSEMBLY_RECEIPT.json`. Override with `--source-root`,
`--output-dir`, or `--receipt`. Source bytes are rechecked after all writes.

## Runner and freeze sequence

1. Accept each checker's complete normal and optimized goldens, dependency
   inventory, named checks and exact intended failure reasons. Record the
   tested exact Python and dependency versions in the approved build-time
   `publication/RUNNER_CONFIG.json` described by `RUNNER_CONFIG_API.md`.
2. Build tools without executing scientific code:

   ```sh
   python -I -B publication/build_replay_runner.py --config publication/RUNNER_CONFIG.json --source-root . --output-dir publication/replay-tools
   ```

   `--replace` is required to replace accepted prefreeze runner outputs.
   The builder authenticates source-input pins and derives actual named counts.
3. Assemble the accepted source roster, then run the separate path/exact-body
   scope audit. It checks only explicit or available private corpora; an empty
   corpus is accurately reported as zero body checks:

   ```sh
   python -I -B publication/verify_publication_scope.py --source assembled-source --private-reference literature/_internal --output publication/PUBLICATION_SCOPE_CHECK.json
   python -I -B publication/freeze_source.py --source assembled-source --output-dir frozen-source --capsule HDBLAST_COORDINATE_SHELL_SOURCES_20261015 --created-utc ACTUAL_FIXED_UTC_TIMESTAMP
   python -I -S -B publication/authenticate_fresh_zip.py
   ```

4. Execute the authenticated fresh extracted runner with the freeze receipt's
   actual manifest SHA256 and a new output directory:

   ```sh
   python -I -B fresh-zip-source/HDBLAST_COORDINATE_SHELL_SOURCES_20261015/reproduce.py --manifest-sha256 ACTUAL_MANIFEST_SHA256 --output-dir fresh-zip-replay
   python -I -B publication/record_replay_tool_validation.py
   ```

   Optional repeated `--review FILE` records actual independent tooling reviews.
   Every complete process result/log/additional output and all source inputs
   are rechecked. Counts are derived from the accepted roster, never inherited.
5. Root supplies an isolated research worktree based on the reviewed current
   remote main. Install the source/archive and dedicated workflow:

   ```sh
   python -I -B publication/install_checkpoint.py --worktree publication/hdblast-worktree
   ```

   This authenticates every frozen source/ZIP member before copying, pins the
   actual runtime, and requires exactly 44 unchanged existing workflow blobs.
   Use `SITE_INTEGRATION_API.md` for the two website summaries/optional assets.

## Initial evidence roster

`install_initial_evidence.py --worktree publication/hdblast-worktree --config
publication/INITIAL_EVIDENCE_CONFIG.json` accepts:

```json
{
  "schema": "hdblast-initial-evidence-config-v1",
  "status": "APPROVED_FOR_INITIAL_EVIDENCE",
  "capsule": "HDBLAST_COORDINATE_SHELL_SOURCES_20261015",
  "readme": "Approved Markdown describing completed local replay and pending remote publication.\n",
  "files": [
    {"source": "frozen-source/SOURCE_FREEZE_RECEIPT.json", "path": "SOURCE_FREEZE_RECEIPT.json", "bytes": 123, "sha256": "complete actual SHA256"}
  ]
}
```

Every configured input needs its actual accepted pin. Required target receipts
are `SOURCE_FREEZE_RECEIPT.json`, `FRESH_ZIP_AUTHENTICATION.json`,
`FIGURE_INTEGRATION_RECEIPT.json` and `REPLAY_TOOL_VALIDATION.json`. Include the
scope, assembly, CI, site, synthetic tooling and independent review receipts
that actually exist. The helper authenticates and copies every replay log,
complete JSON result and additional CSV/output automatically; do not list those
again. `README.md`, `REPRODUCTION_RECEIPT.json`, `CONFIG_APPROVAL.json` and
`INITIAL_REPLAY_MANIFEST.json` are reserved helper outputs.

The new dated delivery directory must be absent. The helper records the actual
complete initial manifest and never claims remote CI, merge or live facts.
After root stages all source, **force-add ignored manifest members explicitly**,
then run `verify_git_inventory.py` against the Git index and again against the
reviewed immutable commit. A successful local ZIP replay cannot establish Git
membership. Continue with `MONITOR_READBACK_USAGE.md` for real publication.
