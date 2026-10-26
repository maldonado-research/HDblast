# Portable capsule helpers

Freeze `reproduce.py` at the root of `HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026`.
The manifest schema is `HDBLAST_SCREEN_MAGNETIC_CURVATURE_MANIFEST_V1` and its
`files` list contains every regular file except `MANIFEST.json`, with `path`,
integer `bytes`, and lowercase `sha256` fields. The manifest itself is pinned
externally. Empty directories are immaterial. Symlinks and special members are
rejected; the full control suite uses a filesystem that supports symlinks.

Build a fresh external archive:

```sh
python -I -S -B package_capsule.py --capsule-root CAPSULE \
  --manifest-sha256 EXTERNAL_MANIFEST_SHA --archive NEW_EXTERNAL_ZIP
```

Verify existing bytes without project imports or output writes:

```sh
python -I -S -B package_capsule.py --capsule-root CAPSULE \
  --manifest-sha256 EXTERNAL_MANIFEST_SHA --archive ZIP \
  --verify --archive-sha256 EXTERNAL_ZIP_SHA
python -I -S -B CAPSULE/reproduce.py --repo-root REPO \
  --manifest-sha256 EXTERNAL_MANIFEST_SHA --archive ZIP \
  --archive-sha256 EXTERNAL_ZIP_SHA --verify-only
```

Run the full new-capsule replay into a nonexistent directory whose parent
already exists. Add `-O` to execute the entry point in optimized mode; each
child's normal/optimized mode is controlled separately by the replay plan.

```sh
python -I -S -B CAPSULE/reproduce.py --repo-root REPO \
  --manifest-sha256 EXTERNAL_MANIFEST_SHA --archive ZIP \
  --archive-sha256 EXTERNAL_ZIP_SHA --output-dir FRESH_EXTERNAL_DIRECTORY
```

`REPO/research/` must contain all ten upstream capsules at the fixed pins in
`reproduce.py`. Complete capsule inventories and all ZIP member bytes are
checked before any project import or execution and checked again afterward.
The ZIP is mandatory. No upstream campaign is launched.

## Replay plan

`REPLAY_PLAN.json` is an authenticated capsule member, supplied by the final
inventory owner. Minimal example:

```json
{
  "schema": "HDBLAST_SCREEN_MAGNETIC_REPLAY_PLAN_V1",
  "pins": {"inputs": "numeric/ACTUAL_APEX_INPUTS_NORMAL.json", "map": "numeric/CONSENSUS_SELECTED_MAGNETIC_MAP.json"},
  "runs": [
    {
      "id": "actual-normal",
      "script": "numeric/export_actual_inputs.py",
      "mode": "normal",
      "argv": ["--repo-root", "{repo}", "--output", "{output}/actual-normal.json"],
      "goldens": [{"output": "actual-normal.json", "reference": "numeric/ACTUAL_APEX_INPUTS_NORMAL.json"}]
    }
  ],
  "expected_counts": {"positive_runs": 1, "negative_runs": 0, "golden_matches": 1},
  "program_roles": {},
  "parity_groups": []
}
```

`runs` execute sequentially in declared order. Every run has a unique safe
`id`, a capsule-relative Python `script`, a `mode` of `normal` or `optimized`,
and a string `argv` list. All children receive `-I -S -B`, plus `-O` in optimized
mode. Argument expansion accepts `{output}`, `{repo}`, `{capsule}`. Each entry
`pins.NAME` is a capsule-relative authenticated file and supplies `{NAME}` as
its full path and `{NAME_sha}` as its SHA-256. Thus `{inputs_sha}` can pin fresh
actual output against the frozen input reference without reading stale input.

Every positive run needs at least one `goldens` row. Its `output` is relative
to the external replay directory; `reference` is capsule relative. Optional
`sha256` must equal the reference digest. All bytes must match, including the
final newline. For a program that prints its scientific JSON, add
`"stdout": "direct-normal.json"` and a golden with that output. Otherwise stdout
and stderr are saved separately under `logs/`.

For an expected failure, omit goldens and set `reject_reason` to the exact
terminal exception marker and `absent` to a nonempty list of output paths that
must remain absent. A nonzero exit alone is insufficient: stderr must end in
the marker or `ValueError: MARKER`. The helper also checks absence of dangling
symlinks. Expected output options must resolve inside the external directory.

`parity_groups` lists groups of positive run IDs whose golden outputs, in
declared order, must agree byte for byte. `expected_counts` must contain exactly
`positive_runs`, `negative_runs`, and `golden_matches` and agree with the roster.

Every `.py` member must be accounted for. Direct scripts are inferred from
`runs`; `reproduce.py` is the entry point. All remaining Python members must
appear in `program_roles`, with either:

```json
{"role": "dependency", "description": "Authenticated helper imported by the evaluator.", "used_by": ["numeric/evaluate_selected_magnetic.py"]}
```

or:

```json
{"role": "provenance", "description": "Retained derivation or packaging source; not executed by the scientific replay."}
```

`used_by` entries name direct script paths, not run IDs. The final receipt lists
all programs with their hashes, all runs, expected rejection reasons, exact
golden hash rows, before/after closure receipts, and output inventory. Optional
`scope` and `limits` fields are carried into the receipt.

## Integrity controls

The replay authenticates a bounded fixture built from two actual Python source
payloads and a complete new fixture manifest. The copied sources are never
imported or executed. Fresh isolated authentication subprocesses check both
modes for wrong external pins, changed primary/dependency bytes, undeclared
files, source symlinks, duplicate entries, parent traversal, and absolute paths.
There are 16 reason-specific source rejection cases. Three ZIP controls check
an incorrect external pin, a correctly repinned changed member with valid CRCs,
and a correctly repinned incomplete member roster. None launches a project
calculation. The original capsule, ten upstream roots, archive, and entry bytes
are rechecked after these controls, before the first project run.
