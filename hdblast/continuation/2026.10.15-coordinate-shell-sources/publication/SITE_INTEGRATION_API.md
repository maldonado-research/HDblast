# Coordinate-source summary integration API

`integrate_coordinate_summaries.py` is prepared tooling. It contains no final
scientific numbers, inherited publication receipts or fixed repository commits.
The root publication agent must supply the reviewed result and actual baseline
pins. Preparing this helper did not integrate or publish any coordinate-source
result.

## Invocation and write scope

Run from this task with the SHA-256 of the exact reviewed config bytes:

```sh
python -I -B publication/integrate_coordinate_summaries.py --config publication/SITE_INTEGRATION_CONFIG.json --approved-config-sha256 ACTUAL_CONFIG_SHA256
```

The default invocation is read-only and prints `PASS_READ_ONLY_PLAN` only after
all local checks pass. It writes no receipts or summaries. Add `--apply` to
write the nine fixed summaries; add `--review-only` to require their exact
expected integrated bytes without writing summaries. These modes are mutually
exclusive. Both emit these four local publication receipts:

- `SITE_INTEGRATION_REVIEW.json`
- `SUMMARY_INTEGRATION_PREFLIGHT.json` (identical content)
- `FIGURE_INTEGRATION_RECEIPT.json`
- `SITE_BASELINES.json` (the actual configured immutable baseline commits)

The helper does not install an asset, copy a capsule, create a commit, change a
branch, fetch, push, create or merge a PR, run CI, deploy, or request a URL.
All target files are preflighted before the first summary write. Writes occur
one file at a time; this is not a multi-file atomic transaction. A rerun accepts
only exact baseline or exact planned bytes, so an interrupted application can
be completed without overwriting an independent edit.

## Config contract

The JSON must have exactly these top-level keys. Unknown keys and duplicate JSON
object keys are rejected. The following is a shape illustration; uppercase
placeholders must be replaced with actual reviewed facts and pins.

```json
{
  "schema_version": 1,
  "checkpoint": "coordinate-shell-sources-20261015",
  "capsule": "HDBLAST_COORDINATE_SHELL_SOURCES_20261015",
  "delivery": "hdblast/continuation/2026.10.15-coordinate-shell-sources",
  "approval": {
    "approved": true,
    "approved_by": "ACTUAL_REVIEWER_OR_ROOT_AGENT",
    "basis": "ACTUAL_REVIEWED_RESULT_AND_WORDING_BASIS"
  },
  "scientific_scope": {
    "boundary_trace_and_proper_time_derivative": true,
    "one_copy_coordinate_scalar_and_isotropic_israel_sources": true,
    "scalar_degrees": [0, 2, 4],
    "full_u2_quantitatively_unresolved": true
  },
  "repositories": {
    "hdblast-worktree": {
      "repository": "maldonado-research/HDblast",
      "worktree": "/workspace/hdblast-research-work/coordinate-shell-sources-20261015/publication/hdblast-worktree",
      "baseline_commit": "ACTUAL_40_HEX_BASE_COMMIT"
    },
    "shared-worktree": {
      "repository": "maldonado-research/maldonado-research.github.io",
      "worktree": "/workspace/hdblast-research-work/coordinate-shell-sources-20261015/publication/shared-worktree",
      "baseline_commit": "ACTUAL_40_HEX_BASE_COMMIT"
    }
  },
  "result_note": {"path": "ACTUAL_RESULT_NOTE.md", "sha256": "ACTUAL_64_HEX_SHA256"},
  "result_inputs": [],
  "words": {"en": {}, "es": {}},
  "figure": null,
  "summaries": []
}
```

Approval is a root-reviewed config assertion, not an additional user permission
flow or evidence that scientific checking was executed. The explicit CLI digest
binds that assertion to the exact config bytes. A missing config, a differing
digest, `approved=false`, missing approval attribution/basis, or an unsupported
scientific scope blocks every mode.

Paths in `result_note` and `result_inputs` are canonical task-relative paths.
Every member has exactly `path` and `sha256`; duplicates are rejected. The note
must be Markdown. Supply the note and additional actual source/result inputs
that justify the approved wording. Their bytes and hashes are checked and
reported. The note and any figure must also exist at the same relative paths
in the published source capsule; capsule membership is checked by the separate
source-graph/publication tooling, not this helper.

`words.en` must have exactly `title`, `result`, `limits`, `description`,
`proof_label`, `checks_label`, `evidence_label`, `card_result`, `card_limits`, and
`current_summary`. `words.es` must have the first seven keys only. When a figure
is selected, both also require `alt` and `caption`. Values are nonempty, trimmed,
plain single-paragraph text; numeric comparisons are allowed. HTML tags,
newlines, the new checkpoint marker, and Markdown delimiters in labels/title
are rejected. Titles should omit the trailing colon added by the renderer.

The result wording should describe the approved continuous metric boundary trace
and its proper-time derivative, and the complete one-copy coordinate scalar and
isotropic Israel sources, including the approved Sp/SK degree-0/2/4 projections,
norms and sign certificates. Only actual reviewed numerical results belong in
these fields. The helper does not infer coefficients or validate the mathematics.
Its fixed scope sentence is appended to both body-limit paragraphs and the
shared current-version paragraph: “Full U2 remains quantitatively unresolved.”
The Spanish body uses its Spanish equivalent. Avoid duplicating that sentence
in the configured limits. The approved descriptions must retain an accurate
current scope as well.

Every new entry links to the configured note, the new capsule `README.md`, and
the new dated `PUBLICATION_RECEIPT.json`. The rendered entries direct readers to separate publication evidence for
replay, CI and live-readback status, without a public sentence that becomes
stale after publication. The internal integration receipt remains marked
PENDING_FINAL_ADDITIVE_RECEIPT until the final additive evidence stage.

## Nine summary rows and metadata pins

Supply exactly one row per path below. Each row has exactly `repository`
(the slot name), `path`, `baseline_sha256`, `anchor`, and `metadata`.
`baseline_sha256` is the hash of the actual immutable `git show BASE:PATH`
bytes, including their original line endings. Do not derive it from an edited
working file.

| Repository slot | Path | Kind | Language |
| --- | --- | --- | --- |
| `hdblast-worktree` | `README.md` | Markdown | English |
| `hdblast-worktree` | `hdblast/README.md` | Markdown | English |
| `hdblast-worktree` | `hdblast/CURRENT_STATUS.md` | Markdown | English |
| `hdblast-worktree` | `hdblast/PROGRAM_MAP.md` | Markdown | English |
| `hdblast-worktree` | `site/index.html` | Detailed HTML | English |
| `hdblast-worktree` | `site/es/index.html` | Detailed HTML | Spanish |
| `shared-worktree` | `README.md` | Markdown | English |
| `shared-worktree` | `index.html` | HDBLAST project card | English |
| `shared-worktree` | `projects/hdblast/index.html` | Detailed HTML | English |

The literal insertion anchor must occur exactly once in the pinned baseline.
Use its complete historical marker/opening tag. Existing baseline patterns are:

```text
<!-- quadratic-shell-forcing-20261015 -->
<div data-hdb-quadratic-forcing="quadratic-shell-forcing-20261015">
<p data-hdb-quadratic-forcing="quadratic-shell-forcing-20261015">
```

These are historical insertion anchors, not claims carried into the new result.
The tool also permits another exact prior HDBLAST checkpoint marker/opening
with the same literal shape when the actual baseline requires it. Markdown uses
the comment, detailed HTML uses a `div`, and the shared card uses a `p`. An
existing new marker or an ambiguous/missing old anchor is rejected.

Markdown and card rows have `metadata: null`. Each detailed HTML row has exactly:

```json
{
  "old_description": "COPY_EXACT_BASELINE_DESCRIPTION",
  "old_version": "COPY_EXACT_BASELINE_JSONLD_VERSION",
  "current_line": null,
  "old_is_based_on": "COPY_EXACT_BASELINE_ISBASEDON_IF_UPDATING"
}
```

`old_description` must match all three actual description meta fields and the
single described Dataset/CreativeWork JSON-LD node. The three meta fields and
that node's description are replaced with the approved language's `description`.
The node's `version` changes to
`2026.10.15-coordinate-shell-sources (GitHub checkpoint)`.

For the shared detailed page only, `current_line` must be the exact full existing
HDBLAST `<p class="version-line" ...>...</p>` containing `old_version`. Its
opening/closing tags and attributes are preserved; its text becomes the new
version, `words.en.current_summary`, and the fixed unresolved-U2 sentence. The
two dedicated-site rows require `current_line: null`. Other dated version lines
remain unchanged.

Set `old_is_based_on` to the exact baseline JSON-LD string to advance the current
source pointer to the new capsule tree URL. Set it to JSON `null` to preserve
that pointer. No other JSON-LD field, historical `subjectOf` member or modified
date is rewritten. Current metadata pins are intentionally explicit; the tool
does not guess an old version or copy prior scientific facts.

## Optional figure and readback compatibility

For no figure, use `figure: null`. The emitted figure receipt is exactly:

```json
{
  "status": "PASS",
  "mode": "no_figure_assets",
  "source": null,
  "site_assets": [],
  "all_asset_bytes_equal_to_source": true
}
```

For a figure, use exactly:

```json
{
  "path": "computed/ACTUAL_FIGURE.png",
  "sha256": "ACTUAL_64_HEX_SHA256",
  "asset_name": "coordinate_shell_sources_20261015.png"
}
```

The asset filename must match `[a-z0-9_-]+\.png`. Root must first install the
exact authoritative PNG bytes in both worktrees at `site/assets/ASSET_NAME` and
`projects/hdblast/assets/ASSET_NAME`. Every mode requires those bytes to be
present and identical to the pinned task source. The helper renders a figure
only on the three detailed pages; the Spanish dedicated page uses
`../assets/ASSET_NAME`.

The selected-figure receipt has `status: PASS`, `mode: figure_assets`,
`all_asset_bytes_equal_to_source: true`, a `source` object with capsule-relative
`path`, `bytes`, and `sha256`, and exactly two `site_assets` objects. Each asset
has full `repository`, repository-relative `path`, public Pages `url`, `bytes`,
and `sha256`. These fields match the general published-readback contract.
They are local byte observations, not HTTP or deployment observations.

## Preservation and receipt meaning

Each worktree must be the exact declared isolated slot, have the pinned full
commit object locally, and have that baseline as an ancestor of its actual
HEAD. The helper reports the actual HEAD and baseline tree. It compares tracked
existing files against that baseline and rejects any change outside the nine
summary paths. Additional new source/evidence files are outside this helper's
inventory scope and are explicitly not authenticated here.

Each summary must be a regular, nonsymlinked UTF-8 file with its baseline
executable mode and an unconflicted, matching-mode index entry. An independently
edited summary is rejected. Every rendered change must reverse to all original
baseline bytes; the receipt includes baseline, observed, expected-updated and
reversed SHA-256 values. `updated_sha256` is populated only when the observed
file actually equals the expected integration. Plan mode can therefore describe
expected bytes without asserting that a write occurred.

Prior checkpoint text is preserved. Detailed-page changes are confined to the
new additive block, declared description/current-version metadata, optional
current source pointer, and the explicitly pinned shared current paragraph.
Navigation and every historical HTML anchor-link occurrence must remain.
The shared home-page insertion must be within the unique HDBLAST project card;
all its original bytes and every other original project card are retained.

The review receipt binds its checks to the config hash, helper hash and result
input hashes and records the actual baseline pins. It explicitly reports that
no network request, Git mutation, CI, PR merge, deployment, or live readback was
performed. Approval in the receipt remains the supplied assertion. A local
integration `PASS` is not a completed-publication receipt or scientific proof.

## Preparation checks actually performed

The helper's CLI help executed successfully; missing required arguments and a
missing config were rejected. Python AST parsing passed without writing a
bytecode cache. Read-only, in-memory synthetic rendering against all nine prior
summary templates passed exact reversal both without a figure and with an
optional figure block. The shared-card check retained all six other project
cards. Descriptions containing quotes and ampersands round-tripped correctly.

Synthetic negative checks rejected duplicate/missing insertion anchors, an
existing new checkpoint marker, duplicate JSON config keys, a mismatched
approved-config digest, and `approved=false`. The no-figure receipt matched the
readback schema exactly. No scientific result, actual approved integration,
figure installation, capsule installation, public worktree mutation, CI, merge,
deployment, or live readback was executed as part of these preparation checks.
