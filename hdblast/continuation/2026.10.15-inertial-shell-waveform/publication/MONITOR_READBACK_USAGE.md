# Monitor, readback and additive evidence tools

Preparation only until the real receipts exist. Run these newly authored tools
from `/workspace/hdblast-research-work/inertial-waveform-20261015`. They do not
create PRs, push, merge, fetch or change a checkout. Monitors observe GitHub via
`gh api`; normal readbacks use HTTPS. Local preflight modes make no HTTP request.

## Input schemas

`PREPARED_COMMITS.json` is a two-row list, one for each full repository name:
`maldonado-research/HDblast` and
`maldonado-research/maldonado-research.github.io`. Each row has `repository`,
`worktree`, `branch`, `base`, `head_sha`, `reviewed_tree`, and `pr_number`.
Commit/tree values are full immutable 40-character SHA values. `pr_number`
is the actual positive PR number; it may be omitted during local source
preflight, before a PR exists. Normal monitors/readback require it.

`MERGE_RESULTS.json` is a two-row list with `repository`, `number`, `url`,
`merged=true`, `head_sha`, `merge_commit_sha`, `merged_at`, `reviewed_tree`,
`merged_tree`, and `merged_tree_equals_reviewed_tree=true`.

Other source readback inputs are:

- `frozen-source/SOURCE_FREEZE_RECEIPT.json`: current external manifest/ZIP
  pins and actual source counts from the freeze tool.
- `publication/FIGURE_INTEGRATION_RECEIPT.json`: either two committed asset paths,
  live Pages URLs, hashes and shared source-figure identity, or an explicit
  no-figure receipt (`status=PASS`, `mode=no_figure_assets`, `source=null`,
  `site_assets=[]`, `all_asset_bytes_equal_to_source=true`).
- `PRE_MERGE_CHECKS.json` and `POST_MERGE_WORKFLOW_RUNS.json`: produced by the
  monitors below. These bind their observations to the input receipt hashes.

The readback records all six input-file hashes from the same bytes it parses.
Counts come from the complete committed source and initial-delivery manifests.
Any selected site figure files receive both raw-GitHub and live-Pages readbacks.

## Source/site observation

After reviewed commits and initial delivery are committed, local preflight is:

```sh
python -I -B publication/readback_published.py --preflight-only --output publication/READBACK_HELPER_PREFLIGHT.json
```

After actual PR numbers are recorded, observe their three public workflows
and the shared PR's clean state:

```sh
python -I -B publication/monitor_premerge.py
```

Both monitors accept `--once`, which exits with code 2 when the one observation
is still pending. Default polling is every 20 seconds, bounded by 90 observations.
They reject moved reviewed heads/bases, failed checks and ambiguous duplicate
workflow runs rather than silently choosing one.

After the authorized merges, fill `MERGE_RESULTS.json` and observe all five
workflows, including both deployments:

```sh
python -I -B publication/monitor_postmerge.py
```

Root must obtain each immutable merge commit object in the corresponding local
reviewed worktree before readback. The helpers deliberately perform no Git fetch.
Keep each worktree at its reviewed head with tracked files unchanged.

Only after all five exact-merge workflows pass:

```sh
python -I -B publication/readback_published.py --output FINAL_LIVE_READBACK.json
```

Normal readback checks the three successful reviewed-head premerge workflows,
the five successful exact-merge workflows, deployment identity, actual Git merge
trees, source/ZIP/initial evidence inventories, figure bytes and unique URLs
before starting HTTP. A pending or missing workflow blocks HTTP entirely.

## Stage the final evidence

The complete immutable-commit source audit must be saved as
`publication/COMMITTED_SOURCE_GRAPH_AUDIT.json` using `verify_git_inventory.py`
at the public reviewed head. Its schema is the verifier's normal output.
The fresh replay roster is not inherited from an earlier checkpoint. It must
be finalized in source `REPLAY_ROSTER.json` before freezing. The builder checks
its authenticated counts against the actual normal/optimized process records.
Source, replay, initial evidence and live target counts are dynamic. See
`REPLAY_ROSTER_CONTRACT.md` for the required schema.

When the real merge, CI and readback receipts pass:

```sh
python -I -B publication/build_final_evidence.py
```

The fresh `evidence-staging/` copies every initial-delivery file from the
reviewed immutable Git commit. It adds completed facts, source-graph audit,
readback/preflight, support tool source and a support manifest, a new
`PUBLICATION_COMPLETION.md` addendum, `PUBLICATION_RECEIPT.json`, and a complete
`PUBLICATION_EVIDENCE_MANIFEST.json`. It preserves the initial `README.md` and
all initial bytes. Use repeated `--include task/relative/file` arguments for
additional planning or audit documents; they must be regular files inside the
task. It refuses to replace an initial file with different bytes.

Root reviews this stage, copies only its new files into the dated delivery
directory on a separate evidence branch, and records
`PREPARED_EVIDENCE_COMMIT.json`: `base`, `head_sha`, `reviewed_tree`, `branch`,
`added_files`, `total_evidence_files`, `delivery`, and `manifest_sha256`.
`delivery` must be exactly
`hdblast/continuation/2026.10.15-inertial-shell-waveform`.

Before publishing that evidence branch:

```sh
python -I -B publication/readback_evidence.py --preflight-only --output EVIDENCE_READBACK_PREFLIGHT.json
```

After its PR is merged unchanged, obtain that immutable merge object locally
and provide the actual evidence PR number:

```sh
python -I -B publication/readback_evidence.py --pr ACTUAL_PR_NUMBER
```

This verifies an addition-only diff confined to the dated delivery, the exact
reviewed/merged tree, every committed evidence member and its raw-main bytes.
`FINAL_EVIDENCE_READBACK.json` stays local outside the evidence worktree. It is
not added to the package whose upload it verifies.

## Preparation checks

```sh
python -I -B publication/validate_monitor_readback_tools.py
```

This exercises only synthetic observations, CLI help and compilation. Its
receipt explicitly distinguishes these checks from actual CI, merged PR,
source inventory and website observations. Final stage/readback scripts are
prepared and help/compile checked before real prerequisites; they are not run
against nonexistent publication facts.
