# Bounded reproducibility and literature monitor: inactive candidate

This package proposes a GitHub-hosted mechanical routine: two public literature
metadata queries daily, and full reproduction of one immutable completed
checkpoint weekly or by explicit manual dispatch. It invokes no AI service,
proposes no new physical result and performs no automatic GitHub publication,
issue/comment, merge, release, website update or Zenodo operation. It is not a
24/7 research service. No schedule or workflow has been installed or activated
by this package's preparation.

**Recommendation: retain this concrete candidate until the active-source
checkpoint is completed and reviewed, then perform one hosted manual
validation before enabling the schedule.** `TARGET.json` intentionally contains
no guessed science, freeze, registration or result pins. In this pending state
the active gate returns `ready=false` and launches no metadata or science job.
The files are staged outside the working repository; `DEPLOYMENT_MAP.json`
provides the future public paths for a separate reviewable automation PR.

## What the proposed routine does

- A daily cron at `03:23 UTC` selects metadata watching. Monday invocations in
  UTC also select a full completed-checkpoint replay. The schedule is expressed
  in UTC, rather than the user's changing daylight-saving offset.
- An explicit UTC campaign start/end limits scheduled eligibility to at most
  30 days. Expired or future windows fail closed. The GitHub enable variable
  must be exactly `true`; the default is inactive. Re-running the same GitHub
  invocation is rejected; a deliberate new manual dispatch is distinct.
- Both fixed INSPIRE queries read their first page once: at most 20 hits each,
  4 MiB per response, 35-second socket timeout and no redirects or retries.
  TLS validation and normal proxy handling remain enabled. No API key, token,
  authentication header, abstract/full-text download or remote write is used.
- Records are normalized into a small bibliographic schema and compared with
  the hash-pinned reviewed `SOURCE_BASELINE.json`. Added, modified and absent
  first-page records are retained separately. Absence is not a retraction;
  metadata changes are candidates for later reading, not evidence for HDBLAST.
- Replay uses a separate checkout at the full completed science Git SHA. Every
  frozen file, registration, driver and requirements hash is authenticated
  before dependency installation or checkpoint-code execution. Unregistered
  executable/module files and symlinks are rejected. Git ancestry and the
  checkpoint's own freeze verifier are then checked.
- The existing complete active-source driver runs once with Python3.12.14,
  NumPy2.2.6, SciPy1.15.3, SymPy1.14.0, mpmath1.3.0 and Matplotlib3.10.1;
  all numeric thread variables are one. Its nonroot Linux isolation, original
  per-entire-route `900 s /262144 KiB` budgets, full matrices and both precision
  contexts remain unchanged. The monitor neither creates smaller substitute
  matrices nor resets budgets per case.
- All registered scientific outcomes, including `CONSISTENCY_FAILURE`, remain
  possible completed outcomes. A reproduction is accepted only when both
  routes, both validators, source checks, recorded classification and pinned
  retained scientific bytes match the reviewed completed reference. Earlier metric
  `FAIL` stays explicit. Resource measurements are retained without requiring
  runtime/RSS to equal earlier observations.
- Complete logs, progress, outputs, hashes/diffs and failure receipts become
  artifacts retained for seven days. Job summaries contain counts and fixed
  status labels; untrusted titles/URLs are not rendered into summary Markdown.
  The package creates no issues, comments or repeated public failure messages.

## Prerequisites and deployment sequence

1. Complete the separately registered active-source experiment and fresh replay
   first. Retain every original outcome, even a negative outcome. Record final
   scientific files and resources; do not rerun uncomputed physical science to
   populate this monitor's target.
2. Populate `TARGET.json` from that completed evidence, changing its state to
   `COMPLETED_REVIEWED_CHECKPOINT`. Bind the full science commit, ancestor freeze
   commit, full-registration SHA256, requirements and driver SHA256, original
   classification and every relevant deterministic scientific output hash.
   Both `primary/diagnostic.json` and `independent/diagnostic.json` are mandatory;
   include all deterministic flow files used by the scientific result. Exclude
   timestamps, elapsed times, memory measurements and runtime logs from the
   expected scientific-hash map. For those two diagnostics the expected hash is
   `scientific_projection.project(parsed_diagnostic, route)[0]`, under the exact
   `HDBLAST_ACTIVE_SOURCE_SCIENTIFIC_PROJECTION_V1` policy bound in `TARGET.json`.
   Other declared scientific `.json`/`.jsonl` files use raw SHA256. The hash keys
   are relative to `fresh/` in the outer replay directory. Never insert a raw
   whole-diagnostic SHA into a projection-hash slot.
3. Review these controller assets separately and merge them onto the public
   default branch through the repository's existing author-review workflow.
   The completed science can remain at its reachable immutable research
   commit; the monitor's scientific checkout does not assume default `main`
   contains that research. GitHub cron requires the workflow on the default
   branch.
4. Configure the three **repository variables**, not secrets:
   `HDBLAST_FOLLOWUP_ENABLED=false`,
   `HDBLAST_FOLLOWUP_START_UTC=<actual start, YYYY-MM-DDTHH:MM:SSZ>` and
   `HDBLAST_FOLLOWUP_UNTIL_UTC=<actual end, no more than30 days later>`.
   Timestamp values are validated data, never interpolated into shell code.
   Check GitHub Actions is available and enabled for the public repository.
5. Dispatch `guards`; then dispatch `all` once on the default branch while
   scheduled enable remains `false`. Manual dispatch selects the requested
   bounded jobs independently of the schedule flag, but still requires the
   reviewed completed target, default branch and live campaign window.
   Inspect all job outcomes and artifact contents. A green upload step or job
   metadata alone does not verify numerical rows.
6. Enable `HDBLAST_FOLLOWUP_ENABLED=true` only after the complete hosted manual
   validation succeeds and its artifact is inspectable. If hosted artifact
   downloads remain inaccessible, retain the candidate pending verification;
   don't infer content from a green workflow. Disable or let the window expire
   after a replay mismatch, persistent API failure or the campaign's useful
   purpose is complete. Review changed sources or baseline bytes before a new
   campaign; the watcher never updates its own baseline or target.

The workflow requests only `contents: read`; there are no user-configured
secret references or write permissions. GitHub's checkout/artifact actions
still use the platform's own ephemeral runner authentication internally.
`persist-credentials: false` keeps checkout credentials out of local Git
configuration. No ChatGPT login, OpenAI key or Zenodo secret is requested,
copied or exported. Third-party actions are pinned to full commits already
used in the reviewed inactive scheduler package.

## Completeness, data and schema risks

The versioned diagnostic projection authenticates exact top-level schemas and
the resource/provenance objects that contain its exclusions. It keeps all
scientific fields, decimal strings, profile samples, list order, classification,
producer/freeze/registration/input pins, fixed resource limits and scope in a
canonical JSON hash. It excludes only `resources.seconds` for the primary route,
`resources.elapsed_seconds` for the independent route, both
`resources.peak_rss_kib` values, and independent
`provenance.full_frozen_verification.package_manifest_sha256`. Elapsed time/RSS
must remain finite, positive and within the fixed complete-route limits. The
package receipt must be null or a SHA256; actual ZIP membership/provenance is
verified by the frozen driver and validators. The original raw diagnostic hash,
all excluded observations and all raw files remain inspectable in artifacts.
The package receipt can differ between an unpackaged original run and a sealed
ZIP replay; it does not alter the frozen scientific inputs. No recursive removal
of keys, numerical tolerance or decimal-string normalization is used. The
complete frozen normal/optimized validators remain responsible for the detailed
scientific schemas. A new unknown retained scientific field changes the hash;
unknown fields in top-level or exclusion-containing objects are rejected.

This corrects the earlier staged candidate's raw-diagnostic hash mismatch, which
incorrectly treated elapsed-time differences as scientific differences.

The checkpoint must contain all four authenticated nineteen-member input
capsules, complete registered sources and final schema/validators at the
pinned commit. The proposed sparse checkout retains the entire checkpoint
directory and full Git history needed for freeze verification. It does not
substitute cached arrays, last-run results, release assets or a mutable latest
artifact. Neither dependency data nor the science/bibliography baseline uses
Actions caches. Package installation still depends on availability of the
exact pinned wheels; successful setup on this current machine does not prove
the future runner can restore them.

This candidate was designed against the currently staged active driver
interface (`FULL_REGISTRATION.json`, `REPLAY.json`, four complete commands and
`fresh/{primary,independent}`). The active scientific schema is still being
reviewed. If the completed interface changes, adjust this controller before
activation and run a new hosted manual validation. Never reinterpret missing
files, changed schemas or classifications as a pass. INSPIRE may return
unexpected schemas or rate limits; one failed query retains a failure receipt
and causes the watch job to fail without a retry. The two fixed searches and
first-page window do not cover all relevant publications.

GitHub may delay or omit cron runs, and inactivity may disable public-repo
schedules. A 30-day window gives at most30 nominal daily invocations and five
nominal Monday replay invocations, aside from deliberately dispatched manual
runs. Job ceilings are3 minutes for the gate,5 for metadata,75 for replay and5
for offline guards; these are ceilings, not measured usage or a monetary cap.
No guaranteed recurrence, new discovery, paid-model autonomy or 24/7 uptime is
claimed. Real research progression remains a separately reviewed queue of
newly justified registered experiments and primary-source reading.

## Verification performed during preparation

The offline guard suite uses fabricated files and responses only. It tests
pending/guessed targets, branch/repository/retry/campaign restrictions, negative
scientific outcomes, missing/mutated outputs, pre-execution source mutation and
module shadowing, symlinks, malformed metadata, redirects, limits, exactly two
requests after failures, and summary injection. It performs no network request,
physical evaluation or publication. The expanded suite also accepts only the
declared resource/package variations, retains their raw hashes, and rejects
changed scientific profiles, resource limits/scope, producer/freeze receipts,
unknown exclusion fields, malformed membership and changed comparison policy.
Run normal and optimized variants:

```bash
PYTHONDONTWRITEBYTECODE=1 python test_monitor.py
PYTHONDONTWRITEBYTECODE=1 python -O test_monitor.py
```

Consult `PREPARATION_VALIDATION.json` for actual local counts, syntax/actionlint
checks and the explicit unperformed hosted/deployment tests. Pending pins are
an intentional readiness condition, not a request for a new credential.
