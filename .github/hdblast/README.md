# Bounded recurring HDBLAST research — prepared, inactive

No recurring research job was created or activated. No cron daemon was started in this ephemeral environment. The files here are a reviewable scheduler package, not evidence of 24/7 execution.

The current machine has Codex CLI 0.159.0-alpha.3 and reports a ChatGPT login. It can support a current-session round through that existing host. Its login does not establish authentication on a new GitHub runner, and no login files or credentials have been read, copied or exported. OPENAI_API_KEY and CODEX_API_KEY were absent at inspection. Existing scoped GitHub authentication supports repository reads and the current authorized publication workflow, but the Actions settings, secret-name and variable-name endpoints returned HTTP 403 for both repositories. Those settings and bindings remain unknown.

## What is prepared

- `NATIVE_AUTOMATION_PROMPT.md`: a paste-ready bounded research prompt for a supported native scheduler. This session has no callable scheduler tool; the prompt alone schedules nothing.
- `hdblast_research_round.yml`: a GitHub-hosted schedule every six hours, guarded by explicit configuration; an independently invoked exact-model smoke check; one concurrent round; separate model and publication jobs; no automatic main merge, release or DOI.
- `queue.json`: de Sitter sources, stationary closure and fixed-geometry variance calibration are externally completed in `SCIENTIFIC_PROGRESS.md`. Next publicly register matched fixed-geometry linear density/pressure/current calibration from its analytic proposal. Coupled initial-data/evolution still requires metric, shifted-root, bulk/boundary and quantum-state inputs. The automated ledger remains empty.
- `round_control.py` and `runtime_gates.py`: trusted selection, scope/provenance/hash gates, negative-result preservation, idempotence, finite campaign/day guards and duplicate-task suppression. Offline tests exercise the control plane, not scientific correctness.

The model job gets repository read permission. The publication job runs on a fresh runner with the trusted scheduler checkout and obtains write/PR permissions only after the model job completes. The publication gate reconstructs its expected context independently, checks hashes and limits, and copies only a fresh checkpoint plus the generated round ledger. It does not execute submitted research code. An independent internal review and actual recorded checks are required, but the structural gate cannot establish mathematical truth or peer review.

## Verified model and runner configuration

The [pinned official Codex 0.160.0 catalog](https://github.com/openai/codex/blob/a956835d020762cb2b570053af06f643a11c0ecc/codex-rs/models-manager/models.json) lists `gpt-6.1-sol`, reasoning level `ultra`, API support, and minimum client version 0.153.0. The [official Action](https://github.com/openai/codex-action/tree/86365089eb2b84e0a8fb0717b304f8bdcb13b20e) accepts separate `model` and `effort` inputs. The template pins the Action to v1.12's immutable commit and the CLI/proxy to published npm version 0.160.0. This verifies a supported request configuration; no fresh-runner account/model request has succeeded yet.

The GitHub alternative requires an OpenAI API account/key with access to that exact model and effort, API billing/usage availability, GitHub-hosted runner availability/minutes, and repository settings permitting the required workflow and PR operations. ChatGPT subscription authentication in this cloud machine is not substituted for the runner's API credential. Enter API credentials securely in repository Actions settings; never put them in chat, a committed file, a CLI auth export or the workflow body.

## Deployment and activation requirements

`DEPLOYMENT_MAP.json` maps these staged files into a separate automation PR based on the public default branch. Scheduler controls must be on the default branch for GitHub cron. A separate exact `HDBLAST_SCIENCE_SHA` checkout provides the reviewed scientific baseline, so activation need not merge the entire historical scientific PR chain. The source must be an existing pushed 40-character commit containing the completed causal-response checkpoint and inherited stationary, de Sitter, junction and prior baselines. Recommended reviewed source: `e784b18128012825f126c5bd96a4bf2033ab63d9` (`research/20261002-causal-response`). Earlier source and registration-only commits cannot satisfy this refreshed queue. Preflight checks all required source digests before paid work. Main-branch round outputs are copied as additional sources, and source/scheduler commits are both retained in provenance.

Required secure binding: `OPENAI_API_KEY`, entered in GitHub Actions settings. Its name/presence is currently unverified because the scoped GitHub API denies secret listing. No new secret binding was declared or activated here.

Required repository variables:

| Name | Purpose / fail-closed setting |
| --- | --- |
| `HDBLAST_CODEX_MODEL` | Exactly `gpt-6.1-sol` |
| `HDBLAST_CODEX_REASONING_EFFORT` | Exactly `ultra` |
| `HDBLAST_SCIENCE_SHA` | Exact pushed commit for the reviewed science input |
| `HDBLAST_SMOKE_ENABLED` | Explicitly enable the bounded account smoke invocation |
| `HDBLAST_MODEL_SMOKE_PASSED` | `true` only after the exact-model smoke invocation succeeds |
| `HDBLAST_RESEARCH_ENABLED` | Leave absent/`false` until all activation prerequisites are verified |
| `HDBLAST_CAMPAIGN_START_UTC` | An explicit timezone-aware ISO UTC campaign start |
| `HDBLAST_RESEARCH_UNTIL_UTC` | Explicit UTC expiry, at most seven days after start |
| `HDBLAST_MAX_PAID_ROUNDS_PER_DAY` | Integer 1–4; every workflow invocation conservatively consumes allowance |
| `HDBLAST_MAX_CAMPAIGN_INVOCATIONS` | Integer 1–12; finite campaign bound including smoke/failed/skipped runs |

With credentials/configuration supplied securely and the disabled workflow merged onto the default branch, a bounded manual `smoke` dispatch can test the exact provider/model request. Enable research only after the account test succeeds and GitHub permissions/runner availability are verified. No silent model or effort fallback is permitted. Do not mark the smoke flag from a catalog lookup alone.

Each research round has a 60-minute model-step timeout, a 75-minute research-job timeout, at most three documented meaningful evaluations, at most two concurrent agent threads per session, and a 20 MiB publication artifact bound. The official Action has no token or dollar-cap input; elapsed time and invocation limits are real bounds, not an exact spending guarantee. Configure provider-side account controls separately. Re-runs of a prior workflow are blocked from paid work. All attempts must leave honest invocation/failure accounting; an aborted model call is not a completed scientific result.

An open automation PR pauses new paid research. A closed matching task/source PR also prevents rerunning unchanged work. Passing merged ledger records unlock the next queued task; FAIL/BLOCKED/INCONCLUSIVE records stop dependent work. A failed/cancelled/timed-out prior research invocation also pauses paid research for the same campaign, even when no scientific PR could publish. Correct and diagnose the failure before explicitly starting a new bounded campaign. Further work after a scientific outcome needs a new evidence-based task/diagnosis, not a forced retry. The schedule can therefore pause well before expiry. GitHub can delay or drop scheduled invocations and may disable inactive public schedules. This is recurring bounded cloud execution when activated, not an uninterrupted runtime or a guaranteed breakthrough.

## Native scheduling evidence and limitation

After the network allowlist update, the [official app introduction](https://openai.com/index/introducing-the-codex-app/) was reachable. It confirms scheduled Automations with a review queue; its future-looking section discusses adding cloud triggers so runs do not depend on an open computer. The [current ChatGPT-plan help page](https://help.openai.com/en/articles/11369540-using-codex-with-your-chatgpt-plan) mentions Scheduled Tasks and distinguishes local device workflows from OpenAI-managed cloud tasks. The current developer Automations/Scheduled Tasks pages still returned proxy HTTP 403. Current cloud-vs-local/app-awake prerequisites could not be established from the reachable sources. Do not promise that a local desktop Automation will continue while the user's computer sleeps. A supported native cloud schedule would retain plan authentication on its supported host; it must not be constructed by exporting this environment's auth files.

## Validation and scientific boundaries

Run `PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s . -p 'test_*.py' -v` from this directory. The 36 offline tests passed using fixtures and make no model, GitHub or Zenodo calls. Both staged workflows passed checksum-verified actionlint 1.7.12; shell block and embedded Python syntax are checked separately. No paid model invocation, remote schedule, fresh-runner execution, automated PR creation or DOI publication is claimed as tested.

The original smooth-FRW curved cutoff outcome remains FAIL. The separately registered K384 follow-up passed 13 groups for a prescribed background. The reviewed action/junction specification and static de Sitter bridge remain prerequisites, with any corrected derivation and open blockers preserved. The de Sitter checkpoint now supplies the common-action vacuum source and resolved classical susceptibility. The stationary checkpoint now declares benchmark normalization, reference scale, quadratic mass law and source amplitudes and solves the paired stationary equations. These are research choices, not empirical measurements. Its original 48-TypeError validator FAIL and separately repaired 1584-check PASS remain retained. Only positive gamma=.01 passes the declared sampled tenfold gravitational screen; the resolved nonlinear gamma=1e6 case fails it. The publicly registered fixed-geometry variance calibration now covers twelve response points, thirty-six finite-K comparisons, seventeen rejected controls and seventy-two independent archive reconstructions. Only its omitted UV-tail component has the stated analytic enclosure; quadrature and integration errors remain estimates. The new matched-stress proposal is analytic preparation only: direct density and pressure, complete contacts, stress/derivative error bounds and numerical gates must be publicly frozen before stress calibration. Metric/bulk/boundary/actual-root/state inputs remain missing for coupled evolution or stability. None of these results establishes a coupled hot epoch, thermalization, uniform infinite-cutoff theorem, or a breakthrough.

Zenodo remains outside this scheduler: publication credentials have not been validated, and no DOI is created by these scripts.
