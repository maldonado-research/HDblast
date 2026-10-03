You are carrying out one bounded, useful HDBLAST research round. The user authorizes ongoing research and useful narrow GitHub publication, preserving negative evidence. Use the configured model gpt-6.1-sol at ultra effort; do not change the model or reduce effort silently. This invocation is one round, not proof of a persistent background process.

Read AGENTS.md if present, .github/hdblast/queue.json, .github/hdblast/round_state.json, and the JSON context at $HDBLAST_ROUND_CONTEXT. Read the newest relevant checkpoints and their independent review, not just the frozen benchmark. The context identifies exactly one task and a fresh output directory at $HDBLAST_ROUND_DIR. Follow its acceptance requirements. Never recreate an unchanged completed calculation merely to fill a schedule slot. If current source evidence already resolves the queued question, cite and audit that evidence rather than pretending it is new execution.

Before calculation, write REGISTRATION.md in the fresh round directory with pinned input/source hashes, the actual model action and conventions, one scientific question, fixed acceptance/stop rules and planned checks. Preserve registration bytes when producing the result. Use no more than three meaningful evaluations and at most the scheduler's 60-minute model step. Set OPENBLAS_NUM_THREADS, OMP_NUM_THREADS and MKL_NUM_THREADS to 1. Investigate actual failures only when a new diagnosis changes the next attempt. Do not repeat unchanged failing or passing tests.

Perform a genuinely separate review of the derivation or implementation using supported delegation if available within the finite round budget. Retain the check's actual evidence and label it independent internal review, not external peer review. A claim of independent review is insufficient without its reasoning or executed comparison. A BLOCKED or INCONCLUSIVE result can be useful if it identifies a precise mathematical or environmental obstruction from evidence.

Maintain the known scientific boundaries. The original smooth-FRW experiment remains FAIL for curved pressure at its registered K192 comparison. The separate K384 follow-up passed all 13 groups for a prescribed smooth background. Neither result is a coupled shell evolution, a uniform infinite-cutoff theorem, decay, thermal radiation, sustained expansion, a hot Big Bang, or a breakthrough. The existing model's earlier negative evidence is retained. Shared-action stress and scalar sources, all three junctions, finite local terms counted once, constraint/energy ledgers, UV state, and the higher-derivative/EFT treatment must be consistent before any coupled claim. Never insert phenomenological radiation, alter the physical model, relax a tolerance after seeing a failure, or select parameters to manufacture an expected cosmology.

Write only within $HDBLAST_ROUND_DIR. Read frozen checkpoints without changing their code, data, manifests, outcomes or archive. Do not change .github workflows, trusted scheduler helpers, queue/state, main, credentials or deployment settings. Do not push, send messages to others, merge, create a release or mint a DOI. The trusted publication job handles a narrow PR after validation. The scheduler itself does not have validated Zenodo publication credentials.

Retain REGISTRATION.md, RESULTS.md, INDEPENDENT_REVIEW.md, source/code, finite output evidence and ROUND_RESULT.json. ROUND_RESULT.json must have:

    schema_version: 1
    run_id, task_id, source_sha, scheduler_sha, model, reasoning_effort: exactly from context
    status: PASS, FAIL, BLOCKED or INCONCLUSIVE
    historical_failures_preserved: true
    summary: concrete question, result and meaning
    claims: nonempty array of strictly supported claims
    limitations: nonempty array of material limits
    next_step: concrete next calculation or prerequisite
    evaluations: 1 to 3 objects containing purpose, command as an argv array,
      inputs_sha256, actual exit_code, outcome, and evidence_paths
    artifacts: map of EVERY retained relative file path except ROUND_RESULT.json
      to its actual SHA-256 digest

Do not invent execution evidence, hashes, exit codes, API model authorization or publication success. Do not claim symbolic verification proves numerical convergence. A failed scientific gate remains failed even if the packaging check succeeds. The trusted publication gate checks shape, hashes, provenance and scope; your scientific claims still require the independent review. Stop when this round's result is concrete and reviewable.
