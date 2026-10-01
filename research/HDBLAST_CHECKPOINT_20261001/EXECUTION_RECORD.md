# Actual execution record

## Before publication

C1 and C2 were exploratory before registration; C3's exact-pulse cases were registered before integration.

The root independently executed code/cohort_verify.js and code/duration_verify.js in isolated V8 on 1 October 2026. Both returned PASS. The C1 source integrates Gaussian moments with Simpson quadrature from 1024 to 8192 bins, tests five positive decay histories and two assumption-violation controls, reproduces four archived endpoint budgets, and checks the exploratory coupling corridor. The C2 source independently integrates the quadratic-density continuation with RK4 at 500, 1000 and 2000 steps and reproduces saved-point fractions and endpoint targets.

The unchanged C3 source supplied by the implementing agent was independently replayed by root at 2026-10-01T21:39:13.488Z. For execution inside V8 only, the module export keyword was stripped and the exported function was called; no integration, parameter or acceptance logic changed. The saved .mjs source preserves the supplied executable module. It returned PASS on all 90 occupation controls, 30 finest Wronskian checks and the independent coefficient round-trip.

Root C3 replay exactly reproduced the reported maximum occupation error 2.438444030028464e-10, Wronskian drift 1.5013545962005992e-10 and domain change 2.2970839119729192e-10. The full root replay output, including all rows, is saved.

## GitHub runner

The read-only GitHub Actions workflow completed successfully: run [36930743823](https://github.com/maldonado-research/HDblast-archive/actions/runs/36930743823), source commit 0594f3159031bce095d646668df8c6adff47e5a0, job 110599088143. Checkout, Node setup, all three script replays, the explicit registered-mode pass assertion and output-artifact upload completed successfully. The fetched run concluded at 2026-10-01T21:45:35Z. This runner reads the checkpoint only; it does not modify source or publish the website. Its private artifact is hdblast-exact-controls. Raw-trajectory status remains separately recorded in TRAJECTORY_AUDIT.md.

No new full five-dimensional evolution, Python field-equation simulation or self-consistent renormalized quantum calculation was performed by the V8 tests.

The separate raw-trajectory runner 36931967040 completed successfully with nine requested B3 histories and zero load/audit errors. The full stdout report was fetched from job logs and saved under outputs/trajectory_audit.json. Original arrays were read with pickle disabled. This is a data audit, not a new field-equation evolution.

## Curated input verification before final replay

Root fetched the eleven required scientific NPZ binaries and their matching already-public summaries. Their SHA-256 values match the original trajectory report's recorded inputs. The final read-only runner also replays this bundled data and verifies the curated ZIP's CRCs, inventory and per-file SHA-256 values. Its final check status is available from the pull request; the earlier completed run IDs above remain pinned evidence and are not reassigned to later commits.
