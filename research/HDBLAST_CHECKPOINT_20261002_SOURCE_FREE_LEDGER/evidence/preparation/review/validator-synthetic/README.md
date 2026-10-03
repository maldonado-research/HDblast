# Synthetic validator mutation review

Final review: **89/89 cases pass in normal Python and `python -O`**, with identical case decisions and error messages. Thirteen valid or threshold-boundary fixtures are accepted; 76 invalid mutations are rejected. Six aggregate source/mode/scope checks pass.

Validator SHA256: `cd1bdaddd4f530ece599587cc4f487d62312c97113a7c00c0badda36208e14d0`. Both runs pin the validator before and after import/testing; the aggregate receipt also compares the current source. The script SHA and each fabricated input fixture SHA are recorded.

`test_cross_validate.py` imports only the validator's pure functions and fabricates four ordered source/resolution records, two precision levels and three cutoffs, including every full profile field. No saved physical arrays/results, registered physical inputs or physical routes are loaded or executed. No checkout file or remote state is changed. Imports disable bytecode writes inside the checkout.

Coverage includes complete fixtures for all three classifications; exact and overrun arithmetic, scientific, attribution, profile and precision thresholds; shared biased recurrence with an unbiased direct primitive; missing/reordered/duplicate cases, cutoffs, precision levels and profiles; lying full-profile maxima; understated precision gaps; false classifications, witnesses and failure evidence; route identities; separate route/outer time and RSS budgets; execution status, exits, timeout, memory failure, and single-process containment; nondecimal and nonfinite scientific scalars/profile entries.

The tests also reject four tolerance-laundering mutations: a true profile maximum above the science gate reported at the gate, a true direct defect above the gate reported at the gate, actual recurrence differences above their gate reported at the gate, and an actual signed decomposition residual above its gate reported at the gate. Both primary and independent recurrence variants are covered. False reported failure values, thresholds and nonfinite values are rejected.

The review identified failure-value and signed-decomposition integrity gaps, which the root editor repaired before the final run. Initial and intermediate JSON receipts remain as historical evidence. An initial exact-boundary fixture omitted its own precision-gap report and was corrected in the synthetic harness. Inner and outer resource measurements have different scopes; within-budget cross-inequalities are accepted, while any actual route or outer overrun is rejected.

Reproduce outside the checkpoint with the pinned interpreter:

```sh
/workspace/hdblast-cloud-setup/venv-frw/bin/python /workspace/hdblast-research-work/ledger-review-20261002/validator-synthetic/test_cross_validate.py --output /tmp/ledger-validator-synthetic-normal.json
/workspace/hdblast-cloud-setup/venv-frw/bin/python -O /workspace/hdblast-research-work/ledger-review-20261002/validator-synthetic/test_cross_validate.py --output /tmp/ledger-validator-synthetic-optimized.json
```

Receipts: `FINAL_NORMAL.json`, `FINAL_OPTIMIZED.json`, and `SUMMARY.json`. This is a schema and validation-logic review. It does not establish input-data truth, scientific outcomes or successful physical calculations.
