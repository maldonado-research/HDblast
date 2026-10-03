# Post-result design review and next-study draft

This directory is outside the sealed active-source checkpoint. It contains an output-only independent internal review and a draft for a future registration. No new physical calculation, source callback, checkpoint-array decode, or edit of the sealed package was performed.

- `RESULT_DESIGN_REVIEW.md`: supported interpretation, cancellation and numerical scope, publishing GO with the original limitations.
- `SAVED_RESULT_AUDIT.json`: exact-rational audit of the already saved twelve canonical cases and six sealed report/result files, including before/after byte hashes.
- `audit_saved_result.py`: reusable read-only audit, using only Python's standard library.
- `NEXT_REGISTERED_STUDY.md`: concrete fixed-target integral verification, full ledger replacement, independent convergence comparisons, and prerequisites for coupled backreaction.
- `NEXT_STUDY_PROPOSAL.json`: machine-readable draft status and carryforward thresholds; it is not a prospective registration.
- `REVIEW_RECEIPT.json` and `MANIFEST.json`: review scope and exact staged-file inventory.

The sealed result report and audit recorded the fresh ZIP replay as pending. Later replay evidence belongs outside the immutable ZIP. Parent publication may add that evidence without rewriting the sealed checkpoint or implying it was already completed at this review.

To repeat the saved-output audit without executing a physics producer:

```sh
python audit_saved_result.py /path/to/HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER /fresh/path/SAVED_RESULT_AUDIT.json
```

The future-study proposal uses the known completed result to select a useful next question. Its source bytes, resource limits and complete execution schema still require preparation and a separate verified public freeze before any new physical evaluation.
