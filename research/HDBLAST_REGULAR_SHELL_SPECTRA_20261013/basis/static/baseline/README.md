# Preserved 8 October theorem and review provenance

These files preserve the accepted baseline rather than replacing it with the
11 October improvement. The immutable original checkpoint is
[`30c2ef623a1b2a322ebf245402433eca57ce24f9`](https://github.com/maldonado-research/HDblast/tree/30c2ef623a1b2a322ebf245402433eca57ce24f9/research/HDBLAST_SCALAR_BOUNDARY_CURVATURE_20261008).

The older theorem and exact calculation retain historical labels such as
“pending independent review” and “pending at computation.” Those labels describe
the state when the documents or calculation were first produced. The subsequent
8 October internal analytical and arithmetic acceptance is recorded in the
unchanged [REVIEW_ACCEPTANCE.json](REVIEW_ACCEPTANCE.json). It authenticates the
proof, exact report and two independent review receipts by SHA256. It does not
claim external peer review, priority, APS acceptance or a cosmological origin.

The preserved [proof review](independent-math/UNIFORM_PROOF_REVIEW.md) checks 22
implications. Its [identity receipt](independent-math/UNIFORM_PROOF_REVIEW.json)
pins the theorem bytes. The separate
[arithmetic receipt](independent-math/INDEPENDENT_EXACT_UNIFORM_CONSTANTS.json)
checks the accepted exact constant report.

The unchanged receipts use paths from their original research checkpoint.
Their corresponding locations in this portable capsule are:

| Original checkpoint path | Capsule location |
|---|---|
| `branch-theorem/UNIFORM_BRANCH_THEOREM.md` | [UNIFORM_BRANCH_THEOREM.md](UNIFORM_BRANCH_THEOREM.md) |
| `independent-math/UNIFORM_PROOF_REVIEW.md` and `.json` | `baseline/independent-math/` |
| `independent-math/INDEPENDENT_EXACT_UNIFORM_CONSTANTS.json` | `baseline/independent-math/` |
| `numerical-audit/uniform-constants/EXACT_UNIFORM_CONSTANTS.json` | [public exact baseline fixture](../bounds/public_baseline/EXACT_UNIFORM_CONSTANTS.json) |
| `numerical-audit/uniform-constants/REVIEW_ACCEPTANCE.json` | [REVIEW_ACCEPTANCE.json](REVIEW_ACCEPTANCE.json) |

The historical evaluator source remains available at the immutable original
checkpoint. This capsule uses its stored public exact report as a regression
fixture; it does not require or execute the old evaluator. The new self-contained
producer reproduces all stored constants and polynomial coefficients for both
accepted baseline models when the original norm, root bound and endpoint estimate
are selected.

The older registered exact-decimal interval ends at approximately
`1.0915542913337286e-7`. The new theorem and its separate review receipts support
every `0 < delta <= 0.002` for the specified model and local tube. This extension
is a new mathematical bound; it does not retroactively validate archived
floating-point ODE trajectories.
