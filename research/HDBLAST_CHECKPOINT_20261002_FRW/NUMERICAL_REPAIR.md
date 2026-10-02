# Preserved numerical failures and representation repair

The first automated run, [36957517031](https://github.com/maldonado-research/HDblast-archive/actions/runs/36957517031), stopped at compilation because of a copied dictionary-comprehension syntax error. No new physical calculation ran. Commit 0ee14880ffa04b938bda62cb54fad26ba2c70b93 repaired that mechanical error without changing the protocol or thresholds.

The next run, [36957629852](https://github.com/maldonado-research/HDblast-archive/actions/runs/36957629852), passed input/protocol hashes, 72 exact-rational graded Ward checks, 36 flat-matching checks, omission-detection controls and all 18 oscillator evolutions. It computed all eight archive reconstructions, but the registered 1e-9 continuity test failed after conversion of the higher-degree B-splines into separate floating-point monomial pieces:

| Reconstruction | Fine | Coarse |
|---|---:|---:|
| C4 quintic converted pieces | 4.65e-9 | 5.86e-9 |
| C6 septic converted pieces | 1.09e-7 | 4.00e-8 |

Those complete results are preserved in outputs/geometry_audit_float_conversion.json. Hermite and cubic continuity checks passed.

The intended native B-splines with simple interior knots are mathematically C4 and C6. Rounding converted polynomial coefficients and then subtracting high derivatives loses that consistency numerically. The repair evaluates the intended native basis consistently and audits one-sided limits with higher precision. It retains the original 1e-9 threshold and all underlying samples. It is a representation repair, not a relaxation of the test or a newly chosen physical smoothing fit.

A recovered numerical pass cannot establish that the high derivatives are physical. The separate reconstruction, state, matching and Hadamard limitations remain. Descriptive comparisons added after the first outcomes are labeled post hoc and do not add or change acceptance criteria.

The repaired source at commit 2bad8db9b1531732b917a7783e45e8257b0f7478 passed [run 36959579410](https://github.com/maldonado-research/HDblast-archive/actions/runs/36959579410). All eight registered reconstruction continuity checks passed. C4 and C6 report zero residual at the working precision because both sides are evaluated from the same native knot/coefficient representation; no theoretical zero replaces an evaluated residual. PPoly is retained only for candidate roots, which are refined and checked on the native function. The highest-derivative jumps are retained. An additional post hoc diagnostic compares native SciPy evaluation with the 80-digit de Boor evaluator and explicitly measures interpolation fidelity at stored values.
