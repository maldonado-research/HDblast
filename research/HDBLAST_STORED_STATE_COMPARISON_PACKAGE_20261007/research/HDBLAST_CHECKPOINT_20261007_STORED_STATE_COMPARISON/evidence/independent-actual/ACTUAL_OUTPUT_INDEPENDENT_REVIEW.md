# Independent review of the actual stored-state comparison

**PASS.** No blocking finding remains in this review scope. Both actual runs have successful authenticated entry and bounded-execution receipts. All eight scientific artifacts are byte-identical between normal and optimized Python.

The reviewed public registration is `05e9943f0b4c8134252a2fecef7631ddba8bd398d18553d6e126fbf50fc3aacc` at freeze `bedcca7e86995da1230c12a20fae3755b31f4a94`. The externally pinned GO is `27784552d9d7e10167066927d54f1757351b0f925b499e951b79849a44b62769`. All 117 registered files authenticate, and all 17 previously reviewed production Python files are unchanged.

The frozen strict readback recomputed all 49,152 exported nodes and all 12 prefixes. A second, independent standard-library computation, importing no production numerical module, exactly reproduced the U/W maxima, weighted error sums, canonical suprema, first-order Wronskian defects, stable kA quantities, and finite measure sums. It also recorded exact maximizing node witnesses. Neither review opened or decoded an original retained array, constructed a source, or evaluated a target.

## Quantified incoming-state errors

The table encloses the true maximum Cartesian complex L1 error over each finite prefix. Endpoints are rounded outward. It uses the rigorous interval `[max(0,S−2R),S]`, where S is the registered rectangle-norm upper maximum and R is the complete target L1 radius for that coordinate. This is an L1 statement, not a complex-modulus or continuum statement.

| Source / grid / K | Maximum normalized U error | Maximum normalized W error |
| --- | --- | --- |
| positive_B/coarse/64 | [1.568250397E-16, 1.568250638E-16] | [5.695761712E-16, 5.695762340E-16] |
| positive_B/coarse/128 | [1.568250397E-16, 1.568250638E-16] | [5.695761712E-16, 5.695762340E-16] |
| positive_B/coarse/256 | [1.568250397E-16, 1.568250638E-16] | [1.913175710E-14, 1.913175717E-14] |
| positive_B/fine/64 | [1.570727708E-16, 1.570727948E-16] | [5.774951219E-16, 5.774951848E-16] |
| positive_B/fine/128 | [1.570727708E-16, 1.570727948E-16] | [5.774951219E-16, 5.774951848E-16] |
| positive_B/fine/256 | [1.570727708E-16, 1.570727948E-16] | [5.774951219E-16, 5.774951848E-16] |
| signed_uB/coarse/64 | [8.280069530E-17, 8.280071924E-17] | [5.251323659E-16, 5.251324288E-16] |
| signed_uB/coarse/128 | [8.280069530E-17, 8.280071924E-17] | [5.251323659E-16, 5.251324288E-16] |
| signed_uB/coarse/256 | [8.280069530E-17, 8.280071924E-17] | [1.908431753E-14, 1.908431760E-14] |
| signed_uB/fine/64 | [8.306256020E-17, 8.306258415E-17] | [5.264466551E-16, 5.264467179E-16] |
| signed_uB/fine/128 | [8.306256020E-17, 8.306258415E-17] | [5.264466551E-16, 5.264467179E-16] |
| signed_uB/fine/256 | [8.306256020E-17, 8.306258415E-17] | [5.264466551E-16, 5.264467179E-16] |

The maximum complete exported target radii are approximately `1.19667170096e−23` for U and `3.13774347544e−23` for W, below the preregistered `1e−18` gate. These contain source-representation, cap, arithmetic, Taylor truncation, and export uncertainty. The two coarse-grid W maxima occur at the following retained nodes; indices are zero-based.

| Case | Index | Exact represented k | Approximate k |
| --- | ---: | --- | ---: |
| positive_B/coarse/256 | 8186 | `9218490285978407865/36028797018963968` | 255.8645041944 |
| signed_uB/coarse/256 | 8187 | `1152491271506017559/4503599627370496` | 255.9044690611 |

These are observed locations of maxima, not an explanation of their cause. Normalized U/W errors divide the stored u_1/w_1 perturbations by the exact represented epsilon once; multiplying by epsilon recovers raw incoming perturbation errors.

## Source certificate and normalization

All 22 actual source vectors, comprising 2,486 exact real coefficients, match their prior authenticated coefficient digests. Source uncertainty uses `max(rebuilt_error, prior_uniform_error, prior_analytic_tail + prior_coefficient_error)`. The source L1 model-error totals are approximately `4.66967072853e−28` and `4.63601799593e−28`; these are components of the complete target uncertainty, not substitutes for it.

Real forcing conserves `c=Re U−Im W/(2k)`, and the prescribed zero-initial-data target has exactly c=0. The saved c is retained without projection. The largest absolute first-order Wronskian defect across these cases is approximately `2.75365155786e−22`; it is `2 epsilon |c|`, excluding the quadratic finite-epsilon term. Exact values and witnesses are in the norm receipt.

## Finite stress scope and replay

The inherited positive weights are used exactly with `mu=w k²/(2 Pi²)` and the represented Pi. Stress bounds use the inherited `a_0^4/epsilon` linear-response convention. For K=256, the normalized finite pressure-difference upper bounds are:

| Source / grid | Uniform finite pressure upper bound |
| --- | ---: |
| positive_B/coarse | 3.268406612E-10 |
| positive_B/fine | 2.321551142E-13 |
| signed_uB/coarse | 3.181395735E-10 |
| signed_uB/fine | 3.058293984E-13 |

These bounds compare two exact linear evolutions initialized at the compared incoming states and using identical subsequent forcing and complete contacts. They do not certify later saved numerical trajectories, source/contact implementation error, continuum quadrature, or ultraviolet completion. The original full pressure/contact gate remains **UNRESOLVED**, historical calibration remains **FAIL**, higher-dimensional origin remains **NOT_ESTABLISHED**, and external novelty remains **NOT_ASSESSED**.

Each actual run records 22 physical source constructions, 20 selected retained-array decodes, 49,152 target-node evaluations, and 18 additional diagnostic evaluations. The diagnostic probes do not certify other nodes by interpolation. The reviewer added zero source constructions, target evaluations, or original-array reads/decodes.

Independent byte comparison reproduced the root replay result for eight complete scientific artifacts. Entry receipts agree outside their derived RESOURCE file size/hash; both custodian receipts correctly bind their mode, registration, GO, successful exit, entry bytes, and child log. Custodian wall times were 125.082 seconds normal and 133.682 seconds optimized; peak resident sizes were 89,856 and 90,272 KiB.

The exact review remains conditional on the registered analytic enclosure proof and the pinned worker/decoder for original-input identity and target truth. Serialized readback independently verifies the exported intervals, budgets, finite sums, and custody; it is not a fresh physical-target construction.
