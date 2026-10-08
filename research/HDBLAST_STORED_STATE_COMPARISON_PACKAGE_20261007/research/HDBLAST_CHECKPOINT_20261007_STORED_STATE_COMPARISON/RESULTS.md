# Exact comparison of the original HDBLAST incoming states with the prescribed BD target

Scientific study dated 7 October 2026; actual executions and this report completed 8 October 2026 UTC.

The complete original incoming state is now **ENCLOSED on every retained node** in the registered finite experiment. The comparison resolves a nonzero saved-minus-target mismatch in both normalized complex coordinates, including the imaginary component of U. The most prominent difference is in W on the two coarse grids near k=256: its maximum Cartesian complex L1 error is about 1.91e-14, while the fine-grid maxima stay near 5.3e-16 to 5.8e-16. These are comparisons with the prescribed scalar-model BD target, with fixed H=1. Their physical or numerical cause has not been established.

The result covers all 49,152 capsule-node occurrences, all twelve source/grid/cutoff cases, and the original retained positive momentum weights. It also gives continuous-time bounds for the finite weighted density and pressure contribution transported from this incoming mismatch under identical subsequent forcing and complete contacts. This completes the finite incoming-state comparison; the full continuous pressure/contact certificate remains **UNRESOLVED** and the historical metric calibration remains **FAIL**.

## Compared states and exact evidence

The two source prescriptions are `positive_B` and `signed_uB`, each with a coarse capsule of 8,192 nodes and a fine capsule of 16,384 nodes. Each capsule has three nested positional prefixes, exactly equivalent to k<K for K=64,128,256. Coarse counts are 2,048,4,096,8,192; fine counts are 4,096,8,192,16,384. The prefixes therefore contain 86,016 node appearances, while each of the 49,152 capsule-node occurrences receives its own target enclosure. Nodes at different sources/grids remain distinct occurrences even where represented coordinates agree. The actual run separately established equality of the two sources' represented nodes and weights within each resolution; no cross-resolution equality or quadrature accuracy is inferred.

The exact represented constants and compared coordinates are

```text
epsilon = 3777893186295716171 / 37778931862957161709568
Pi      = 14488038916154245685 / 4611686018427387904
U_saved = u_1 / epsilon,   W_saved = w_1 / epsilon
delta U = U_saved - U_BD,  delta W = W_saved - W_BD.
```

Division by epsilon occurs exactly once. Multiplication by the positive represented epsilon recovers the corresponding raw u_1/w_1 incoming perturbation errors. The original binary80 stored values and all their complex components are retained; no state projection, recalibration, new initialization, or later physical trajectory is introduced in this comparison.

The authoritative numerical result is [DATA.json](evidence/actual/normal/DATA.json); [SCIENCE_SUMMARY.json](evidence/actual/normal/SCIENCE_SUMMARY.json) is byte-identical. Their SHA256 is `2618232d6a2f713ce335d9becea301c9d6f30a314b088834398dde6a9a0ff98e`. Every number in the tables and figure below is **DISPLAY ONLY**. Table bounds use outward decimal rounding, and all sign, ordering, and identity decisions use the exact rational exports. Exact derived maximum-error endpoints and node witnesses are supplied in [the independent norm receipt](evidence/independent-actual/INDEPENDENT_SERIALIZED_NORMS_NORMAL_001.json).

## Incoming errors are resolved against complete target uncertainty

For X=U or W, the Cartesian complex L1 norm is ||X||_1=|Re X|+|Im X|. Let S_X be the reported maximum of the saved-minus-target rectangle upper norm in a prefix, and let R_X be the largest complete exported target L1 radius for that coordinate over the registered experiment. The exact true finite-prefix maximum obeys

```text
max(0, S_X - 2 R_X) <= max_j ||delta X_j||_1 <= S_X.
```

The rectangle upper at each node equals its midpoint L1 norm plus its L1 radius. A true target within that rectangle can differ from this upper by at most twice the radius; choosing a node attaining S_X proves the lower endpoint. Thus a large upper here is accompanied by a certified positive lower. For the complex modulus, the corresponding lower is max(0,S_X-2R_X)/sqrt(2), with S_X still a valid upper. The table is specifically an L1 enclosure, with no continuum interpolation between nodes.

The maximum complete target radii are exactly

```text
R_U = 948101 / 79228162514264337593543950336
R_W = 4971953 / 158456325028528675187087900672.
```

Their approximate magnitudes are 1.19667170096e-23 and 3.13774347544e-23. Both satisfy the registered target-radius gate of 1e-18. These are complete source-representation, cap, arithmetic, Taylor-truncation, and export uncertainties. Every retained node was evaluated; the eighteen additional diagnostic target evaluations are crosschecks and do not certify other nodes by interpolation.

**Table 1. True finite-prefix maximum incoming L1 errors, in U=u_1/epsilon and W=w_1/epsilon units.** All endpoints are outward-rounded display values. The global coordinate-specific target radius gives a conservative bracket for each prefix.

| Source / grid / K | Nodes | Maximum U error enclosure | Maximum W error enclosure |
| --- | ---: | --- | --- |
| positive_B/coarse/64 | 2,048 | [1.568250397e-16, 1.568250638e-16] | [5.695761712e-16, 5.695762340e-16] |
| positive_B/coarse/128 | 4,096 | [1.568250397e-16, 1.568250638e-16] | [5.695761712e-16, 5.695762340e-16] |
| positive_B/coarse/256 | 8,192 | [1.568250397e-16, 1.568250638e-16] | [1.913175710e-14, 1.913175717e-14] |
| positive_B/fine/64 | 4,096 | [1.570727708e-16, 1.570727948e-16] | [5.774951219e-16, 5.774951848e-16] |
| positive_B/fine/128 | 8,192 | [1.570727708e-16, 1.570727948e-16] | [5.774951219e-16, 5.774951848e-16] |
| positive_B/fine/256 | 16,384 | [1.570727708e-16, 1.570727948e-16] | [5.774951219e-16, 5.774951848e-16] |
| signed_uB/coarse/64 | 2,048 | [8.280069530e-17, 8.280071924e-17] | [5.251323659e-16, 5.251324288e-16] |
| signed_uB/coarse/128 | 4,096 | [8.280069530e-17, 8.280071924e-17] | [5.251323659e-16, 5.251324288e-16] |
| signed_uB/coarse/256 | 8,192 | [8.280069530e-17, 8.280071924e-17] | [1.908431753e-14, 1.908431760e-14] |
| signed_uB/fine/64 | 4,096 | [8.306256020e-17, 8.306258415e-17] | [5.264466551e-16, 5.264467179e-16] |
| signed_uB/fine/128 | 8,192 | [8.306256020e-17, 8.306258415e-17] | [5.264466551e-16, 5.264467179e-16] |
| signed_uB/fine/256 | 16,384 | [8.306256020e-17, 8.306258415e-17] | [5.264466551e-16, 5.264467179e-16] |

All twelve U brackets and all twelve W brackets have positive lower endpoints. Within each capsule the U maximum is unchanged across the three prefixes; the pronounced high-k increase appears in coarse W only. At K=256 the coarse W maxima have exact maximizing-node witnesses:

| Capsule | Zero-based index | Exact represented k | Approximate k, display only |
| --- | ---: | --- | ---: |
| positive_B/coarse | 8186 | 9218490285978407865/36028797018963968 | 255.8645041944 |
| signed_uB/coarse | 8187 | 1152491271506017559/4503599627370496 | 255.9044690611 |

The exact positive lower endpoints establish a genuine discrepancy relative to the prescribed target at the retained nodes. The locations and coarse/fine contrast identify where further preparation or solver-provenance investigation is useful. They do not establish why the discrepancy occurred, a convergence order, or a bound for momenta between nodes.

![Certified finite stored-state comparison bounds](figures/stored_state_bounds.png)

**Figure 1. Display of the exact reported upper bounds.** The top panels show S_U and S_W, the Cartesian complex L1 incoming rectangle-norm supremum uppers in U=u_1/epsilon and W=w_1/epsilon units. The bottom panels show the uniform finite weighted density and pressure contribution uppers over eta in [-9/2,-7/2], in the inherited a_0^4/epsilon scaled first-order R/P response convention with H=1. All four vertical axes are logarithmic; K takes exactly 64,128,256 in inherited momentum units. Blue denotes positive_B and orange signed_uB; solid curves denote coarse and dashed curves fine. Joining values guides the eye across nested prefixes. The drawing uses binary64 coordinates only for display; exact DATA and Table 1 determine the enclosures. These panels show bounds rather than sampled maxima of later trajectories. Coarse/fine curves represent different finite node sets and retained weights.

## Canonical normalization is retained as an exact diagnostic

For the real-forcing system U'=W and W'=2ikW-g, the conserved quantity is

```text
c_saved = (2 k Re U_saved - Im W_saved) / (2 k)
c' = 0,     c_BD = 0,     c_error = c_saved.
```

The prescribed target has zero initial data, so its invariant is exactly zero. The saved invariant is evaluated from the original state and preserved. The stable phase amplitude satisfies kA=-i delta W/2. The constant imaginary component of delta U-A still contributes to the full incoming U norm even though it does not enter the displayed first-order stress operators.

For f=exp(-ik eta)(1+epsilon U)/sqrt(2k), the normalization identity is

```text
(f f*' - f' f*)/i
  = 1 + 2 epsilon c
      + epsilon^2 [|U|^2 - Im(U* W)/k].
```

The reported 2epsilon|c| is the absolute **first-order** Wronskian defect. It excludes the quadratic finite-epsilon term; cancellation with that term is not asserted. For each capsule these maxima already occur inside K=64 and remain identical in its larger prefixes.

**Table 2. Canonical and weighted incoming diagnostics at K=256.** The c and first-order defect columns show outward upper displays of exact represented maxima. The weighted columns are positive retained-dk-weighted L1 uppers, sum_j w_j B_j, with inherited momentum-weight units; their measure differs from a supremum and from the stress measure mu.

| Capsule | max abs(c) | max 2epsilon abs(c) | Weighted U upper | Weighted W upper |
| --- | ---: | ---: | ---: | ---: |
| positive_B/coarse | 6.140849718e-19 | 1.228169944e-22 | 2.158889122e-12 | 9.805191621e-10 |
| positive_B/fine | 1.376825779e-18 | 2.753651558e-22 | 1.994582103e-13 | 6.963858168e-13 |
| signed_uB/coarse | 4.076916142e-19 | 8.153832283e-23 | 2.023907985e-12 | 9.544160225e-10 |
| signed_uB/fine | 8.918190818e-19 | 1.783638164e-22 | 1.102777003e-14 | 9.174081719e-13 |

The largest first-order Wronskian defect is approximately 2.75365155786e-22 on positive_B/fine. This is an exact canonical diagnostic with its stated perturbative scope; it neither repairs the state nor accounts by itself for the coarse high-k W discrepancy.

## Continuous-time transport of the finite incoming contribution

The transport interval is a=-9/2 to b=-7/2, with L=-1/eta, L_a=2/9 and L_b=2/7. Two exact linear evolutions start respectively from the saved state and prescribed target at a, and use the **same subsequent real forcing and complete contacts**. Their inhomogeneous difference cancels, leaving

```text
A = delta W(a)/(2 i k),   E = exp(2 i k(eta-a))
delta W(eta) = 2 i k E A
delta U(eta) = delta U(a) + A(E-1).
```

These statements describe transport of the incoming error component. They provide no residual certificate for later saved numerical trajectories and no proof that the implemented subsequent forcing or contacts equal those prescribed.

The finite positive measure is exactly mu_j=w_j k_j^2/(2Pi^2), with q_j=w_j/(4Pi^2). The direct operators used for density R and pressure P are

```text
delta R = [(2k^2+3L^2) Re delta U - k Im delta W - L Re delta W]/(2k)
delta P = [(2k^2/3-L^2) Re delta U - k Im delta W - L Re delta W]/(2k).
```

All reported R/P quantities are finite sums of these normalized first-order responses under the inherited a_0^4/epsilon scaling, with H=1. Epsilon is not inserted again inside the operators. Multiplication by epsilon/a_0^4 converts a normalized stress difference to the corresponding raw linear-stress convention at that time. These are model quantities, with no SI or observational interpretation assigned.

The uniform bound preserves cancellation in the exact signed canonical sum: for X=R or P it equals the larger absolute canonical endpoint plus a positive phase-envelope sum,

```text
upper_X = max(|C_X(L_a)|, |C_X(L_b)|) + phase_upper_X.
```

The canonical sums are affine in L^2, which makes their endpoint maximum exact. The phase term is an analytic majorant over continuous time; unknown oscillatory cancellation is not assumed. The separately exported canonical triangle sums remain available for comparison. At K=256 the phase majorant accounts for more than 99.99% of each coarse pressure upper. This identifies the dominant **reported bound component**, not a measured phase stress maximum.

**Table 3. Uniform finite stress and integrated-contribution uppers for all twelve cases.** R and P columns use the inherited a_0^4/epsilon scaled first-order response units. Work is |integral L(delta R-3delta P) d eta|=|delta R(b)-delta R(a)| in those density units. The last column bounds |integral delta P d eta|, in scaled pressure times inherited conformal-time units. Every display rounds upward.

| Source / grid / K | Uniform abs(R) upper | Uniform abs(P) upper | Work upper | Signed pressure integral abs upper |
| --- | ---: | ---: | ---: | ---: |
| positive_B/coarse/64 | 6.675155099e-16 | 3.370792044e-14 | 1.182030907e-15 | 1.550814190e-15 |
| positive_B/coarse/128 | 1.663318022e-15 | 2.848993197e-13 | 2.926645858e-15 | 3.845187261e-15 |
| positive_B/coarse/256 | 5.949313693e-13 | 3.268406612e-10 | 1.057377825e-12 | 1.387788615e-12 |
| positive_B/fine/64 | 6.645398279e-16 | 3.312296111e-14 | 1.177335326e-15 | 1.544539317e-15 |
| positive_B/fine/128 | 8.223809298e-16 | 6.541101571e-14 | 1.436703266e-15 | 1.888897736e-15 |
| positive_B/fine/256 | 1.192581116e-15 | 2.321551142e-13 | 2.082424730e-15 | 2.738681350e-15 |
| signed_uB/coarse/64 | 6.351393109e-16 | 3.389494782e-14 | 1.122602761e-15 | 1.473352772e-15 |
| signed_uB/coarse/128 | 1.626256238e-15 | 2.851996871e-13 | 2.878984669e-15 | 3.779399411e-15 |
| signed_uB/coarse/256 | 5.791460378e-13 | 3.181395735e-10 | 1.029400935e-12 | 1.351054868e-12 |
| signed_uB/fine/64 | 6.311077363e-16 | 3.308932393e-14 | 1.114393629e-15 | 1.462776306e-15 |
| signed_uB/fine/128 | 8.278068713e-16 | 7.650190027e-14 | 1.461525109e-15 | 1.918806839e-15 |
| signed_uB/fine/256 | 1.349200001e-15 | 3.058293984e-13 | 2.349734456e-15 | 3.091766641e-15 |

The absolute value of the signed pressure integral differs from the integral of absolute pressure. The latter has its own upper, duration times the uniform pressure upper; the exact interval duration is one, so its numerical upper equals the P column. The signed work and pressure-integral bounds retain exact signed canonical cancellation before adding the phase envelope. Pressure is computed directly from its operator, with continuity identities used as checks.

The coarse K=256 pressure uppers are about 3.27e-10 and 3.18e-10, while the fine uppers are about 2.32e-13 and 3.06e-13. These finite upper-bound differences accompany the coarse high-k W feature. A ratio of such uppers would be a ratio of bounds, rather than a certified ratio of realized pressure errors or a continuum convergence estimate.

Signed anchor intervals provide stronger evidence than a positive phase upper alone. Every one of the twelve R anchor intervals and twelve P anchor intervals excludes zero exactly. The full-prefix examples are:

**Table 4. Signed finite anchor differences at eta=-9/2 and K=256.** Units are the same scaled first-order finite density/pressure units as Table 3. Each endpoint is rounded outward; exact endpoint signs are checked in DATA.

| Capsule | Signed density difference interval | Signed pressure difference interval |
| --- | --- | --- |
| positive_B/coarse | [-1.774450715e-14, -1.774450198e-14] | [-6.678297407e-12, -6.678294699e-12] |
| positive_B/fine | [-2.120901192e-17, -2.120385524e-17] | [-2.033755236e-13, -2.033728165e-13] |
| signed_uB/coarse | [1.852460912e-14, 1.852461428e-14] | [1.192159981e-11, 1.192160253e-11] |
| signed_uB/fine | [-3.308393063e-17, -3.307877341e-17] | [2.805702342e-13, 2.805729416e-13] |

The finite anchor pressure difference is negative for positive_B and positive for signed_uB at all three cutoffs on both grids. The density signs vary with source, grid and prefix. These nonzero anchor intervals certify a realized finite direct-operator contribution at the incoming time. The uniform envelopes still need not be attained at any later time, and no phase-sensitive all-time maximum was measured.

## Source reuse, exact replay, and independent checking

All 22 source coefficient vectors were reused from the prior authenticated representation in the sense that their actual coefficient digests match exactly. They contain 2,486 exact real coefficients. The registered runs nevertheless each performed 22 physical source constructions to obtain and authenticate the current certificates. The conservative source-error rule is max(rebuilt_error, prior_uniform_error, prior_analytic_tail+prior_coefficient_error). Its source L1 model-error totals are approximately 4.66967072853e-28 for positive_B and 4.63601799593e-28 for signed_uB. These terms are already included in the complete target enclosures and are not added again to the stored-state difference.

Historical continuous prehistory source/coefficient stress bounds use a different momentum measure from this finite retained-weight comparison. Their scale may supply context, but these results do not establish a rigorous dominance ratio across those measures or an additive full error budget. A compatible source-component upper and state-contribution lower under the same measure would be needed for such a claim.

The normal and `python -O` actual executions are successful, and [the exact replay receipt](evidence/ACTUAL_EXACT_REPLAY.json) records byte identity for eight complete scientific artifacts: DATA, SCIENCE_SUMMARY, NODE_TARGETS, SOURCE_CERTIFICATE, DIAGNOSTIC_CROSSCHECK, OUTPUT_READBACK, DECODE_ATTEMPTS and SOURCE_ATTEMPTS. Timing, memory and derived resource-receipt fields have their documented exclusions from this scientific-byte identity. Each run recorded 20 selected retained-array decodes, 294,936 selected scalar binary80 slots, 49,152 retained-node target evaluations, and 18 diagnostic target evaluations. There were no fabricated nodes, unselected member decodes, physical trajectories, or likelihood evaluations.

The [independent actual-output review](evidence/independent-actual/ACTUAL_OUTPUT_INDEPENDENT_REVIEW.md) reports PASS with no blocking findings. Frozen strict readback checked all exported nodes and all twelve prefixes. A separate standard-library computation, importing no production numerical modules, reproduced incoming maxima, weighted L1 sums, canonical/first-order Wronskian diagnostics, stable kA values, finite measure sums, and exact maximizing-node witnesses. It used zero original retained-array reads or decodes, zero source constructions, and zero target evaluations. The present report likewise uses only serialized results and the frozen mathematical documentation, with exact summary-arithmetic checks.

The interpretation rests on the registered analytic target enclosure proof, authenticated source construction, and pinned decoder/worker for original-state identity and physical target truth. Independent serialized readback verifies the exported intervals, budgets, finite arithmetic and custody; it does not independently reconstruct physical targets. Internal execution and mathematical review are not external peer review or proof-assistant formalization.

## Registration and remaining scientific questions

The frozen study is bound to commit `bedcca7e86995da1230c12a20fae3755b31f4a94`, registration SHA256 `05e9943f0b4c8134252a2fecef7631ddba8bd398d18553d6e126fbf50fc3aacc`, and externally pinned GO SHA256 `27784552d9d7e10167066927d54f1757351b0f925b499e951b79849a44b62769`. The completed independent review authenticated all 117 registered files and verified the 17 previously reviewed production Python files unchanged. The frozen theory, transport and preparation documents describe prospective or fabricated validation where marked. This RESULTS report, actual execution artifacts, independent actual review and figures are clearly separated **post-freeze delivery supplements**; they do not rewrite the registration or replace its evidence.

The independent serialized-norm review script is archived in the outer delivery package at `preparation/root-independent-review/review_serialized_norms.py`, outside the research checkpoint C. Its serialized review results are included under `evidence/independent-actual/`. It is an archival independent-review implementation, not a registered production target evaluator.

The resolved result is the complete original incoming-state error on the finite retained nodes, together with its conditional finite-weight continuous-time stress transport. Remaining prerequisites include continuous-momentum incoming-state control; later solver/arithmetic residuals; actual source and complete-contact evaluation errors; time and momentum quadrature errors; and ultraviolet completion. All-node coverage does not certify these terms. The original full twelve-case pressure/contact gate remains **UNRESOLVED**, historical metric calibration **FAIL**, a higher-dimensional Big Bang cause **NOT_ESTABLISHED**, and external mathematical novelty **NOT_ASSESSED**. No numerical comparison with the historical full gate is made without its exact units and full-budget definition.

The next useful scientific step is a separately registered investigation of incoming-state preparation and the coarse high-k feature, preserving these original states and the failed calibration. A prescribed BD initialization with independent phase and normalization control, followed by a solver-residual and complete-contact analysis, would test the origin of the discrepancy. If tighter transported stress bounds are needed, a proved phase-sensitive estimate is required; sampled extrema or a projection of the saved states would not provide that certificate.

