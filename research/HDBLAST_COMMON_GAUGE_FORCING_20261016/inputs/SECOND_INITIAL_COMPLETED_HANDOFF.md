# HDBLAST quantitative initial-data checkpoint — 15 October 2026

Published, merged, deployed and independently read back. This supersedes CURRENT_CONTINUATION_STATUS.md/json in this task. The full second-order solution remains unfinished.

This continues the already merged coordinate-shell-source, interior-master and explicit compact initial-constraint checkpoints. Those frozen sources, archives and the earlier manuscript are preserved. The selected state remains the registered aligned, admissible quadrupole. No empirical cosmological conclusion follows.

## New mathematics and exact scope

The previous checkpoint constructed a particular physical second Cauchy jet cancelling the Hamiltonian and all four momentum sources. The present result quantitatively bounds that same jet, including nonzero second scalar and unit-normal scalar-velocity corrections, on the entire original packet annulus ell<=z<=2ell with ell=1e-6, all of S3, and its smooth flat zero extension.

| Physical second initial field | Absolute component bound | Full first spatial covariant derivative component bound |
|---|---:|---:|
| gamma2 | 1.06e24 | 4.25e29 |
| K2 | 7.96e7 | 1.96e13 |
| phi2 | 5.13e33 | 3.06e39 |
| Pi2 | 2.00e11 | 5.02e16 |

Bounds concern individual adapted background frame components, including the derivative slot. K=+Lie_n(gamma)/2 is covariant physical extrinsic curvature; Pi=n(phi) is an ordinary normal velocity, not a momentum density. Expansions use epsilon²X2/2. The numbers are coarse upper bounds in model units, not actual response sizes, lower bounds, arbitrary boosted-frame norms, Sobolev energies or finite-amplitude intervals.

Two exact reductions improve the interval calculation: the first constraint eliminates explicit u_y,u_yy,a_y from the Hamiltonian source, and the registered velocity relation cancels v*v_y in the scalar momentum correction. The interval evaluations preserve these correlations before differentiating/bounding. The fourth physical bump derivative is certified, including ell^-4 scaling; background coefficients through order four are enclosed throughout the packet annulus.

The five initial lapse/shift time velocities are explicit through both parameter orders. They enforce the relative harmonic coefficients initially while preserving the physical Cauchy data. The fixed reference connection is contracted with the perturbed inverse metric. The second spacetime conversion retains g2_tt=-2rho²U² and (dot g_ty)_2=dot beta2_y-8U dot beta1_y. Its background-frame connection cancels a further mixed term: (nabla0 g2)_0y=-26U*U_y-(3/2)f_y-3a*f.

For the complete second metric/scalar INITIAL spacetime one-jet, |g2|<1.06e24, |nabla[g0]g2|<6.4e29, |phi2|<5.13e33 and |dphi2|<3.06e39 in that fixed frame on t=0. This is neither a positive-time estimate nor a spatial C1 bound for every gauge time velocity. The initial normal gauge coefficients vanish through order two conditional on the corresponding reduced acceleration equations and verified physical constraints. Those accelerations are not independently prescribed Cauchy data.

## Sources and verification

Scientific capsule: https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_SECOND_INITIAL_NORMS_GAUGE_20261015

Local reproduction evidence: https://github.com/maldonado-research/HDblast/tree/main/hdblast/continuation/2026.10.15-second-initial-norms-gauge

The frozen capsule has 66 files including its manifest. Manifest SHA256: 4c0dae7de0b63acca15e5204681a859f49e5b5028340457dea5c29246bf0626c. ZIP SHA256: d5cde2f6340e4379afced4464e0061e229a887bd4774adc903487456a3445efb (499673 bytes). The shared replay runner is externally pinned too.

All 21 registered fresh extracted processes matched complete expected reports; all 70 intended-error subprocesses were rejected for their specified reasons and left no output report. The immutable source remained unchanged. A distinct implementation uses unnormalized derivatives and a different outward rounding grid to reproduce all 17 rounded physical caps. Direct five-dimensional Christoffel checks are separate from all-point formal ADM and source identities. The fresh referee independently verified the required 8ell first-exit collar enclosure, avoiding reliance on an older prose note's shorter collar.

The reviews are independent internal AI reviews, not external specialist peer review. Five targeted queries and selected passages from five primary papers support the established harmonic initialization/constraint architecture; they establish no priority or general new theorem. Mathematical proofs and interval coverage support the field identities; check counts alone do not.

- Research PR75: https://github.com/maldonado-research/HDblast/pull/75 — merge `bc24776801b74fde85e299ed3ff4be27649338f9`, reviewed head `a9e23f3713fd8b18d87fd569ea4965674d873136`.
- Shared-site PR64: https://github.com/maldonado-research/maldonado-research.github.io/pull/64 — merge `1b98351e4f8ff6bf88faf58bf83c16b29126cfee`, reviewed head `29c4f7f3146f1ca764489907481e96664e2e5e00`.
- The shared branch was rebased onto newer main `95d1bcad7b082e0666482f565641ddc1ee065d44`, keeping the unrelated TMD changes. Our shared commit changes only the HDBLAST project page.
- All three research premerge checks passed. The new dedicated postmerge replay and both website deployments passed. Two automatic postmerge reruns of unchanged legacy workflows remain running in the final saved snapshot; their premerge counterparts passed. Do not describe all postmerge checks as complete.
- Pinned remote readback passed for all141 changed files and all3 live pages with their full new blocks and prior markers. Both merge trees equal the reviewed trees. Fresh before/after checks preserve all5 saved original/prior-publication workspaces, all116 earlier frozen source files and all3 prior ZIP hashes.
- The final receipt is publication/INDEPENDENT_FINAL_READBACK.json, frozen read-only. SHA256: `df1c4af476df913764f6b1e5ccba0bdde943ee895ad221ae36c403d4ea8d7107`.


## Next scientific calculation

The full U2 evolution remains unfinished. No first boundary-preserving gauge wave, consistently transformed evolving bulk/shell source, coupled second evolution, high-Sobolev existence estimate, exact nonlinear family, remainder bound or empirical cosmological origin has been established here. The old longitudinal forcing cannot be relabeled as harmonic forcing after changing the first representative.

Continue with NEXT_COMMON_GAUGE_EVOLUTION.md. A separate unpublished note, future_work/gauge-wave/FIRST_GAUGE_VECTOR_SYSTEM_UNPUBLISHED.md, now specifies the three-amplitude gauge-vector equations and candidate shell-preserving rows. Root and another agent derived the same system; its initial acceleration agrees with the published gauge completion. It permits tangential shell reparametrization/clock transport. It supplies no gauge-wave solution, well-posedness theorem or numerical certificate and is deliberately outside the frozen public capsule.

Next construct and bound that actual gauge wave and the needed derivatives, transform both parameter orders with all Lie, embedding and clock terms, evaluate the common-gauge quadratic source and all physical/gauge/subsidiary boundary rows, then solve and validate the coupled second-order problem. A bulk/energy residual bound alone does not control required shell time/normal derivatives: a suitable observation estimate remains necessary.

No new Zenodo version or journal submission was made. Private editorial history remains outside public supplements. The earlier manuscript and its separate journal preparation were untouched. Machine-readable details are in COMPLETED_HANDOFF.json; exact source/readback/CI receipts are linked there.
