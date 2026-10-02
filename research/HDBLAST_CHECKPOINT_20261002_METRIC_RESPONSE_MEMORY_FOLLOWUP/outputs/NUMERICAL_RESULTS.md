# HDBLAST registered metric response: completed data, scientific FAIL

Ricardo Maldonado · 2 October 2026 · ORCID 0009-0009-3937-6527

The memory-only third registration completed all 12 primary observations and all four independent mode evolutions within the unchanged resource limits. It **fails** the registered conservation gates. The separately frozen reproduction wrapper correctly confirms this expected scientific failure; wrapper PASS does not establish calibration PASS or support a cosmological-origin claim.

| Registration | Completed physical data | Actual outcome |
|---|---|---|
| Original metric | One pre-source primary row; independent modes not reached | FAIL: fatal SciPy weighted-log continuum roundoff warning |
| High-precision followup | Twelve primary rows; positive-pulse coarse/fine archives | FAIL: peak 306088 KiB exceeds 262144 KiB |
| Memory-only followup | Twelve primary rows; both pulses coarse/fine, four full raw archives | FAIL:59 Ward endpoint/refinement gates |

The original freeze is 57668b8fadd75df8738565e0bbd1eb852c1ebae8 and its published failure is 5d4e92c90331808e2baf6d26d69717bf2d0a5590. The separate high-precision freeze is 5ff571cdb7729398834ff3c025b50191c4ee29f2 and its failure payload is published at c53287e51e34193575eb19d11a3f66025a52dbc8. Both prior failures, complete available outputs and original frozen sources remain preserved separately.

The executed third public freeze is **19fde76912af6e2f30d6f55e066b27d88cc7e34b**. Registration SHA256 is **1c5bc9b21b34b1d036b7a7d12ccdb6766ac2ba7835ea04182dba868843866607**; independent manifest SHA256 is 71830eb8aaf31bccd2cb628efeb9735490e41a7d7e5b5bffc527177f7c70e075. All 371 prospective files were checked against exact GitHub blob identities before any third physical calculation. The prior third-round preparation commit 9619921 was superseded only to correct introductory prose and include independent review before any physical evaluation; its immutable registration and verification remain in Git history and preparation evidence.

## What completed and what failed

The unchanged 32-command driver attempted 28 commands. Twenty-six pure/preflight commands and the primary producer passed. The independent producer completed all four evolutions, then deliberately exited 1 at its scientific gates. Normal/optimized main validators, accepted-calibration summary and figure commands were never reached. All stdout/stderr, exit records, source checks and raw data are retained.

The mode evolutions used two compact metric profiles at six observation times, fixed K=64,128,256, coarse dt=1/128 and fine dt=1/256, incoming Bunch–Davies state and extended precision. Each archive preserves 631 members. Independent elapsed time was **322.835271128 seconds**, with peak **103400 KiB**, below 900 seconds and 262144 KiB. No resource or timeout failure occurred in this round.

All declared observable refinement differences satisfy their unchanged thresholds. The 59 actual internal failures are 29 fine Ward endpoint checks and 30 Ward-ledger refinement checks. Maximum fine endpoint residual is **1.9552444892255006e-5**, above its **2e-6** gate. Maximum coarse/fine ledger difference is **0.02105427079282522**, above its **1e-6** gate. These are failed numerical conservation checks, not evidence of a physical instability.

## Separate saved-data diagnosis

A post-run checker review found an incorrect background derivative metadata formula in the frozen raw auditor: it used (n+2)!a^(n+2), while C=a''/a=2a² gives **C_n=2(n+1)!a^(n+2)**. The original checker and its failures remain unchanged. A separately labeled V2 changes that metadata reconstruction only, with exact identities checked under normal and optimized Python. It does not alter physical modes, operators, grids, outputs or scientific thresholds.

The strict V2 reconstructs both fine archives. All 180 primary-versus-reconstructed response comparisons and 144 baseline/current comparisons meet unchanged gates; the largest core gap is 1.752154e-12. However 29/36 fine conservation endpoints fail, with the independently reconstructed maximum matching the producer. All 48 saved primary differentiated-Ward residuals are at most 1.11e-16; these reuse the primary values and are not an independent second physical method.

Both coarse histories stop at the original raw-baseline tolerance: gap 1.283694e-12 exceeds its approximately 1.000056e-12 allowance. The completed record-only diagnostic reconstructs all 72 saved points and 72 conservation endpoints under normal and optimized Python, retaining exactly two coarse-history baseline failures and marking archive consistency UNVERIFIED. All 180 response refinement and 36 current refinement comparisons meet unchanged gates; all 29 endpoint and 30 ledger-refinement failures are independently recovered. It never claims validator or scientific PASS. Detailed reports, sources, exact proof and diagnostic plots are under review/saved_data_diagnostics.

Independent archival comparison checks all 1262 positive-pulse member pairs and 14151140 elements against the preceding failed round: dtypes, shapes and values agree exactly. Raw ZIP/storage hashes differ from unused long-double padding; the storage repair changes neither physical values nor numerical acceptance. The detailed byte-lane audit is retained separately.

## Mathematical next step and limits

The sampled conservation history contains frequencies up to 2K=512. On the coarse grid, the maximum phase per step is 4, exceeding the Nyquist phase π; even the fine maximum phase 2 is outside a uniform small-phase Simpson approximation. The standard single-harmonic Simpson transfer is θ(2+cosθ)/(3sinθ). This motivates a separately registered ledger diagnosis; it does not certify that every residual is explained by aliasing. The source-free terminal-mode ledger proposal in theory/next_ledger_diagnostic derives an exact phase-aware antiderivative with arbitrary unprojected normalization defect. Twenty-two identities and six mutations pass normally and with optimization. Twelve fixed proposed tests can separate signed Simpson error, free-flow residual and stress projection without defining density through conservation. No new physical run of that proposal is included.

This is one homogeneous conformal direction at H=1,x=r=2,xi=0 with fixed physical mass/reference/scalar and prescribed geometry. Actual shifted-root propagators, general lapse response, bulk/shell matching, state response and coupled initial data remain incomplete; the full metric prerequisite remains BLOCKED. Omitted-band UV bounds do not certify total finite numerical errors or continuum pressure/signs. No coupled stability, heat, particle yield, hot Big Bang, observation-based confirmation or novel fundamental discovery is established.

The deterministic package includes all curated inputs and raw evidence. Larger complete ZIPs are published as 64-MiB parts with authenticated external reassembly. Final ZIP replay receipts live outside this payload to avoid self-reference. Zenodo metadata describes failure evidence; a prepared record alone is not a deposit or DOI.
