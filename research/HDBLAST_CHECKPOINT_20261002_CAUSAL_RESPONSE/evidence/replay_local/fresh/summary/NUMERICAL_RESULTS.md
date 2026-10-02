# Registered fixed-geometry causal-response calibration

The bounded calibration passed all 36 finite-cutoff comparisons and all 17 wrong-formula controls, both normally and under Python `-O`. The independent forced modes agree with the exact finite-cutoff expression at the twelve registered source/observation points. Both pulses leave a negative response at the two registered observations after the source has vanished.

These are responses of the scalar variance on a prescribed de Sitter geometry. The test establishes one matched retarded susceptibility; it does not establish stress response, shell stability, relaxation, particle production, radiation transfer, thermalization, or heating.

## Frozen experiment and normalization

H=1, r=2, a(eta)=-1/eta, epsilon=1e-4. With u=eta+4 and B=exp(1-1/(1-u²)) on |u|<1, the two sources are s=epsilon B and s=epsilon uB. Both vanish outside -5<eta<-3. The incoming BD state and geometry are fixed. The observation grid is {-5.5,-4.5,-4,-3.5,-2.5,-1.5}; cutoffs are K={64,128,256}.

Every response and response discrepancy below uses **y=a² deltaQ/epsilon**. To recover the physical response, deltaQ=epsilon y/a². The small changes between finite K and removed cutoff are regulator errors, not changes in geometry or source amplitude.

## Numerical evidence

| Quantity | Recorded value | Interpretation |
|---|---:|---|
| Maximum modes / finite-K expression difference | 6.938893904e-17 | Cross-route numerical agreement |
| Maximum coarse / fine mode difference | 1.283695372e-16 | Empirical refinement estimate |
| Maximum two-memory-form difference | 6.938893904e-18 | Internal representation check |
| Maximum finite-K / removed-cutoff difference, all K | 1.838939684e-06 | Actual regulator discrepancy |
| Maximum K=256 / removed-cutoff difference | 1.144687630e-07 | Actual regulator discrepancy |
| K=256 nonzero analytic tail bounds | 1.153855612e-06 to 2.340779456e-06 | Constructive UV remainder enclosure |
| Maximum primary quadrature estimate | 1.375462960e-14 | QUADPACK estimate, not a proof |
| Maximum linear Wronskian residual / epsilon | 9.757819552e-14 | All independent mode steps/nodes |
| Maximum reconstructed Wronskian residual / epsilon | 9.692347929e-14 | Observation-time canonical modes |

The analytic tail bound encloses only the omitted combined momentum integral. Quadrature estimates, mode refinement differences, and cross-route differences remain numerical evidence; their small observed values do not certify all displayed digits or make the entire response a rigorous interval. Before support the response and tail are exactly zero. Post-pulse regulator differences are much smaller than their conservative analytic bounds.

## All registered observations

| Source | eta | s/epsilon | Removed-cutoff y | abs(y256-y) | Analytic K=256 tail |
|---|---:|---:|---:|---:|---:|
| positive_B | -5.5 | 0 | 0 | 0.00000e+00 | 0.00000e+00 |
| positive_B | -4.5 | 0.716531311 | 0.0127511894936 | 9.91541e-08 | 1.22896e-06 |
| positive_B | -4 | 1 | 0.000624222328826 | 6.64308e-08 | 1.23682e-06 |
| positive_B | -3.5 | 0.716531311 | -0.0152323526152 | 1.05853e-07 | 1.31682e-06 |
| positive_B | -2.5 | 0 | -0.0110398935227 | 8.48314e-14 | 2.34078e-06 |
| positive_B | -1.5 | 0 | -0.00627782869388 | 2.25306e-14 | 2.34078e-06 |
| signed_uB | -5.5 | -0 | 0 | 0.00000e+00 | 0.00000e+00 |
| signed_uB | -4.5 | -0.358265655 | -0.00449596958902 | 1.11120e-07 | 1.15898e-06 |
| signed_uB | -4 | 0 | 0.00764278557548 | 1.79292e-14 | 1.15386e-06 |
| signed_uB | -3.5 | 0.358265655 | 0.0057897696856 | 1.14469e-07 | 1.37432e-06 |
| signed_uB | -2.5 | 0 | -0.00127426913307 | 1.22976e-13 | 2.30771e-06 |
| signed_uB | -1.5 | 0 | -0.000409000583736 | 4.92385e-14 | 2.30771e-06 |

The zero source at eta=-2.5 and -1.5 coexists with a negative variance response for both histories. For the positive pulse this is the signed retarded-memory test. For the signed pulse, pairing u and -u makes the post-pulse convolution positive before the common overall minus sign. The finite-part local contact vanishes at these post-pulse observations; it cannot account for their nonzero response.

## References and controls

The independent analytic stationary identity gives Q_x=-(2 gamma_E+ln2)/(16 pi²)=-0.011699927596383. This is a separate past-infinite-source identity, not a substitution of an instantaneous equilibrium response for the pulse histories. Normalization, matching contacts, fixed reference, strict causality, both conformal factors, occupation-state sensitivity, and the initial beta boundary term were tested. The frozen symbolic suite supplies the additional exact mass-law and coordinate-rescaling identities.

The deliberately occupied state and initial-beta diagnostics are separate controls; they do not alter the BD state of the registered primary experiment. Every numerical wrong-formula control was rejected under normal and optimized Python execution.

## Figures and reproducibility

`causal_memory.svg/.png/.pdf` shows only the six registered times per source. Source markers are raw recorded values; shaded support comes from registration. Response lines guide the eye and are not additional evaluated histories. `causal_agreement.svg/.png/.pdf` compares the registered K=256 results with the removed-cutoff memory response and its analytic UV bound. Exact pre-pulse zeros are omitted on the logarithmic axes. `causal_cutoff.svg/.png/.pdf` reduces the existing points to maxima at each of the three registered cutoffs.

The figures and this report are explicitly post-run presentation. No additional source, response, mode, or denser time-grid evaluations were performed. `SUMMARY.json` contains full-precision values, source hashes, raw archive hashes, and the registration digest. The renderer has separate presentation provenance and does not modify any frozen scientific input.

Public prospective freeze: `a01d17e015f070be2ad2018745e877148fe2f4e4`. Full registration SHA256: `f59abda221a6a0296e839bc3022298bdf81fac913feb3877ac8f1485f4d71da6`. Primary runtime: 0.134226 s. Independent mode runtime: 2.318429 s; peak RSS: 49908 KiB. These timings are execution metadata, not portable performance guarantees.

## Exact raw and presentation hashes

| Artifact | SHA256 |
|---|---|
| FULL_REGISTRATION.json | `f59abda221a6a0296e839bc3022298bdf81fac913feb3877ac8f1485f4d71da6` |
| outputs/primary/results.json | `713d0081b18948ffd51f8fa3f7c288fd3173295e8d90a7bcc052e6efab3dbf50` |
| outputs/independent/results.json | `a6443d96258e92b0f7869f0edad4ca2bb6ef50e8db735761d7b9319191402bd8` |
| outputs/CHECKS.json | `c1001dcc487648584bfc027786ecdf00edec0e3042d9ede167374c77e436230f` |
| outputs/CHECKS_OPTIMIZED.json | `32a8d5c9c46a792921647b18fc873b8e884d3c06e6bd2b1679c49ede7eaa68c8` |
| postrun/summarize_causal.py | `664b33b62a5357165c8d83a3cd95b251681b589ae892486bd42e1fdd2eca192e` |
| outputs/independent/modes_positive_B_coarse.npz | `c40b2838511a369d0d5fd28d08c204a819bb179b8539a2512e16e8db9c76de52` |
| outputs/independent/modes_positive_B_fine.npz | `648a1bba299be7d1d2e3401e0aa9058c636642d4162a2fc6d224cf237ec0cf7d` |
| outputs/independent/modes_signed_uB_coarse.npz | `deb4073258c1940a06626009219c1962da9aeeba209f92bf6a38d3d104251398` |
| outputs/independent/modes_signed_uB_fine.npz | `6bd3cfe0bbb2b9223597adc3ca21814f74976e72985c7b341d3ba6bf4e320d09` |
| outputs/primary/started.json | `2a259fd00316ab3d4d630c6f722b6d1e0eea1443247a552aa046137acc9ca3bc` |
| outputs/primary/EXECUTION.json | `128805d7627dc1c1205247128648a304ffe5a88d0fa0d3828c014dc6fee4390d` |
| outputs/primary/partial_results.json | `524919227f6365ae1869f886a3475e4337e4c6561cd267114eeacd5c107d3c87` |

All frozen source hashes are preserved in `SUMMARY.json` under `frozen_source_sha256` and in `FULL_REGISTRATION.json`.
