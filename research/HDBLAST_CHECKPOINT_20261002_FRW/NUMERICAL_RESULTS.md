# Executed numerical results

Primary repaired run: [36959579410](https://github.com/maldonado-research/HDblast-archive/actions/runs/36959579410), source 2bad8db9b1531732b917a7783e45e8257b0f7478. Original conversion failure: [36957629852](https://github.com/maldonado-research/HDblast-archive/actions/runs/36957629852). See the audit for later supplementary diagnostics and packaging replay.

## Archived background

The fine and coarse tuned source-free histories use delta=0.1, Y=0, G=100 and phi_star=0.5, restricted to s in [0,6.9]. These are dimensionless model quantities with a(0)=1, not SI measurements. All four input byte hashes and both protocol hashes passed. Each history has 300 retained stored samples in the window. Diagnostics use their union with 2001 uniform times; reported maxima are sampled maxima, not bounds.

| Formal Hermite diagnostic | Fine | Coarse |
|---|---:|---:|
| Interior knots | 299 | 299 |
| Maximum absolute Delta U from -a² Delta Hdot | 0.0030185104 | 0.0059493905 |
| sum (Delta U)² | 2.82512699e-5 | 9.85115820e-5 |
| Interior-knot coefficient of ln K | 2.27795500e-11 | 7.94625405e-11 |
| Initial amplitude-WKB coefficient of ln K | 3.20130485e-6 | 3.20211722e-6 |
| Initial/knot coefficient ratio | 140534.16 | 40297.19 |
| a(6.9) | 7.916090993 | 7.915325537 |
| phi=0.5 crossing time | 1.486121632 | 1.486015123 |

The initial and knot coefficients diagnose different formal extensions; they are not observed radiation energies. The initial-state issue is not removed by smoothing. The factor 3.4883 coarse/fine knot coefficient shows reconstruction dependence, not an extrapolated continuum convergence rate.

An independent closed-form Hermite calculation agrees with the polynomial-limit jumps within 2.48e-10 in the registered scaling. An independently written JavaScript calculation agrees with the fine sum of squared jumps to numerical precision appropriate to the cancellation. Its result is retained in outputs/independent_geometry_check.json.

## Smoothness and derivative sensitivity

All eight reconstruction continuity tests pass at the unchanged 1e-9 threshold after the native-basis repair. Hermite residuals are at most 2.22e-16, cubic at most 2.66e-13, and C4/C6 one-sided residuals round to zero at the 80-digit working precision. These are evaluations of the specified interpolants, not estimates of physical interpolation error.

Fine C4 versus C6 differences on all 2001 uniform times:

| Quantity | Maximum absolute difference | Fraction of sampled C4 maximum |
|---|---:|---:|
| phi | 9.1358e-9 | 8.8175e-9 |
| H | 1.9967e-7 | 1.9960e-7 |
| Hdot | 9.4963e-5 | 7.7036e-5 |
| Hdddot | 4.48235 | 0.20533 |
| U | 3.85512e-4 | 2.1405e-9 |
| Rddot | 26.28140 | 0.15037 |
| regular unit-R² pressure | 103.92074 | 0.15329 |

The last two maxima occur at s=0. The C4/C6 fits impose stored values rather than stored derivatives: fine saved-time phi_dot discrepancies are approximately 0.0040811 and 0.0040878. These matters are physical input/reconstruction uncertainty, separate from the repaired evaluation error. Curvature distributions at knots are omitted from regular-part plots and expressly retained as a limitation.

Generic interior-knot beta powers are k^-2 (Hermite), k^-3 (cubic), k^-5 (C4), and k^-7 (C6). Higher powers do not constitute a certified finite-cutoff bound: high derivative coefficients are large and the asymptotic onset has not been bounded.

## Algebra and independent evolution controls

- Positive-reference FRW subtraction: 24 exact rational local jets, 72 graded exchange assertions; 12 flat jets, 36 flat-matching assertions. Both pressure-omission and current-omission controls are detected.
- DOP853: all 9 cases pass; maximum occupation absolute error 2.220446e-16.
- Independently written real Radau: all 9 cases pass; maximum occupation absolute error 2.872702e-15.
- Maximum DOP853/Radau occupation difference 3.053113e-15. Agreement at tiny occupation is an absolute-error result, not meaningful relative accuracy.
- Abrupt-step final k^4 occupation ratio differs from its asymptote by 7.63e-5 (registered tolerance 1e-3).
- Separate opposite-jump leading-tail control yields logarithmic coefficient ratio 1.00006738 (tolerance 1e-3). This analytically integrates leading amplitudes; it does not evolve archived high-k modes.

## Post hoc smooth-width calculation

For the separate exactly solvable tanh step, the late static-region particle energy is finite at each positive width. Across epsilon=0.1 to 0.0003 it increases from 0.02608656 to 0.18954231. The last adjacent-width logarithmic slope is 0.0284962622 versus the analytic coefficient 9/(32 pi²)=0.0284965829, a relative difference of -1.1254e-5. The largest tight-tolerance/extended-band change is 1.39e-16. This refinement agreement is not an interval bound.

This calculation confirms the nonuniform sharp-transition limit in the toy model. It selects neither a physical HDBLAST smoothing width nor a shell energy budget. It was added after the primary outcomes and has no retroactive registered acceptance criterion.

## Scientific verdict

CONTROL PASS; ARCHIVED ABSOLUTE CURVED SOURCE NOT ESTABLISHED. A fixed local vacuum subtraction cannot cure a history-dependent excitation tail from an inadmissible unchanged state/interpolant. A smooth controlled background, admissible state and explicitly matched curved action are required before coupling the source to the shell.

## Supplemental fidelity and archive interference checks

The explicit nodal-fidelity diagnostic passed execution in [36959791142](https://github.com/maldonado-research/HDblast-archive/actions/runs/36959791142), source 84d33b785edf351426adc083510c4c4c923445cd. Maximum stored-value discrepancies over all eight reconstructions are 6.66e-16 for ln a and 4.44e-16 for phi. Native SciPy and independently implemented 80-digit de Boor evaluations were compared at 51 probes per smooth reconstruction, including both endpoints. Fine ln a fourth-derivative differences between evaluators are at most 2.37e-8 for C4 and 1.27e-8 for C6, far below the 4.48235 difference between the two reconstruction choices. Higher derivatives are more cancellation-sensitive: the largest seventh-derivative scaled discrepancy is 1.81e-7. Native evaluation is not being claimed to provide 80-digit field accuracy.

The smooth value-only fits do not preserve the archived phi_dot(0)=0; their roughly 0.0041 maximum derivative discrepancy is about 0.31% of the full-history velocity scale. This and the initial-boundary maxima motivate a physically constrained prehistory/reconstruction, not automatic preference for the highest spline degree.

A supplementary leading-amplitude archive interference calculation uses the fine 299 jumps, minimum conformal knot spacing 0.003742357676, and K=[2.672112306e5,2.672112306e8]. The coefficient ratio to the diagonal logarithm is 1.00000295920, within the original 1e-3 target. The JavaScript calculation evaluates the large-argument Ci expansion (minimum argument 2000); an independent SciPy sici cross-check is supplied and executed in the final replay. The scale multipliers 1e3 and 1e6 were selected after the primary audit and are explicitly supplementary. This completes the archive-spacing form of the proposed check without identifying it with exact high-k mode evolution.

An independently derived phase-independent cross-term inequality also yields an absolute interference bound of 1.21136e-7, or 6.20724e-4 of the diagonal logarithm, below the 1e-3 target. It uses |Ci(z)|<=2/z and the saved coefficients/times; it bounds interference within the leading-amplitude model, not physical UV corrections. The explicit oscillatory integral and this bound were independently checked by the reviewer.
