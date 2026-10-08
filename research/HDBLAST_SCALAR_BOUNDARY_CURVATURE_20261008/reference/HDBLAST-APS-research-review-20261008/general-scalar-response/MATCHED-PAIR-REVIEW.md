# Matched bulk/tension variations

8 October 2026. Twelve additional radial shooting integrations completed; all recorded solver flags and numerical residual checks pass. The original family_benchmark.py and family-benchmark-results.json were not changed. matched_pair_benchmark.py imports the original solver, so this is a controlled model-variation test, not a second independent software implementation. Files and source bindings are in matched-pair-results.json; recalculated diagnostics are in matched-pair-review-summary.json.

Starting from registered k,w,f0,f1, the first variation changes only v=W''' from 2 to -3; the second changes only f2=f'' from 0 to 1.2. Each uses delta=.002,.001,.0005 with both original numerical settings.

The shared predicted second-order scalar coefficient is -0.0009299992353414982. Both variations approach it with coefficient-error halving ratios:

| Variation | .002→.001 | .001→.0005 | Relative error at .0005 |
|---|---:|---:|---:|
| Bulk cubic only | 2.000092 | 2.000029 | 0.12085% |
| Tension curvature only | 1.999907 | 1.999909 | 0.11652% |

Their differences from the baseline enter at the next order. The measured (H²_variant-H²_baseline)/delta³ values are:

| delta | Bulk cubic only | Tension curvature only |
|---:|---:|---:|
| .002 | 7.68517635e-5 | 1.57712207e-4 |
| .001 | 7.66673441e-5 | 1.57323902e-4 |
| .0005 | 7.63550193e-5 | 1.56880991e-4 |

This supports absence of v/f2 from second order while resolving nonzero higher-order changes. It cannot establish exact independence or a uniform remainder theorem; the exact derivation supplies the independence statement. Subtraction of nearby floating-point solutions and division by delta³ amplify the numerical floor, particularly for the smallest delta. Do not treat the last digits of these differences as certified third-order coefficients.

A separate local formal expansion checks the expected higher-order sensitivity. With lambda=w-4k, write p=lambda*eta+A*eta²+B*z*eta+..., z=1/rho². The radial equations give A=v*(3w/2-2k)/(3lambda+4k), B=-lambda/[k(lambda+k)], and a=-f1/[2(lambda+w)]. If eta_b=a*delta+b*delta²+..., then b=[-(A+v/2)*a²-B*(k*f0/3)*a-f2*a/2]/(lambda+w). Substituting this into the exact boundary identity gives:

    d h3/d v  = -k*f1³/[96*(3w-8k)*(w-2k)²],
    d h3/d f2 =  k*f1²/[96*(w-2k)²].

The corresponding formal differences are 7.64551451e-5 and 1.56937371e-4, consistent with the finite-detuning table. These derivatives assume the same local growing-branch expansion and do not add an existence claim. No external priority claim is made.

The original read-only review remains valid as historical evidence. Its unmatched-variation limitation is now addressed by this separate supplement. Remaining limitations still apply: exploratory parameter selection, positive detunings only, sampled residual checks, shared binary64/SciPy implementation, and no interval proof, stability, attraction, observational confirmation or external novelty assessment. The combined inventory is 48 radial integrations, of which 42 fit shooting parameters; the six original zero-coupling controls use the exact constant branch.
