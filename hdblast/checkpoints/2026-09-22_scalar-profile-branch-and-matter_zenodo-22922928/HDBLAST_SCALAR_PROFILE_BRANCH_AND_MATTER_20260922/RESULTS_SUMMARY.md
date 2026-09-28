# Completed calculations and their limits

Seven new finite-amplitude PDE evolutions are retained, including failed accuracy controls. The static branch has six detunings, two producer tolerances, and a separate 18-integration review. The matter extension is derived and screened, not evolved.

| Run | ε | h at shell | Stretch | L | Final t | Final φ_b | Max C_H/background scale over stored times | Max C_M/background scale over stored times |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| balanced_epsp01_h2 | 0.01 | 0.0002 | 0.05 | 6.0 | 2.5 | -0.03947126743 | 1.4053 | 0.69409 |
| balanced_epsp01_h1 | 0.01 | 0.0001 | 0.05 | 6.0 | 2.5 | -0.03953402834 | 2.228 | 1.1059 |
| balanced_epsp01_wide_h2 | 0.01 | 0.0002 | 2.0 | 3.0 | 1.0 | -0.002771327592 | 0.30286 | 0.15011 |
| balanced_epsp01_wide_h1 | 0.01 | 0.0001 | 2.0 | 3.0 | 1.0 | -0.002771030394 | 0.034033 | 0.016915 |
| balanced_epsp001_wide_h2 | 0.001 | 0.0002 | 2.0 | 3.0 | 1.0 | -3.615146718e-05 | 0.025596 | 0.012667 |
| balanced_epsm001_wide_h2 | -0.001 | 0.0002 | 2.0 | 3.0 | 1.0 | 3.306855042e-05 | 0.024588 | 0.012164 |
| balanced_epsp01_wide_h05 | 0.01 | 5e-05 | 2.0 | 3.0 | 0.5 | -0.00120805468 | 0.0065507 | 0.0032669 |

The denominator is 1+6Hc²+φs,z². The independent residual review also gives actual-term cancellation ratios. These are constraints, not error bars on φ_b. The monitor domain is z>−0.8L with the first six grid nodes excluded; a separate near-shell domain is z>−0.2. Maxima in this table cover stored times and can exceed the final-time values quoted in comparisons.

The corrected registered static rate differs from the old constant-scalar metric benchmark by -7.868927 parts per million. This small correction cannot account for the old late-time evolution gap of several percent.

All 26 completion/accounting checks pass. The global nonlinear-accuracy gate remains open. The JSON retains exact values, source hashes, first/final rows and maximum-over-time diagnostics.
