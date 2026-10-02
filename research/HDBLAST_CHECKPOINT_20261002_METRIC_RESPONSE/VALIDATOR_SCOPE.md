# Prospective metric evidence validation

The four source files in `code/` are preparation artifacts for root integration
and public registration. They import no physical producer and perform no
import-time response calculation. No physical mode evolution, registered source
sample, response quadrature or tail evaluation was run in this stage.

`validate_metric.py` accepts explicit primary and independent result paths,
the registration SHA-256 and public freeze commit, and a fresh external output
directory. Normal and optimized Python use the same explicit guards. It checks
the unchanged prospective gates and requires the same frozen configuration,
source files, independent manifest, state, runtime, arrays and resource budgets.

The comparison coverage is 180 core quantities, 144 physical baseline/current
quantities, 180 core refinements, 36 current refinements, 72 raw observation/K
points, 72 reconstructed Ward endpoints (36 fine endpoints gated), 36 Ward
refinements, 432 saved-memory/local-contact reductions, 48 differentiated primary
Ward identities, and 324 removed-band comparisons. Density and pressure are
separately evaluated from physical operators; the Ward ledger defines neither.

`raw_metric_audit.py` independently rebuilds every saved observation from complex
unconstrained modes. It checks full metric D-operator terms, opposite explicit
mass contacts, variance/stress prefactors, and complete grade-2/4 subtractions
with the separately authored generic Taylor module in `theory/metric_wkb.py`.
Its stable baseline implementation evaluates explicit polynomial power sums;
the producer uses a separately implemented Horner grouping. Three exact
symbolic identities and six mutations compare those rationalizations against
the complete generic WKB expressions. Raw literal bare/subtraction pieces are
retained and checked with a separately declared leading-term ulp allowance.
Every baseline history is also checked against displayed finite-band primitives.

Raw LD-to-LD comparisons use `1e-12 + 512*LD_epsilon*scale`. A separately named
JSON comparator uses `1e-12 + 16*binary64_epsilon*scale`, allowing float
serialization and short metadata arithmetic. Neither changes the physical
cross-route, refinement or Ward acceptance gates. The exact synthetic large
binary64 roundtrip test demonstrates why these comparators are separate.

The directed metric-tail module encloses only omitted UV bands. Its conservative
derivative caps are loose. Primary quadrature estimates, independent refinement
estimates, cross-route allowances and arithmetic allowances remain separately
reported empirical diagnostics; their sum is not a certified total error.
No precise continuum pressure/sign, actual shifted-root propagator, general
lapse response, bulk closure, coupled evolution, stability or heating follows.

Ten saved-mode/history sensitivity controls and five synthetic metadata controls
are reported with their actual residuals and operational-gate status. A nonzero
witness is not automatically an operational failure and is not an alternate
physical-state experiment. The separate guard suite passes 66 synthetic algebra,
format, provenance and static checks in normal and optimized Python.

Actual stdout, exit codes and source hashes are in `evidence/`. Earlier passing
63-check captures are preserved in `development/`; their full source snapshots
were not preserved, and they are not the authoritative final guard evidence.
