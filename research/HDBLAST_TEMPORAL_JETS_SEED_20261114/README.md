# Nonlinear initial temporal jets and an evaluated compatible seed

14 November 2026 UTC. Conditional mathematical checkpoint for HDBLAST,
continuing the [completed Cauchy-data bound](../HDBLAST_UNIFORM_CAUCHY_DATA_20261113/README.md)
at immutable base `2592e90650cb440cab2ce74f726e3dcf9b9d1247`.
The same actual background, coefficient chart, selected canonical family,
original parameter and full reflection-invariant angular space are retained.

The exact reduced PDE now gives a triangular nonlinear initial hierarchy
`J_0,...,J_14`. An explicit regularized heat extension matches all fifteen
jets in the strong mixed Sobolev topology and supplies a compatible auxiliary
starting field with fixed zero endpoint collars. This closes the numerical
initialization interface; the seed does not solve the nonlinear PDE away
from the initial section.

| Map or restriction | Evaluated bound |
|---|---|
| Initial pair size | `d(a)<=min(10^-2200,10^4200 a)` |
| Nonlinear initial jets | `I14=sum(j=0..14) norm(J_j,H^(14-j))<=C_J d(a)`, `C_J=10^6503000` |
| Linear jet-table lift | `sup_s X14(seed;s)<=C_lift I14`, `C_lift=10^255` |
| Completed initial pair to seed | `C_seed=C_J*C_lift=10^6503255` |
| Selected original parameter interval | `abs(epsilon)<=a_seed=10^-6507485` |
| On that interval | `I14<=10^-285`, `sup_s X14(seed;s)<=10^-30`, `norm_star(seed)<=10^-23` |

The literal mixed norm is

```
X14(W;s)=sum_(k=0)^14 ||partial_s^k W(s)||H^(14-k).
```

It uses `Sigma=[0,6ell] x S3`, `ell=10^-6`, and the full
`Gamma=(Z2)^4` invariant space. The result covers one copy and the retained
identical-two-copy Hilbert union. In the latter case I14 already uses the
global union norms, so the linear lift adds no copy factor.
The prescribed past elapsed slab has `0<=s<=13ell/8`. The extension exists
for all auxiliary real time, but this asserts no nonlinear solution lifetime.
The interval is deliberately smaller than the preceding initial interval;
the exponents are conservative proof bounds in the registered normalization,
without a measured physical amplitude or energy interpretation.

[Nonlinear recursion and initial-jet proof](NONLINEAR_INITIAL_JETS.md),
[all-state coefficient/product ledger](MIXED_INITIAL_JET_PRODUCTS.md), and
[compatible heat-lift proof](QUANTITATIVE_COMPATIBLE_SEED.md) supply the
analytic arguments. They retain the original lapse, half-action matter
normalization, fixed density, all mixed metric velocities, same-state time
reflection and all symmetry-allowed modes. The nonzero-state horizontal
coefficient derivative includes its state-transport term.

The lift uses compatible DD/ND spectral domains, evaluated covariant/graph
norm comparisons, finite positive heat counterterms, strong traces, a flat
cutoff and natural symmetry averaging after the lift. The positive
counterterms are finite differential expressions, never a backward heat
operator applied to rough data. Smooth compatible approximants smooth the
initial pair and recompute its nonlinear jets before lifting.

[Reproduction and fault controls](REPRODUCIBILITY.md), [audit scope](AUDIT.md),
[primary literature and search limits](LITERATURE.md), and [next calculation](CONTINUATION.md)
separate finite replay from analytic proof. No nonlinear evolved solution,
observer tube, evolved C2 family, calibrated observable remainder or empirical
higher-dimensional event follows. This working GitHub checkpoint creates
no new Zenodo edition.
