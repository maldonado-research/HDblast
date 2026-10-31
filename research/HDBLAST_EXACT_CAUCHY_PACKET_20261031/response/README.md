# Same-function exact response certificate refinement

The candidate is still `W1_exact_packet`. Its mathematical function and true
error against the original actual target `U` are unchanged. This experiment
tightens absolute displacement allowances and propagates them through the
existing conditional physical-side coordinate-error estimate. It produces no
new waveform, point/quadrature values, natural or constraint residuals,
invariant observable, or nonlinear shell-contact result.

The exact target-error caps come from the frozen published `INITIAL.json`.
Both absolute velocity caps, both velocity-error caps, `F`, original unfiltered
`F_x`, `G`, `G_x`, `M`, all operators, the domain, and the full 23 `H` norm remain
fixed. The adapter first reconstructs all 35 published nested C1/C2 ledger
nodes and every constant, with `fractions.Fraction`. It imports no scientific
implementation. Assertions are not used.

The complete differentiated ledger retains `f00`, `r00`, `M*Bs+G_x`, the
coefficient derivative commutator, both feedback inverses, the original `F_x`,
and the cancellation-preserving `r0x`. The correction has zero first Cauchy
jet, so `r1(0)=0`; no fictitious acceleration charge is introduced there.
Separate absolute `W1` prescribed-forcing allowances are also reconstructed.

All 35 refined ledger nodes are no worse; 32 shrink strictly and 3 stay equal.
The exact arithmetic proves each second-error upper certificate is at most
5/8 of its published value and the C1 upper certificate is at most 6/7.

| Conditional upper certificate | Published | Refined | Old/new |
| --- | ---: | ---: | ---: |
| C1 X | 1.60214883e23 | 1.36453813e23 | 1.17413 |
| xx and xr H | 2.13407313e25 | 1.32984325e25 | 1.60476 |
| rr H | 2.14322475e25 | 1.33757082e25 | 1.60233 |
| Absolute W1 rr H | 3.17812330e35 | 9.40271305e34 | 3.38001 |

Exact rationals and comparison inequalities are in
`reference/RESPONSE.json` and `reference/IMPROVEMENT.json`. These are upper
certificate improvements, not a new true-error comparison. Physical-side C2
estimates retain the published source/coefficient regularity, compatible
corner, and one-sided shell-trace premises.

`checks/FREEZE.json` pins the original sources and the arithmetic recipe.
Every selected source, imported proof source, and absolute-cap report is
hashed before and after. `reference/CONTROLS.json` records 62 rejected fresh
normal/optimized subprocesses: 8 effective formula corruptions in each mode,
plus frozen-error and velocity changes, worse/missing caps, source changes,
report-node changes, table-integrity changes, recipe/freeze changes, and a
symlinked source ancestor. Table-integrity controls reject changes against the
imported report's frozen hash; interval-table validity belongs to its separate
producer and independent proof checks.

Run with Python3 standard library only:

```sh
python3 checks/run_replay.py \
  --repo-root /path/to/published/repository \
  --absolute-dir /path/to/absolute-refinement \
  --output-dir /path/to/fresh/replay
```

The absolute directory contains the frozen `REFINEMENT.json`, `refine.py`,
`check_refinement.py`, `auth.py`, and `FREEZE.json`. The repository contains
the reports named by the source pins. The replay requires exact canonical
report bytes in normal and optimized Python, and reruns all 62 corruption
subprocesses. `reference/SENSITIVITY.json` separately records the hypothetical
zero-absolute-rr ledger calculation; it is neither a candidate nor a lower
bound on true error.
