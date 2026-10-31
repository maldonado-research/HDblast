# Complete absolute displacement cap refinement for the same exact packet

The completed `W1_exact_packet` function is unchanged. Refining its interval
certificate at 128 bits reduces the complete absolute initial `u_rr` cap by
an exact factor greater than 3.38. The `u` cap is identical and the `u_r` cap
is smaller. Every paired initial target-error charge, the smoothed P1
velocity, formal target `U`, original `F/F_x/G/M`, physical-side domain and
positive all23 angular norm remain frozen.

| Complete H upper cap | Published initial report | Refined certificate |
| --- | ---: | ---: |
| `u` | 5.343420368515216e31 | 5.343420368515216e31 |
| `u_r` | 3.1953900575389437e33 | 2.6981603984020674e33 |
| `u_rr` | 3.1768732556099413e35 | 9.398955573842927e34 |

The displayed decimals describe exact rational entries in `REFINEMENT.json`;
acceptance and every comparison use the exact rationals. The exact ratio of
old to new `u_rr` caps is approximately 3.380027951670615. This is a reduction
of an upper certificate. The same candidate's unknown true error is unchanged.
The result does not assert that the evolved error certificate halves.

## Frozen experiment and authentication

`FREEZE.json` was written before any new cell evaluation and is pinned by
`auth.py`. It fixes the same candidate, a first 128-bit campaign, at most
4096 total leaves, depth at most 16 relative to the original seeds, no worse
`u/u_r/u_rr` caps, a strict `u_rr` reduction and a desired factor of at least
two. The total leaf budget includes the collar and tail and is stricter than
a separate allowance per seed. The 192-bit fallback was unnecessary and was
not executed.

Before importing any inherited evaluator, the checker authenticates the
closed published exact-packet capsule against its pinned manifest, including
its exact file set, and authenticates every selected inherited manifest and
member named by the published initial freeze. The initial report, producer,
original freeze and full stored velocity profile have explicit pins. The same
source closure is authenticated again after calculation. No source file or
Git state was modified.

## Complete interval proof

The published partition has 49 cells: the zero collar `[0,1]`, analytic tail
`[1,65/64]`, and 47 consecutive active cells through `7/4`. The authenticated
published component tables supply valid parent caps. Only three cells need
bisection to meet the frozen goal:

* `[17/16,69/64]`, midpoint `137/128`;
* `[69/64,35/32]`, midpoint `139/128`;
* `[35/32,71/64]`, midpoint `141/128`.

The result has 52 leaves, three splits and depth one. All 55 tree nodes,
including all 49 original roots and six new children, are numerically
re-evaluated by the checker. The final maximum `u_rr` cap occurs on unchanged
cell `[71/64,9/8]`. Both collar and exact analytic tail remain in the complete
coverage. The tail uses the inherited validated exponential/polynomial
inequalities through the required packet derivatives, never an omitted cell
or a raw denominator at the flat support edge.

At each node the evaluator reconstructs all 23 completed registered nominal
displacement components, including the five mean overrides, using the exact
analytic `q` and `E` with outward exponential enclosures. It retains the
nominal background coefficients and radial jet scaling from the pinned
producer. No exponential midpoint defines the candidate; no unknown interval
remainder is differentiated or discarded.

For each derivative and component, a child's raw cap and its parent's retained
cap both bound the same function on the child's exact subset. Their minimum
is therefore a valid child cap. The checker independently reconstructs every
binary midpoint, parent/child interval, depth, componentwise minimum and
connected tree. It verifies the ordered leaf endpoints form exactly `[0,7/4]`.
No uncovered interval, overlap of interiors, disconnected node or omitted tail
can be accepted.

For every derivative table, the all23 H cap is the outward dyadic square root
of the sum of the 23 positive fixed Gram weights times the squared component
caps. The checker implements this rational aggregation directly, independently
of the producer's norm helper. It verifies 165 such node aggregations and
then takes the maximum over the complete leaf cover. The completed quadratic
extension and its exact value/r/rr matching at `7/4` are unchanged; this initial
physical footprint ends at `7/4`. Global C3 regularity is not asserted.

## Reproduction and controls

The assertion-free checker uses explicit runtime rejection and standard-library
arithmetic. Both the producer and full numerical checker pass under ordinary
and optimized Python, with byte-identical corresponding reports.
`IDENTITY.json` records the digests. `CHECK_NORMAL.json` records the full
55-node replay; `controls/CONTROLS.json` records 14 genuine corruption cases
executed in both modes, giving 28 rejected subprocesses.

The controls mutate exact interval endpoints to create gaps or overlaps,
remove a component or the tail, alter a label or Gram weight, violate parent
inheritance, understate local or global norms, change frozen initial-error or
velocity charges, and corrupt source formula bytes or the closed source tree.
One control halves a nonzero raw completed `w:0` rr component and recomputes
its dependent local minima and H norm; it passes the structural relationships
but fails the actual numerical re-evaluation. Every corrupt run exits nonzero
with an explicit failure witness and produces no successful checker report.

Run the portable integration replay from a complete checkout with a fresh
output directory:

```sh
python -I -S -B run_refinement.py --repo-root /path/to/HDblast \
  --output-dir /tmp/hdblast-absolute-refinement-fresh
```

The seven adjacent files needed by this runner are `run_refinement.py`,
`auth.py`, `refine.py`, `check_refinement.py`, `run_controls.py`, `FREEZE.json`
and `REFINEMENT.json`. Relocation changes no authoritative path or report
content. The native source-root argument may point at the authoritative
checkout; none of these scripts writes there.

## Scope

This calculation tightens complete absolute initial displacement allowances
for a fixed auxiliary finite approximation. It retains the full published
initial error ledger verbatim; it does not recalculate or reduce any paired
error channel. A separate response adapter can use these smaller absolute
charges while proving every retained response identity and correction term.

No new natural-equation or constraint residual, shell point value, invariant
observable, physical relative error, cosmological realization or observational
fit follows from this certificate. Existing conditional analytic initial-data
and precontact results retain all their hypotheses. Full finite-amplitude
nonlinear evolution through shell contact remains unestablished.
