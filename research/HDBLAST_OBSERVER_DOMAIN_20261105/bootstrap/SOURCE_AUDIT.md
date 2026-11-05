# Public source audit and remaining input boundary

5 November 2026, UTC. This audit states the public source interface for the
conditional finite-dimensional observer bootstrap. The continuous proof
and exact checker are unchanged by this audit. The source identities below
are byte-integrity pins and immutable source locations, not signatures or
a new verification of the retained analytic theorems.

The public source checkpoint is commit
[`21b652a69f5b2c49b24cec17971a6d6b44d4afdd`](https://github.com/maldonado-research/HDblast/tree/21b652a69f5b2c49b24cec17971a6d6b44d4afdd/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105).
All links below identify files at that immutable commit.

| Public source | Interface used in this note |
| --- | --- |
| [Observer-flow budget](https://github.com/maldonado-research/HDblast/blob/21b652a69f5b2c49b24cec17971a6d6b44d4afdd/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/observer/OBSERVER_FLOW_BUDGET.md), §§1–4 | Assumes a common compact proper-clock tube, distinguishes complete physical shell pullbacks from a coordinate restriction, fixes the past signs, charges scaled norms and chart Hessians, and proves the conditional V/W budgets. |
| [Source dependencies](https://github.com/maldonado-research/HDblast/blob/21b652a69f5b2c49b24cec17971a6d6b44d4afdd/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/SOURCE_DEPENDENCIES.md) | Explicitly records genuine uniform first B13 and second B12 variation bounds and C^2 dependence into X10, while leaving their numerical constants, the reflected past-domain application, full-field norm comparisons and observer-domain inputs unresolved. |
| [Variation completion](https://github.com/maldonado-research/HDblast/blob/21b652a69f5b2c49b24cec17971a6d6b44d4afdd/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/inputs/FINITE_VARIATION_COMPLETION.md), §§2–4 | Uniform B13/B12 bounds, actual variation limits and C^2 dependence into X10 with derivative-loss accounting. Constants depend on the common nonlinear X14 bound, fixed interval, hyperbolic margins and initial-curve derivatives; observer coefficient constants are not evaluated. |
| [Exact quasilinear family](https://github.com/maldonado-research/HDblast/blob/21b652a69f5b2c49b24cec17971a6d6b44d4afdd/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/inputs/EXACT_QUASILINEAR_FAMILY.md), hypotheses and theorem | Conditional finite-regularity exact evolution on a prescribed slab, regular state and hyperbolic margins, and an unevaluated sufficiently small parameter interval. |
| [Nonlinear bridge](https://github.com/maldonado-research/HDblast/blob/21b652a69f5b2c49b24cec17971a6d6b44d4afdd/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/inputs/NONLINEAR_BRIDGE.md), §§1–4 | Requires the registered past-domain application; fixes the proper-clock actual-area observable with complete observer pullbacks, the old coarse endpoint 6145 ell/4096 and coordinate slab 13 ell/8. The old 106 guard belongs to its retained coarse bracket. |

The new step is conditional finite-dimensional continuation: bounds over
every convex fiber ball in a fixed coefficient neighborhood of a known
background imply nonexit and continuation of perturbed trajectories over
a prescribed duration. Only after that step are C^2 parameter-flow
estimates applied. This does not reduce the retained uniform analytic
variation theorem to tangent data and does not add a nonlinear PDE theorem.
Full-parameter bounds, compact regular-domain margins, complete initial
normal and compatible Jacobi data, scaled norms and coefficient extensions
remain explicit premises.

The public source-dependency note supplies the source-context summary used
here. It distinguishes the existing conditional analytic theorems from
the missing evaluated joins. The observer coefficient constants, their
common parameter radius, and the full pullback/norm comparisons for the
registered physical past family remain unevaluated. A parameter-zero
Taylor coefficient or a first-field gradient estimate cannot replace those
uniform inputs.

Physical embedding and injectivity remain separate from this ODE result.
The fixed auxiliary ambient representation does not prove that the actual
observer map is globally embedded. A coefficient-only extension near t=0
does not assert a PDE solution beyond the physical past slab: the time
guard and strict monotonicity must show that the constructed trajectory
stays on that slab. No effective numerical amplitude is claimed.

`check_flow_bootstrap.py` binds this audit, the unchanged proof and its own
source in each canonical report. Its 383 finite Fraction checks cover the
algebraic companion identities. They do not establish analytic Gronwall,
local ODE existence, compact continuation, extension hypotheses or the
physical/PDE inputs. The continuous argument uses the standard
finite-dimensional facts stated in the proof.
