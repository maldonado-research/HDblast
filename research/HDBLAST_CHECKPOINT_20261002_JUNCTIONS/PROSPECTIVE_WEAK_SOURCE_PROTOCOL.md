# Prospective first weak-source experiment: choices to freeze before execution

Status: design and readiness gates only. This is **not a completed registration**
and contains no executed physical source or shell experiment. No new numerical
parameter, physical prior, state or cutoff is invented here. Its purpose is to
make a bounded next calculation registerable after the missing model data are
selected and the associated constrained initial data are constructed.

The question is: does a short weak-source addition of the common-action
quantum field to a specified existing HDBLAST branch preserve the bulk
constraints, all three junctions and the signed work/bulk/Weyl budgets at the
expected convergence and perturbative order? A PASS would validate that
coupled implementation on its declared interval. It would not establish a
radiation epoch, substantial energy conversion, decay, thermalization or an
observational fit.

## 1. Exact inputs still missing

| Item | Available inherited information | Must be fixed before a physical run |
| --- | --- | --- |
| Classical branch | Same W and U; original d=0 and A1/B1 tuned d are distinguished | Select the branch, delta,c,d and renormalized tension meaning explicitly; do not merge the two historical models |
| Dimensionful mass and coupling | x=m0²+G²(phi−phi_star)²; B1 used m0=0 and a prior G_hat/phi_star grid | Select physical m0,G,phi_star,H0/M5 and multiplicity within that model; old grid entries are historical choices, not measured priors |
| Finite action | Frozen reference prescription specifies one finite matching convention | Declare whether its V/F/alpha/beta/gamma conditions define this model or supply alternate matching; pin the retained EFT operator basis, and give each retained term one location and coefficient, including any induced Einstein or scalar-derivative term |
| Running/reference scale | Constant positive r belongs to subtraction; physical x may touch zero | Select r as a fixed matching input and define any scale-running rule by the action; r must not track x or the numerical cutoff |
| Physical EFT domain | B1 used max(m_chi,sqrt(q)) relative to M5 | Declare the physical EFT scale and relevant dynamical mass, curvature, scalar-gradient and production-momentum screens before the run |
| Initial quantum state | Prescribed control had exactly static smooth past; inherited de Sitter shell expands | Construct a full Gaussian state with admissible UV asymptotics and stated infrared data on the actual shell; positive-reference subtraction is not the state |
| Initial bulk/shell data | B1 used a neighboring-tension seed that is off the final junction by O(dc) | Solve the constraints and sourced junctions for the chosen state and action before the starting slice; specify initial quantum-state and geometric perturbations together |
| Higher derivative evolution | Pressure/subtraction needs a''''; source depends on acceleration; differentiated p can need a''''' | Specify a constrained enlarged system or justified order reduction, its independent initial data, causal response, retained order and remainder; apply it to rho,p,J_phi together |
| Boundary treatment | Bulk doubled/GHY signs are fixed | Pin the far-boundary causal separation, initial/final variational data and any quantum state/surface counterterms |
| Discretization and tolerances | Existing classical and prescribed-mode code provides separate controls | Freeze physical interval, base grids, source iteration, refinements, absolute scales and gates in a new hashed registration before execution |

Conservation does not fill any finite matching entry in this table. Numerically
small bare/renormalized vacuum stress does not determine an unknown induced
Einstein coefficient. A fixed reference convention is an explicit model
definition, not a fitted or empirical measurement.

## 2. Minimal useful order of work

1. **Stationary initialization bridge.** Use an existing exact static bulk
   branch with its induced de Sitter shell. For constant x_b>0, a declared
   Euclidean/Bunch–Davies Gaussian state is a possible initialization route,
   because its quantum source is de Sitter invariant, p_Q=−rho_Q. Compute that
   source with the same reference action and solve the sourced static
   boundary-value conditions, rather than imposing the old vacuum shell.
   A massive de Sitter state is a proposed choice, not an already evaluated
   source. A minimally coupled x_b=0 state requires additional infrared
   information; do not extend the massive Euclidean state through that point
   without examining its zero mode. Another state/prehistory is acceptable
   only when its UV and IR properties and boundary data are supplied.

2. **Constrained smooth disturbance.** For the selected initialized state,
   construct a short smooth family of initial data satisfying both bulk
   constraints and all sourced junctions. Its disturbance must derive from
   that pinned action/state and the selected branch; this note supplies no new
   physical shell-mode spectrum. The quantum state perturbation must be
   consistent with the same initial slice. The old neighboring-tension seed,
   which fails the target junction at O(dc), cannot by itself establish this
   prerequisite. Do not use the old nonsmooth time-series interpolation or
   the prescribed control's a=exp(A B) history as the evolved geometry.

3. **Short first-order response.** Evaluate one interval with a nonzero,
   resolved J_phi v signal, away from an unresolved chart endpoint or EFT
   boundary. A formal constant epsilon multiplying the entire matched quantum
   action is useful to derive/check first-order response:

       S(epsilon)=S_classical+epsilon[Gamma_Q,ref+DeltaGamma_loc],
       bulk/shell = classical + epsilon response + O(epsilon²).

   Epsilon here is a mathematical loop/source bookkeeping variable. It is not
   a new physical multiplicity, coupling, parameter prior or reheating model.
   It must multiply rho,p,J_phi and every moved local action term together.
   A physical run has the selected action with epsilon=1 and is weak only when
   its actual kappa5²-scaled source is small. For a controlled first-order
   calculation, varying the leading-background quantum source on the O(epsilon)
   corrected geometry first affects the equation at O(epsilon²); this may
   justify a one-pass linear-response experiment. It does not justify claiming
   a fully nonlinear self-consistent source. If the retained order includes
   the response of the quantum state, register its causal response operator.

4. **Compare with the selected coupled/EFT treatment.** A short iterated or
   enlarged-system calculation may be compared with first-order response only
   after its derivative-order/state prerequisites are met. Use the discrepancy
   expected from the stated truncation and independent discretization error.
   A stable result after deleting terms is not evidence for the original
   action. Extending beyond weak response is a later registration.

This ordering allows the static bridge and algebra/code adapter work to proceed
without inventing a physical perturbation, scanning unselected parameters or
launching a new shell trajectory. A stationary correction by itself is not
evidence of energy transfer; the time-dependent stage must resolve nonzero work.

## 3. Data to preserve and prospective acceptance structure

Save all three renormalized sources from the same mode/state prescription,
including coherent pressure/current, raw normalized Wronskians, exact cutoff
label, local action terms, subtraction jets, bulk fields/constraints and all
three signed junction residuals. Retain the matching/placement ledger,
physical EFT screens, source-update convergence and all endpoints. Record
the subtraction reference separately from the physical EFT scale and numerical
mode integration limit. A large numerical subtraction cutoff is not itself a
physical gravity scale; assess actual resolved dynamics and the declared EFT
matching domain separately.

Prospectively freeze the tolerances and comparison scales after the constrained
initialization construction and resource estimate, before physical execution.
The prescribed benchmark's finite-cutoff tolerances do not automatically apply
to the dynamic bulk source. A finite refinement comparison is not an
infinite-cutoff error bound.

The registration should require all of the following, with distinct outcomes:

- The source-free control reproduces the chosen classical initial-value
  problem and its bulk/junction budgets within the registered numerical gate.
- Initial Hamiltonian and momentum constraints and the three junctions satisfy
  their own initial gates; a wrong initial slice must fail before evolution.
- Renormalized quantum work and full tension-plus-matter/bulk flux agree over
  the same interval and improve under source/time and momentum refinement.
  Record bare work and subtraction work as separate diagnostics.
- The interior bulk constraints, each of A/B/phi junction residuals and Weyl
  balance improve under independent bulk-grid, timestep, mode-band,
  momentum-quadrature and source-update refinement. A correct algebraic shell
  Ward identity does not excuse interior constraint loss.
- The response is weak in absolute kappa5² source scales relative to the
  selected background curvature/extrinsic scale; avoid dividing by a vanishing
  sigma_phi or H. If formal epsilon refinement is used, separate its expected
  truncation order from discretization and physical weak-source validity.
- Deliberate wrong J_phi sign, omitted J_phi, incorrect scalar junction half
  factor and tension-only A=B data fail their respective budgets on the same
  nonzero-work diagnostic interval. Preserve those failures without evolving
  them into physical claims. Include the complete action-counting ledger,
  because a coherently doubled local source can still pass Ward checks.
- Compare the stated EFT/order-reduced system with controlled first-order
  response at the retained order. A numerical mismatch or missing causal/state
  input is an unresolved outcome, not grounds to change the registered gates.

Stop and report a failed constraint/junction gate, unresolved state singularity
or turning-point behavior, source iteration failure, chart boundary,
unresolved UV pressure/current, or departure from the declared EFT domain.
Keep completed raw data, failure timing and reason. A negative or inconclusive
result is a valid next outcome. No occupation rescaling, fitted tail, patched
radiation fluid, phenomenological Yv term or post hoc threshold change may
turn this experiment into a PASS.

## 4. Scope and reproducibility

Only the prerequisite algebra was executed in this package:

    python code/verify_action_junctions.py --output evidence/ACTION_JUNCTION_CHECKS.json
    python -O code/verify_action_junctions.py --output evidence/ACTION_JUNCTION_CHECKS_OPTIMIZED.json

The physical registration and all physical stages above remain UNRUN until
their exact inputs are supplied and pinned. No physical source magnitude,
conversion efficiency, radiation ratio, matching measurement or stable
coupled solution is inferred from the algebra PASS. Earlier B1 INCONCLUSIVE
and smooth-FRW original FAIL/separate follow-up PASS remain part of the
scientific record.
