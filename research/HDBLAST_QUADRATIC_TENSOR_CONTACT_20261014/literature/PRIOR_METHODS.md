# Prior local methods for the shell-contact calculation

Inventory date: 14 October 2026 UTC. This is a read-only survey of useful
local sources in `/workspace/HDblast`, `/workspace/HDblast-archive` and
`/workspace/hdblast-research-work`. The inventory itself is new; the sources
were preserved. No archived scientific script, pickle, opaque cache, saved
trajectory or numerical solver was executed. No network retrieval or publishing
was performed for this subtask.

## Controlling current target

The current target is
`/workspace/hdblast-research-work/publish-live-readback-20261014/research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/NEXT_SCIENTIFIC_TARGET.md`.
It requires the total quadratic geometric-plus-matter bulk source for
`A=B=diag(1,-1,0,0)`, material-to-evolution-gauge interface conversion,
constraint/gauge-compatible boundary estimates, and an independently verified
second-order contact solution. A previously nonzero matter-only tensor
projection does not determine the total tensor source.

The corresponding material supplement, pinned as `C2` below, is the controlling
interface algebra. It retains two outward normal derivatives and **one** common
shell tension. With the derivative/factorial Taylor convention it gives

\[
\sum_i K_{i,2,ab}-\frac{\sigma_0}{3}h_{2,ab}
-\frac{\sigma'_0}{3}\psi_2h_{0,ab}
=\frac{2\sigma'_0}{3}\psi_1h_{1,ab}
+\frac{\sigma''_0}{3}\psi_1^2h_{0,ab},
\qquad
\sum_i p_{i,2}+\sigma''_0\psi_2=-\sigma'''_0\psi_1^2.
\]

For the registered tension, `sigma'''=4`, so the summed scalar source is
`-4 psi_1^2`. In the exactly reflection-symmetric even-master class, the
material first-order induced metric is conformal to `h0`. Consequently the
displayed metric shell source is pure trace. If
`h2_TF=rho^2 t2_TF`, the necessary trace-free boundary law is
`partial_r t2_TF=0` at the shell, rather than `partial_r h2_TF=0`.
This distinction must survive any variable or gauge change. These observations
do not provide the total tensor bulk forcing or a contact existence theorem.

## Methods that can be reused after fresh derivation

| Method and local sources | Actual reusable content | Limitation for the current calculation |
|---|---|---|
| Tensor linearization and warped Neumann law: `T1`–`T4`; closer archival variants in the companion archive note | Derive the full Einstein residual from a metric perturbation `delta g_TT=rho^2 t_TT`, then compute the Israel residual directly. The Sept27 script constructs its off-diagonal TT test component and obtains a factor `rho^2/2` times the normal derivative of the dimensionless tensor. This reproduces the rescaling logic in `C2`. | The older calculation is linear and largely separated. It does not compute `E''[h1,h1]`, the total second-order bulk forcing, general gauge boundary data, or the material/Eulerian conversion. Its test-component derivation is a useful check, not a substitute for covariant TT/angular-sector reconstruction. |
| Tensor Sturm–Liouville identity: `T1` | The separated equation is written `-(rho^4 t')'=mu^2 rho^2 t`. Integration by parts with regular cone data and shell Neumann data gives `int rho^4 |t'|^2 = mu^2 int rho^2 |t|^2`. The endpoint Green term must actually vanish. | This is a radial spectral identity. It is not an energy estimate for forced generalized-harmonic Einstein–scalar evolution, nor an all-sector stability theorem. No old eigenvalue, spectral gap or normalization is imported. |
| Independent tensor spectral diagnostics: `T3`, `T4` | Use a scale-invariant shell mismatch and a separate Prüfer-angle evolution for `(t,rho t')`. The latter avoids amplitude overflow and does not import the primary `gi_core.py`. Track the background first integral instead of integrating an unconstrained free radial variable near the cone. | Useful for checking a newly derived radial operator or resolvent. The numerical regular-cone cutoff, tolerances, branch data and spectrum belong to the old run. They do not certify the present contact solution. |
| Fixed-domain wave formulation and source energy: `S1`, `S2` | The Oct13 scalar calculation transforms the moving shell to a hemisphere, derives the transformed normal and boundary operators, retains the coefficient derivative in a dynamic boundary law, and obtains a compact-slab Gronwall estimate from time-dependent mass/stiffness forms. A source enters the energy identity through its actual bulk/boundary pairing. | Reuse the procedure and bookkeeping, deriving the tensor/Einstein–scalar coefficients afresh. The scalar master Wentzell theorem is explicitly insufficient for the full forced metric/scalar system. The auxiliary control energy is not a gravitational Hamiltonian or a ghost criterion. |
| Robin coercivity and trace pitfall: `S1`, `S2` | For a Robin problem keep the actual Robin potential in the stiffness form; use the epsilon trace inequality and a sufficiently large bulk mass shift to obtain coercivity. `S2` identifies why adding an arbitrary positive `boundary z^2` term and estimating its derivative from bulk `H1 x L2` data fails: that energy does not control `trace(u_t)`. | Boundary conditions derived from the action must be retained. A convenient positive artificial boundary term cannot establish the required Einstein–scalar boundary estimate. |
| Corner compatibility and characteristic forcing: `S1` | Recursive bulk time jets plus differentiated boundary equations expose the complete smooth corner hierarchy. A fixed finite spherical-harmonic class has an explicit characteristic/Volterra system coupled to a boundary ODE. Interior data vanishing in a shell collar satisfy every homogeneous initial corner relation. | The current inhomogeneous sources require corresponding boundary and source jets. Vanishing initial data alone cannot establish later source compatibility. Finite-mode scalar construction does not imply a full tensor/all-angular Sobolev theorem. |
| Matched two-copy graph geometry and shell current: `J1`, `J2`, `C2` | Compute unit-normal normalization, displacement evaluation, induced metric/scalar matching and graph curvature before reducing sectors. `J1` checks all linear tensor junction components, including trace and exceptional Hessian kernels. `J2` keeps Einstein divergence, GHY, moving-domain and single-shell terms together; dropping moving-domain terms reverses a boundary kinetic sign. | `J1` is a sufficient linear scalar reconstruction family; `J2` is a local quadratic kinetic identity. Neither determines an unrestricted second-order Einstein–scalar boundary domain. The current material supplement, rather than its older graph formulas alone, controls quadratic interface sources. |
| Source-ledger and Duhamel verification: `G1` | Propagate with a phase-aware forced Green/Duhamel integral, then verify an independently formed integral of the direct operator/current, retaining contact and work terms. `expm1` avoids loss in a short-step phase difference. Define the forcing and all derivatives from the same analytic source so a source-definition defect is visible. | `G1` concerns shell quantum modes on a fixed four-dimensional de Sitter/FRW experiment. Its masses, times, subtractions and finite-momentum formulas do not belong to the present five-dimensional contact problem. Only the verification pattern is reusable. |
| Radial Green-function/Wronskian construction: companion archive note | Select independent solutions satisfying the two endpoint conditions; form a weighted Wronskian; impose the derivative jump; then distinguish the resolvent from a causal retarded propagator. Conjugating a warped radial field can turn a physical Neumann condition into Robin data for the conjugated Schrödinger field. | The January RS2 notes use an AdS/Poincare branch, different coordinates/endpoints, and their own parameters. They also have sign and undo-rescaling cautions detailed in the companion note. Fresh derivation is required. Do not transplant `I2/K2`, continuum weights, suppression exponents or Newton constants. |

For an operator `A_lambda u=-(p u')'+(q-lambda w)u`, a useful convention is
`W_p=p(u_L' u_R-u_L u_R')`. The endpoint-matched radial resolvent kernel is
`G_lambda(y,y')=u_L(y_<)u_R(y_>)/W_p`, with `A_lambda G=delta(y-y')`;
the sign follows from the normal-derivative jump. For constant `lambda`,
`W_p` is constant. At an eigenvalue, its zero marks a pole and requires
explicit treatment of the kernel/compatibility. This generic construction
does not assert a current HDBLAST Green function has been solved. A causal
contact evolution additionally needs its correct gauge/constraint boundary
domain and retarded time prescription.

## Files and exact provenance

`PRIOR_METHODS_FILES.json` records absolute source filenames, SHA256 values,
byte lengths and read/line hints for the entries below. It also records
identical hashes where the September/October sources were found in both the
active local tree and the archive. Matching local copies are not independent
scientific replications.

| ID | Absolute source filename | Read hint |
|---|---|---|
| C1 | `/workspace/hdblast-research-work/publish-live-readback-20261014/research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/NEXT_SCIENTIFIC_TARGET.md` | Full target, especially steps 1–4. |
| C2 | `/workspace/hdblast-research-work/publish-live-readback-20261014/research/HDBLAST_MATERIAL_JUNCTIONS_20261014/MATERIAL_SECOND_ORDER_JUNCTIONS.md` | Sections 1–5, equations (2), (9)–(17). |
| T1 | `/workspace/HDblast/research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant/README.md` | Lines 45–98: linear tensor equation, junction and radial integration identity. |
| T2 | `/workspace/HDblast/research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant/derive_linearized.py` | Lines 357–378: actual TT metric test, bulk residual and shell curvature. |
| T3 | `/workspace/HDblast/research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant/verify/indep_spectrum.py` | Docstring and `tensor_theta`; independent Prüfer method and constraint-form background. |
| T4 | `/workspace/HDblast/research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant/gi_core.py` | Lines 1–17 and 168–187: operator conventions and tensor mismatch. |
| S1 | `/workspace/hdblast-research-work/publish-HDblast-20261014/research/HDBLAST_REGULAR_SHELL_CONTACT_20261013/existence/FINITE_SLAB_SHELL_CONTACT_EXISTENCE.md` | Sections 2–5: fixed hemisphere, corner jets, energy, characteristic/Volterra construction. |
| S2 | `/workspace/hdblast-research-work/publish-HDblast-20261014/research/HDBLAST_REGULAR_SHELL_CONTACT_20261013/reviews/ROOT_MATHEMATICAL_AUDIT.md` | Section 3: correct Robin coercivity and uncontrolled `trace(u_t)` pitfall. |
| J1 | `/workspace/hdblast-research-work/publish-HDblast-20261014/research/HDBLAST_REGULAR_SHELL_CONTACT_20261013/junctions/FULL_PACKET_JUNCTION_SUFFICIENCY_AND_CURRENT.md` | Sections 2–4: displaced graph and full residual; later current/action section. |
| J2 | `/workspace/hdblast-research-work/publish-HDblast-20261014/research/HDBLAST_REGULAR_SHELL_CONTACT_20261013/basis/prior-checkpoint/operator/ODD_KINETIC_COMPLETION_20261013.md` | Sections 3–5, especially boundary action and integration-by-parts endpoint. |
| J3 | `/workspace/hdblast-research-work/publish-HDblast-20261014/research/HDBLAST_REGULAR_SHELL_CONTACT_20261013/basis/prior-checkpoint/operator/PHYSICAL_OPERATOR_AUDIT_20261013.md` | Sections 1–3: action normalization, unreduced constraints and graph residuals before harmonic division. |
| G1 | `/workspace/HDblast/research/HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER/theory/ACTIVE_SOURCE_DERIVATION.md` | “Exact forcing and source work” and “Phase-aware forced propagation and a genuine integral check”. |
| X1 | `/workspace/HDblast/research/HDBLAST_CHECKPOINT_20261002_JUNCTIONS/ACTION_AND_JUNCTION_SPECIFICATION.md` | Sections 1 and 4: branch distinction and action-derived two-copy sourced junctions. This is an excluded model extension, included to document the mismatch. |

The companion `ARCHIVE_GREEN_METHODS.md` and `ARCHIVE_GREEN_METHOD_FILES.json`
pin the older archive paths, including files with spaces in their names, and
contain the additional September tensor and January RS2 read hints.

## Explicit exclusions

No earlier numerical background, radial endpoint, detuning, floating-point
`c` value, eigenvalue, tolerance, KK spectral density, RS2 curvature radius,
Planck normalization or quantum-subtraction parameter is adopted here. The
Sept27 scalar/tensor scripts discuss more than one branch and a finite cone
cutoff; the present calculation must use the exact registered background and
its proven local domain. In `X1`, the B1 quantum/matter adapter explicitly uses
a changed tension with a quadratic detuning term and a tuned branch, unlike
the controlling registered tension. It is useful for action/normalization
bookkeeping only.

An older linear scalar-master theorem, a radial nonnegative spectrum, a
homogeneous TT shell Neumann relation, a conserved current, or a valid
Wronskian formula does not settle the remaining coupled second-order contact
problem. The reusable material is the derivation and verification method,
with every current coefficient, source, domain and gauge condition recomputed.
