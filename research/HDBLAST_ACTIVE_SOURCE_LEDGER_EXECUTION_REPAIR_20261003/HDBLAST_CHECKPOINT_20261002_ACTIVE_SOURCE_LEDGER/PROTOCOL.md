# Prospective active-source saved-trajectory diagnostic

This is the fixed prospective protocol. The retained arrays and earlier
active-source failures were already published; this diagnostic is motivated by
those known failures and is not blind confirmation. At preparation, its public
freeze and result report have not yet been made, and the new reference-flow,
contact, direct-integral and attribution evaluations are **UNCOMPUTED**.
Physical evaluation requires the
verified public registration described below. `EXPERIMENT.json`,
`OUTPUT_SCHEMA.json` and the two implementation notes define the exact
registered numerical and output interfaces.

The motivating result is the registered source-free diagnostic at public
science commit `a8394d0127e58200cd4b9a87c14ae63a9dd69f02`: it demonstrated a
stored-history integration defect on `[-2.5,-1.5]`. The original metric
experiments remain scientific **FAIL**, including the retained 29 endpoint and
30 ledger-refinement failures. This new experiment tests the local numerical
pipeline while the same prescribed source is active.

## Fixed question, interval and input cases

Use both `positive_B` and `signed_uB`, both saved `coarse` and `fine` settings,
and every `K` in `{64,128,256}`: exactly twelve cases. The single diagnostic
interval is `eta_a=-4.5` to `eta_b=-3.5`, with fixed midpoint `eta_m=-4`.
Anchor unrestricted forced flow on the actual saved `u_1,w_1`. Compare predicted
flow with saved `u_2,w_2` and `u_3,w_3`; do not project a normalization,
Wronskian or canonical invariant. Compare every inclusive retained density,
pressure and Ward-integrand history sample: 129 coarse and 257 fine samples.
No interval, source, cutoff, control, tolerance or outcome subset is selected
after inspection of physical results.

The four deterministic input capsules copy nineteen NPY members verbatim from
the original 631-member archives. Their manifest binds full original archive
sizes, SHA256 and published Git objects at original science lineage commit
`71d00cc423e9049c8166b7ee7afbefbac28d8a18`, selected member bytes and headers,
and deterministic capsule reconstruction. The input manifest's preparation
SHA256 is `b70dea97bc7223a77f32c5e0964ee2a43812215d84bdac6393bedbb6601a291d`;
the final full registration binds these same bytes. Provenance checks and
header inspection must complete before an array-value reader is called.

Retain the existing finite comoving momentum nodes, weights, incoming state,
`H=1`, `x=r=2`, `xi=0`, fixed physical mass, reference scale and subtraction
inventory. Represented long-double input values convert through exact integer
ratios. Split complex long-double values into long-double real and imaginary
components; conversion through binary64 `float` or `complex` is prohibited.
Exact ratios preserve saved rounding; they do not recover unsaved exact data.

The producer normalization constants are fixed:

| Constant | Exact represented ratio |
| --- | --- |
| epsilon | `3777893186295716171 / 37778931862957161709568` |
| pi | `14488038916154245685 / 4611686018427387904` |

## Actual analytic source and independently constructed operators

For `z=eta+4`, the source is `B(z)=exp(-z^2/(1-z^2))` inside `abs(z)<1`,
with its smooth zero extension outside. The signed source is `z*B(z)`.
Directional source derivatives are evaluated from the exact integer bump
polynomial inventory; the algorithms use native long-double arithmetic at
their fixed quadrature nodes. They use no polynomial or sampled source
surrogate. The primary generates Gaussian nodes and weights with a native
binary80 Newton iteration. The independent route uses NumPy `leggauss` nodes
and weights computed in binary64 and promoted to long double; promotion does
not recover more accurate Gaussian constants. These route-specific quadrature
constants are separate from the preserved native long-double physical input
nodes and weights, which are never downcast.
Set `L=-1/eta`, `g=4*L^2*h-2*L*h_prime-h_second` and

```text
u_prime = w
w_prime = 2*i*k*w - epsilon*g.
```

The physical direct density and pressure are derived independently from the
minimal-scalar operators and complete retained fixed-reference metric
subtraction. Include mode, operator, mass, scale/prefactor and baseline work
contacts. Neither density nor pressure may be defined by a Ward equation.
For scaled response stresses `R,P`, the finite-band Ward integrand is

```text
F = L*(R-3*P) - 3*h_prime*(R0K+P0K).
```

Here `R0K,P0K` are the matching conformal baseline integrals. The source work
term cannot be dropped or replaced with a continuum anomaly. Moving a momentum
cutoff changes the target and is excluded from this calculation.

## Primary primitive and independent direct-integral routes

The primary route uses exact phase variation of constants on each inherited
time cell, with fixed GL24 and GL32 forcing controls that evaluate the actual
analytic source. Native binary80 source and complex160 phase arithmetic is
explicit. With `U=u/epsilon`, `W=w/epsilon` and
`mu_j=momentum_weight_j*k_j^2/(2*pi^2)`, define the mode moments

```text
C = sum(mu_j*Re(U_j)/(2*k_j))
B = sum(mu_j*Re(W_j)/(2*k_j))
E = sum(mu_j*k_j*Re(U_j))
J = sum(mu_j*Im(W_j)/2)
G = 3*L^2*C-L*B.
```

The separately expanded bare mode stresses obey
`Rm=E+3*L^2*C-L*B-J` and `Pm=E/3-L^2*C-L*B-J`.
For real forcing, `E-J` is an invariant and
`G_prime=L*(Rm-3*Pm)+M_d*L*g`, where
`M_d=sum(mu_j/(2*k_j))` uses the saved momentum quadrature.
The primary canonical integral uses independently propagated Green endpoints,
the separately integrated forcing term and direct full contact work:

```text
I_primary = delta(G_d) - M_d*integral(L*g) + integral(F_contact).
```

This endpoint reduction is permitted only after the derivative identity is
proved against the independently expanded direct integrand for unrestricted
complex amplitudes and all contacts, and after independent direct integration
is implemented. **Measured final endpoint density is never an input to I.**

The independent route explicitly reconstructs the direct mode/operator
density and pressure, differentiates the published W0/W2/W4 subtraction
inventory with separate directional algebra, and integrates direct `F` at
dense local quadrature nodes. It uses its own Duhamel kernels, source jets and
contact evaluation. The three fixed controls are GL16 source/GL16 ledger,
GL16 source/GL24 ledger and GL24 source/GL24 ledger. The first two share the
same GL16 endpoint trajectory and isolate direct ledger quadrature order.
The last two hold ledger order at GL24 and isolate forcing quadrature order.
Shared geometry is cached. These are empirical numerical controls; they are
not rigorous error bounds. Also compare GL16 source/GL16 ledger against GL24
source/GL24 ledger across the full declared numerical universe at `2e-7`;
passing the two adjacent comparisons does not replace this joint control.
No extra control is introduced after physical evaluation.

Both routes retain all declared controls at both MP80 and MP100 reduction
contexts. These precision contexts apply to final accumulation, and to any
explicitly identified primary analytic contact evaluation. Native source,
phase and propagated mode arithmetic stays native long double. Precision
agreement does **not** certify phase roundoff, source quadrature, dense output,
inherited state, continuum momentum or physical model accuracy.

## Saved-discrete and analytic finite-band targets

Raw results from different contact momentum rules must not be compared as if
they were the same numerical target. If the primary route uses analytic
finite-band contacts, it must label them `analytic_finite_band`; the independent
direct-stress target is `saved_discrete_momentum`. Report raw analytic versus
discrete contact, baseline, source-work, profile and integral gaps and their
signed `E_momentum` corrections at both precision contexts and controls.

In particular, `M_A=K^2/(8*pi^2)` is the analytic bare coefficient. A primary
analytic-contact primitive check must retain
`delta(G_d)+delta(C_R^A)+(M_A-M_d)*integral(L*g)`, with matching signs.
The exact signed momentum correction and matched primary scalar are

```text
E_momentum = delta(C_R^d-C_R^A) + (M_d-M_A)*integral(L*g)
I_primary_discrete = I_primary_analytic + E_momentum.
```

These corrections are fixed analytically before the outcomes; they are not a
fit. The primary independently evaluates the per-node retained direct
subtraction and stable baseline at the fixed initial, midpoint and final knots
through the pinned reference inventory. It must retain separate density and
pressure contact and baseline values for discrete and analytic targets there.
Match those three knot profiles and the integral scalar to the discrete target
prospectively. Preserve full-history raw analytic-versus-discrete profile gaps
and identify their target difference; three knot audits do not certify contact
matching over the entire history. A momentum correction diagnoses the
specified finite-band target difference, not continuum-momentum convergence.

Cross-route full-profile comparison uses mode-only density, pressure and Ward
integrands after separately removing each route's own contact contributions,
at every retained sample. Compare all two primary controls against all three
independent controls at both reduction contexts. Cross-route scalar comparison
uses the matched discrete integral. Compare direct discrete density/pressure
contacts and conformal baselines at the three fixed knots. Primary physical
baseline values `rho0,p0` are multiplied by `L^4` to match the independent
conformal `r0,p0` units. Full raw profile/contact/baseline differences are
retained; mode-only agreement plus three knot audits does not establish full
contact-target equivalence over the entire interval.

## Unchanged Simpson ledger and signed attribution

Compute original full-history Simpson global prefixes in NumPy long double,
then subtract endpoint prefixes in long double before exact-ratio conversion.
This is `S_ab`; do not replace it with a reset local Simpson ledger or a
freshly evolved trajectory. Report reset local Simpson and the difference
separately. For the fine same-trajectory control, subsample the complete stored
`F_fine[::2]`, double its spacing, use corresponding even prefix indices and
again subtract global prefixes in long double. This changes only Simpson
sampling, not modes, stresses or source evolution.

For separately converted saved endpoint density values and the declared
matched-target integral form

```text
DeltaR = R_b-R_a
D_S = DeltaR-S_ab
D_cont = DeltaR-I_ab
E_Q = I_ab-S_ab
D_S = D_cont+E_Q.
```

Keep signed quantities and raw mismatched-target variants. Report the signed
mode-density projection of measured midpoint and endpoint forced-flow
differences and separate coefficient triangle envelopes. Initial, midpoint and
endpoint mode-contact profiles are also checked. A projection/envelope concerns
the measured anchored differences only; it does not bound errors inherited
before `eta_a`, source quadrature or total trajectory error.

The mandatory signed direct-operator increment mismatch is

```text
E_operator = [R_saved-R_mode_on_measured-C_R^d]_b
             - [R_saved-R_mode_on_measured-C_R^d]_a.
```

Report its two bracketed endpoint values separately, and the corresponding
midpoint audit. It is not permissible to suppress the initial bracket, fold
the contact difference into a fitted operator correction, or regard this
increment check as a certificate of the incoming anchor state. The final
flow/momentum/operator decomposition must also retain the weighted canonical
invariant drift and independent source-primitive residual with their signs.
The primary's fourth signed reconstruction term is

```text
E_reconstruction = delta(E-J) + delta(C_R^A)
                   + M_A*integral(L*g) - integral(F_contact_A).
```

The uncorrected analytic residual and matched discrete residual must carry
different field names. Their fixed signed decompositions are

```text
D_cont_analytic = E_flow + E_momentum + E_operator + E_reconstruction
D_cont_matched = E_flow + E_operator + E_reconstruction
D_cont_matched = DeltaR - (I_A + E_momentum).
```

Check these exact signed decompositions at `1e-12`, while testing the
independently computed numerical reconstruction term itself at `2e-7`.
Do not confuse these two checks or infer numerical accuracy from exact
definitional closure.

## Fixed gates and failure semantics

All twelve cases at both declared precision contexts must meet:

- Direct density, pressure and `F` profile matching at absolute `2e-7`.
- Baseline source-work contact profile matching at absolute `2e-7`, retaining
  the raw analytic-versus-saved-discrete target gap explicitly.
- Matched-target continuous-ledger residual and signed measured-flow density
  projections, including the midpoint projection, at absolute `2e-7`.
- Initial/midpoint/contact matching and forcing/integral control differences at
  absolute `2e-7`, using the final prospectively fixed field inventory.
- Initial, midpoint and final direct operator, contact and baseline audits at
  absolute `2e-7`; operator endpoint increments alone cannot replace the
  separate knot mismatches.
- `E_operator` and `E_reconstruction` at absolute `2e-7`. The raw analytic
  momentum correction is also reported and gated at `2e-7`; the principal
  ledger comparison uses the prospectively matched discrete integral.
- Canonical invariant drift and the independent source/contact primitive
  residual individually at absolute `2e-7`, in addition to their summed
  `E_reconstruction` gate; compensating errors cannot suppress these controls.
- Common matched-target numerical agreement across the two routes at absolute
  `2e-7`, which reflects their native long-double numerical scope.
- Final MP80/100 accumulation gaps, exact definitional signed closure and decimal
  serialization roundtrip checks at absolute `1e-12`.

Every nonmetadata reported arithmetic scalar/profile participates in its
declared precision or numerical control universe. `OUTPUT_SCHEMA.json` records
the exact literal field inventories under `source_constants`, including
`PRIMARY_SCIENCE_SCALARS`, `PRIMARY_SCIENCE_ANCHORS` and
`INDEPENDENT_SCIENCE_SCALARS`. Its `additional_evaluator_gates` also specifies
the independently recomputed scalar, profile, three-knot and complete control
universe checks; the literal science-field lists alone are not the full gate
inventory. The authoritative validator applies both. Independent native
profiles are retained once as
exact represented values; they are not separately recomputed with MP80/100.
Identical stored-reader scalar values must agree at `1e-12`; a discrepancy is
an integrity/arithmetic failure, distinct from the native numerical controls.
This includes the original global-prefix Simpson result, reset result and
fine-only double-step result from the two readers. The validator separately
recomputes `DeltaR` from serialized stored endpoint profiles. Its JSON inputs
contain only the local `F` interval, not the complete prefix from `eta=-6`;
reader agreement is therefore not a third independent validation of full
native global-prefix Simpson arithmetic. The registered reader implementations
must still preserve the complete-prefix expression and long-double subtraction
specified above.
Signed closure is a necessary accounting check, not independent evidence of
conservation or trajectory accuracy. Decimal output values remain strings and
must not roundtrip through binary64.

The canonical attribution rows are the primary GL32/MP100 matched-discrete
rows and the independent GL24 source/GL24 ledger/MP100 rows. Every declared
control at both precision contexts must pass all consistency criteria before
these canonical rows can establish a witness. Only classify an active-interval
ledger error as demonstrated when at least one same registered case has
canonical support from both routes satisfying all of

```text
abs(D_S) > 2e-6
abs(D_cont) <= 0.1*abs(D_S)
abs(E_Q) >= 0.9*abs(D_S).
```

No witness and otherwise passing consistency gives
`NO_GATE_SCALE_ATTRIBUTION`. A failed scientific consistency gate gives
`CONSISTENCY_FAILURE` as a retained completed calculation. Neither alters
`old_metric_status="FAIL"`. Every native long-double numerical GL,
cross-route, contact-reconstruction, profile, flow or operator control above
`2e-7` is a retained **scientific consistency failure**, with execution exit
zero. It is not a fatal arithmetic error. Integrity, resource, schema,
MP80/100 arithmetic-gap, exact definitional closure, serialization or route
execution failures prevent an accepted classification and are retained as
execution failures (fatal result gates exit two; other nonzero exits also
fail execution). Do not erase failed
trials, retune gates or edit a frozen implementation to obtain a desired result.
A revised method requires a new prospective registration before new evaluation.

## Resource, reproduction and publication boundary

Each complete route has one **900 second / 262144 KiB** allowance for all twelve
cases, all controls and both MP contexts, beginning before provenance/member
decoding and ending after output serialization. Do not reset it per case.
Outer receipts capture actual child exits, timeout, peak RSS, stdout, stderr
and partial progress. Symbolic identity checks run during preparation; the
registration binds their proof code and receipts. The physical route budget
covers byte authentication, input decoding, numerical evaluation and output
serialization. It does not claim a fresh execution of the SymPy proofs inside
each physical route. Run physical routes in normal Python 3.12.14 under the
pinned dependency versions and x87 binary80 ABI, with one numeric thread and
no child processes. The independent complete execution has its own equal
budget; two individual passes do not establish a combined runtime claim.

Recorded fabricated feasibility executions used 219.155 seconds / 183944 KiB
for the primary and 473.506 seconds / 110420 KiB for the independent numerical
route. Both include all twelve cases, fixed numerical controls, reduction
contexts and serialized numerical outputs. They do not evaluate the study
source or original array values. The independent numerical benchmark used the
unchanged final numerical engines before later production entry/schema guard
additions; those additions were checked separately in normal and optimized
Python. Its benchmark did not execute every final production entry byte as
one whole application. These preparation checks support feasibility; actual
production exits, wall time and peak RSS still come from the authoritative
outer receipts and must satisfy the frozen complete-route budget.

Before any physical source, phase or array evaluation, publish every frozen
protocol, schema, implementation, proof, input manifest and capsule in one
full registration. Independently read the registered public commit and verify
all blob SHA256 values. Record the external full commit and registration hash
in a protected `FREEZE_RECEIPT.json`. CLI bootstrap authenticates those pins
before physical source calls or numerical array decoding. The primary imports
ordinary numerical libraries before this authentication but does not evaluate
physical source values or decode arrays first. The independent has a stronger
stdlib-only entry wrapper, authenticating registered code/input bytes before
importing its numerical core. Subsequent replay rechecks frozen bytes before,
between and after separately budgeted commands. Use fresh external output
directories and preserve original and fresh standalone-package replay.

The round identifier retains `20261002`; the prospective registration occurs
on 3 October 2026 UTC. The recorded public freeze receipt, rather than the
round identifier, establishes the actual publication and verification times.

The report will state what this local anchored experiment can establish. It
cannot establish a new physical law, repair the full metric calibration,
bound the incoming state, validate lapse/bulk/state response, demonstrate
coupled Einstein evolution or heating, or provide evidence for a higher
dimensional origin of the hot Big Bang.

Automatic Zenodo GitHub archiving is OFF for all eight public repositories,
superseding earlier integration-on notes. Publication, if warranted after
review, uses the existing manual token/API workflow and existing HDBLAST
concept DOI `10.5281/zenodo.17088132`. Latest published main record is
`22347452` (v24); linked unpublished main draft is `23112891`. Its unchanged
metadata HTTP500 blocker remains unresolved. Do not retry the same failing
request, re-enable autoarchiving, create a test release, duplicate the family
or treat a reserved draft DOI as publication. Saved metadata plus filename,
byte size and checksum verification is required before any publication.
