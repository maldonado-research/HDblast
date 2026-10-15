# Next target: exact quasilinear contact evolution for the aligned seed

Planning note, not a proved nonlinear existence theorem. Artifact identifier:
15 October 2026 UTC; this work began 14 October Pacific time.

The next target is to carry the already constructed exact localized initial
family through the physical shell-contact problem. The current checkpoint
provides all-angular **linear** estimates at the registered background and
an **exact algebraic** transformation of the nonlinear boundary laws. These
two results supply preparation for the next theorem; they do not establish it.

## Established inputs to retain

The [exact initial-data theorem](../../../research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md)
already proves, under its background, support, reflection and weighted-inverse
hypotheses, a parameter-smooth exact constraint family for each fixed integer
`k ≥ 12`. The metric is in `H^{k+2}` and the contravariant metric momentum
density is in `H^{k+1}`. The scalar and scalar momentum sectors are prescribed
as the background plus their chosen first parameter derivatives. Corrections
are supported in a fixed closed annulus away from the shell and cone vertex;
the fields equal the background on open neighborhoods of both.

Use the selected `Γ`-invariant aligned seed, with
`A = B = diag(1,−1,0,0)` and the already stated smooth annular profiles.
The same checkpoint supplies a conditional exact short precontact evolution
with `C²` parameter dependence in a lower piecewise spacetime `C⁴` norm.
The family is at one fixed high finite spatial regularity. It does not give
a common spatial `C∞` family on one parameter interval.

The [linear contact estimate](../../../research/HDBLAST_ANGULAR_UNIFORM_CONTACT_20261015/ANGULAR_UNIFORM_LINEAR_CONTACT.md)
controls the complete sixteen-field system uniformly over the full tensorial
SO(4) cutoff. Its positive norm is auxiliary analytic control. The complete
nonsymmetric twelve-row Robin matrix is removed by a bounded collar
conjugation; smooth inhomogeneous data have an explicit lifting cost. This
does not supply arbitrary `L²` Neumann boundary control.

The [exact nonlinear boundary reduction](../../../research/HDBLAST_ANGULAR_UNIFORM_CONTACT_20261015/NONLINEAR_BOUNDARY_STRAIGHTENING.md)
uses a local pointwise state flow to turn all twelve physical and normal-gauge
rows into homogeneous Neumann conditions. The four coordinate Dirichlet rows
remain homogeneous. The reduced equations retain their common scalar
quasilinear wave principal symbol. The transformation is defined only in
the stated state neighborhood and fixed collar; its exactness alone is not
a quasilinear well-posedness theorem.

## Precise theorem to attempt

Fix one sufficiently high finite `k`, a regular positive-cone collar and a
background time interval extending beyond the intended first shell-contact
time. Construct a sufficiently small common parameter interval on which the
aligned exact initial family has a unique solution of the complete nonlinear
Einstein–scalar equations, both-copy matching and the physical shell laws
on that controlled interval. State the solution's finite Sobolev class,
parameter regularity and all compatibility assumptions explicitly.

The theorem must cover a prescribed interval through contact. Short local
existence at the initial section, while the shell collar is still background,
does not establish this target. The background contact time
`t_c = log(b/R₀)` is a reference time for the selected packet; an exact
perturbed characteristic arrival need not equal that background value.
Choose and prove a controlled lifespan beyond the intended contact region,
using the actual perturbed metric cone bounds.

## Required arguments

1. **Initial harmonic completion and corners.** Retain the exact canonical
   constraint data. Choose lapse and shift data and their initial time
   derivatives so that the relative harmonic vector and its normal initial
   derivative vanish. Verify this using the trace-reversed physical
   Hamiltonian and momentum constraints. Preserve the open background shell
   collar, then transport the exact data through the boundary state flow.
   Check every corner identity required at the chosen finite regularity.
   Finite-order jets cannot be declared compatible merely because each
   frozen iterate has homogeneous boundary traces.

2. **Variable-principal-coefficient control.** Derive mixed boundary energy
   estimates for the actual transformed metric wave operator. Preserve
   `g_{za}=g^{za}=0` at the shell in every iterate, together with the fixed
   Dirichlet/Neumann split. Control the timelike time covector, spacelike
   shell normal and positive spatial energy margins on a common state
   neighborhood. The current frozen-background estimate does not prove
   these variable-coefficient bounds automatically.

3. **Compatible quasilinear iteration.** Build a compatible initial jet
   extension and a sequence of iterates satisfying the mixed traces and
   required corner jets. Prove tame high-order estimates using time/angular
   commutation and recovery of normal derivatives from the equations.
   Track coefficient derivatives multiplying second derivatives of the
   unknowns. Establish a uniform high norm and convergence in an explicitly
   lower norm, followed by uniqueness and continuous dependence. If the
   intended iteration loses derivatives, supply a justified replacement
   rather than infer closure from the common principal symbol.

4. **Parameter derivative budget.** Start from the previously proved
   parameter-smooth curve in fixed `H^{k+2} × H^{k+1}` canonical spaces.
   Reserve enough finite Sobolev regularity for two solution variations,
   coefficient differences and finite-difference quotient convergence.
   Write the first and second variation bulk and boundary problems and
   prove the claimed parameter regularity in their declared lower norms.
   The earlier precontact argument reserved four Sobolev losses; its budget
   is a useful starting point, not a proved through-contact budget. A
   top-norm smooth solution-map assertion would require a separate proof.

5. **Prescribed lifespan.** Prove a continuation criterion and a background
   comparison bootstrap on the chosen interval. Show that a sufficiently
   small common parameter neighborhood keeps the necessary coefficient
   bounds, signature, boundary timelikeness and support inside the controlled
   positive-cone collar until after the intended contact. Identify every
   constant depending on the fixed background interval and collar.

6. **Exact subsidiary constraints and both copies.** Use the actual
   Einstein–scalar Noether identity, scalar equation and physical shell
   laws to derive the homogeneous subsidiary bulk equation and normal
   Dirichlet/tangential Robin conditions. Prove variable-coefficient
   subsidiary uniqueness with zero exact Cauchy data; then conclude that
   the reduced solution solves the physical equations. Verify induced-field
   matching and both-copy parity. Derive any normal-force consequence from
   the established exact identities without adding an extra primal boundary
   condition.

## Reviewable completion criteria

The next checkpoint should contain a complete theorem with the actual model,
boundary split, variable-coefficient estimates, compatible iteration,
parameter-loss budget, lifespan bound and subsidiary proof. Exact algebraic
checks should test the model-specific transformations and signs, with
wrong-formula controls. Independent internal review should inspect the
analytic proof separately from those finite checks. Any computed example
must declare its initial-data regularity, support, constraints and boundary
errors; numerical evolution alone does not certify exact existence.

The existing exact initial family need not have the particular second initial
correction selected in the earlier frozen Taylor construction. Equality of
those coefficients, if desired, requires an additional argument. No such
identification, nonlinear contact theorem, spatial `C∞` family, stability
result, external novelty or higher-dimensional Big Bang conclusion is claimed
in this planning note. The target follows from the current Einstein–scalar
model; adding quantum or string sectors would require a separate physical
model and justification. No new Zenodo version or APS submission is needed
for this planning checkpoint.
