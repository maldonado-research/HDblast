# Claims, limits and the next calculation

## Scientific scope

This checkpoint combines the inherited massive de Sitter action and corrected
classical shell sensitivity in a stationary semiclassical boundary-value
calculation. Both action-derived quantum sources are recomputed at the moving
shell endpoint. The selected regular bulk, positive metric branch, and two
shell junctions are solved together. This is a local numerical family in an
explicitly declared model, with empirical convergence evidence.

The model has one canonical free, minimally coupled shell scalar in the
Euclidean/Bunch–Davies de Sitter state. At fixed positive reference
`r=0.00011848029588657912`, its dimensionless mass law is

```text
x = r [1 + b(eta − eta0)/2]²,
x_phi = b r [1 + b(eta − eta0)/2],
x_phiphi = b² r/2,
b in {−1,0,+1},  eta0 = −0.00008405268308670538.
```

For nonzero `b` this is the original quadratic interaction with `m0=0`,
`Ghat=sqrt(r)/2`, and `phi_chi=1+eta0−2/b`. The `b=0` member has constant
positive mass squared `r`. These choices specify a research family; the
calculation does not infer them from observations. The reference relation
`r=2H0²` is imposed once. It is held fixed during variation and continuation.

The action density `W` supplies `rho=W−zW_z/2`, `Q=2W_x`, and
`j=x_phi Q/2`, with `z=H²`. The invariant stress has `p=−rho`. For positive
`gamma`, the chosen normalization is `N=1` and `gamma=kappa5²/L³`.
The classical control at `gamma=0` disables the paired sources. In every
nonzero case, the sources entering the two junctions are `gamma*rho` and
`gamma*j`. No counterterm is retuned as gamma or the endpoint changes.

## What the derivation establishes

The full Newton Jacobian retains mass and curvature feedback from the common
action Hessian, including `x_phiphi*Q/2` in `j_phi`. Scalar variation is
independent of metric variation: in general `j` is not `rho_phi`. The stable
metric residual is algebraically equivalent to the original unsquared
junction only on the enforced positive-q, positive-total-tension branch. The
radial constraint and moving-endpoint derivatives are checked separately.

The archived symbolic verifier passes 51 exact identities and detects 20
deliberately wrong formulas, both normally and with Python optimization.
These include the mass-zero current obstruction. At fixed positive H and
for `b=±1`, putting `d=phi−phi_chi` gives

```text
x = r b² d²/4,
Q ~ 3 H⁴/(8 pi² x),
j = x_phi Q/2 ~ 3 H⁴/[8 pi² (phi−phi_chi)].
```

The zero of `x_phi` therefore does not remove the singular current. For finite
positive gamma and finite classical junction data, smooth continuation
through this mass zero is obstructed within the same invariant-state
prescription at fixed positive H. The registered local positive-mass domain
excludes that point. This result does not cover a simultaneous H-to-zero
limit or a different infrared state prescription.

## How to interpret numerical evidence

The registered 18 model points use `gamma={0,.01,1,100,10000,1000000}` and
the three slopes. Primary and refined integrations cover the entire grid;
Radau and a separate cone-start refinement cover both gamma endpoints. The
independent acceleration integrator solves its own six-point subset at two
resolutions, starting from inherited classical coordinates and its own
continuation rather than producer roots. Direct proper-time sources are
checked at four nonzero independent endpoints at two fixed settings.

Root acceptance, observable-shift detection, nonlinear-feedback detection,
and gravitational hierarchy are separate questions. Shifts are measured
against each resolution's own solved classical control. A shift is resolved
only above ten times the observed cross-resolution spread plus the registered
floating-point floor. Nonlinear departure must also exceed one part in a
thousand of the observed shift. An accepted root alone does not establish
either kind of detection. The complete classifications and signed values
are retained in the output records summarized in [00_READ_FIRST.md](00_READ_FIRST.md).

Analytic determinant-tail bounds cover specified source and Hessian series
remainders. Separate proper-time spectral and infrared bounds cover those
contributions only. Arithmetic, quadrature, ultraviolet truncation and radial
integration have no combined rigorous enclosure. The convergence envelopes
and condition numbers supply numerical evidence; they do not prove exact
existence, uniqueness, absence of folds, or the maximal extent of a branch.

For each positive-gamma root, the separately declared screen is
`Ehat*gamma^(1/3) <= 0.1`, with `M5 L=gamma^(−1/3)`. `Ehat` includes the shell
mass and Hubble scales, endpoint `abs(q)`, and sampled intrinsic bulk curvature
scales. It excludes the coordinate divergence of cone-slicing q. The samples
do not certify a continuum supremum. Passing this convention establishes
neither a universal cutoff nor full EFT control; failures are retained and
limit any physical interpretation. A small H alone cannot replace the bulk
screen. The gamma-zero control has no inferred physical Planck ratio.

Solving the sourced equations at finite gamma includes repeated feedback of
the specified one-loop determinant. It does not calculate all two-loop
effects, graviton loops, source-field loops, or omitted higher operators.
For this N=1 family, increasing gamma provides no large-N suppression of
metric fluctuations. Resolved feedback in a formal high-gamma root does not
establish a controlled physical quantum correction.

## Provenance and review limits

The initial checkpoint contents, including both numerical implementations,
their registrations and the symbolic files, were publicly frozen at commit
`b4f77f5826e0a603400513832aa53b0673df7533` before any new quantum-source or
radial/root evaluation in this round. The full manifest pins 17 files; the
manifest itself is the eighteenth initial checkpoint file. The inherited
scientific baseline is `e24643d579f7560db719775b2998c0ffa22898b0`.

The original primary validator is retained byte for byte. It initially
reported `FAIL` because the diagnostic keyword `condition=` collided with
the same named positional argument of its `gate` helper, producing one
TypeError per root. The separately named v2 validator changes that diagnostic
keyword only. It is a disclosed post-run implementation correction: original
evidence and frozen source hashes remain preserved, and registered scientific
tolerances, model choices and grids remain the basis of acceptance.
The repair was published as `d8a0b28fa7440215108ddfff51355816c932b541`
before v2 validation. Later record audits, summaries and replay infrastructure
are also identified as post-run work rather than prospective registration.

The independent post-run comparison audit had a separate initial JSON
serialization failure involving a NumPy boolean. Its failed script, log and
exit record are retained under `independent/`. The passing
`audit_primary_v2.py` converts that metadata to ordinary JSON scalar types
and supplies portable input/output arguments; it changes no scientific
parameters or gates. The prospectively frozen independent root and
proper-time calculations completed with their original sources.

Separate derivations and implementations in this AI-assisted workflow are
internal cross-checks. They are not external peer review or independent human
replication. The targeted literature search supplies source and dimensional
context; it does not establish novelty or priority. The inherited DESITTER
checkpoint and [PR #14](https://github.com/maldonado-research/HDblast/pull/14)
remain explicit dependencies.

## Next work supported by this result

1. **Establish a physical hierarchy and matching model.** Specify the physical
   length, gravitational scale, source-field normalization, finite local
   coefficients and any lower cutoff. Assess bulk curvature and higher
   operators as well as shell mass and H. A new model or enlarged parameter
   scan requires a new prospective registration.
2. **Determine whether nonlinear feedback can be resolved in a controlled
   regime.** Improve the error analysis or supply an independently justified
   hierarchy before assigning physical meaning to a finite-amplitude effect.
   Preserve the present detection thresholds and unresolved cases as the
   record of this experiment.
3. **Certify a local stationary family if a theorem is needed.** Bound cone
   initialization, source derivatives and radial sensitivities with validated
   arithmetic, then apply an explicit interval inclusion argument. Floating
   point residuals and a finite determinant do not supply that proof.
4. **Derive causal response and stability separately.** Fix an initial state,
   derive the in-in metric/scalar response with its memory terms, retain the
   same finite local matching, and enforce the coupled conservation identity.
   A bounded next calibration can use the inherited exact-mode point
   `x=r=2H²` to test causal linear response before general stability or
   time-evolution calculations. This proposed calibration is not computed in
   the stationary checkpoint. The static action Hessian and stationary
   shooting Jacobian are insufficient for a dynamical stability or attraction
   conclusion.
5. **Evaluate heating only after a dynamical energy-transfer model exists.**
   Define particle or detector observables, track the full bulk/shell energy
   budget, and provide a thermalization mechanism before making a radiation
   or hot-universe claim.

The next bounded task is specified in the [prospective causal-response protocol](theory/PROSPECTIVE_CAUSAL_RESPONSE_PROTOCOL.md): test the derived retarded variance kernel, local subtraction, state dependence and memory on a prescribed geometry before attempting metric response or shell stability. Its analytic calibration is separate from the completed stationary numerical experiment.
