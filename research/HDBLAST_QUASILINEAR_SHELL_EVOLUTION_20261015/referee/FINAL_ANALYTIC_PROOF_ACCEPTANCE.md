# Independent final analytic proof acceptance

15 October 2026 UTC. Internal critical mathematical review of the actual
proofs listed below. I independently audited the constructions and estimates;
this is not an endorsement inferred from matching conclusions or finite
algebra checks.

## Exact reviewed sources

| Source | SHA-256 |
| --- | --- |
| `astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md` | `d6f80e6158226b253c56f2f60c8d9c6cfa2ee2940e79bde42dd8fe8e48a5d6c4` |
| `sol/INDEPENDENT_QUASILINEAR_MIXED_PROOF.md` | `76b83bf96b1f6cb82f63c8b2ea3acab9e9c9394e3f961db350a05d6ea68bb1f8` |
| `sol/linear_audit/LINEAR_VARIABLE_WAVE_AUDIT.md` | `fcd7dc55b9459f0ebc56b1e6a77584be6a101878126efb692614378b4884af86` |
| `gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md` | `74bce09080d4015b2eccf8c183b33f47b80516216a9547837e3e0e34251d7b03` |
| `gauge-parameter/FINITE_REGULARITY_VARIATION_COMPLETION.md` | `b965f5d698a7e2382737376ad904130f3773f5658374bfc04ded5ea44e482e32` |

All five hashes were read directly from the reviewed filesystem bytes.

## Verdict and mathematical reasons

I accept the analytic conclusions in these sources, with their stated
hypotheses, finite regularity and scope. The actual proof now closes the
obstacles identified in the preliminary audit.

1. The variable-coefficient energy has positive time and spatial blocks.
   Both time-normal and spatial-normal cross coefficients vanish at the
   shell. The orthogonal homogeneous D/N split therefore annihilates the
   complete principal flux. Lower-order matrices need no symmetry.

2. The Astra frozen solver supplies existence in addition to an estimate.
   A boundary-preserving smooth Borel lift removes initial projection
   incompatibility. Its fixed elliptic basis is used only with L2 projection
   bounds. The principal C0 difference is absorbed in a discrete H2 estimate
   on short coefficient intervals; coefficient derivative differences need
   only be bounded because they multiply H1 terms. Differentiated time
   energies and the discrete H2 bounds close simultaneously, before taking
   a limit. Full elliptic spatial regularity is obtained from the limiting
   equation, rather than assumed for projected equations.

3. The independent angular/radial solver is also sufficient. Projection
   coefficient commutators are retained. The projected normal stiffness has
   a uniform positive inverse in L2, which suffices for differentiated
   normal recovery. At fixed angular dimension, radial variational energies
   corrected by lower-level cross terms give time regularity before normal
   regularity. These construction constants may depend on angular dimension;
   the separate geometric energy/normal-recovery estimates used for the
   angular limit do not. The specially selected global angular frame in the
   child audit avoids the invalid arbitrary coordinate-component cutoff.

4. The high-order nonlinear estimates are genuinely tame. Astra's finite
   derivative count puts the lower factor in L-infinity with three remaining
   spatial derivatives. Sol instead substitutes the equations for mixed
   time derivatives and applies ordinary spatial Moser estimates to a finite
   set of recent iterates. Neither argument interpolates independent time
   jets at a fixed time. Normal derivatives are recovered by induction in
   normal count, without commuting normal derivatives through the boundary.

5. Exact nonlinear initial jets are retained by the protected Cauchy seed
   and every frozen successor. Instantaneous lower source terms in normal
   recovery are handled by the fundamental theorem of calculus in the jet
   norm (or Sol's explicit lower spatial interpolation). This closes a
   uniform high norm on the whole prescribed background interval for small
   enough initial data. Low-energy differences have a Volterra estimate,
   yielding factorial convergence or an equivalent weighted contraction.

6. The strong top Sobolev conclusion includes a real smoothing comparison.
   Its higher factor is the coarser smoother solution; its lower difference
   is controlled one order below the top. Smoothing gives o(h) at that lower
   order against an O(h^-1) extra high norm. The high estimate at the extra
   order depends only on the already controlled base low norm, so it does
   not impose extra high smallness or shorten the prescribed lifespan.

7. The gauge/subsidiary companion has the correct exact ADM signs, selected
   initial lapse tangent, trace-reversed constraints, Codazzi cancellation
   and coordinate Robin coefficient N sigma/3. The finite-coefficient
   addendum constructs variations using actual smoothed nonlinear solutions
   rather than independently smoothing boundary coefficients. Bounded first
   and second variation jets at orders 13 and 12, with strong comparison at
   orders 12 and 11, suffice for the asserted C2 parameter curve at order10.
   No top-order differentiability is presumed.

## Retained conditions and limits

This acceptance concerns a fixed finite background interval and controlled
collar; sufficiently small data; scalar principal coefficients depending on
the state, not its derivatives; preserved fixed D/N traces; positive energy
margins; and the existing exact initial-family/background hypotheses when
used physically. It applies to the whole-interval construction from data
flat on an open initial shell collar. A general finite-Sobolev restart at
an arbitrary nonzero corner is not proved and is not needed for this result.

Application to the physical contact family must explicitly combine these
sources with the retained exact nonlinear boundary flow, initial canonical
constraint family, harmonic completion and subsidiary uniqueness. The
consolidated physical theorem should be reviewed against its own exact
bytes before publication. Existence through a reference interaction window
alone does not prove a nonzero invariant shell response or an exact arrival
time. No common spatial C-infinity family, global stability, external
novelty, external peer review, or cosmological origin conclusion follows.
