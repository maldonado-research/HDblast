# Independent domain and spectral-gap review

13 October 2026. Read-only mathematical audit of
`../SCALAR_SELF_ADJOINT_DOMAIN_AND_GAP.md`; no network use, numerical scan, or
edits to the reviewed note.

The corrected version was rechecked and pinned to SHA-256
`ad3b23eedad94e47c97952db3b273ff8204f4e898bcb86221be64eb11ac357fd`.

**Verdict:** no blocking mathematical defect identified in the stated abstract
fixed-background theorem. The hypotheses are sufficient for the form
closure, operator domain, ground-state identity, strictly positive gap,
essential threshold, and conditional ordinary-Robin companion.

## Review findings resolved in the pinned version

1. **Essential-spectrum comparison clarified.** The Dirichlet operator used
   in the finite-rank resolvent comparison in section 6 must have the same
   potential `V_L` as (24). The preceding compact-perturbation argument first
   transfers the essential spectrum from the constant-potential Dirichlet
   operator to this `V_L` Dirichlet operator. The phrase "that Dirichlet
   operator" was ambiguous. With the same potential, the homogeneous-equation
   and finite-rank proof is correct; with the original constant potential,
   that particular proof would be false. This is a precision repair, not an
   obstruction to the conclusion. The corrected note explicitly defines
   `D_V` using `V_L` and distinguishes this finite-rank comparison from the
   compact potential comparison. The issue is resolved.
2. **Displayed-equation rendering repaired.** The original domain formula had
   a single backslash before another alignment marker, and the optional Hardy
   calculation contained two form-feed characters in place of fraction
   commands. The corrected version combines the domain conditions properly
   and restores both fraction commands. A byte-level scan found no unexpected
   control characters. These issues are resolved.

No outstanding corrective finding remains within this review's scope.

## Arguments independently checked

- **Closed form and cone core:** local `H^1` convergence identifies the regular
  endpoint trace. The cutoff estimate follows from
  `alpha/w = rho^2 = O(y^2)` and controls the only extra derivative term. The
  independent endpoint component does not obstruct Hilbert-space density.
- **Operator realization:** compact interior tests give the bulk equation;
  regular-endpoint tests give the stated positive flux component. Conversely,
  integration by parts against the cone core and form continuity establishes
  the complete domain without evaluating a cone boundary product. The
  proposed core operator is dense, symmetric, and semibounded, so the literal
  Friedrichs statement is valid.
- **Cone classification:** reduction of order gives the second solution of
  order `y`; its weighted square is nonintegrable, whereas `P_0` is square
  integrable. The endpoint is limit-point without differentiated asymptotics.
- **Picone identity on the closure:** the core identity equates the form norms
  `e + N` and `D + 5N`. Local convergence identifies the transformed weak
  derivative. The reverse inclusion follows from the analogous transformed
  cutoff approximation. No endpoint expansion of arbitrary form vectors is
  needed.
- **Gap:** the power-weight integration by parts is valid for complex
  functions on each positive truncation. Its nonpositive lower boundary term
  permits passage to zero without a cone trace assumption. The resulting
  Hilbert-space triangle and Cauchy-Schwarz estimate proves the displayed
  coercivity on the full form domain. The bound is a background functional,
  not a certified numerical constant or a uniform family bound.
- **Liouville realization:** the orientation, endpoint-mass sign,
  coefficient-only potential (22), and threshold `25/4` check out.
  Differentiated asymptotics make `s` and `V_L` bounded; thus the stated
  `H^1` form and `H^2` operator domains follow. With the comparison clarified
  above, the essential-spectrum and simplicity arguments are valid.
- **Ordinary Robin companion:** the Robin sign, transformed endpoint penalty,
  positive gap with zero endpoint mass, and unchanged essential threshold are
  correct under its separately stated hypotheses.

## Scope and limits

This audit checks the abstract operator analysis. It does not independently
derive the perturbation reduction, physical kinetic norm, branch sign
certificates, candidate odd-sector Robin condition, or family limit of the
physical scale. Equation (23), which invokes background equations not stated
in the theorem, was not independently verified; the coefficient-only formula
(22) suffices for all audited spectral conclusions. No claim about a uniform
dimensionless gap, exceptional harmonics, or nonlinear stability follows.
