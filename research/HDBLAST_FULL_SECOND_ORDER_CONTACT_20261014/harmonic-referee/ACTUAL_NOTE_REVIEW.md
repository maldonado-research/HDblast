# Review of the actual complete-contact note

Reviewed file: `../FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md`.

Reviewed SHA-256:
`59c58b560b257a7311d4d00a72a1d4facb95c23377b88cf0323a4d57edc81949`.

**Decision: accept the stated conditional finite-angular second-coefficient
theorem.** The entire file with this hash was read. The conclusion agrees
with the independent proof in `INDEPENDENT_CONTACT_REFEREE.md`. There is
no substantive correction required to the equations or construction.

The accepted scope is a selected smooth second Taylor coefficient in the
actual five-dimensional Einstein–scalar model on a controlled finite
positive-cone collar through shell contact. The proof depends explicitly
on the inherited smooth first contact solution and smooth localized sourced
second initial-constraint lemma. It does not derive those inputs from the
ten charge tests alone. It does not use a finite-regularity nonlinear
implicit-function theorem to manufacture a smooth formal coefficient.

The following points were independently reviewed:

1. The background-reference gauge is compatible with the Gaussian radial
   background. The normal lapse law in (8) and the four-plus-twelve
   independent boundary rows are correct.
2. The exact second boundary source expressions in (9)–(10), including
   the `-16 a_b nu_1^2` lapse term and `k_yy=2 nu_2+2 nu_1^2`, follow from
   direct second differentiation of (8).
3. The first wave-gauge conversion preserves the boundary, and equation
   (6) includes both quadratic Lie terms. The complete reduced Hessian
   source and `Q_C` prescription is essential and correctly retained.
4. SO(4) isotypic projection of the inherited smooth localized correction
   preserves its sourced constraint equation and radial support. The
   complete representation includes the exceptional scalar monopole.
5. Initial harmonic completion has the stated invertible lapse/shift
   derivative coefficients. The trace-reversed physical constraints,
   rather than an isolated normal Ricci row, give the zero initial time
   derivative of the subsidiary field.
6. The curved Codazzi cancellation uses the scalar junction. It yields
   the covariant subsidiary Robin condition, whose coordinate form has
   the required factor two. There is no extra tangential harmonic
   prescription on the primal problem.
7. The characteristic construction solves the Robin trace ODE with the
   newly integrated outgoing field before forming the incoming field.
   This order supplies the contraction factor and handles a general
   bounded nonsymmetric matrix. The independent proof expands the
   smoothness/uniqueness argument by differentiated systems and difference
   quotients.
8. Finite propagation keeps the coefficients away from the cone at finite
   global de Sitter time. No uniform all-angular, late-time, or exact
   nonlinear conclusion is inferred.

An earlier complete version with SHA-256
`158290f79f3159581ea0473003d78de16cc4bc1741fef50965a11f3e40698189`
was also accepted. The current version was read in full after its update.
It incorporates both optional wording refinements sent to the author:
the chosen second lapse/shift values are explicitly smooth, localized,
and in representation (11), and the exact normal-force identity makes
the second force row a consequence of the solved residuals. It also
displays the correct complete reduced source in (6a) and the correct
covariant source-conversion identity. There are no outstanding requested
corrections.

The frozen final version adds only the explicit requirement that the
chosen second lapse/shift values vanish near the initial shell corner,
with zero values admissible. Replacing that exact added text, including
its line break, by the previously reviewed text reconstructs SHA-256
`69833a130d3e234c3ed304e1aa04671d869e5679254df4ee4aae1b306691dd8f`.
`FINAL_NOTE_DELTA.json` records this exact comparison. The clarification
is correct and preserves acceptance; it changes neither checker code nor
any tested identity, so no scientific rerun is required.

The companion checker evaluates direct connection contractions, exact
second boundary derivatives, initial-gauge and subsidiary identities, the
SO(4) character decomposition, and a nonsymmetric Robin trace ODE. Its
ordinary and optimized runs and actual corruption controls are recorded
in `VERIFICATION_RECEIPT.json`. Algebra checking does not replace the
analytic existence argument or constitute external peer review.
