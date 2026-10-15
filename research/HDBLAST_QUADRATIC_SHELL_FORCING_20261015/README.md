# Quadratic material tension contribution at the HDBLAST shell

For one registered aligned quadrupole, this checkpoint bounds the explicit
quadratic tension term in the second-order material shell laws using the
previously certified actual first variation. It supplies continuous harmonic
coefficient intervals, a nonzero monopole tail, and proper-time source norms.

- [Derivation, conventions, bounds and unresolved terms](QUADRATIC_MATERIAL_TENSION_FORCING.md)
- [Exact enclosure and source norms](computed/QUADRATIC_FORCING_ENCLOSURE.json)
- [All 4353 continuous coefficient cells](computed/QUADRATIC_FORCING_CELLS.csv)
- [Source figure](computed/QUADRATIC_MATERIAL_FORCING.png) · [PDF figure](computed/QUADRATIC_MATERIAL_FORCING.pdf)
- Independent [Astra derivation](astra/ASTRA_QUADRATIC_MATERIAL_BOUNDARY_FORCING.md),
  [Sol derivation and complete-cell audit](sol/QUADRATIC_MATERIAL_FORCING.md),
  [registered model audit](archive/REGISTERED_SHELL_COEFFICIENT_AUDIT.md),
  and [referee review](referee/FINAL_QUADRATIC_SOURCE_REVIEW.md)
- [Bounded primary-literature check](literature/LITERATURE_CHECK.md)
- [Authenticated prior public inputs](inputs/INPUT_PINS.json)
- [Fresh reproduction instructions](REPRODUCE.md) · [Exact process/control roster](REPLAY_ROSTER.json)

In the outward-normal sum convention and derivative-coefficient Taylor
convention, the material scalar row is
`sum(p_i,2) + sigma'' psi2 = -4 psi1^2`. For identical copies the right side
is `-2 psi1^2`. Its scalar angular degrees are 0, 2 and 4. The squared
proper-time L2 norm of this explicit term lies strictly between
`5.778e21` and `8.644e21` in the registered model units.

This is the explicit constitutive contribution. The complete source written
in bulk-coordinate unknowns also contains normal, metric and embedding terms.
The full second-order response, total second-order stress and a nonlinear
remainder are not evaluated here. The earlier conditional construction of the
complete second Taylor coefficient is preserved. Finite arithmetic and internal
AI-assisted review do not establish the analytic PDE assumptions or replace
external peer review. No observational calibration, novelty claim or evidence
for a higher-dimensional origin of the Big Bang is supplied.

The separate dated publication delivery records real replay, GitHub CI, merge
and live-readback facts when completed. No Zenodo version or APS submission is
part of this checkpoint.
