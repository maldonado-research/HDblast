# Continuous coordinate shell sources for HDBLAST

This checkpoint reconstructs the actual first metric boundary trace for the
registered admissible quadrupole and encloses the complete one-copy coordinate
scalar/isotropic Israel source coefficients over 4353 continuous time cells.
It also evaluates the retained degree-four Robin forcing, proper-shell source
norms and the distinction between integrated smallness and signs at scalar
zeros. All statements are conditional on the registered analytic inputs.

- [Derivation, bounds, conventions and unresolved work](COORDINATE_SHELL_SOURCES.md)
- [Exact report](computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json) · [Whole-cell table](computed/COORDINATE_SHELL_SOURCE_CELLS.csv)
- [Astra independent derivation and checks](astra/) · [Sol all-cell reconstruction](sol/) · [Independent scope/bulk-source review](referee/)
- [Interior value contraction and full tensor source recipe](referee/FULL_BULK_SOURCE_AND_INTERIOR_NEXT_STEP.md)
- [Primary-literature check](literature/LITERATURE_CHECK.md) · [Input identities](inputs/INPUT_PINS.json)
- [Reproduction instructions](REPRODUCE.md) · [Registered process/control roster](REPLAY_ROSTER.json)

In the retained longitudinal representative,
`Sp=-4p u^2+2sigma''uv-2v^2` and `SK=sigma''v^2/6-12a u^2`.
The scalar material summand dominates the stated L2 norm, but both complete
coordinate coefficients have at least two distinct zeros inside the certified
central scalar-zero bracket. These are coordinate source coefficients; they
do not determine a physical force or the sign of the second-order response.

The earlier conditional construction is preserved. The full numerical U2,
second-order stress, nonlinear remainder and cosmological interpretation remain
open. Internal AI-assisted checks are not external peer review. No Zenodo
version or APS submission is made for this checkpoint.
