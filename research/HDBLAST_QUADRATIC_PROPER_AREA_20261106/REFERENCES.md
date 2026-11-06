# Sources and comparison

## Immutable HDBLAST inputs

The four copied sources are pinned at commit
`a715e8ee2963467b8988b7984a1cd0d402468f3d`; exact URLs, byte lengths,
SHA256 and Git blob identities are in [SOURCE_PINS.json](SOURCE_PINS.json).

1. [First-cut bound](inputs/FIRST_CUT_BOUND.md): actual first physical cut
   G1, j1, first-tangent support, reconstruction and joint observer/anchor
   gauge identification. It supplies no uniform nonlinear amplitude bound.
2. [Higher master-jet proof](inputs/HIGHER_JETS.md) and
   [retained report](inputs/HIGHER_JETS_NORMAL.json): actual future y and
   scaled derivative yprime on 0<=x<=3, with the same old actual-state
   input hash. These are consumed, not recomputed by the present runner.
3. [Past–future proof](inputs/PAST_FUTURE_PROOF.md), sections 2–3:
   positive dynamic-endpoint operator, common extension, nonzero velocity
   reversal, regularity and complete energy. Its basic energy gives the
   derivative transfer used here.

The [retained nonlinear quiet-window proof](https://github.com/maldonado-research/HDblast/blob/a715e8ee2963467b8988b7984a1cd0d402468f3d/research/HDBLAST_QUIET_SHELL_INTERVAL_20261105/EXACT_QUIET_WINDOW.md)
is a prior premise, distinct from the first tangent's longer quiet interval.

## Closest primary literature read

- Rohit Dhormare, *Scale-Invariant Non-Collapse from Raychaudhuri Comparison*,
  [arXiv:2610.07879v1](https://arxiv.org/html/2610.07879v1), sections 2–7,
  equations (9), (22)–(25), (28), (34), (37), Theorem 1. This recent
  comparison uses full intrinsic curvature, expansion and shear bounds
  and metric distortion to pass from transported volumes to geodesic balls.
  It does not provide the actual HDBLAST curvature projection or uniform
  amplitude derivatives. No noncollapse constants are imported here.
- Madhukrishna Chakraborty and Subenoy Chakraborty, *Curvature form of
  Raychaudhuri equation and its consequences: A geometric approach*,
  [arXiv:2409.19900v1](https://arxiv.org/html/2409.19900v1), section II,
  equations (3), (5). Its scalar shear convention uses twice sigma² for
  the tensor contraction; this note uses the contraction directly.

The determinant identity, Raychaudhuri identity, sphere moments and
endpoint energy method are established mathematics. The evaluated
registered-field application and the explicit weighted-curvature gap
are the present contribution. This focused comparison establishes neither
general mathematical originality nor journal suitability. A particular
negative source cannot be promoted to net focusing or a singularity theorem.
