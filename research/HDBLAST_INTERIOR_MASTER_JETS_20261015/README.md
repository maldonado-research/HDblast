# Certified interior master derivatives and longitudinal bulk-source bounds

15 October 2026. Conditional mathematics for the registered aligned quadrupole
in the five-dimensional doubled Einstein–scalar shell model.

This continuation supplies the missing first-master derivative certificate
through order three and uses it to bound the complete quadratic bulk residuals
in the retained first longitudinal representative. It follows the published
[coordinate shell-source checkpoint](https://github.com/maldonado-research/HDblast/tree/3dd61d552bd2e1808ab950b9afa3e5856e157b07/research/HDBLAST_COORDINATE_SHELL_SOURCES_20261015).
Earlier research and the frozen rotational-paper manuscript are preserved.

For `ell=10^-6`, `U(x,r)=w(ell*x,ell*r)`, on
`0<=x<=3, r>=0, x+r<=8`, every mixed partial of each total order obeys:

| Total scaled derivative order | Strict uniform upper bound |
| --- | ---: |
| 0 | 1.000193 |
| 1 | 16.000488 |
| 2 | 256.002081 |
| 3 | 6665.267480 |

Physical derivatives of order `j` require the factor `ell^-j`. Exact fractions,
reflected boundary-ray terms, the incoming trace derivative and initial-corner
compatibility are retained in the derivation and reports. The actual potential
now has the additional certified bound `|W_zz|<3e10` on the full audited collar.

The reconstruction and complete tensor formulas then give coarse continuous
bounds over the triangle times all of S3:

- Every Einstein quadratic-source component has absolute value below `1.2e21`
  in the stated adapted background orthonormal frames.
- The complete scalar quadratic-source coefficient has absolute value below
  `1e13`; its angular dependence is `Y^2`.

These are model-coordinate source bounds for the actual registered first
solution, conditional on the inherited analytic assumptions. They are neither
source trajectories nor physical forces. The Einstein source needs master
derivatives only through order two; the scalar source still needs order three.
An exact scalar cancellation avoids bounding terms that cancel analytically.

Read the three derivations in dependency order:

1. [Actual potential, Wzz and reconstruction coefficient jets](potential/ACTUAL_POTENTIAL_WZZ_DERIVATION.md).
2. [Reflected C2/C3 master derivatives and boundary trace closure](jets/HIGHER_JETS.md).
3. [All longitudinal Einstein components and the scalar source](bulk/LONGITUDINAL_BULK_SOURCES_AND_DERIVATIVE_ORDER.md).

[Primary-literature update](literature/TARGETED_PRIMARY_UPDATE.md) records a
bounded review and theorem applicability limits. [Internal independent review](review/REVIEW.md)
and exact-check reports accompany the source. Internal AI-assisted checks and
finite algebra tests do not replace external peer review or establish novelty.

## Reproduce

The executable calculations use Python's standard library, tested with
Python 3.12.14. No network or historical scientific program is needed. From a
fresh extraction, run:

```sh
python -I -B reproduce.py --manifest-sha256 ACTUAL_EXTERNAL_MANIFEST_HASH --output-dir /tmp/hdblast-interior-replay-new
```

The output directory must be new and outside the capsule. The accompanying
publication record and CI workflow pin the manifest externally. The runner
authenticates every source/input/golden byte, executes the explicit roster in
ordinary and optimized Python, compares complete reports, records process logs
and rechecks source integrity. Use the registered control jobs to reproduce
the intentional-error rejections. All controls raise explicit exceptions under
optimized Python. `REPLAY_ROSTER.json` defines the actual processes.

## Remaining scientific target

The first longitudinal bulk forcing is bounded. A numerical second-order
solution still requires an explicit compatible second initial constraint
correction, bulk and physical/gauge/constraint boundary data in one common
representative, and a validated coupled evolution with continuous residual and
boundary-observation bounds. The earlier conditional coefficient construction
retains its stated scope; this checkpoint does not supply its quantitative U2.
No effective nonlinear amplitude, remainder bound, observational prediction or
higher-dimensional Big Bang origin is established.

This is a GitHub research checkpoint. It does not create a new Zenodo version
or a journal submission. Code and text licensing follow the included repository
LICENSE. Public inputs are pinned by the individual checkers and capsule manifest;
private editorial correspondence and journal submission files are excluded.
