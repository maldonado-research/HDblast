# Complete quadratic tensor forcing and a shell-contact sector

14 October 2026 · conditional Einstein–scalar mathematics · GitHub checkpoint

The geometric part of the quadratic source has now been calculated. In the
inherited first-order longitudinal gauge, for the aligned quadrupole
`Y=n₁²−n₂²`, `Ψ=uY`, `χ=vY`, and the convention
`g=g₀+εh₁+ε²h₂/2`, the complete spatial sphere-TT residual is

\[
\boxed{\Pi_{\mathrm{TT}}D^2\mathscr E_{ij}
=\left(16u^2+\frac83v^2\right)\mathcal T_{ij}.}
\]

Here the `16u²` contribution is geometric and `8v²/3` is matter.
They reinforce each other wherever `(u,v)≠(0,0)`. The actual forcing
of the second-order linear equation has the opposite sign. Pure velocity
master data can have `u=v=0` on the initial slice and generate this source
later. A scalar-only second-order correction in this fixed gauge cannot
solve the nonzero tensor equation.

For an even-copy first-order master with zero first shell displacement,
the rescaled tensor amplitude has the shell condition `H_y=0`. The
companion theorem constructs its smooth forced evolution on a controlled
finite positive-cone collar extending through shell contact. It derives
the incoming characteristic, compatible corners, gauge-source conversion
and a coercive auxiliary energy estimate. The tensor initial data can be
chosen zero without changing the sourced initial constraints.

- [Complete bulk source and explicit tensor wave equation](bulk/FULL_QUADRATIC_TT_BULK_SOURCE.md)
- [Independent direct Ricci derivation](independent/FULL_QUADRATIC_TT_SOURCE.md)
- [Forced tensor contact construction and precise domain](boundary/FORCED_TENSOR_CONTACT_IBVP.md)
- [Independent boundary review](reviews/BOUNDARY_INDEPENDENT_REVIEW.md)
- [Source/boundary cross-review](reviews/SOURCE_BOUNDARY_CROSS_REVIEW.md)
- [Recent primary literature and applicability](literature/CONTACT_LITERATURE_NOTE.md)

The projection is spatial TT on the closed unit `S³`, not four-dimensional
Lorentzian TT. The source formula fixes the first gauge. A change of first
gauge requires transforming the complete bulk source as well as the boundary
data. The homogeneous Eulerian Neumann law is proved for the specified even
master; generic odd graph states require quadratic embedding corrections.

This solves a required **tensor coefficient** and its forced boundary
problem. The coupled second-order scalar/vector/gauge boundary system,
constraint propagation through contact, and an exact nonlinear contact
family remain open. The auxiliary norm is not a physical gravitational
energy or a global stability theorem. No Big Bang origin, external novelty
or submission readiness is established. Reviews are internal AI-assisted
checks, not external peer review. No new Zenodo version or APS submission.

## Reproduction

Use Python 3.12 and the pinned dependency in `requirements.txt`. Authenticate
the ZIP/source with the [external source receipt](../../hdblast/continuation/2026.10.14-quadratic-tensor-contact/FROZEN_SOURCE_RECEIPT.json)
before executing project code. Then run:

```sh
python -m pip install -r requirements.txt
python reproduce.py --manifest-sha256 EXTERNAL_MANIFEST_SHA256 --output-dir /tmp/hdblast-tensor-fresh
```

The output directory must be new or empty. The runner copies only the six
scientific checkers into fresh directories, runs each normally and with
`python -O`, compares all scientific output to the saved results, and
requires real corrupted-source/operator executions to fail. It does not
numerically evolve a spacetime or certify the analytic theorem from counts.

Four independently written bulk checkers pass 42 exact identities per mode;
the separate direct Ricci/sphere calculation passes its saved component
checks and wrong-formula controls. The boundary verifier passes 28 positive
identities and 11 wrong-formula controls per mode. Historical agent receipts
retain their original workspace paths; the package manifest pins the
distributed files. Full downloaded primary-paper texts and unaccepted
odd-graph candidates are excluded.
