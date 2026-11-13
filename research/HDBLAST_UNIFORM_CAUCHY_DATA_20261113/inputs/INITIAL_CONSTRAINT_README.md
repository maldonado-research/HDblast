# Quantified nonlinear initial-constraint family

13 November 2026. This checkpoint continues the [all-mode weighted inverse](../HDBLAST_ALL_MODE_INVERSE_20261113/README.md) merged in PR #139. Under the retained smooth-background and selected weighted-isomorphism premises, it evaluates a uniform nonlinear constraint estimate and an explicit interval for the registered initial family:

```
||D^2 G|| <= K = 10^6400,
|epsilon| <= a = 10^-6402,
||u(epsilon)||_X <= r = 10^-6404.
```

Here `G=L^-1 F` uses the same fixed background inverse, `F=w^-2 Chat`, and the product norm is `|delta epsilon|+||delta u||_X`. The potential space is the reflection-invariant Hilbert direct sum `X=(dot H16(q,w) ⊕ dot H14(q,w))^Gamma`, with `Gamma=(Z2)^4`. The constraint source uses the fixed background density: `Chat=(H,M_density/rho0)`. The scalar and scalar-momentum-density data remain affine in epsilon; the selected corrections act on the metric and gravitational momentum density. The original seed includes its nonzero velocity.

The resulting constraint family is unique in the stated correction ball and is C2 on a neighborhood of this parameter interval, with `u(0)=u'(0)=0`, `||u'||<=1` and `||u''||<=8K`. The metric remains positive. Every symmetry-allowed angular mode is included; the nonlinear completion is not truncated to the seed's finite angular support.

**The interval is a deliberately conservative mathematical proof threshold.** It is not a measured or fitted cosmological amplitude, a useful physical scale, evidence for a higher-dimensional event, or a new fundamental law. The construction remains conditional on the retained background premises. No quantitative nonlinear evolution tube, observer continuation or observable remainder is established here.

## What closes the initial-family estimate

| Input | Evaluated bound | Proof |
|---|---:|---|
| All-mode background inverse | `10^5000` | [Prior inverse](../HDBLAST_ALL_MODE_INVERSE_20261113/README.md) |
| Physical weighted correction operators | E0 metric `<10^30`; E1 momentum tensor `<10^27` | [Correction maps](PHYSICAL_CORRECTION_MAPS.md) |
| Weighted pointwise embedding and products | `10^7`, `10^12`, `10^17` | [Embedding and products](WEIGHTED_EMBEDDING_PRODUCTS.md) |
| Derivatives of full source normalization | `10^18` | [Normalization proof](WEIGHTED_EMBEDDING_PRODUCTS.md) |
| Four-sector actual seed | combined graded norm `<10^180`, charged as `10^200` | [Seed envelope](SEED_ENVELOPE.md) |
| Whole-ball nonlinear Hessian and contraction | `K=10^6400` | [Nonlinear constraint ball](NONLINEAR_CONSTRAINT_BALL.md) |

The physical graded norms include the additional decay of the actual correction maps. The potential space with decaying weight is not itself a same-weight algebra. The proof uses unitization for the identity metric, finite graded products, explicit background multipliers, and uniform inverse/determinant composition on the positive ball. The density momentum constraint is bilinear in the canonical fields. The Hamiltonian proof retains every inverse-metric, curvature, determinant, scalar-gradient, momentum and potential term.

## Reproduction and interpretation

[Manifest and replay protocol](REPRODUCIBILITY.md) explains the exact source pins, closed file roster, external manifest hash, integer/rational arithmetic, symbolic component checks and refreshed fault controls. The executable checks audit finite identities and bounds used by the written analytic proof; they are not a formal proof assistant verification of the continuum theorem. An independent mathematical audit is recorded in [AUDIT.md](AUDIT.md).

[LITERATURE.md](LITERATURE.md) identifies the primary methods and collection limits. This is a quantified application of established geometric-analysis and contraction methods, with model-specific estimates. No claim of mathematical novelty or exhaustive literature coverage is made. [CONTINUATION.md](CONTINUATION.md) identifies the next missing map before a numerical nonlinear evolution calculation. GitHub publication receipts and live-page verification belong in the external handoff, after merge. This checkpoint creates no new Zenodo edition.
