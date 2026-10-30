# Independent arithmetic audit of the local forcing defect

At `(x,r)=(9/8,0)`, the forcing defect of the prior boundary-compatible
reconstruction is

```
r0 = La W0 + F - L0 B_s.
```

The main evaluator uses the direct difference and the decomposition
`(La-L0)B_s+La(W0-B_s)+F`. The independent audit imports neither that evaluator
nor its arithmetic helpers. It forms a third exact identity,

```
r0 = (La-L0)W0 + L0(W0-B_s) + F,
```

with signed rational interval arithmetic and independently indexed input
tables. The earlier free-wave and reconstruction enclosures are authenticated.
The coefficient driver is also authenticated before import; its constructor
authenticates local dependencies and all seven complete scientific roots.
The audit independently reproduces every component of `L0 B_s` and `La W0`,
then compares its complete residual intervals with the main report. All 23
intervals overlap; both calculations certify the same 17 nonzero rows.

The scalar-mean interval is approximately

```
[-9.416751844939703e22, -9.416550797141028e22],
```

and in particular lies strictly between `-9.417e22` and `-9.416e22`.
The remaining six rows have unresolved signs under these enclosures. They
have not been set to zero.

## Source convention

The source coefficients are from the authenticated complete actual
common-gauge geometric-source report, with its complete source-map premise.
Its global cover contains the sample. Its source values are already in the
reduced wave equation normalization:

```
metric_forcing = -2 ell^2 rho^2 QEred,
scalar_forcing = ell^2 rho^2 Qphi.
```

These are the conventions in the retained operator's `bulk_forcing`. The
component caps are already scaled forcing caps; multiplying by these factors
again would change the equation. The older operator's matter-only caller
warning does not replace the complete geometric source-map premise of the
actual source report.

The independent audit reconstructs the positive comparison Gram weights from
the full23 labels. For seven source blocks, it checks all 93 Gram entries,
including multiplicity two for mixed metric and vector blocks, and rebuilds
the summed source norm square. All 23 source coefficients are mapped exactly
once. This proves agreement of the source comparison norm and the registered
`H` convention. It is not a new proof of the inherited source geometry.

The genuine uniform pointwise `F` cap is used here, not the earlier response
bound `X(KF+LG)`. The affine source is magnitude-enclosed on a uniform cover;
its exact value and sign have not been newly resolved at this point.

## Reproduction and scope

Run `check_point_residual_audit.py --repo-root PATH --residual-report FILE` in
normal and optimized Python; JSON outputs must match byte for byte. By
default the residual report is `../reference/POINT_RESIDUAL.json`. Optional
`--robin-capsule` and `--contact-capsule` arguments select separately located
pinned earlier tables. The supplied primary residual report is authenticated
before comparison. The same coefficient assembler and earlier scientific
certificates are shared: the independence is in arithmetic decomposition,
indexing and source-norm verification, not independent field evolution.

This is a conditional point enclosure of a reduced-evolution forcing defect.
The principal forcing of `W0` is known from its analytic construction, so its
defect can be computed from first jets without numerically differentiating
the field twice. It does not supply quantitative second-derivative errors, a
full pointwise Einstein/constraint residual, an evaluated bulk correction,
an exact nonlinear spacetime or a physical observable.
