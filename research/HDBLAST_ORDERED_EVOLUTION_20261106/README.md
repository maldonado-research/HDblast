# Ordered background evolution and quadratic forcing

6 November 2026. Conditional mathematical checkpoint following the merged
initial-transport checkpoint (PR #134).

This evaluates the complete lower-order maps of the zero-state background
linearized sixteen-field wave, matches the coordinate and bundle equations,
and gives a finite-Sobolev H12 bound for its full quadratic second source.
The selected initial data and retained background are the same as before.

With C=10^150000, the complete ordered background maps satisfy
`l_t<=2C`, `l_w<=9C`, `l_0<=40C`. For `0<=t<=13/8000000`, the propagator
satisfies `G_bg(t,0)<exp(10^150000)`. Thus the specified first linear wave has

```
||V1(t)||H1 + ||V1,t(t)||L2 < 10^5203 exp(10^150000).
```

This is an explicit, very conservative bound for a low-order evolved norm,
conditional on the retained background linear solution and boundary setup.
It supplies no evolved H14 bound or finite-amplitude radius.

A second proof charges every term in the quadratic source, including
`V1 D²V1` and its factor two. Using the fixed product covariant norms, put
`N14(t)=||V1(t)||H14+||V1,t(t)||H13`. Then

```
||F2(t)||H12 <= 10^240 (10^152000)^2 N14(t)^2,
||F2(0)||H12 < 10^315000.
```

The map estimate applies when the displayed high norms exist; their uniform
evolved numerical bound remains open. The source estimate costs two spatial
derivatives of the first wave. It does not close the high-order second-wave
or proper-time actual-area quadrupole remainder.

Read [ORDERED_BACKGROUND.md](ORDERED_BACKGROUND.md) and
[QUADRATIC_FORCING.md](QUADRATIC_FORCING.md) for exact norms, reconstruction,
operator order, hypotheses and constants. [REPRODUCE.md](REPRODUCE.md)
describes the exact finite replay and corruption controls.

Eight immutable public source snapshots are authenticated in
[SOURCE_PINS.json](SOURCE_PINS.json). [LITERATURE.md](LITERATURE.md) and
[PHYSICAL_NEIGHBORS.md](PHYSICAL_NEIGHBORS.md) compare eight inspected primary
full texts, including recent theorem revisions and distinct bubble/brane
models. This is targeted research, not an exhaustive literature census.

The analytic background, canonical construction and retained evolution
premises remain conditional inputs. Internal review and exact arithmetic
checks are not external peer review or physical confirmation. Constants use
the fixed model normalization without SI calibration. [CONTINUATION.md](CONTINUATION.md)
records the remaining mathematical work.
