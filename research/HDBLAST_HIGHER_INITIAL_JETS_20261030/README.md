# Higher initial jets and a matched outer extension

This checkpoint identifies and repairs an artifact in the initial profiles
used by the finite reconstruction. Here “curvature” means the coordinate
second radial derivative of an initial amplitude, `u_rr`; it is not the
spacetime curvature tensor. The outer endpoint `R=7/4` is the end of the
initial radial footprint, not the physical shell at `r=0`.

The archived candidate `W1` is preserved. A distinct candidate `W1_quadratic`
matches the nominal profile's second derivative when extending it beyond R
before convolution. The actual target, model, first field, sources and
boundary conditions remain the same.
The transfer report calls this same new candidate `W1_dagger`.

## A real endpoint error, and its repair

The old nominal displacement is extended by its tangent line. Its second
derivative therefore jumps from `kappa=u_nom,rr(R-)` inside to zero outside.
The even smoothing kernel averages the two sides. As its width tends to zero,
the smoothed second derivative at R tends to `kappa/2`, rather than kappa.
The actual-minus-smoothed scalar-mean error at the fixed width `2^-30` obeys

```
|e0_rr(R)_w:0| > 1.6245e33.
```

This is an actual data-error lower bound proved by two independent jet
calculations. It is not merely a large upper error allowance. Decreasing the
kernel width alone cannot remove the mismatch with the same affine extension.

The new continuation is defined using the same exact analytic nominal
endpoint values and left derivative:

```
u_ext_q(R+h)=u_nom(R)+u_nom,r(R) h+(1/2)u_nom,rr(R-) h^2,  h>=0.
```

The nominal function on `[0,R]`, reflection parity, smoothing widths and
smoothed P1 velocity are unchanged. This gives

```
|e0_rr(R)_w:0| < 1.820e27   for the new candidate.
```

The new upper bound is below the old proved lower bound divided by 800,000.
This compares actual errors for one initial component at one endpoint; it
does not establish a corresponding improvement of the evolved solution or
of the global second-derivative error.

## First-derivative control and protected calculations

The same global C1 upper certificate can be re-established for the new
candidate:

```
X(U-W1_quadratic) < 3.541e25.
```

The proof checks all 23 endpoint curvature coefficients against the existing
global second-derivative caps. The new extension obeys the same Lipschitz
bound for the first derivative; consequently the earlier displacement and
first-derivative smoothing allowances still apply. The complete initial,
free-wave, boundary-history, bulk-residual and one-update majorants are
recomputed. This is the same numerical upper certificate for a different
candidate, not a claim that their unknown true errors are equal.

An independent comparison also bounds the change of W1 in the C1 norm.
The changed initial data occupy only `[R-delta_u,R]`. Finite propagation
therefore leaves all events with `x+r<R-delta_u` unchanged. In particular,
the selected point `(9/8,0)` and the retained shell history keep the same
fields and first jets as the archived reconstruction. The published point
and harmonic-value calculations in the preceding checkpoints remain valid.

## Other higher-derivative inputs

Continuous higher-initial-error calculations give the following H-norm
upper bounds over the initial footprint:

| Initial quantity | Old affine continuation | New quadratic continuation |
| --- | ---: | ---: |
| Actual-minus-smoothed displacement `u_rr` | `<8.064e33` | `<6.219e33` |
| Actual-minus-smoothed P1 velocity `v_r` | `<9.672e24` | unchanged |

These are conservative initial-data bounds, not spacetime C2 error bounds.
The two inner cutoff joins at `65/64` and `66/64` still have second-derivative
jumps and are charged explicitly. Rounded P1 nodal values are also charged
when bounding velocity gradients.

The complete actual homogeneous coefficient time derivatives are now bounded
continuously in the same scaled-coordinate norm:

```
||A_x|| <= 3e-12,    B_x=M_x=0,    ||D_x|| < 2.354e-16.
```

A separate audit differentiates the complete collected bulk-source maps.
It finds a nonzero coefficient multiplying a fourth characteristic gauge
jet, including on a retained actual first-field cell. Principal wave-trace
relations alone do not cancel that coefficient. This is not a bound or a
nonzero value for the complete actual `F_x`, and does not rule out additional
data-specific identities or a new higher-jet estimate.

## Reproduction and limits

The `checks/`, `reference/` and `theory/` directories contain the paired
order-five initial calculation, an independent ODE-jet endpoint calculation,
the extension repair and C1 transfer, coefficient derivative checks, and
the complete-source principal-symbol audit. The workflow runs each check in
normal and optimized Python and compares the reference reports byte for byte.
Only the standard library is needed. Source functions and upstream complete
closures are authenticated before project imports and read back afterward.

Start with the [repair report](reference/REPAIR.json),
[global C1 reproof](reference/TRANSFER.json),
[higher initial bounds](reference/INITIAL.json),
[independent endpoint proof](theory/ENDPOINT_CURVATURE_NOTE.md), and
[matched-extension derivation](theory/QUADRATIC_EXTENSION_NOTE.md).
The [complete-source audit](theory/BULK_DERIVATIVE_ROUTE.md) explains the
remaining differentiated-source obligation.

The full differentiated bulk-source bound remains open. No complete C2
evolution-error gate, constraint-derivative or full Einstein residual,
self-consistent second-order solution, exact nonlinear family or evidence
of a higher-dimensional Big Bang origin is established. The repair concerns
a numerical initial extension, not an obstruction to the physical hypothesis.

Frozen manuscripts, prior checkpoints, websites and private editorial history
are preserved. No new Zenodo version or journal submission is included.
