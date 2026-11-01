# Conditional intrinsic shell geometry and the required second mean response

This checkpoint moves from an interior bulk curvature projection to intrinsic
geometry of the physical timelike shell. It gives continuous, outward rational
bounds for first-field quadratic contributions to the shell expansion budget,
under the retained analytic and actual-first-field premises.

For background shell time `x in [31/25,63/50]`, with `ell=10^-6`, the completed
calculation gives the following deliberately loose outward bounds:

| Second parameter derivative | Continuous enclosure |
|---|---|
| Intrinsic curvature defect `C=<R3>-6/Rvolume^2` | `[3e-24,2e-22]` |
| Expansion/shear correction `Q_D=(2/3)Var(theta)-<|shear|^2>` | `[1e-9,2e-8]` |
| Mean intrinsic scalar gradient norm `<(D phi)^2>` | `< -4e17` |
| Ratio of spatial to temporal scalar-gradient squares | `< 1e-12` |

These numbers use the model normalization. They have no supplied SI scale.
Every coefficient is a second parameter derivative, rather than the coefficient
after dividing by two in a Taylor expansion. The very different numerical sizes
do not by themselves identify a physical energy or a cosmological effect.

The curvature defect quantifies why the shell's volume radius cannot simply be
substituted into the spatial curvature formula for an exactly round FRW space.
The positive `Q_D` coefficient means the expansion-variance contribution exceeds
the shear-square contribution at this order. Neither fact proves accelerated
expansion, heating, or a hot Big Bang history.

The scalar-gradient bound also implies that its second-order Lorentzian norm is
negative on more than **99.9999%** of background angular volume throughout the
interval. For each admitted exact family, the scalar gradient is timelike on
the same strict fraction of actual proper volume for sufficiently small nonzero
perturbations. The allowed
amplitude interval is uncomputed and need not be uniform over nonlinear
completions. Both angular signs occur in the coefficient; no all-angle claim is
made. This concerns the bulk scalar restricted to the shell, not created matter.

Read the [assembled statement and proof](DERIVATION.md) and the component
derivations linked there. The scalar-gradient identity itself was already
published in the October 16 work. Standard Gauss, averaging and normalized
scalar-curvature machinery are attributed in [the references](REFERENCES.md).
The contribution here is their application to the retained exact-family tangent,
the new continuous actual-field enclosures, and the explicit remaining mean data.

## How the stronger shell bounds were obtained

The old master comparison vanishes at the shell. Its global error bound alone
leaves the actual shell metric amplitude in an interval containing zero. A
different valid consequence of the same physical equations supplies a forced
intrinsic trace ODE. Integrating its actual scalar source on a frozen whole-cell
cover gives the new asymmetric metric and clock-history bounds. The calculation
retains all inherited source errors and the complete lower-order feedback.
It adds no boundary condition and assumes no exact precontact-zero shortcut.

The strict coefficients depend only on the common first field. Unspecified
second scalar/metric data cannot be silently replaced by an older formal
completion. The exact volume equation contains a second scalar mean with a
provably nonzero coefficient, so the full second volume response remains open.

An independently reviewed initial-constraint calculation additionally proves a
strict negative, fixed-gauge weighted monopole moment for this packet. It holds
for all admissible weight choices, but does not select a unique radial mean
profile or an energy sign. The defining function and weight parameters in the
earlier abstract construction were not numerically instantiated.

## Reproduce and assess the evidence

Use Python 3.12.14 and the pinned requirements from a complete source checkout:

```sh
python -m pip install -r research/HDBLAST_INTRINSIC_SHELL_BUDGET_20261101/requirements.txt
python -I -S -B research/HDBLAST_INTRINSIC_SHELL_BUDGET_20261101/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-shell-budget-fresh
```

The output path must be fresh and outside the source tree. The root runner
authenticates the closed new capsule and retained dependencies before scientific
imports and afterward. It distinguishes canonical report pairs, actual rejected
subprocesses, in-process formula checks and successful relocation controls.
Source authentication is not a new proof of the retained analytic PDE results.
Independent rational reconstruction shares the recorded actual-jet inputs; it
does not independently solve that first PDE. Internal AI-assisted review is not
external specialist or proof-assistant verification.

No frozen manuscript or older capsule is changed. No journal submission or
Zenodo version accompanies this checkpoint. The higher-dimensional origin of
the Big Bang remains unestablished; [the next calculation](NEXT_STEP.md) is stated
around the specific data still missing.
