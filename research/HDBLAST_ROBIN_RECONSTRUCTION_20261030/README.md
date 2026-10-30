# Robin reconstruction and a sharper conditional error bound

The [previous checkpoint](../HDBLAST_SMOOTH_SHELL_CONTACT_20261030/) found a
nonzero boundary defect in the first-feedback candidate. This checkpoint makes
two separate improvements:

* It constructs an integral-defined reconstruction `W` satisfying the actual
  affine Robin and homogeneous Dirichlet laws of an auxiliary wave problem,
  including the uncertain constant Robin matrix and the actual shell source.
  Its initial data are smoothed nominal data and its bulk forcing is still
  frozen at the free wave. It is not the full coupled second-order solution.
* A direct residual calculation gives the whole-cone conditional bound
  `X(U-W)<1.735e27`, over twelvefold tighter than the earlier transferred
  certificate. This is an improvement in an error **upper bound**, not a
  measurement of the true error. It avoids the old piecewise-linear bulk
  forcing interpolation charges; the boundary correction alone is not claimed
  to cause a twelvefold accuracy gain.

All bounds concern the registered formal 23-mode problem and inherit its
background, initial-field, operator, affine-source, zero-wedge and formal
completion premises. No physical observable, nonlinear family or cosmological
origin is established.

## Exact boundary reconstruction

Let `B0` be the free wave from the tighter smoothed initial data, and freeze
the bulk source at

```
f0=A0 B0_x+B0mat B0_r+D0 B0.
h=PN(B0+K f0)|shell.
Z=h-integral_0^x [M Z+G](s) ds.
W=B0+K f0+L(MZ+G).
```

`B0mat` denotes the frozen bulk matrix, distinct from the free wave `B0`.
`M` and `G` are the same fixed actual Robin matrix and actual affine shell
source as in the retained formal target `U`; they are not replaced by their
midpoints or set to zero. The complete physical bulk source is coupled before
its output is reflected. The [derivation](theory/ROBIN_RECONSTRUCTION.md) uses a
convergent Volterra resolvent, with no inverse of `M` or commutation of the bulk
matrices with parity. It proves

```
PN W_r=M PN W+G,     PD W=0,
W_xx-W_rr=f0.
```

The last identity leaves the bulk lower-order coupling frozen. Its replacement
by the actual coupling acting on `W`, and the affine bulk forcing, remain in
the error estimate. These identities concern the reconstructed auxiliary
problem, rather than a completed Einstein/constraint solve.

The numerical code encloses the integral-defined reconstruction with a degree-3
homogeneous trace series and a proved factorial tail. The affine source
convolution is retained and enclosed using uniform inherited source bounds;
its actual pointwise sign is not resolved. A finite partial sum by itself does
not satisfy the boundary law exactly. Machine-report zero residual entries
denote the analytic identity above, not interval subtraction yielding a zero
numerical residual. Continuity, differentiability and flat-corner compatibility
of the actual source are separate inherited premises; a magnitude cap alone
does not prove them. Qualitative physical-side C2 regularity under those
premises is distinct from the missing numerical C2 error certificate.

## Direct error estimate and one shell sample

For the registered diagonal Gram norm `H`, the comparison norm is

```
X(e)=max(sup ||e||_H/3, sup ||e_x||_H, sup ||e_r||_H).
```

The [direct proof](theory/DIRECT_ERROR_NOTE.md) accounts for actual initial-data
differences, smoothing, velocity interpolation, actual-minus-frozen
coefficients, remaining bulk feedback and affine forcing, then inverts a
verified contraction. The conservative calculation using the inherited
contraction gives

```
X(U-W) <= 1.734335783354e27 < 1.735e27.
```

The [separate triangle-transfer check](theory/TRANSFER_NOTE.md) gives the looser
valid bound `X(U-W)<2.185e28` and bounds the global boundary correction by
`X(W-V)<6.344e24`. It is a cross-check on the reconstruction, not the source of
the tighter certificate. The affine forcing bound inherited from earlier work
is an `X(KF+LG)` response bound; the direct proof explicitly uses
`X(KF)<=X(KF+LG)+X(LG)`. It does not treat that number as a pointwise bulk
source cap.

About 90% of the numerator of this conservative direct error budget comes
from the initial-data term. Tightening the initial radial-derivative enclosure
is therefore another concrete target alongside completing the bulk evolution.

At `(x,r)=(9/8,0)`, the [full interval evaluation](theory/NUMERIC_NOTE.md) and an
[independent residual-resolvent enclosure](theory/COARSE_NOTE.md) agree in all
23 rows. The new scalar-mean trace correction is strictly positive,
`7.45e20 < (W-V)_w:0 < 1.526e21`. Combining the reconstruction with the direct
global certificate gives these outward-rounded strict conditional target bounds:

| Registered formal scalar-mean quantity | Lower | Upper |
| --- | ---: | ---: |
| `U_w:0` | `8.2901e30` | `8.2907e30` |
| `partial_x U_w:0` | `1.97460e32` | `1.97465e32` |
| `partial_r U_w:0`, shared Robin law | `2.15401e27` | `2.15414e27` |

The direct C1 derivative interval already fixes the radial sign; the shared
Robin law narrows it further. Because `U` and `W` share the same actual `M,G`,
`PN(U-W)_r=M PN(U-W)` on the shell. `G` cancels in that error identity. Zero
Robin rows retain their sharper source-only radial bounds, including the two
tensor slots. Broad value/error tables do not supersede earlier specialized
tensor-mode certificates.

These are second-order modal coefficients at one event, in a fixed gauge and
normalization, with `x=t/ell`, `ell=1e-6`, and inward scaled coordinate `r`.
They do not determine the full angular scalar's sign, a physical force,
proper-time observable or finite-amplitude signal.

## Reproduction

Python 3 standard library is sufficient. From the repository root:

```sh
C=research/HDBLAST_ROBIN_RECONSTRUCTION_20261030
python "$C/checks/check_robin_resolvent.py"
python "$C/checks/check_robin_transfer.py" --repo-root . > /tmp/robin-transfer.json
python "$C/checks/check_robin_direct_error.py" --repo-root . --transfer-report /tmp/robin-transfer.json > /tmp/robin-direct.json
python "$C/checks/check_robin_resummed.py" --repo-root . > /tmp/robin-resummed.json
python "$C/checks/check_robin_coarse.py" --repo-root . > /tmp/robin-coarse.json
python "$C/checks/check_posterior.py" --repo-root . --resummed /tmp/robin-resummed.json --coarse /tmp/robin-coarse.json --direct /tmp/robin-direct.json --transfer /tmp/robin-transfer.json
```

Normal and optimized Python outputs must match byte for byte and reproduce the
five [reference reports](reference/). The workflow performs both checks. The
polynomial audit has 44 named exact groups; the posterior comparison checks
156 cross-method intervals. Input hashes, block Gram operator bounds, source
units and actual-versus-frozen coefficient uncertainty are checked explicitly.
The nominal interval producer authenticates its upstream roots before scientific
imports and reauthenticates them afterward. Other arithmetic programs use
authenticated JSON without importing scientific code. Underlying field
certificates remain inherited assumptions. These internal AI-assisted reviews
are not external specialist or proof-assistant verification.

[Classical primary-source context](RELATED.md) makes no novelty claim. The next
task is to replace the frozen bulk source with compatible coupled evolution and
certify differentiated bulk and subsidiary residuals. A complete second-order
solution, quantitative C2 control, exact nonlinear family and higher-dimensional
Big Bang origin remain unestablished.
