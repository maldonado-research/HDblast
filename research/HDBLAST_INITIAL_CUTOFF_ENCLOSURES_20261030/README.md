# Tighter initial-data enclosures for the unchanged Robin reconstruction

This checkpoint tightens the conditional whole-cone C1 error bound for the
[existing Robin reconstruction](../HDBLAST_ROBIN_RECONSTRUCTION_20261030/)
from about `1.73434e27` to

```
X(U-W) <= 2.021348047763e26 < 2.0214e26.
```

The new upper certificate is over 8.57 times tighter. The mathematical target
`U`, candidate `W`, initial functions, smoothing widths, velocity table,
operators and forcing are unchanged. This improves an enclosure of the
unknown error; it does not measure a change in the true error.

## One problematic cell, ten rigorous subintervals

The previous initial radial-derivative certificate was dominated by the last
cell of a packet cutoff transition, `[263/256,33/32]`. The same authenticated
paired interval evaluator now runs on ten dyadic subintervals of that cell.
The other 191 cell certificates remain intact. An independent audit verifies
the exact cover, all three channels and all 23 components, clipping against
the inherited parent bounds, component maxima, positive Gram norms and channel
sums. The zero collar covers the rest of the initial footprint.

The [continuous refinement](theory/REFINEMENT_NOTE.md) gives these strict
initial-profile H-norm bounds:

| Actual-minus-nominal initial quantity | New upper bound |
| --- | ---: |
| Displacement | `1.121e21` |
| Radial derivative | `1.890e24` |
| Scaled-time derivative | `2.683e13` |

The initial radial upper certificate is over 800 times smaller than its earlier
value of about `1.53344e27`. It is a continuous interval bound, rather than an
estimate from sample points. All smoothing and P1 velocity interpolation errors
remain charged separately in the response estimate.

[Finite point diagnostics](theory/POINT_DIAGNOSTIC_NOTE.md) identified the cell
worth refining. They are separate evidence: nine points with small enclosures
do not prove a small continuous maximum. A nonzero difference at one point
concerns the actual-minus-**unsmoothed nominal** initial profile, not a lower
bound for the smoothed reconstruction's response error. The
[aggregation proof](theory/AGGREGATION.md) states why the new cover preserves
the bound. No new interval method or physical discovery is claimed.

## Same candidate, tighter response certificate

The [transfer calculation](theory/TRANSFER_NOTE.md) independently reconstructs
the direct C1 ledger. Only the three initial-difference caps change. Twelve
bulk/operator quantities are checked for exact equality with the preceding
report, including the whole-cone first feedback, boundary correction,
coefficient uncertainty, remaining bulk residual and affine response bound.

The norm and retained cone remain

```
X(e)=max(sup ||e||_H/3, sup ||e_x||_H, sup ||e_r||_H),
0<=x<=3/2, r>=max(0,x-5/4), x+r<=7/4.
```

The source estimate is still used in its original response units:
`X(KF)<=X(KF+LG)+X(LG)`. It is not treated as a pointwise bulk source cap.
The target and reconstruction use the same fixed actual Robin matrix and
affine shell source, so the source cancels in their boundary error equation.

All 23 candidate value/first-derivative boxes at `(x,r)=(9/8,0)` are reused
unchanged. The new target intervals lie inside the earlier ones; 92
containment checks cover the three generic quantities and shared Robin radial
quantity in every slot. They are checks on updated enclosures, not new
independent field evaluations. For the registered scalar mean, strict
outward-rounded bounds are

| Conditional formal target quantity | Lower | Upper |
| --- | ---: | ---: |
| `U_w:0` | `8.29035e30` | `8.29044e30` |
| `partial_x U_w:0` | `1.974619e32` | `1.974625e32` |
| `partial_r U_w:0`, shared Robin law | `2.154062e27` | `2.154085e27` |

These are parameter-derivative modal coefficients in a registered gauge, with
scaled time `x=t/ell`, `ell=1e-6`, and inward scaled coordinate `r`. They are
not the full angular scalar field, a physical force or an invariant observable.
Zero-M tensor rows retain the source-only radial enclosures already present
in the earlier reconstruction; no tensor-source zero is inferred.

The initial-data term now contributes about 17.5% of the conservative direct
error numerator, down from about 90.4%. The unchanged bulk-residual term now
contributes about 82.5%. The next useful target is coupled bulk feedback and
its differentiated/subsidiary residuals; further subdivision of this initial
cell alone has diminishing returns.

## Reproduction

Python 3 standard library is sufficient. From the repository root:

```sh
C=research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030
python "$C/checks/check_initial_refinement.py" --repo-root . > /tmp/initial-refined.json
python "$C/checks/check_refined_response.py" --repo-root . --refined-report /tmp/initial-refined.json > /tmp/initial-response.json
python "$C/checks/check_refinement_guards.py" --repo-root . --refined-report /tmp/initial-refined.json
python "$C/checks/check_initial_points.py" --repo-root .
```

Run with `-O` as well. Normal/optimized output must agree byte for byte and
match the four [reference reports](reference/). The workflow does this from a
fresh checkout. The unchanged scientific producer authenticates seven complete
upstream roots before scientific imports and reauthenticates afterward.
The JSON-only response checker independently audits coverage and the error
ledger. Sixteen deliberately corrupted certificates are rejected by the
coverage/aggregation auditor in both modes. These controls do not independently
re-prove the scientific interval evaluator or its inherited field assumptions.

All conclusions remain conditional on the retained first-field, background,
operator, source, zero-collar and formal-response premises. The reconstruction
still has frozen bulk forcing. This work does not complete the coupled
second-order solution, provide a numerical C2 error or pointwise Einstein
residual, establish an exact nonlinear family, or provide evidence for a
higher-dimensional cause of the Big Bang. Internal AI-assisted mathematical
reviews are not external specialist review or proof-assistant verification.
